/*
  Expression Evaluator Library (NS-EEL) v2
  Copyright (C) 2004-2013 Cockos Incorporated
  Copyright (C) 1999-2003 Nullsoft, Inc.
  nseel-compiler.c
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include "eelCommon.h"
#include "glue_port.h"
// C string manipulation utilities -- [v]snprintf for Win32, also snprintf_append, lstrcatn, etc
#if defined(_WIN32) && defined(_MSC_VER)
// provide snprintf()/vsnprintf() for win32 -- note that these have no way of knowing
// what the amount written was, code should(must) be written to not depend on this.
#ifdef vsnprintf
#undef vsnprintf
#endif
#define vsnprintf WDL_vsnprintf
#endif
// use wdlcstring.h's lstrcpyn_safe rather than the real lstrcpyn.
#ifdef _WIN32
  #ifdef lstrcpyn
  #undef lstrcpyn
  #endif
  #define lstrcpyn lstrcpyn_safe
#endif
  #if defined(_WIN32) && defined(_MSC_VER)
    void WDL_vsnprintf(char *o, size_t count, const char *format, va_list args)
    {
      if (count>0)
      {
        int rv;
        o[0]=0;
        rv=_vsnprintf(o,count,format,args); // returns -1  if over, and does not null terminate, ugh
        if (rv < 0 || rv>=(int)count-1) o[count-1]=0;
      }
    }
    void WDL_VARARG_WARN(printf,3,4) WDL_snprintf(char *o, size_t count, const char *format, ...)
    {
      if (count>0)
      {
        int rv;
        va_list va;
        va_start(va,format);
        o[0]=0;
        rv=_vsnprintf(o,count,format,va); // returns -1  if over, and does not null terminate, ugh
        va_end(va);
        if (rv < 0 || rv>=(int)count-1) o[count-1]=0;
      }
    }
  #endif
#if defined(__ANDROID__) || defined(__APPLE__) || defined(__linux__)
#define WDL_snprintf snprintf
#endif
  static void lstrcpyn_safe(char *o, const char *in, int count)
  {
    if (count>0)
    {
      while (--count>0 && *in) *o++ = *in++;
      *o=0;
    }
  }
  static void lstrcatn(char *o, const char *in, int count)
  {
    if (count>0)
    {
      while (*o) { if (--count < 1) return; o++; }
      while (--count>0 && *in) *o++ = *in++;
      *o=0;
    }
  }
  static void WDL_VARARG_WARN(printf,3,4) snprintf_append(char *o, int count, const char *format, ...)
  {
    if (count>0)
    {
      va_list va;
      while (*o) { if (--count < 1) return; o++; }
      va_start(va,format);
      vsnprintf(o,count,format,va);
      va_end(va);
    }
  }
#define NSEEL_VARS_MALLOC_CHUNKSIZE 8
#define RET_MINUS1_FAIL(x) return -1;
#define MIN_COMPUTABLE_SIZE 32 // always use at least this big of a temp storage table (and reset the temp ptr when it goes past this boundary)
#define COMPUTABLE_EXTRA_SPACE 16 // safety buffer, if EEL_VALIDATE_WORKTABLE_USE set, used for magic-value-checking
/*
  P1 is rightmost parameter
  P2 is second rightmost, if any
  P3 is third rightmost, if any
  registers on x86 are  (RAX etc on x86-64)
    P1(ret) EAX
    P2 EDI
    P3 ECX
    WTP RSI
    x86_64: r12 is a pointer to ram_state.blocks
    x86_64: r13 is a pointer to closenessfactor
  registers on PPC are:
    P1(ret) r3
    P2 r14 
    P3 r15
    WTP r16 (r17 has the original value)
    r13 is a pointer to ram_state.blocks
    ppc uses f31 and f30 and others for certain constants
  */
// used by //#eel-no-optimize:xxx, in 0
#define OPTFLAG_NO_FPSTACK 2
#define OPTFLAG_NO_INLINEFUNC 4
#define MAX_SUB_NAMESPACES 32
typedef struct
{
  const char *namespacePathToThis;
  const char *subParmInfo[MAX_SUB_NAMESPACES];
} namespaceInformation;
static int nseel_evallib_stats[5]; // source bytes, static code bytes, call code bytes, data bytes, segments
int *NSEEL_getstats()
{
  return nseel_evallib_stats;
}
static int findLineNumber(const char *exp, int byteoffs)
{
  int lc=0;
  while (byteoffs-->0 && *exp) if (*exp++ =='\n') lc++;
  return lc;
}
static void *__newBlock(llBlock **start,int size);
#define OPCODE_IS_TRIVIAL(x) ((x)->opcodeType <= OPCODETYPE_VARPTRPTR)
enum {
  OPCODETYPE_DIRECTVALUE=0,
  OPCODETYPE_DIRECTVALUE_TEMPSTRING, // like directvalue, but will generate a new tempstring value on generate
  OPCODETYPE_VALUE_FROM_NAMESPACENAME, // this.* or namespace.* are encoded this way
  OPCODETYPE_VARPTR,
  OPCODETYPE_VARPTRPTR,
  OPCODETYPE_FUNC1,
  OPCODETYPE_FUNC2,
  OPCODETYPE_FUNC3,
  OPCODETYPE_FUNCX,
  OPCODETYPE_MOREPARAMS,
  OPCODETYPE_INVALID,
};
struct opcodeRec
{
 int opcodeType; 
 int fntype;
 void *fn;
 union {
   struct opcodeRec *parms[3];
   struct {
     double directValue;
     EEL_F *valuePtr; // if direct value, valuePtr can be cached
   } dv;
 } parms;
 int namespaceidx;
 // OPCODETYPE_VALUE_FROM_NAMESPACENAME (relname is either empty or blah)
 // OPCODETYPE_VARPTR if it represents a global variable, will be nonempty
 // OPCODETYPE_FUNC* with fntype=FUNCTYPE_EELFUNC
 const char *relname;
};
static void *newTmpBlock(compileContext *ctx, int size)
{
  const int align = 8;
  const int a1=align-1;
  char *p=(char*)__newBlock(&ctx->tmpblocks_head,size+a1);
  return p+((align-(((INT_PTR)p)&a1))&a1);
}
static void *__newBlock_align(compileContext *ctx, int size, int align, int isForCode) 
{
  const int a1=align-1;
  char *p=(char*)__newBlock(
                            (                            
                             isForCode < 0 ? (isForCode == -2 ? &ctx->pblocks : &ctx->tmpblocks_head) : 
                             isForCode > 0 ? &ctx->blocks_head : 
                             &ctx->blocks_head_data) ,size+a1);
  return p+((align-(((INT_PTR)p)&a1))&a1);
}
static opcodeRec *newOpCode(compileContext *ctx, const char *str, int opType)
{
  const size_t strszfull = str ? strlen(str) : 0;
  const size_t str_sz = wdl_min(NSEEL_MAX_VARIABLE_NAMELEN, strszfull);
  opcodeRec *rec = (opcodeRec*)__newBlock_align(ctx,
                         (int) (sizeof(opcodeRec) + (str_sz>0 ? str_sz+1 : 0)),
                         8, ctx->isSharedFunctions ? 0 : -1); 
  if (rec)
  {
    memset(rec,0,sizeof(*rec));
    rec->opcodeType = opType;
    if (str_sz > 0) 
    {
      char *p = (char *)(rec+1);
      memcpy(p,str,str_sz);
      p[str_sz]=0;
      rec->relname = p;
    }
    else
    {
      rec->relname = "";
    }
  }
  return rec;
}
#define newCodeBlock(x,a) __newBlock_align(ctx,x,a,1)
#define newDataBlock(x,a) __newBlock_align(ctx,x,a,0)
#define newCtxDataBlock(x,a) __newBlock_align(ctx,x,a,-2)
static void freeBlocks(llBlock **start);
#ifndef DECL_ASMFUNC
#define DECL_ASMFUNC(x)         \
  void nseel_asm_##x(void);        \
  void nseel_asm_##x##_end(void);
void _asm_megabuf(void);
void _asm_megabuf_end(void);
#endif
  DECL_ASMFUNC(booltofp)
  DECL_ASMFUNC(fptobool)
  DECL_ASMFUNC(fptobool_rev)
  DECL_ASMFUNC(sin)
  DECL_ASMFUNC(cos)
  DECL_ASMFUNC(tan)
  DECL_ASMFUNC(1pdd)
  DECL_ASMFUNC(2pdd)
  DECL_ASMFUNC(2pdds)
  DECL_ASMFUNC(1pp)
  DECL_ASMFUNC(2pp)
  DECL_ASMFUNC(sqr)
  DECL_ASMFUNC(sqrt)
  DECL_ASMFUNC(log)
  DECL_ASMFUNC(log10)
  DECL_ASMFUNC(abs)
  DECL_ASMFUNC(min)
  DECL_ASMFUNC(max)
  DECL_ASMFUNC(min_fp)
  DECL_ASMFUNC(max_fp)
  DECL_ASMFUNC(sig)
  DECL_ASMFUNC(sign)
  DECL_ASMFUNC(band)
  DECL_ASMFUNC(bor)
  DECL_ASMFUNC(bnot)
  DECL_ASMFUNC(bnotnot)
  DECL_ASMFUNC(if)
  DECL_ASMFUNC(fcall)
  DECL_ASMFUNC(repeat)
  DECL_ASMFUNC(repeatwhile)
  DECL_ASMFUNC(equal)
  DECL_ASMFUNC(equal_exact)
  DECL_ASMFUNC(notequal_exact)
  DECL_ASMFUNC(notequal)
  DECL_ASMFUNC(below)
  DECL_ASMFUNC(above)
  DECL_ASMFUNC(beloweq)
  DECL_ASMFUNC(aboveeq)
  DECL_ASMFUNC(assign)
  DECL_ASMFUNC(assign_fromfp)
  DECL_ASMFUNC(assign_fast)
  DECL_ASMFUNC(assign_fast_fromfp)
  DECL_ASMFUNC(add)
  DECL_ASMFUNC(sub)
  DECL_ASMFUNC(add_op)
  DECL_ASMFUNC(sub_op)
  DECL_ASMFUNC(add_op_fast)
  DECL_ASMFUNC(sub_op_fast)
  DECL_ASMFUNC(mul)
  DECL_ASMFUNC(div)
  DECL_ASMFUNC(mul_op)
  DECL_ASMFUNC(div_op)
  DECL_ASMFUNC(mul_op_fast)
  DECL_ASMFUNC(div_op_fast)
  DECL_ASMFUNC(mod)
  DECL_ASMFUNC(shl)
  DECL_ASMFUNC(shr)
  DECL_ASMFUNC(mod_op)
  DECL_ASMFUNC(or)
  DECL_ASMFUNC(or0)
  DECL_ASMFUNC(xor)
  DECL_ASMFUNC(xor_op)
  DECL_ASMFUNC(and)
  DECL_ASMFUNC(or_op)
  DECL_ASMFUNC(and_op)
  DECL_ASMFUNC(uplus)
  DECL_ASMFUNC(uminus)
  DECL_ASMFUNC(invsqrt)
  DECL_ASMFUNC(dbg_getstackptr)
  DECL_ASMFUNC(stack_push)
  DECL_ASMFUNC(stack_pop)
  DECL_ASMFUNC(stack_pop_fast) // just returns value, doesn't mod param
  DECL_ASMFUNC(stack_peek)
  DECL_ASMFUNC(stack_peek_int)
  DECL_ASMFUNC(stack_peek_top)
  DECL_ASMFUNC(stack_exch)
static void *NSEEL_PProc_Stack(void *data, int data_size, compileContext *ctx)
{
  codeHandleType *ch=ctx->tmpCodeHandle;
  if (data_size>0) 
  {
    UINT_PTR m1=(UINT_PTR)(NSEEL_STACK_SIZE * sizeof(EEL_F) - 1);
    UINT_PTR stackptr = ((UINT_PTR) (&ch->stack));
    ch->want_stack=1;
    if (!ch->stack) ch->stack = newDataBlock(NSEEL_STACK_SIZE*sizeof(EEL_F),NSEEL_STACK_SIZE*sizeof(EEL_F));
    data=EEL_GLUE_set_immediate(data, stackptr);
    data=EEL_GLUE_set_immediate(data, m1); // and
    data=EEL_GLUE_set_immediate(data, ((UINT_PTR)ch->stack&~m1)); //or
  }
  return data;
}
static void *NSEEL_PProc_Stack_PeekInt(void *data, int data_size, compileContext *ctx, INT_PTR offs)
{
  codeHandleType *ch=ctx->tmpCodeHandle;
  if (data_size>0) 
  {
    UINT_PTR m1=(UINT_PTR)(NSEEL_STACK_SIZE * sizeof(EEL_F) - 1);
    UINT_PTR stackptr = ((UINT_PTR) (&ch->stack));
    ch->want_stack=1;
    if (!ch->stack) ch->stack = newDataBlock(NSEEL_STACK_SIZE*sizeof(EEL_F),NSEEL_STACK_SIZE*sizeof(EEL_F));
    data=EEL_GLUE_set_immediate(data, stackptr);
    data=EEL_GLUE_set_immediate(data, offs);
    data=EEL_GLUE_set_immediate(data, m1); // and
    data=EEL_GLUE_set_immediate(data, ((UINT_PTR)ch->stack&~m1)); //or
  }
  return data;
}
static void *NSEEL_PProc_Stack_PeekTop(void *data, int data_size, compileContext *ctx)
{
  codeHandleType *ch=ctx->tmpCodeHandle;
  if (data_size>0) 
  {
    UINT_PTR stackptr = ((UINT_PTR) (&ch->stack));
    ch->want_stack=1;
    if (!ch->stack) ch->stack = newDataBlock(NSEEL_STACK_SIZE*sizeof(EEL_F),NSEEL_STACK_SIZE*sizeof(EEL_F));
    data=EEL_GLUE_set_immediate(data, stackptr);
  }
  return data;
}
#if defined(_MSC_VER) && _MSC_VER >= 1400
static double __floor(double a) { return floor(a); }
static double __ceil(double a) { return ceil(a); }
#define floor __floor
#define ceil __ceil
#endif
#define FUNCTIONTYPE_PARAMETERCOUNTMASK 0xff
#define BIF_NPARAMS_MASK       0x7ffff00
#define BIF_RETURNSONSTACK     0x0000100
#define BIF_LASTPARMONSTACK    0x0000200
#define BIF_RETURNSBOOL        0x0000400
#define BIF_LASTPARM_ASBOOL    0x0000800
//                             0x00?0000 -- taken by FP stack flags
#define BIF_TAKES_VARPARM      0x0400000
#define BIF_TAKES_VARPARM_EX   0x0C00000 // this is like varparm but check count exactly
#define BIF_WONTMAKEDENORMAL   0x0100000
#define BIF_CLEARDENORMAL      0x0200000
#if defined(GLUE_HAS_FXCH) && GLUE_MAX_FPSTACK_SIZE > 0
  #define BIF_SECONDLASTPARMST 0x0001000 // use with BIF_LASTPARMONSTACK only (last two parameters get passed on fp stack)
  #define BIF_LAZYPARMORDERING 0x0002000 // allow optimizer to avoid fxch when using BIF_TWOPARMSONFPSTACK_LAZY etc
  #define BIF_REVERSEFPORDER   0x0004000 // force a fxch (reverse order of last two parameters on fp stack, used by comparison functions)
  #ifndef BIF_FPSTACKUSE
    #define BIF_FPSTACKUSE(x) (((x)>=0&&(x)<8) ? ((7-(x))<<16):0)
  #endif
  #ifndef BIF_GETFPSTACKUSE
    #define BIF_GETFPSTACKUSE(x) (7 - (((x)>>16)&7))
  #endif
#else
  // do not support fp stack use unless GLUE_HAS_FXCH and GLUE_MAX_FPSTACK_SIZE>0
  #define BIF_SECONDLASTPARMST 0
  #define BIF_LAZYPARMORDERING 0
  #define BIF_REVERSEFPORDER   0
  #define BIF_FPSTACKUSE(x) 0
  #define BIF_GETFPSTACKUSE(x) 0
#endif
#define BIF_TWOPARMSONFPSTACK (BIF_SECONDLASTPARMST|BIF_LASTPARMONSTACK)
#define BIF_TWOPARMSONFPSTACK_LAZY (BIF_LAZYPARMORDERING|BIF_SECONDLASTPARMST|BIF_LASTPARMONSTACK)
EEL_F NSEEL_CGEN_CALL nseel_int_rand(EEL_F f);
#define FNPTR_HAS_CONDITIONAL_EXEC(op)  \
  (op->fntype == FN_LOGICAL_AND || \
   op->fntype == FN_LOGICAL_OR ||  \
   op->fntype == FN_IF_ELSE || \
   op->fntype == FN_WHILE || \
   op->fntype == FN_LOOP)
// Add custom functions
#include <time.h>
#ifndef _WIN32
#include <unistd.h>
#include <sys/time.h>
#endif
static EEL_F NSEEL_CGEN_CALL _eel_sleep(void *opaque, EEL_F *amt)
{
	if (*amt >= 0.0)
	{
#ifdef _WIN32
		if (*amt > 30000000.0) Sleep(30000000);
		else Sleep((DWORD)(*amt + 0.5));
#else
		if (*amt > 30000000.0) usleep(((useconds_t)30000000) * 1000);
		else usleep((useconds_t)(*amt*1000.0 + 0.5));
#endif
	}
	return 0.0;
}
static EEL_F * NSEEL_CGEN_CALL _eel_time(void *opaque, EEL_F *v)
{
	*v = (EEL_F)time(NULL);
	return v;
}
static EEL_F * NSEEL_CGEN_CALL _eel_time_precise(void *opaque, EEL_F *v)
{
#ifdef _WIN32
	LARGE_INTEGER freq, now;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&now);
	*v = (double)now.QuadPart / (double)freq.QuadPart;
	// *v = (EEL_F)timeGetTime() * 0.001;
#else
	struct timeval tm = { 0, };
	gettimeofday(&tm, NULL);
	*v = tm.tv_sec + tm.tv_usec*0.000001;
#endif
	return v;
}
#ifdef EEL_WANT_DOCUMENTATION
static const char *eel_misc_function_reference =
"sleep\tms\tYields the CPU for the millisecond count specified, calling Sleep() on Windows or usleep() on other platforms.\0"
"time\t[&val]\tSets the parameter (or a temporary buffer if omitted) to the number of seconds since January 1, 1970, and returns a reference to that value. "
"The granularity of the value returned is 1 second.\0"
"time_precise\t[&val]\tSets the parameter (or a temporary buffer if omitted) to a system-local timestamp in seconds, and returns a reference to that value. "
"The granularity of the value returned is system defined (but generally significantly smaller than one second).\0"
;
#endif
void discreteHartleyTransform(EEL_F *A, const int nPoints, const EEL_F *sinTab)
{
	int i, j, n, n2, theta_inc, nptDiv2;
	EEL_F alpha, beta;
	for (i = 0; i < nPoints; i += 4)
	{
		const EEL_F	x0 = A[i];
		const EEL_F	x1 = A[i + 1];
		const EEL_F	x2 = A[i + 2];
		const EEL_F	x3 = A[i + 3];
		const EEL_F	y0 = x0 + x1;
		const EEL_F	y1 = x0 - x1;
		const EEL_F	y2 = x2 + x3;
		const EEL_F	y3 = x2 - x3;
		A[i] = y0 + y2;
		A[i + 2] = y0 - y2;
		A[i + 1] = y1 + y3;
		A[i + 3] = y1 - y3;
	}
	for (i = 0; i < nPoints; i += 8)
	{
		alpha = A[i];
		beta = A[i + 4];
		A[i] = alpha + beta;
		A[i + 4] = alpha - beta;
		alpha = A[i + 2];
		beta = A[i + 6];
		A[i + 2] = alpha + beta;
		A[i + 6] = alpha - beta;
		alpha = A[i + 1];
		const EEL_F beta1 = 0.70710678118654752440084436210485*(A[i + 5] + A[i + 7]);
		const EEL_F beta2 = 0.70710678118654752440084436210485*(A[i + 5] - A[i + 7]);
		A[i + 1] = alpha + beta1;
		A[i + 5] = alpha - beta1;
		alpha = A[i + 3];
		A[i + 3] = alpha + beta2;
		A[i + 7] = alpha - beta2;
	}
	n = 16;
	n2 = 8;
	theta_inc = nPoints >> 4;
	nptDiv2 = nPoints >> 2;
	while (n <= nPoints)
	{
		for (i = 0; i < nPoints; i += n)
		{
			int theta = theta_inc;
			const int n4 = n2 >> 1;
			alpha = A[i];
			beta = A[i + n2];
			A[i] = alpha + beta;
			A[i + n2] = alpha - beta;
			alpha = A[i + n4];
			beta = A[i + n2 + n4];
			A[i + n4] = alpha + beta;
			A[i + n2 + n4] = alpha - beta;
			for (j = 1; j < n4; j++)
			{
				EEL_F sinval = sinTab[theta];
				EEL_F cosval = sinTab[theta + nptDiv2];
				EEL_F alpha1 = A[i + j];
				EEL_F alpha2 = A[i - j + n2];
				EEL_F beta1 = A[i + j + n2] * cosval + A[i - j + n] * sinval;
				EEL_F beta2 = A[i + j + n2] * sinval - A[i - j + n] * cosval;
				theta += theta_inc;
				A[i + j] = alpha1 + beta1;
				A[i + j + n2] = alpha1 - beta1;
				A[i - j + n2] = alpha2 + beta2;
				A[i - j + n] = alpha2 - beta2;
			}
		}
		n <<= 1;
		n2 <<= 1;
		theta_inc >>= 1;
	}
}
#define M_PIDouble 3.1415926535897932384626433832795
void getAsymmetricWindow(EEL_F *analysisWnd, EEL_F *synthesisWnd, int k, int m, double freq_temporal)
{
	int i;
	memset(synthesisWnd, 0, k * sizeof(EEL_F));
	if (freq_temporal < 0.4)
		freq_temporal = 0.4;
	if (freq_temporal > 1.8)
		freq_temporal = 1.8;
	int n = ((k - m) << 1) + 2;
	for (i = 0; i < k - m; ++i)
		analysisWnd[i] = (EEL_F)pow(sqrt(0.5 * (1.0 - cos(2.0 * M_PIDouble * (i + 1.0) / (double)n))), freq_temporal);
	n = (m << 1) + 2;
	if (freq_temporal > 1.5)
		freq_temporal = 1.5;
	for (i = k - m; i < k; ++i)
		analysisWnd[i] = (EEL_F)pow(sqrt(0.5 * (1.0 - cos(2.0 * M_PIDouble * ((m + i - (k - m)) + 1.0) / (double)n))), freq_temporal);
	n = m << 1;
	for (i = k - (m << 1); i < k; ++i)
		synthesisWnd[i] = (EEL_F)(0.5 * (1.0 - cos(2.0 * M_PIDouble * (double)(i - (k - (m << 1))) / (double)n))) / analysisWnd[i];
}
void STFT_DynInit(int *indexFw, EEL_F *analysisWnd)
{
	int i;
	int ovpSmps = indexFw[0] / indexFw[1];
	int bufferSize = (indexFw[0] * 6) + indexFw[3] + indexFw[3] + (int)((double)(indexFw[0] * sizeof(unsigned int)) / (double)(sizeof(EEL_F) / sizeof(unsigned int)));
	memset(analysisWnd, 0, bufferSize * sizeof(EEL_F));
	EEL_F *synthesisWnd = analysisWnd + indexFw[0];
	unsigned int *bitRevTbl = (unsigned int*)(analysisWnd + (indexFw[0] * 6) + indexFw[3] + indexFw[3]);
	unsigned int bitsConst = 0;
	unsigned int v = indexFw[0];
	while (v > 1)
	{
		++bitsConst;
		v >>= 1;
	}
	for (i = 0; i < indexFw[0]; ++i)
	{
		unsigned int bits = bitsConst;
		unsigned int x = i;
		bitRevTbl[i] = 0;
		while (bits--)
		{
			bitRevTbl[i] = (bitRevTbl[i] + bitRevTbl[i]) + (x & 1);
			x >>= 1;
		}
	}
	EEL_F *sineTbl = synthesisWnd + indexFw[0];
	double pi2dN = (M_PIDouble * 2.0) / indexFw[0];
	for (i = 0; i < indexFw[0]; ++i)
		sineTbl[i] = (EEL_F)sin(pi2dN * i);
	getAsymmetricWindow(analysisWnd, synthesisWnd, indexFw[0], ovpSmps, indexFw[5] / (EEL_F)32767);
	// Pre-shift window function
	for (i = 0; i < indexFw[0] - indexFw[2]; i++)
		synthesisWnd[i] = synthesisWnd[i + indexFw[2]];
	for (i = 0; i < indexFw[0]; i++)
		analysisWnd[i] *= (1.0 / indexFw[0]) * 0.5;
}
int STFTCartesian(EEL_F *indexer, EEL_F *analysisWnd, EEL_F *ptr)
{
	int *indexFw = (int*)indexer;
	EEL_F *mSineTab = analysisWnd + indexFw[0] * 2;
	EEL_F *mInput = mSineTab + indexFw[0];
	EEL_F *mTempBuffer = mInput + indexFw[0];
	unsigned int *bitRevTbl = (unsigned int*)(analysisWnd + (indexFw[0] * 6) + indexFw[3] + indexFw[3]);
	int i, symIdx;
	for (i = 0; i < indexFw[0]; ++i)
		mTempBuffer[bitRevTbl[i]] = mInput[(i + indexFw[4]) & (indexFw[0] - 1)] * analysisWnd[i];
	discreteHartleyTransform(mTempBuffer, indexFw[0], mSineTab);
	ptr[0] = mTempBuffer[0] * 2.0f;
	ptr[1] = 0.0;
	EEL_F lR, lI;
	for (i = 1; i < ((indexFw[0] >> 1) + 1); i++)
	{
		symIdx = indexFw[0] - i;
		lR = mTempBuffer[i] + mTempBuffer[symIdx];
		lI = mTempBuffer[i] - mTempBuffer[symIdx];
		ptr[i << 1] = lR;
		ptr[(i << 1) + 1] = lI;
	}
	return indexFw[0] + 2;
}
int STFTPolar(EEL_F *indexer, EEL_F *analysisWnd, EEL_F *ptr)
{
	int *indexFw = (int*)indexer;
	EEL_F *mSineTab = analysisWnd + indexFw[0] * 2;
	EEL_F *mInput = mSineTab + indexFw[0];
	EEL_F *mTempBuffer = mInput + indexFw[0];
	unsigned int *bitRevTbl = (unsigned int*)(analysisWnd + (indexFw[0] * 6) + indexFw[3] + indexFw[3]);
	int i, symIdx;
	for (i = 0; i < indexFw[0]; ++i)
		mTempBuffer[bitRevTbl[i]] = mInput[(i + indexFw[4]) & (indexFw[0] - 1)] * analysisWnd[i];
	discreteHartleyTransform(mTempBuffer, indexFw[0], mSineTab);
	ptr[0] = fabs(mTempBuffer[0] * 2.0);
	ptr[1] = ((mTempBuffer[0] < 0.0) ? M_PIDouble : 0.0);
	EEL_F lR, lI;
	for (i = 1; i < ((indexFw[0] >> 1) + 1); i++)
	{
		symIdx = indexFw[0] - i;
		lR = mTempBuffer[i] + mTempBuffer[symIdx];
		lI = mTempBuffer[i] - mTempBuffer[symIdx];
		ptr[i << 1] = hypot(lR, lI);
		ptr[(i << 1) + 1] = atan2(lI, lR);
	}
	return indexFw[0] + 2;
}
int STFTCartesianInverse(EEL_F *indexer, EEL_F *analysisWnd, EEL_F *ptr)
{
	int *indexFw = (int*)indexer;
	EEL_F *synthesisWnd = analysisWnd + indexFw[0];
	EEL_F *mSineTab = synthesisWnd + indexFw[0];
	EEL_F *timeDomainOut = mSineTab + indexFw[0] * 3;
	EEL_F *mOutputBuffer = timeDomainOut + indexFw[0];
	EEL_F *mOverlapStage2Ldash = mOutputBuffer + indexFw[3];
	unsigned int *bitRevTbl = (unsigned int*)(analysisWnd + (indexFw[0] * 6) + indexFw[3] + indexFw[3]);
	int i;
	timeDomainOut[0] = ptr[0];
	EEL_F lR, lI;
	for (i = 1; i < ((indexFw[0] >> 1) + 1); i++)
	{
		lR = ptr[i << 1];
		lI = ptr[(i << 1) + 1];
		timeDomainOut[bitRevTbl[i]] = (lR + lI);
		timeDomainOut[bitRevTbl[indexFw[0] - i]] = (lR - lI);
	}
	discreteHartleyTransform(timeDomainOut, indexFw[0], mSineTab);
	for (i = 0; i < indexFw[0] - indexFw[2]; i++)
		timeDomainOut[i] = timeDomainOut[i + indexFw[2]] * synthesisWnd[i];
	for (i = 0; i < indexFw[3]; ++i)
	{
		mOutputBuffer[i] = mOverlapStage2Ldash[i] + timeDomainOut[i];
		mOverlapStage2Ldash[i] = timeDomainOut[indexFw[3] + i];
	}
	return indexFw[3];
}
int STFTPolarInverse(EEL_F *indexer, EEL_F *analysisWnd, EEL_F *ptr)
{
	int *indexFw = (int*)indexer;
	EEL_F *synthesisWnd = analysisWnd + indexFw[0];
	EEL_F *mSineTab = synthesisWnd + indexFw[0];
	EEL_F *timeDomainOut = mSineTab + indexFw[0] * 3;
	EEL_F *mOutputBuffer = timeDomainOut + indexFw[0];
	EEL_F *mOverlapStage2Ldash = mOutputBuffer + indexFw[3];
	unsigned int *bitRevTbl = (unsigned int*)(analysisWnd + (indexFw[0] * 6) + indexFw[3] + indexFw[3]);
	int i;
	EEL_F magnitude = ptr[0];
	EEL_F phase = ptr[1];
	EEL_F lR, lI;
	timeDomainOut[0] = magnitude * cos(phase);
	for (i = 1; i < ((indexFw[0] >> 1) + 1); i++)
	{
		magnitude = ptr[i << 1];
		phase = ptr[(i << 1) + 1];
		lR = magnitude * cos(phase);
		lI = magnitude * sin(phase);
		timeDomainOut[bitRevTbl[i]] = (lR + lI);
		timeDomainOut[bitRevTbl[indexFw[0] - i]] = (lR - lI);
	}
	discreteHartleyTransform(timeDomainOut, indexFw[0], mSineTab);
	for (i = 0; i < indexFw[0] - indexFw[2]; i++)
		timeDomainOut[i] = timeDomainOut[i + indexFw[2]] * synthesisWnd[i];
	for (i = 0; i < indexFw[3]; ++i)
	{
		mOutputBuffer[i] = mOverlapStage2Ldash[i] + timeDomainOut[i];
		mOverlapStage2Ldash[i] = timeDomainOut[indexFw[3] + i];
	}
	return indexFw[3];
}
static EEL_F NSEEL_CGEN_CALL stftInit(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	compileContext *c = (compileContext*)opaque;
	EEL_F **blocks = c->ram_state.blocks;
	EEL_F *start1 = parms[0];
	int offs1 = (int)(*start1 + NSEEL_CLOSEFACTOR);
	EEL_F *start2 = parms[1];
	int offs2 = (int)(*start2 + NSEEL_CLOSEFACTOR);
	EEL_F *indexer = __NSEEL_RAMAlloc(blocks, offs1);
	if (!indexer || indexer == &nseel_ramalloc_onfail)
		return 0;
	EEL_F *stftFloatStruct = __NSEEL_RAMAlloc(blocks, offs2);
	if (!stftFloatStruct || stftFloatStruct == &nseel_ramalloc_onfail)
		return 0;
	int *indexFw = (int*)indexer;
	STFT_DynInit(indexFw, stftFloatStruct);
	return indexFw[3];
}
static EEL_F NSEEL_CGEN_CALL stftForward(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	compileContext *c = (compileContext*)opaque;
	EEL_F **blocks = c->ram_state.blocks;
	EEL_F *start = parms[0];
	int offs = (int)(*start + NSEEL_CLOSEFACTOR);
	EEL_F *ptr = __NSEEL_RAMAlloc(blocks, offs);
	if (!ptr || ptr == &nseel_ramalloc_onfail)
		return 0;
	EEL_F *start1 = parms[1];
	int offs1 = (int)(*start1 + NSEEL_CLOSEFACTOR);
	EEL_F *start2 = parms[2];
	int offs2 = (int)(*start2 + NSEEL_CLOSEFACTOR);
	EEL_F *indexer = __NSEEL_RAMAlloc(blocks, offs1);
	if (!indexer || indexer == &nseel_ramalloc_onfail)
		return 0;
	EEL_F *stftFloatStruct = __NSEEL_RAMAlloc(blocks, offs2);
	if (!stftFloatStruct || stftFloatStruct == &nseel_ramalloc_onfail)
		return 0;
	int cartesian = (int)(*parms[3]);
	int *indexFw = (int*)indexer;
	EEL_F *mInput = stftFloatStruct + indexFw[0] * 3;
	memcpy(&mInput[indexFw[4]], ptr, indexFw[3] * sizeof(EEL_F));
	indexFw[4] = (indexFw[4] + indexFw[3]) & (indexFw[0] - 1);
	if (cartesian)
		return (EEL_F)STFTCartesian(indexer, stftFloatStruct, ptr);
	else
		return (EEL_F)STFTPolar(indexer, stftFloatStruct, ptr);
}
static EEL_F NSEEL_CGEN_CALL stftBackward(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	compileContext *c = (compileContext*)opaque;
	EEL_F **blocks = c->ram_state.blocks;
	EEL_F *start = parms[0];
	int offs = (int)(*start + NSEEL_CLOSEFACTOR);
	EEL_F *ptr = __NSEEL_RAMAlloc(blocks, offs);
	if (!ptr || ptr == &nseel_ramalloc_onfail)
		return 0;
	EEL_F *start1 = parms[1];
	int offs1 = (int)(*start1 + NSEEL_CLOSEFACTOR);
	EEL_F *start2 = parms[2];
	int offs2 = (int)(*start2 + NSEEL_CLOSEFACTOR);
	EEL_F *indexer = __NSEEL_RAMAlloc(blocks, offs1);
	if (!indexer || indexer == &nseel_ramalloc_onfail)
		return 0;
	EEL_F *stftFloatStruct = __NSEEL_RAMAlloc(blocks, offs2);
	if (!stftFloatStruct || stftFloatStruct == &nseel_ramalloc_onfail)
		return 0;
	int cartesian = (int)(*parms[3]);
	int *indexFw = (int*)indexer;
	int ret;
	if (cartesian)
		ret = STFTCartesianInverse(indexer, stftFloatStruct, ptr);
	else
		ret = STFTPolarInverse(indexer, stftFloatStruct, ptr);
	memcpy(ptr, stftFloatStruct + indexFw[0] * 6, indexFw[3] * sizeof(EEL_F));
	return (EEL_F)ret;
}
int STFT_DynConstructor(EEL_F *indexer, int fftLen, int analysisOvp, EEL_F tf_res)
{
	int ovpSmps = fftLen / analysisOvp;
	int sampleShift = (fftLen - (ovpSmps << 1));
	int *indexFw = (int*)indexer;
	indexFw[0] = fftLen;
	indexFw[1] = analysisOvp;
	indexFw[2] = sampleShift;
	indexFw[3] = ovpSmps;
	indexFw[4] = 0;
	indexFw[5] = (int)(tf_res * (EEL_F)32767);
	return ((fftLen * 6) + ovpSmps + ovpSmps + (int)((EEL_F)(fftLen * sizeof(unsigned int)) / (EEL_F)(sizeof(EEL_F) / sizeof(unsigned int))));
}
static EEL_F NSEEL_CGEN_CALL stftCheckMemoryRequirement(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	compileContext *c = (compileContext*)opaque;
	EEL_F **blocks = c->ram_state.blocks;
	EEL_F *start1 = parms[0];
	EEL_F *length = parms[1];
	EEL_F *analysisOvp = parms[2];
	int offs1 = (int)(*start1 + NSEEL_CLOSEFACTOR);
	EEL_F *indexer = __NSEEL_RAMAlloc(blocks, offs1);
	if (!indexer || indexer == &nseel_ramalloc_onfail)
		return 0;
	int fftlen = (int)(*length + NSEEL_CLOSEFACTOR);
	int analyOv = (int)(*analysisOvp + NSEEL_CLOSEFACTOR);
	return (EEL_F)STFT_DynConstructor(indexer, fftlen, analyOv, *parms[3]);
}
/*void EEL_stft_register(NSEEL_VMCTX vm)
{
	compileContext *ctx = (compileContext*)vm;
	NSEEL_addfunc_varparm(ctx, "stftCheckMemoryRequirement", 1, NSEEL_PProc_THIS, &stftCheckMemoryRequirement);
	NSEEL_addfunc_varparm(ctx, "stftInit", 1, NSEEL_PProc_THIS, &stftInit);
	NSEEL_addfunc_varparm(ctx, "stftForward", 1, NSEEL_PProc_THIS, &stftForward);
	NSEEL_addfunc_varparm(ctx, "stftBackward", 1, NSEEL_PProc_THIS, &stftBackward);
}*/
void reverse(EEL_F *arr, int start, int end)
{
	while (start < end)
	{
		EEL_F tmp = arr[start];
		arr[start] = arr[end];
		arr[end] = tmp;
		start++;
		end--;
	}
}
void shift(EEL_F *arr, int k, int n)
{
	k = k % n;
	reverse(arr, 0, n - 1);
	reverse(arr, 0, n - k - 1);
	reverse(arr, n - k, n - 1);
}
EEL_F * NSEEL_CGEN_CALL __NSEEL_circshift(EEL_F **blocks, EEL_F *offptr, EEL_F *shiftptr, EEL_F *lenptr)
{
	int offs = (int)(*offptr + NSEEL_CLOSEFACTOR);
	EEL_F *arr = __NSEEL_RAMAlloc(blocks, offs);
	if (!arr || arr == &nseel_ramalloc_onfail)
		return 0;
	int k = (int)*shiftptr;
	int n = (int)*lenptr;
	k < 0 ? shift(arr, -k, n) : shift(arr, n - k, n);
	return offptr;
}
EEL_F expint(EEL_F x)
{
	if (x <= 1.0)
		return -0.57721566490153286060651209 - log(x) + x * (-x / 4.0 + 1.0);
	EEL_F term = 0.0;
	for (int k = 1 + (int)floor(80.0 / x); k >= 1; k--)
		term = k / (1.0 + k / (x + term));
	return exp(-x) / (x + term);
}
// Domain: 0.001 ~ 0.5
static const float ExpintTable1[500] = { 6.33153936413615f, 5.63939143396494f, 5.23492507691165f, 4.94824125651360f, 4.72609545858444f, 4.54477115683906f, 4.39161723405586f, 4.25908210080260f, 4.14229482717614f, 4.03792957653811f, 3.94361416507444f, 3.85759706007702f, 3.77854812837773f, 3.70543343651052f, 3.63743334995231f, 3.57388711871546f, 3.51425429210118f, 3.45808717909408f, 3.40501076461633f, 3.35470778330971f, 3.30690743883812f, 3.26137674984624f, 3.21791382119158f, 3.17634254828989f, 3.13650840321517f, 3.09827504776313f, 3.06152158606426f, 3.02614031708686f, 2.99203488170511f, 2.95911872402128f, 2.92731380507851f, 2.89654952085833f, 2.86676178682571f, 2.83789225917520f, 2.80988766899124f, 2.78269925022877f, 2.75628224608441f, 2.73059548120986f, 2.70560098950218f, 2.68126368902528f, 2.65755109707773f, 2.63443307960029f, 2.61188163007365f, 2.58987067383731f, 2.56837589440108f, 2.54737457884799f, 2.52684547986474f, 2.50676869229857f, 2.48712554244313f, 2.46789848850997f, 2.44907103095643f, 2.43062763152119f, 2.41255363997231f, 2.39483522770240f, 2.37745932741712f, 2.36041357825803f, 2.34368627578265f, 2.32726632629475f, 2.31114320507867f, 2.29530691814378f, 2.27974796713116f, 2.26445731707363f, 2.24942636673543f, 2.23464692128761f, 2.22011116710190f, 2.20581164846892f, 2.19174124606712f, 2.17789315702680f, 2.16426087644945f, 2.15083818025680f, 2.13761910925637f, 2.12459795432145f, 2.11176924259316f, 2.09912772462139f, 2.08666836236871f, 2.07438631800884f, 2.06227694345737f, 2.05033577057783f, 2.03855850201183f, 2.02694100258574f, 2.01547929125139f, 2.00416953352107f, 1.99300803436109f, 1.98199123151100f, 1.97111568919802f, 1.96037809221915f, 1.94977524036537f, 1.93930404316446f, 1.92896151492086f, 1.91874477003266f, 1.90865101856733f, 1.89867756207924f, 1.88882178965323f, 1.87908117415985f, 1.86945326870865f, 1.85993570328723f, 1.85052618157445f, 1.84122247791703f, 1.83202243445968f, 1.82292395841939f, 1.81392501949540f, 1.80502364740671f, 1.79621792954979f, 1.78750600876935f, 1.77888608123583f, 1.77035639442354f, 1.76191524518360f, 1.75356097790666f, 1.74529198277014f, 1.73710669406567f, 1.72900358860208f, 1.72098118418005f, 1.71303803813460f, 1.70517274594176f, 1.69738393988603f, 1.68967028778564f, 1.68203049177242f, 1.67446328712364f, 1.66696744114317f, 1.65954175208938f, 1.65218504814758f, 1.64489618644475f, 1.63767405210444f, 1.63051755734000f, 1.62342564058417f, 1.61639726565340f, 1.60943142094514f, 1.60252711866667f, 1.59568339409389f, 1.58889930485880f, 1.58217393026415f, 1.57550637062431f, 1.56889574663089f, 1.56234119874219f, 1.55584188659529f, 1.54939698843990f, 1.54300570059289f, 1.53666723691273f, 1.53038082829278f, 1.52414572217289f, 1.51796118206826f, 1.51182648711501f, 1.50574093163171f, 1.49970382469613f, 1.49371448973666f, 1.48777226413783f, 1.48187649885916f, 1.47602655806709f, 1.47022181877919f, 1.46446167052028f, 1.45874551499005f, 1.45307276574158f, 1.44744284787039f, 1.44185519771371f, 1.43630926255935f, 1.43080450036407f, 1.42534037948083f, 1.41991637839478f, 1.41453198546754f, 1.40918669868949f, 1.40388002543983f, 1.39861148225397f, 1.39338059459819f, 1.38818689665114f, 1.38302993109194f, 1.37790924889475f, 1.37282440912947f, 1.36777497876843f, 1.36276053249873f, 1.35778065254024f, 1.35283492846881f, 1.34792295704476f, 1.34304434204620f, 1.33819869410729f, 1.33338563056104f, 1.32860477528666f, 1.32385575856114f, 1.31913821691513f, 1.31445179299275f, 1.30979613541537f, 1.30517089864914f, 1.30057574287617f, 1.29601033386927f, 1.29147434287001f, 1.28696744647023f, 1.28248932649660f, 1.27803966989839f, 1.27361816863814f, 1.26922451958530f, 1.26485842441262f, 1.26051958949533f, 1.25620772581287f, 1.25192254885323f, 1.24766377851975f, 1.24343113904031f, 1.23922435887882f, 1.23504317064902f, 1.23088731103039f, 1.22675652068626f, 1.22265054418389f, 1.21856912991660f, 1.21451203002782f, 1.21047900033696f, 1.20646980026724f, 1.20248419277510f, 1.19852194428150f, 1.19458282460475f, 1.19066660689505f, 1.18677306757053f, 1.18290198625489f, 1.17905314571646f, 1.17522633180872f, 1.17142133341224f, 1.16763794237793f, 1.16387595347169f, 1.16013516432024f, 1.15641537535834f, 1.15271638977706f, 1.14903801347337f, 1.14538005500082f, 1.14174232552131f, 1.13812463875802f, 1.13452681094936f, 1.13094866080393f, 1.12739000945649f, 1.12385068042499f, 1.12033049956838f, 1.11682929504554f, 1.11334689727494f, 1.10988313889530f, 1.10643785472705f, 1.10301088173459f, 1.09960205898947f, 1.09621122763424f, 1.09283823084711f, 1.08948291380742f, 1.08614512366176f, 1.08282470949081f, 1.07952152227697f, 1.07623541487250f, 1.07296624196851f, 1.06971386006443f, 1.06647812743823f, 1.06325890411715f, 1.06005605184912f, 1.05686943407471f, 1.05369891589963f, 1.05054436406786f, 1.04740564693527f, 1.04428263444374f, 1.04117519809585f, 1.03808321093007f, 1.03500654749642f, 1.03194508383258f, 1.02889869744056f, 1.02586726726374f, 1.02285067366439f, 1.01984879840164f, 1.01686152460984f, 1.01388873677738f, 1.01093032072589f, 1.00798616358984f, 1.00505615379652f, 1.00214018104644f, 0.999238136294059f, 0.996349911728865f, 0.993475400756876f, 0.990614497982423f, 0.987767099190299f, 0.984933101328234f, 0.982112402489699f, 0.979304901897022f, 0.976510499884809f, 0.973729097883682f, 0.970960598404302f, 0.968204905021682f, 0.965461922359790f, 0.962731556076425f, 0.960013712848367f, 0.957308300356788f, 0.954615227272931f, 0.951934403244034f, 0.949265738879512f, 0.946609145737375f, 0.943964536310890f, 0.941331824015476f, 0.938710923175825f, 0.936101749013249f, 0.933504217633244f, 0.930918246013275f, 0.928343751990763f, 0.925780654251284f, 0.923228872316966f, 0.920688326535087f, 0.918158938066864f, 0.915640628876432f, 0.913133321720014f, 0.910636940135260f, 0.908151408430781f, 0.905676651675847f, 0.903212595690257f, 0.900759167034383f, 0.898316292999373f, 0.895883901597516f, 0.893461921552770f, 0.891050282291436f, 0.888648913932996f, 0.886257747281089f, 0.883876713814638f, 0.881505745679123f, 0.879144775677990f, 0.876793737264197f, 0.874452564531905f, 0.872121192208285f, 0.869799555645474f, 0.867487590812641f, 0.865185234288193f, 0.862892423252093f, 0.860609095478304f, 0.858335189327350f, 0.856070643738996f, 0.853815398225033f, 0.851569392862186f, 0.849332568285124f, 0.847104865679583f, 0.844886226775585f, 0.842676593840777f, 0.840475909673856f, 0.838284117598100f, 0.836101161455002f, 0.833926985597995f, 0.831761534886266f, 0.829604754678678f, 0.827456590827773f, 0.825316989673861f, 0.823185898039208f, 0.821063263222303f, 0.818949032992211f, 0.816843155583011f, 0.814745579688315f, 0.812656254455868f, 0.810575129482225f, 0.808502154807510f, 0.806437280910244f, 0.804380458702259f, 0.802331639523674f, 0.800290775137951f, 0.798257817727021f, 0.796232719886476f, 0.794215434620836f, 0.792205915338876f, 0.790204115849027f, 0.788209990354839f, 0.786223493450506f, 0.784244580116458f, 0.782273205715011f, 0.780309325986082f, 0.778352897042962f, 0.776403875368148f, 0.774462217809232f, 0.772527881574849f, 0.770600824230679f, 0.768681003695509f, 0.766768378237340f, 0.764862906469556f, 0.762964547347142f, 0.761073260162954f, 0.759189004544038f, 0.757311740448002f, 0.755441428159437f, 0.753578028286381f, 0.751721501756837f, 0.749871809815336f, 0.748028914019544f, 0.746192776236916f, 0.744363358641393f, 0.742540623710147f, 0.740724534220363f, 0.738915053246067f, 0.737112144154998f, 0.735315770605516f, 0.733525896543553f, 0.731742486199606f, 0.729965504085765f, 0.728194914992784f, 0.726430683987183f, 0.724672776408401f, 0.722921157865967f, 0.721175794236724f, 0.719436651662080f, 0.717703696545296f, 0.715976895548810f, 0.714256215591593f, 0.712541623846539f, 0.710833087737889f, 0.709130574938690f, 0.707434053368279f, 0.705743491189805f, 0.704058856807781f, 0.702380118865662f, 0.700707246243463f, 0.699040208055393f, 0.697378973647534f, 0.695723512595532f, 0.694073794702333f, 0.692429789995936f, 0.690791468727175f, 0.689158801367534f, 0.687531758606980f, 0.685910311351833f, 0.684294430722649f, 0.682684088052144f, 0.681079254883129f, 0.679479902966478f, 0.677886004259118f, 0.676297530922046f, 0.674714455318365f, 0.673136750011346f, 0.671564387762516f, 0.669997341529763f, 0.668435584465468f, 0.666879089914658f, 0.665327831413178f, 0.663781782685891f, 0.662240917644895f, 0.660705210387756f, 0.659174635195774f, 0.657649166532257f, 0.656128779040825f, 0.654613447543725f, 0.653103147040173f, 0.651597852704708f, 0.650097539885575f, 0.648602184103115f, 0.647111761048182f, 0.645626246580576f, 0.644145616727490f, 0.642669847681983f, 0.641198915801459f, 0.639732797606177f, 0.638271469777766f, 0.636814909157763f, 0.635363092746164f, 0.633915997699998f, 0.632473601331907f, 0.631035881108752f, 0.629602814650226f, 0.628174379727492f, 0.626750554261822f, 0.625331316323269f, 0.623916644129338f, 0.622506516043683f, 0.621100910574808f, 0.619699806374795f, 0.618303182238032f, 0.616911017099967f, 0.615523290035870f, 0.614139980259605f, 0.612761067122424f, 0.611386530111769f, 0.610016348850085f, 0.608650503093652f, 0.607288972731423f, 0.605931737783878f, 0.604578778401895f, 0.603230074865621f, 0.601885607583366f, 0.600545357090505f, 0.599209304048392f, 0.597877429243282f, 0.596549713585273f, 0.595226138107252f, 0.593906683963851f, 0.592591332430423f, 0.591280064902018f, 0.589972862892376f, 0.588669708032933f, 0.587370582071829f, 0.586075466872932f, 0.584784344414874f, 0.583497196790091f, 0.582214006203880f, 0.580934754973459f, 0.579659425527039f, 0.578388000402911f, 0.577120462248533f, 0.575856793819634f, 0.574596977979323f, 0.573340997697209f, 0.572088836048530f, 0.570840476213289f, 0.569595901475400f, 0.568355095221847f, 0.567118040941839f, 0.565884722225992f, 0.564655122765501f, 0.563429226351332f, 0.562207016873418f, 0.560988478319864f, 0.559773594776161f };
// Domain: 0.5 ~ 1.5
static const float ExpintTable2[251] = { 0.559773594776161f, 0.554950295774654f, 0.550184228566447f, 0.545474419724144f, 0.540819918846431f, 0.536219797845636f, 0.531673150262605f, 0.527179090607617f, 0.522736753726183f, 0.518345294188611f, 0.514003885702249f, 0.509711720545442f, 0.505468009022217f, 0.501271978936804f, 0.497122875087133f, 0.493019958776493f, 0.488962507342567f, 0.484949813703121f, 0.480981185917635f, 0.477055946764210f, 0.473173433331126f, 0.469332996622433f, 0.465534001177009f, 0.461775824700544f, 0.458057857709903f, 0.454379503189402f, 0.450740176258502f, 0.447139303850470f, 0.443576324401589f, 0.440050687550483f, 0.436561853847192f, 0.433109294471587f, 0.429692490960805f, 0.426310934945320f, 0.422964127893360f, 0.419651580863333f, 0.416372814263966f, 0.413127357621880f, 0.409914749356315f, 0.406734536560751f, 0.403586274791166f, 0.400469527860700f, 0.397383867640485f, 0.394328873866420f, 0.391304133951693f, 0.388309242804826f, 0.385343802653065f, 0.382407422870921f, 0.379499719813686f, 0.376620316655740f, 0.373768843233509f, 0.370944935892892f, 0.368148237341009f, 0.365378396502136f, 0.362635068377678f, 0.359917913910035f, 0.357226599850257f, 0.354560798629336f, 0.351920188233033f, 0.349304452080114f, 0.346713278903895f, 0.344146362636973f, 0.341603402299064f, 0.339084101887818f, 0.336588170272547f, 0.334115321090748f, 0.331665272647358f, 0.329237747816627f, 0.326832473946554f, 0.324449182765786f, 0.322087610292923f, 0.319747496748133f, 0.317428586467031f, 0.315130627816727f, 0.312853373114005f, 0.310596578545543f, 0.308360004090134f, 0.306143413442834f, 0.303946573940985f, 0.301769256492064f, 0.299611235503289f, 0.297472288812947f, 0.295352197623386f, 0.293250746435620f, 0.291167722985511f, 0.289102918181468f, 0.287056126043634f, 0.285027143644513f, 0.283015771050992f, 0.281021811267725f, 0.279045070181840f, 0.277085356508928f, 0.275142481740288f, 0.273216260091375f, 0.271306508451441f, 0.269413046334320f, 0.267535695830332f, 0.265674281559272f, 0.263828630624464f, 0.261998572567834f, 0.260183939326000f, 0.258384565187321f, 0.256600286749911f, 0.254830942880567f, 0.253076374674603f, 0.251336425416562f, 0.249610940541772f, 0.247899767598750f, 0.246202756212403f, 0.244519758048025f, 0.242850626776061f, 0.241195218037623f, 0.239553389410736f, 0.237925000377293f, 0.236309912290710f, 0.234707988344255f, 0.233119093540036f, 0.231543094658634f, 0.229979860229366f, 0.228429260501155f, 0.226891167414003f, 0.225365454571042f, 0.223851997211154f, 0.222350672182151f, 0.220861357914491f, 0.219383934395520f, 0.217918283144238f, 0.216464287186559f, 0.215021831031070f, 0.213590800645265f, 0.212171083432249f, 0.210762568207898f, 0.209365145178468f, 0.207978705918639f, 0.206603143349984f, 0.205238351719856f, 0.203884226580681f, 0.202540664769647f, 0.201207564388785f, 0.199884824785426f, 0.198572346533028f, 0.197270031412369f, 0.195977782393090f, 0.194695503615591f, 0.193423100373256f, 0.192160479095018f, 0.190907547328241f, 0.189664213721924f, 0.188430388010206f, 0.187205980996190f, 0.185990904536040f, 0.184785071523391f, 0.183588395874023f, 0.182400792510829f, 0.181222177349038f, 0.180052467281716f, 0.178891580165522f, 0.177739434806718f, 0.176595950947427f, 0.175461049252136f, 0.174334651294438f, 0.173216679544007f, 0.172107057353797f, 0.171005708947475f, 0.169912559407060f, 0.168827534660787f, 0.167750561471176f, 0.166681567423308f, 0.165620480913305f, 0.164567231136999f, 0.163521748078805f, 0.162483962500777f, 0.161453805931852f, 0.160431210657273f, 0.159416109708196f, 0.158408436851462f, 0.157408126579554f, 0.156415114100704f, 0.155429335329182f, 0.154450726875736f, 0.153479226038189f, 0.152514770792199f, 0.151557299782162f, 0.150606752312267f, 0.149663068337703f, 0.148726188455997f, 0.147796053898504f, 0.146872606522027f, 0.145955788800575f, 0.145045543817253f, 0.144141815256283f, 0.143244547395151f, 0.142353685096878f, 0.141469173802417f, 0.140590959523164f, 0.139718988833596f, 0.138853208864020f, 0.137993567293430f, 0.137140012342489f, 0.136292492766605f, 0.135450957849129f, 0.134615357394647f, 0.133785641722382f, 0.132961761659695f, 0.132143668535687f, 0.131331314174900f, 0.130524650891108f, 0.129723631481212f, 0.128928209219219f, 0.128138337850319f, 0.127353971585044f, 0.126575065093524f, 0.125801573499819f, 0.125033452376345f, 0.124270657738376f, 0.123513146038632f, 0.122760874161947f, 0.122013799420013f, 0.121271879546203f, 0.120535072690475f, 0.119803337414337f, 0.119076632685907f, 0.118354917875021f, 0.117638152748432f, 0.116926297465068f, 0.116219312571358f, 0.115517158996635f, 0.114819798048592f, 0.114127191408815f, 0.113439301128373f, 0.112756089623470f, 0.112077519671165f, 0.111403554405148f, 0.110734157311576f, 0.110069292224971f, 0.109408923324170f, 0.108753015128340f, 0.108101532493039f, 0.107454440606342f, 0.106811704985010f, 0.106173291470726f, 0.105539166226369f, 0.104909295732348f, 0.104283646782986f, 0.103662186482950f, 0.103044882243734f, 0.102431701780185f, 0.101822613107082f, 0.101217584535761f, 0.100616584670779f, 0.100019582406633f };
// Domain 1.5 ~ 3.5
static const float ExpintTable3[126] = { 0.100019582406633f, 0.0976709372158198f, 0.0953838384149070f, 0.0931564276459467f, 0.0909869129830044f, 0.0888735660844122f, 0.0868147194907422f, 0.0848087640597390f, 0.0828541465300625f, 0.0809493672062480f, 0.0790929777578066f, 0.0772835791258635f, 0.0755198195311743f, 0.0738003925777604f, 0.0721240354467908f, 0.0704895271756689f, 0.0688956870176292f, 0.0673413728774234f, 0.0658254798189725f, 0.0643469386411095f, 0.0629047145177793f, 0.0614978056992866f, 0.0601252422713888f, 0.0587860849692256f, 0.0574794240432553f, 0.0562043781745348f, 0.0549600934368431f, 0.0537457423032842f, 0.0525605226951504f, 0.0514036570709528f, 0.0502743915536392f, 0.0491719950941388f, 0.0480957586694767f, 0.0470449945137897f, 0.0460190353806841f, 0.0450172338354382f, 0.0440389615756621f, 0.0430836087790745f, 0.0421505834771471f, 0.0412393109534225f, 0.0403492331653865f, 0.0394798081888161f, 0.0386305096836003f, 0.0378008263800681f, 0.0369902615849173f, 0.0361983327058716f, 0.0354245707942568f, 0.0346685201047068f, 0.0339297376712618f, 0.0332077928991558f, 0.0325022671716214f, 0.0318127534710799f, 0.0311388560141009f, 0.0304801898995688f, 0.0298363807694920f, 0.0292070644819443f, 0.0285918867956309f, 0.0279905030656117f, 0.0274025779497217f, 0.0268277851252616f, 0.0262658070155479f, 0.0257163345259216f, 0.0251790667888503f, 0.0246537109177589f, 0.0241399817692477f, 0.0236376017133719f, 0.0231463004116726f, 0.0226658146026547f, 0.0221958878944312f, 0.0217362705642599f, 0.0212867193647089f, 0.0208469973362054f, 0.0204168736257225f, 0.0199961233113807f, 0.0195845272327421f, 0.0191818718265862f, 0.0187879489679693f, 0.0184025558163696f, 0.0180254946667384f, 0.0176565728052777f, 0.0172956023697725f, 0.0169424002143170f, 0.0165967877782784f, 0.0162585909593461f, 0.0159276399905230f, 0.0156037693209245f, 0.0152868175002407f, 0.0149766270667522f, 0.0146730444387568f, 0.0143759198093052f, 0.0140851070441270f, 0.0138004635826326f, 0.0135218503418941f, 0.0132491316235018f, 0.0129821750231990f, 0.0127208513432019f, 0.0124650345071144f, 0.0122146014773545f, 0.0119694321750053f, 0.0117294094020133f, 0.0114944187656559f, 0.0112643486052043f, 0.0110390899207104f, 0.0108185363038491f, 0.0106025838707484f, 0.0103911311967455f, 0.0101840792530053f, 0.00998133134494335f, 0.00978279305239503f, 0.00958837217147620f, 0.00939797865808188f, 0.00921152457297160f, 0.00902892402839181f, 0.00885009313618758f, 0.00867494995735709f, 0.00850341445300474f, 0.00833540843664957f, 0.00817085552784730f, 0.00800968110708622f, 0.00785181227191778f, 0.00769717779428456f, 0.00754570807900947f, 0.00739733512341117f, 0.00725199247801176f, 0.00710961520830429f, 0.00697013985754840f };
// Domain 3.5 ~ 11.436
static const float ExpintTable4[63] = { 0.00697013985754843f, 0.00595165118439979f, 0.00508661641561241f, 0.00435102957490895f, 0.00372481546087525f, 0.00319115004588169f, 0.00273590601139406f, 0.00234719870089483f, 0.00201501303487467f, 0.00173089598574009f, 0.00148770235377602f, 0.00127938403852772f, 0.00110081492605982f, 0.000947645033184157f, 0.000816178756486768f, 0.000703273036244310f, 0.000606252016086186f, 0.000522835399236460f, 0.000451078202821178f, 0.000389320017508240f, 0.000336142209771689f, 0.000290331773353231f, 0.000250850756880151f, 0.000216810375484428f, 0.000187449063132408f, 0.000162113845196945f, 0.000140244512382263f, 0.000121360161310043f, 0.000105047737015578f, 9.09522708167131e-05f, 7.87685555639048e-05f, 6.82340408359011e-05f, 5.91227645848651e-05f, 5.12401661828662e-05f, 4.44186497036028e-05f, 3.85137863511616e-05f, 3.34010618509128e-05f, 2.89730888669951e-05f, 2.51372165380621e-05f, 2.18134793867002e-05f, 1.89328364567230e-05f, 1.64356588152578e-05f, 1.42704297310207e-05f, 1.23926270802452e-05f, 1.07637619830281e-05f, 9.35055145753116e-06f, 8.12420610016641e-06f, 7.05981654290921e-06f, 6.13582477694144e-06f, 5.33356842622777e-06f, 4.63688775713486e-06f, 4.03178666455650e-06f, 3.50614011821845e-06f, 3.04944161624825e-06f, 2.65258510327885e-06f, 2.30767658985817e-06f, 2.00787137790676e-06f, 1.74723336968152e-06f, 1.52061342901041e-06f, 1.32354418525242e-06f, 1.15214903255314e-06f, 1.00306338807548e-06f, 8.73366540296160e-07f };
float linear_value(float Start, float Step, float Input, const float *Space)
{
	int Index = (int)((Input - Start) / Step);
	float X1 = Start + Step * Index;
	return Space[Index] + (Space[Index + 1] - Space[Index]) / (Start + Step * (Index + 1) - X1) * (Input - X1);
}
EEL_F expint_interpolation(EEL_F x)
{
	if (x < (EEL_F)0.001)
		return (EEL_F)6.33153936413615;
	else if (x < (EEL_F)0.5)
		return (EEL_F)linear_value(0.001f, 0.001f, (float)x, ExpintTable1);// 0.001:0.001:0.5
	else if (x < (EEL_F)1.5)
		return (EEL_F)linear_value(0.5f, 0.004f, (float)x, ExpintTable2); // 0.5:0.004:1.5
	else if (x < (EEL_F)3.5)
		return (EEL_F)linear_value(1.5f, 0.016f, (float)x, ExpintTable3);// 1.5:0.016:3.5
	else if (x < (EEL_F)11.436)
		return (EEL_F)linear_value(3.5f, 0.128f, (float)x, ExpintTable4);// 3.5:0.128:11.5
	else
		return (EEL_F)8.73366540296160e-07;
}
EEL_F hypotFast(EEL_F x, EEL_F y)
{
	return sqrt(x * x + y * y);
}
EEL_F invsqrt(EEL_F x)
{
	return 1.0 / sqrt(x);
}
#include "fft.h"
static void fft_reorder_buffer(int bitsz, WDL_FFT_COMPLEX *data, int fwd)
{
	const int *tab = fft_reorder_table_for_bitsize(bitsz);
	if (!fwd)
	{
		while (*tab)
		{
			const int sidx = *tab++;
			WDL_FFT_COMPLEX a = data[sidx];
			for (;;)
			{
				WDL_FFT_COMPLEX ta;
				const int idx = *tab++;
				if (!idx) break;
				ta = data[idx];
				data[idx] = a;
				a = ta;
			}
			data[sidx] = a;
		}
	}
	else
	{
		while (*tab)
		{
			const int sidx = *tab++;
			int lidx = sidx;
			const WDL_FFT_COMPLEX sta = data[lidx];
			for (;;)
			{
				const int idx = *tab++;
				if (!idx) break;
				data[lidx] = data[idx];
				lidx = idx;
			}
			data[lidx] = sta;
		}
	}
}
// 0=fw, 1=iv, 2=fwreal, 3=ireal, 4=permutec, 6=permuter
// low bit: is inverse
// second bit: was isreal, but no longer used
// third bit: is permute
static void FFT(int sizebits, EEL_F *data, int dir)
{
	if (dir >= 4 && dir < 8)
	{
		if (dir == 4 || dir == 5)
			fft_reorder_buffer(sizebits, (WDL_FFT_COMPLEX*)data, dir == 4);
	}
	else if (dir >= 0 && dir < 2)
		WDL_fft((WDL_FFT_COMPLEX*)data, 1 << sizebits, dir & 1);
	else if (dir >= 2 && dir < 4)
		WDL_real_fft((EEL_F*)data, 1 << sizebits, dir & 1);
}
static EEL_F * fft_func(int dir, EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	const int offs = (int)(*start + NSEEL_CLOSEFACTOR);
	const int itemSizeShift = (dir & 2) ? 0 : 1;
	int l = (int)(*length + NSEEL_CLOSEFACTOR);
	int bitl = 0;
	int ilen;
	EEL_F *ptr;
	while (l > 1 && bitl < EEL_FFT_MAXBITLEN)
	{
		bitl++;
		l >>= 1;
	}
	if (bitl < ((dir & 4) ? EEL_FFT_MINBITLEN_REORDER : EEL_FFT_MINBITLEN))  // smallest FFT is 16 item, smallest reorder is 8 item
	{
		return start;
	}
	ilen = 1 << bitl;
	// check to make sure we don't cross a boundary
	if (offs / NSEEL_RAM_ITEMSPERBLOCK != (offs + (ilen << itemSizeShift) - 1) / NSEEL_RAM_ITEMSPERBLOCK)
	{
		return start;
	}
	ptr = __NSEEL_RAMAlloc(blocks, offs);
	if (!ptr || ptr == &nseel_ramalloc_onfail)
	{
		return start;
	}
	FFT(bitl, ptr, dir);
	return start;
}
static EEL_F * NSEEL_CGEN_CALL  eel_fft(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	return fft_func(0, blocks, start, length);
}
static EEL_F * NSEEL_CGEN_CALL  eel_ifft(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	return fft_func(1, blocks, start, length);
}
static EEL_F * NSEEL_CGEN_CALL  eel_fft_real(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	return fft_func(2, blocks, start, length);
}
static EEL_F * NSEEL_CGEN_CALL  eel_ifft_real(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	return fft_func(3, blocks, start, length);
}
static EEL_F * NSEEL_CGEN_CALL  eel_fft_permute(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	return fft_func(4, blocks, start, length);
}
static EEL_F * NSEEL_CGEN_CALL  eel_ifft_permute(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	return fft_func(5, blocks, start, length);
}
static EEL_F * NSEEL_CGEN_CALL eel_convolve_c(EEL_F **blocks, EEL_F *dest, EEL_F *src, EEL_F *lenptr)
{
	const int dest_offs = (int)(*dest + NSEEL_CLOSEFACTOR);
	const int src_offs = (int)(*src + NSEEL_CLOSEFACTOR);
	const int len = ((int)(*lenptr + NSEEL_CLOSEFACTOR)) * 2;
	EEL_F *srcptr, *destptr;
	if (len < 1 || len > NSEEL_RAM_ITEMSPERBLOCK || dest_offs < 0 || src_offs < 0 ||
		dest_offs >= NSEEL_RAM_BLOCKS * NSEEL_RAM_ITEMSPERBLOCK || src_offs >= NSEEL_RAM_BLOCKS * NSEEL_RAM_ITEMSPERBLOCK) return dest;
	if ((dest_offs&(NSEEL_RAM_ITEMSPERBLOCK - 1)) + len > NSEEL_RAM_ITEMSPERBLOCK) return dest;
	if ((src_offs&(NSEEL_RAM_ITEMSPERBLOCK - 1)) + len > NSEEL_RAM_ITEMSPERBLOCK) return dest;
	srcptr = __NSEEL_RAMAlloc(blocks, src_offs);
	if (!srcptr || srcptr == &nseel_ramalloc_onfail) return dest;
	destptr = __NSEEL_RAMAlloc(blocks, dest_offs);
	if (!destptr || destptr == &nseel_ramalloc_onfail) return dest;
	WDL_fft_complexmul((WDL_FFT_COMPLEX*)destptr, (WDL_FFT_COMPLEX*)srcptr, (len / 2)&~1);
	return dest;
}
#ifdef EEL_WANT_DOCUMENTATION
static const char *eel_fft_function_reference =
"convolve_c\tdest,src,size\tMultiplies each of size complex pairs in dest by the complex pairs in src. Often used for convolution.\0"
"fft\tbuffer,size\tPerforms a FFT on the data in the local memory buffer at the offset specified by the first parameter. The size of the FFT is specified "
"by the second parameter, which must be 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, or 32768. The outputs are permuted, so if "
"you plan to use them in-order, call fft_permute(buffer, size) before and fft_ipermute(buffer,size) after your in-order use. Your inputs or "
"outputs will need to be scaled down by 1/size, if used.\n"
"Note that fft()/ifft() require real / imaginary input pairs, so a 256 point FFT actually works with 512 items.\n"
"Note that fft()/ifft() must NOT cross a 65,536 item boundary, so be sure to specify the offset accordingly.\0"
"ifft\tbuffer,size\tPerform an inverse FFT. For more information see fft().\0"
"fft_real\tbuffer,size\tPerforms an FFT, but takes size input samples and produces size/2 complex output pairs. Usually used along with fft_permute(size/2). Inputs/outputs will need to be scaled by 0.5/size.\0"
"ifft_real\tbuffer,size\tPerforms an inverse FFT, but takes size/2 complex input pairs and produces size real output values. Usually used along with fft_ipermute(size/2).\0"
"fft_permute\tbuffer,size\tPermute the output of fft() to have bands in-order. See fft() for more information.\0"
"fft_ipermute\tbuffer,size\tPermute the input for ifft(), taking bands from in-order to the order ifft() requires. See fft() for more information.\0"
;
#endif
int nseel_stringsegments_tobuf(char *bufOut, int bufout_sz, eelStringSegmentRec *list) // call with NULL to calculate size, or non-null to generate to buffer (returning size used)
{
	int pos = 0;
	while (list)
	{
		if (!bufOut)
		{
			pos += list->str_len;
		}
		else if (list->str_len > 1)
		{
			if (pos >= bufout_sz) break;
			pos += nseel_filter_escaped_string(bufOut + pos, bufout_sz - pos, list->str_start + 1, list->str_len - 1, list->str_start[0]);
		}
		list = list->_next;
	}
	return pos;
}
#ifdef CUSTOM_CMD
extern void writeCircularStringBuf(char *cmdCur);
#define EEL_STRING_STDOUT_WRITE(x,len) { writeCircularStringBuf(x); }
#else
#ifndef __ANDROID__
#define EEL_STRING_STDOUT_WRITE(x,len) { fwrite(x,len,1,stdout); fflush(stdout); }
#else
#include <android/log.h>
#define EEL_STRING_STDOUT_WRITE(x,len) { __android_log_print(ANDROID_LOG_INFO, "LiveProg", "%s", x); }
#endif
#endif
void Initeel_string_context_state(eel_string_context_state *st)
{
	st->inuse = 0;
	st->slot = 2;
	st->map = (int*)malloc(st->slot * sizeof(s_str));
	st->m_literal_strings = (s_str*)malloc(st->slot * sizeof(s_str));
	for (int i = 0; i < st->slot; i++)
	{
		st->map[i] = 0;
		st->m_literal_strings[i] = 0;
	}
}
void Freeel_string_context_state(eel_string_context_state *st)
{
	for (int i = 0; i < st->slot; i++)
		s_str_destroy(&st->m_literal_strings[i]);
	free(st->map);
	free(st->m_literal_strings);
	st->slot = 0;
	st->inuse = 0;
}
int arySearch(int *array, int N, int x)
{
	for (int i = 0; i < N; i++)
	{
		if (array[i] == x)
			return i;
	}
	return -1;
}
#define FLOIDX 20000
void* GetStringForIndex(eel_string_context_state *st, EEL_F val, int write)
{
	int castedValue = (int)(val + 0.5);
	if (castedValue < FLOIDX)
		return 0;
	int idx = arySearch(st->map, st->slot, castedValue);
	if (idx < 0)
		return 0;
	if (!write)
	{
		s_str *tmp = &st->m_literal_strings[idx];
		const char *s = s_str_c_str(tmp);
		return (void*)s;
	}
	else
		return (void*)&st->m_literal_strings[idx];
}
int AddString(eel_string_context_state *st, char *ns)
{
	const int l = strlen(ns);
	int x;
	for (x = 0; x < st->inuse; x++)
	{
		s_str *tmp = &st->m_literal_strings[x];
		const char *s = s_str_c_str(tmp);
		if (strlen(s) == l && !strcmp(s, ns))
			break;
	}
	if (x < st->inuse)
		free(ns);
	else
	{
		int currentSlot = st->inuse;
		if (currentSlot > (st->slot - 1))
		{
			st->slot++;
			st->m_literal_strings = (s_str*)realloc(st->m_literal_strings, st->slot * sizeof(s_str));
			st->m_literal_strings[st->inuse] = 0;
			st->map = (int*)realloc(st->map, st->slot * sizeof(int));
			st->map[st->inuse] = 0;
		}
		st->m_literal_strings[st->inuse] = s_str_create_from_c_str(ns);
		st->map[st->inuse] = x + FLOIDX;
		st->inuse++;
		free(ns);
	}
	return x + FLOIDX;
}
static EEL_F addStringCallback(void *opaque, eelStringSegmentRec *list)
{
	if (!opaque) return -1.0;
	compileContext *c = (compileContext *)opaque;
	eel_string_context_state *_this = c->m_string_context;
	if (!_this) return -1.0;
	// could probably do a faster implementation using AddRaw() etc but this should also be OK
	int sz = nseel_stringsegments_tobuf(NULL, 0, list);
	char *ns = (char*)malloc(sz + 32);
	memset(ns, 0, sz + 32);
	sz = nseel_stringsegments_tobuf(ns, sz, list) + 1;
	ns = (char*)realloc(ns, sz);
	int id = AddString(_this, ns);
	return (EEL_F)id;
}
static int eel_string_match(compileContext *c, const char *fmt, const char *msg, int match_fmt_pos, int ignorecase, const char *fmt_endptr, const char *msg_endptr, int num_fmt_parms, EEL_F **fmt_parms)
{
	// %d=12345
	// %f=12345[.678]
	// %c=any nonzero char, ascii value
	// %x=12354ab
	// %*, %?, %+, %% literals
	// * ? +  match minimal groups of 0+,1, or 1+ chars
	for (;;)
	{
		if (fmt >= fmt_endptr)
		{
			if (msg >= msg_endptr) return 1;
			return 0; // format ends before matching string
		}

		// if string ends and format is not on a wildcard, early-out to 0
		if (msg >= msg_endptr && *fmt != '*' && *fmt != '%') return 0;

		switch (*fmt)
		{
		case '*':
		case '+':
			// if last char of search pattern, we're done!
			if (fmt + 1 >= fmt_endptr || (fmt[1] == '?' && fmt + 2 >= fmt_endptr)) return *fmt == '*' || msg < msg_endptr;

			if (fmt[0] == '+')  msg++; // skip a character for + . Note that in this case msg[1] is valid, because of the !*msg && *fmt != '*' check above

			fmt++;
			if (*fmt == '?')
			{
				// *? or +? are lazy matches
				fmt++;

				while (msg < msg_endptr && !eel_string_match(c, fmt, msg, match_fmt_pos, ignorecase, fmt_endptr, msg_endptr, num_fmt_parms, fmt_parms)) msg++;
				return msg < msg_endptr;
			}
			else
			{
				// greedy match
				int len = (int)(msg_endptr - msg);
				while (len >= 0 && !eel_string_match(c, fmt, msg + len, match_fmt_pos, ignorecase, fmt_endptr, msg_endptr, num_fmt_parms, fmt_parms)) len--;
				return len >= 0;
			}
			break;
		case '?':
			fmt++;
			msg++;
			break;
		case '%':
		{
			fmt++;
			unsigned short fmt_minlen = 1, fmt_maxlen = 0;
			if (*fmt >= '0' && *fmt <= '9')
			{
				fmt_minlen = *fmt++ - '0';
				while (*fmt >= '0' && *fmt <= '9') fmt_minlen = fmt_minlen * 10 + (*fmt++ - '0');
				fmt_maxlen = fmt_minlen;
			}
			if (*fmt == '-')
			{
				fmt++;
				fmt_maxlen = 0;
				while (*fmt >= '0' && *fmt <= '9') fmt_maxlen = fmt_maxlen * 10 + (*fmt++ - '0');
			}
			const char *dest_varname = NULL;
			if (*fmt == '{')
			{
				dest_varname = ++fmt;
				while (*fmt && fmt < fmt_endptr && *fmt != '}') fmt++;
				if (fmt >= fmt_endptr - 1 || *fmt != '}') return 0; // malformed %{var}s
				fmt++; // skip '}'
			}
			char fmt_char = *fmt++;
			if (!fmt_char) return 0; // malformed
			if (fmt_char == '*' || fmt_char == '?' || fmt_char == '+' || fmt_char == '%')
			{
				if (*msg++ != fmt_char) return 0;
			}
			else if (fmt_char == 'c')
			{
				EEL_F *varOut = NULL;
				EEL_F vv = 0.0;
				if (!dest_varname)
				{
					if (match_fmt_pos < num_fmt_parms)
						varOut = fmt_parms[match_fmt_pos];
					match_fmt_pos++;
				}
				if (msg >= msg_endptr) return 0; // out of chars
				if (varOut)
				{
					if (varOut == &vv) // %{#foo}c
					{
						s_str *wr = (s_str*)GetStringForIndex(c->m_string_context, vv, 1);
						if (wr)
							s_str_destroy(wr);
						*wr = s_str_create_from_c_str(msg);
					}
					else
					{
						*varOut = (EEL_F)*(unsigned char *)msg;
					}
				}
				msg++;
			}
			else
			{
				int len = 0;
				int lazy = 0;
				if (fmt_char >= 'A'&&fmt_char <= 'Z') { lazy = 1; fmt_char += 'a' - 'A'; }
				if (fmt_char == 's')
				{
					len = (int)(msg_endptr - msg);
				}
				else if (fmt_char == 'x')
				{
					while ((msg[len] >= '0' && msg[len] <= '9') || (msg[len] >= 'A' && msg[len] <= 'F') || (msg[len] >= 'a' && msg[len] <= 'f'))
						len++;
				}
				else if (fmt_char == 'f')
				{
					if (msg[len] == '-') len++;
					while (msg[len] >= '0' && msg[len] <= '9') len++;
					if (msg[len] == '.')
					{
						len++;
						while (msg[len] >= '0' && msg[len] <= '9') len++;
					}
				}
				else if (fmt_char == 'd' || fmt_char == 'u' || fmt_char == 'i')
				{
					if (fmt_char != 'u' && msg[len] == '-') len++;
					while (msg[len] >= '0' && msg[len] <= '9') len++;
				}
				else
				{
					// bad format
					return 0;
				}
				if (fmt_maxlen > 0 && len > fmt_maxlen) len = fmt_maxlen;
				if (!dest_varname) match_fmt_pos++;
				if (lazy)
				{
					if (fmt_maxlen<1 || fmt_maxlen>len) fmt_maxlen = len;
					len = fmt_minlen;
					while (len <= fmt_maxlen && !eel_string_match(c, fmt, msg + len, match_fmt_pos, ignorecase, fmt_endptr, msg_endptr, num_fmt_parms, fmt_parms)) len++;
					if (len > fmt_maxlen) return 0;
				}
				else
				{
					while (len >= fmt_minlen && !eel_string_match(c, fmt, msg + len, match_fmt_pos, ignorecase, fmt_endptr, msg_endptr, num_fmt_parms, fmt_parms)) len--;
					if (len < fmt_minlen) return 0;
				}
				EEL_F vv = 0.0;
				EEL_F *varOut = NULL;
				if (!dest_varname)
				{
					if (match_fmt_pos > 0 && match_fmt_pos - 1 < num_fmt_parms)
						varOut = fmt_parms[match_fmt_pos - 1];
				}
				if (varOut)
				{
					if (fmt_char == 's')
					{
						s_str *wr = (s_str*)GetStringForIndex(c->m_string_context, *varOut, 1);
						if (wr)
							s_str_destroy(wr);
						if (wr)
						{
							const char *strS = s_str_c_str(wr);
							int le = strlen(strS);
							if (msg_endptr >= strS && msg_endptr <= strS + le)
							{
#ifdef EEL_STRING_DEBUGOUT
								EEL_STRING_DEBUGOUT("match: destination specifier passed is also haystack, will not update");
#endif
							}
							else if (fmt_endptr >= strS && fmt_endptr <= strS + le)
							{
#ifdef EEL_STRING_DEBUGOUT
								EEL_STRING_DEBUGOUT("match: destination specifier passed is also format, will not update");
#endif
							}
							else
							{
								*wr = s_str_create_from_c_str(msg);
							}
						}
						else
						{
#ifdef EEL_STRING_DEBUGOUT
							EEL_STRING_DEBUGOUT("match: bad destination specifier passed as %d: %f", match_fmt_pos, *varOut);
#endif
						}
					}
					else
					{
						char tmp[128];
						lstrcpyn_safe(tmp, msg, wdl_min(len + 1, (int)sizeof(tmp)));
						if (varOut == &vv)
						{
							s_str *wr = (s_str*)GetStringForIndex(c->m_string_context, vv, 1);
							if (wr)
								s_str_destroy(wr);
							*wr = s_str_create_from_c_str(tmp);
						}
						else
						{
							char *bl = (char*)msg;
							if (fmt_char == 'u')
								*varOut = (EEL_F)strtoul(tmp, &bl, 10);
							else if (fmt_char == 'x')
								*varOut = (EEL_F)strtoul(msg, &bl, 16);
							else
								*varOut = (EEL_F)atof(tmp);
						}
					}
				}
				return 1;
			}
		}
		break;
		default:
			if (ignorecase ? (toupper(*fmt) != toupper(*msg)) : (*fmt != *msg)) return 0;
			fmt++;
			msg++;
			break;
		}
	}
}
static EEL_F NSEEL_CGEN_CALL _eel_match(void *opaque, INT_PTR num_parms, EEL_F **parms)
{
	if (opaque && num_parms >= 2)
	{
		compileContext *c = (compileContext *)opaque;
		const char *fmt = (const char*)GetStringForIndex(c->m_string_context, *parms[0], 0);
		const char *msg = (const char*)GetStringForIndex(c->m_string_context, *parms[1], 0);
		if (fmt && msg)
			return eel_string_match(c, fmt, msg, 0, 0, fmt + strlen(fmt), msg + strlen(msg), (int)num_parms - 2, parms + 2) ? 1.0 : 0.0;
	}
	return 0.0;
}
static EEL_F NSEEL_CGEN_CALL _eel_matchi(void *opaque, INT_PTR num_parms, EEL_F **parms)
{
	if (opaque && num_parms >= 2)
	{
		compileContext *c = (compileContext *)opaque;
		const char *fmt = (const char*)GetStringForIndex(c->m_string_context, *parms[0], 0);
		const char *msg = (const char*)GetStringForIndex(c->m_string_context, *parms[1], 0);
		if (fmt && msg)
			return eel_string_match(opaque, fmt, msg, 0, 1, fmt + strlen(fmt), msg + strlen(msg), (int)num_parms - 2, parms + 2) ? 1.0 : 0.0;
	}
	return 0.0;
}
static int eel_validate_format_specifier(const char *fmt_in, char *typeOut, char *fmtOut, int fmtOut_sz, char *varOut, int varOut_sz)
{
	const char *fmt = fmt_in + 1;
	int state = 0;
	if (fmt_in[0] != '%') return 0; // ugh passed a non-
	*varOut = 0;
	if (fmtOut_sz-- < 2) return 0;
	*fmtOut++ = '%';
	while (*fmt)
	{
		const char c = *fmt++;
		if (fmtOut_sz < 2) return 0;
		if (c == 'f' || c == 'e' || c == 'E' || c == 'g' || c == 'G' || c == 'd' || c == 'u' || c == 'x' || c == 'X' || c == 'c' || c == 'C' || c == 's' || c == 'S' || c == 'i')
		{
			*typeOut = c;
			fmtOut[0] = c;
			fmtOut[1] = 0;
			return (int)(fmt - fmt_in);
		}
		else if (c == '.')
		{
			*fmtOut++ = c; fmtOut_sz--;
			if (state&(2)) break;
			state |= 2;
		}
		else if (c == '+')
		{
			*fmtOut++ = c; fmtOut_sz--;
			if (state&(32 | 16 | 8 | 4)) break;
			state |= 8;
		}
		else if (c == '-' || c == ' ')
		{
			*fmtOut++ = c; fmtOut_sz--;
			if (state&(32 | 16 | 8 | 4)) break;
			state |= 16;
		}
		else if (c >= '0' && c <= '9')
		{
			*fmtOut++ = c; fmtOut_sz--;
			state |= 4;
		}
		else
			break;
	}
	return 0;
}
int eel_format_strings(void *opaque, const char *fmt, const char *fmt_end, char *buf, int buf_sz, int num_fmt_parms, EEL_F **fmt_parms)
{
	int fmt_parmpos = 0;
	char *op = buf;
	while ((fmt_end ? fmt < fmt_end : *fmt) && op < buf + buf_sz - 128)
	{
		if (fmt[0] == '%' && fmt[1] == '%')
		{
			*op++ = '%';
			fmt += 2;
		}
		else if (fmt[0] == '%')
		{
			char ct = 0;
			char fs[128];
			char varname[128];
			const int l = eel_validate_format_specifier(fmt, &ct, fs, sizeof(fs), varname, sizeof(varname));
			if (!l || !ct)
			{
				*op = 0;
				return -1;
			}
			const EEL_F *varptr = NULL;
			if (fmt_parmpos < num_fmt_parms)
				varptr = fmt_parms[fmt_parmpos];
			fmt_parmpos++;
			double v = varptr ? (double)*varptr : 0.0;
			if (ct == 'x' || ct == 'X' || ct == 'd' || ct == 'u' || ct == 'i')
			{
				WDL_snprintf(op, 64, fs, (int)(v));
			}
			else if (ct == 's' || ct == 'S')
			{
				compileContext *c = (compileContext *)opaque;
				const char *str = (const char*)GetStringForIndex(c->m_string_context, v, 0);
				const int maxl = (int)(buf + buf_sz - 2 - op);
				WDL_snprintf(op, maxl, fs, str ? str : "");
			}
			else if (ct == 'c')
			{
				*op++ = (char)(int)v;
				*op = 0;
			}
			else if (ct == 'C')
			{
				const unsigned int iv = (unsigned int)v;
				int bs = 0;
				if (iv & 0xff000000) bs = 24;
				else if (iv & 0x00ff0000) bs = 16;
				else if (iv & 0x0000ff00) bs = 8;
				while (bs >= 0)
				{
					const char c = (char)(iv >> bs);
					*op++ = c ? c : ' ';
					bs -= 8;
				}
				*op = 0;
			}
			else
				WDL_snprintf(op, 64, fs, v);
			while (*op) op++;
			fmt += l;
		}
		else
			*op++ = *fmt++;
	}
	*op = 0;
	return (int)(op - buf);
}
EEL_F NSEEL_CGEN_CALL _eel_printf(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	if (num_param > 0 && opaque)
	{
		compileContext *c = (compileContext *)opaque;
		const char *fmt = (const char*)GetStringForIndex(c->m_string_context, *(parms[0]), 0);
		if (fmt)
		{
			int stringLength = strlen(fmt);
			const int len = eel_format_strings(opaque, fmt, fmt ? (fmt + stringLength) : NULL, c->printfbuf, (int)sizeof(c->printfbuf), (int)num_param - 1, parms + 1);
			if (len > 0)
			{
				EEL_STRING_STDOUT_WRITE(c->printfbuf, len);
				return 1.0;
			}
			else
			{
				const char *badStr = "printf: bad format string";
				EEL_STRING_STDOUT_WRITE(badStr, strlen(badStr));
			}
		}
		else
		{
			const char *badStr = "printf: bad format specifier passed";
			EEL_STRING_STDOUT_WRITE(badStr, strlen(badStr));
		}
	}
	return 0.0;
}
EEL_F NSEEL_CGEN_CALL _eel_sprintf(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	if (num_param > 0 && opaque)
	{
		compileContext *c = (compileContext *)opaque;
		s_str *wr = (s_str*)GetStringForIndex(c->m_string_context, *(parms[0]), 1);
		if (wr)
			s_str_destroy(wr);
		const char *fmt = (const char*)GetStringForIndex(c->m_string_context, *(parms[1]), 0);
		if (wr && fmt)
		{
			char buf[16384];
			int stringLength = strlen(fmt);
			const int len = eel_format_strings(opaque, fmt, fmt ? (fmt + stringLength) : NULL, buf, (int)sizeof(buf), (int)num_param - 2, parms + 2);
			if (len > 0)
			{
				*wr = s_str_create_from_c_str(buf);
				return 1.0;
			}
			else
			{
				const char *badStr = "printf: bad format string";
				EEL_STRING_STDOUT_WRITE(badStr, strlen(badStr));
			}
		}
		else
		{
			const char *badStr = "printf: bad format specifier passed";
			EEL_STRING_STDOUT_WRITE(badStr, strlen(badStr));
		}
	}
	return 0.0;
}
static EEL_F _eel_strcmp_int(const char *a, int a_len, const char *b, int b_len, int ml, int ignorecase)
{
	// binary-safe comparison (at least if a_len>=0 etc)
	int pos = 0;
	for (;;)
	{
		if (ml > 0 && pos == ml) return 0.0;
		const int a_end = a_len >= 0 ? pos == a_len : !a[pos];
		const int b_end = b_len >= 0 ? pos == b_len : !b[pos];
		if (a_end || b_end)
		{
			if (!b_end) return -1.0; // b[pos] is nonzero, a[pos] is zero
			if (!a_end) return 1.0;
			return 0.0;
		}
		char av = a[pos];
		char bv = b[pos];
		if (ignorecase)
		{
			av = toupper(av);
			bv = toupper(bv);
		}
		if (bv > av) return -1.0;
		if (av > bv) return 1.0;
		pos++;
	}
}
static EEL_F NSEEL_CGEN_CALL _eel_strncmp(void *opaque, EEL_F *aa, EEL_F *bb, EEL_F *maxlen)
{
	if (opaque)
	{
		compileContext *c = (compileContext *)opaque;
		const char *a = (const char*)GetStringForIndex(c->m_string_context, *aa, 0);
		const char *b = (const char*)GetStringForIndex(c->m_string_context, *bb, 0);
		if (!a || !b)
		{
			const char *badStr = "strncmp: bad specifier(s)";
			EEL_STRING_STDOUT_WRITE(badStr, strlen(badStr));
		}
		else
		{
			const int ml = maxlen ? (int)*maxlen : -1;
			if (!ml || a == b) return 0; // strncmp(x,y,0) == 0
			return _eel_strcmp_int(a, a ? strlen(a) : -1, b, b ? strlen(b) : -1, ml, 0);
		}
	}
	return -1.0;
}
static EEL_F NSEEL_CGEN_CALL _eel_strnicmp(void *opaque, EEL_F *aa, EEL_F *bb, EEL_F *maxlen)
{
	if (opaque)
	{
		compileContext *c = (compileContext *)opaque;
		const char *a = (const char*)GetStringForIndex(c->m_string_context, *aa, 0);
		const char *b = (const char*)GetStringForIndex(c->m_string_context, *bb, 0);
		if (!a || !b)
		{
			const char *badStr = "strnicmp: bad specifier(s)";
			EEL_STRING_STDOUT_WRITE(badStr, strlen(badStr));
		}
		else
		{
			const int ml = maxlen ? (int)*maxlen : -1;
			if (!ml || a == b) return 0; // strncmp(x,y,0) == 0
			return _eel_strcmp_int(a, a ? strlen(a) : -1, b, b ? strlen(b) : -1, ml, 1);
		}
	}
	return -1.0;
}
static EEL_F NSEEL_CGEN_CALL _eel_strcmp(void *opaque, EEL_F *strOut, EEL_F *fmt_index)
{
	return _eel_strncmp(opaque, strOut, fmt_index, NULL);
}
static EEL_F NSEEL_CGEN_CALL _eel_stricmp(void *opaque, EEL_F *strOut, EEL_F *fmt_index)
{
	return _eel_strnicmp(opaque, strOut, fmt_index, NULL);
}
static EEL_F NSEEL_CGEN_CALL _eel_strlen(void *opaque, EEL_F *fmt_index)
{
	if (opaque)
	{
		compileContext *c = (compileContext *)opaque;
		const char *fmt = (const char*)GetStringForIndex(c->m_string_context, *fmt_index, 0);
		if (fmt)
			return (EEL_F)strlen(fmt);
	}
	return 0.0;
}
static EEL_F NSEEL_CGEN_CALL _eel_delete_all_strings(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	compileContext *c = (compileContext *)opaque;
	eel_string_context_state *state = c->m_string_context;
	Freeel_string_context_state(state);
	Initeel_string_context_state(state);
	return 0.0;
}
int get_float(char *val, double *F)
{
	char *eptr;
	errno = 0;
	double f = strtod(val, &eptr);
	if (eptr != val && errno != ERANGE)
	{
		*F = f;
		return 1;
	}
	return 0;
}
EEL_F* string2FloatArray(char *frArbitraryEqString, int *elements)
{
	char *p = frArbitraryEqString;
	char *counter = frArbitraryEqString;
	int i = 0, count = 0;
	double number;
	while (*p)
	{
		if (get_float(p, &number))
		{
			strtod(p, &p);
			count++;
		}
		else
			p++;
	}
	*elements = count;
	EEL_F *arrayF = (EEL_F*)malloc(count * sizeof(EEL_F));
	while (*counter)
	{
		if (get_float(counter, &number))
		{
			arrayF[i] = (EEL_F)strtod(counter, &counter);
			i++;
		}
		else
			counter++;
	}
	return arrayF;
}
static EEL_F NSEEL_CGEN_CALL _eel_openText2Sink(void *opaque, EEL_F *fn_index, EEL_F *slot)
{
	compileContext *c = (compileContext *)opaque;
	int idx = arySearch(c->sinksMap, c->numSinks, (int)*slot);
	if (idx >= 0)
	{
		if (c->sinksMap[idx])
			free(c->inputSink[idx]);
	}
	else
	{
		c->numSinks++;
		c->sinksMap = (int*)realloc(c->sinksMap, c->numSinks * sizeof(int));
		c->sinksLength = (int*)realloc(c->sinksLength, c->numSinks * sizeof(int));
		c->inputSink = (EEL_F**)realloc(c->inputSink, c->numSinks * sizeof(EEL_F*));
		idx = c->numSinks - 1;
		c->inputSink[idx] = 0;
	}
	const char *filename = (const char*)GetStringForIndex(c->m_string_context, *fn_index, 0);
	char *buffer = 0;
	long length;
	FILE *textFile = fopen(filename, "rb");
	if (textFile)
	{
		fseek(textFile, 0, SEEK_END);
		length = ftell(textFile);
		fseek(textFile, 0, SEEK_SET);
		buffer = (char*)malloc(length + 1);
		if (buffer)
			fread(buffer, 1, length, textFile);
		fclose(textFile);
		buffer[length] = '\0';
	}
	else
		return 0;
	c->sinksMap[idx] = (int)*slot;
	c->inputSink[idx] = string2FloatArray(buffer, &c->sinksLength[idx]);
	if (buffer)
		free(buffer);
	return (EEL_F)c->sinksLength[idx];
}
static EEL_F NSEEL_CGEN_CALL _eel_sink2Val(void *opaque, EEL_F *start1, EEL_F *slot)
{
	compileContext *c = (compileContext*)opaque;
	EEL_F **blocks = c->ram_state.blocks;
	int idx = arySearch(c->sinksMap, c->numSinks, (int)*slot);
	if (idx < 0)
		return -1;
	int offs1 = (int)(*start1 + NSEEL_CLOSEFACTOR);
	EEL_F *arrayPtr = __NSEEL_RAMAlloc(blocks, offs1);
	memcpy(arrayPtr, c->inputSink[idx], c->sinksLength[idx] * sizeof(EEL_F));
	if (!arrayPtr || arrayPtr == &nseel_ramalloc_onfail)
		return 0;
	return 1;
}
static EEL_F NSEEL_CGEN_CALL _eel_delete_all_sinks(void *opaque, INT_PTR num_param, EEL_F **parms)
{
	compileContext *ctx = (compileContext *)opaque;
	if (ctx->numSinks)
	{
		for (int i = 0; i < ctx->numSinks; i++)
			free(ctx->inputSink[i]);
		free(ctx->sinksMap);
		free(ctx->sinksLength);
		free(ctx->inputSink);
		ctx->numSinks = 0;
		ctx->sinksMap = 0;
		ctx->sinksLength = 0;
		ctx->inputSink = 0;
	}
	return 0.0;
}
static EEL_F NSEEL_CGEN_CALL _eel_importFloatArrayFromString(void *opaque, EEL_F *fn_index, EEL_F *pointer)
{
	compileContext *c = (compileContext *)opaque;
	const char *FLTBuf = (const char*)GetStringForIndex(c->m_string_context, *fn_index, 0);
	int offs1 = (int)(*pointer + NSEEL_CLOSEFACTOR);
	EEL_F *userspaceFLT = __NSEEL_RAMAlloc(c->ram_state.blocks, offs1);
	if (!userspaceFLT || userspaceFLT == &nseel_ramalloc_onfail)
		return 0;
	int elements;
	EEL_F *convertedFLT = string2FloatArray((char*)FLTBuf, &elements);
	memcpy(userspaceFLT, convertedFLT, elements * sizeof(EEL_F));
	free(convertedFLT);
	return (EEL_F)elements;
}
void fractionalDelayLine_clear(EEL_F *fdl)
{
	for (int i = 0; i < (int)fdl[2]; i++)
		fdl[5 + i] = 0.;
}
static EEL_F NSEEL_CGEN_CALL fractionalDelayLineInit(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	int offs1 = (int)(*start + NSEEL_CLOSEFACTOR);
	int max_length = (int)(*length + NSEEL_CLOSEFACTOR);
	EEL_F *fdl = __NSEEL_RAMAlloc(blocks, offs1);
	if (!fdl || fdl == &nseel_ramalloc_onfail)
		return 0;
	fdl[0] = 0;
	fdl[1] = 0;
	fdl[2] = max_length;
	fractionalDelayLine_clear(fdl);
	return (EEL_F)(max_length + 5);
}
static EEL_F NSEEL_CGEN_CALL fractionalDelayLineClear(EEL_F **blocks, EEL_F *start)
{
	int offs1 = (int)(*start + NSEEL_CLOSEFACTOR);
	EEL_F *fdl = __NSEEL_RAMAlloc(blocks, offs1);
	if (!fdl || fdl == &nseel_ramalloc_onfail)
		return 0;
	fractionalDelayLine_clear(fdl);
	return 1;
}
static EEL_F NSEEL_CGEN_CALL fractionalDelayLineSetDelay(EEL_F **blocks, EEL_F *start, EEL_F *lg)
{
	int offs1 = (int)(*start + NSEEL_CLOSEFACTOR);
	EEL_F lag = *lg;
	EEL_F *fdl = __NSEEL_RAMAlloc(blocks, offs1);
	if (!fdl || fdl == &nseel_ramalloc_onfail)
		return 0;
	EEL_F outPointer;
	if (lag > fdl[2] - 1.0)
		outPointer = fdl[0] + 1.0; // force delay to max_length
	else
		outPointer = fdl[0] - lag; // read chases write
	while (outPointer < 0)
		outPointer += fdl[2]; // modulo maximum length
	fdl[1] = (EEL_F)((int)outPointer); // integer part
	fdl[3] = outPointer - fdl[1]; // fractional part
	fdl[4] = 1.0 - fdl[3]; // 1.0 - fractional part (more efficient)
	return 1;
}
static EEL_F NSEEL_CGEN_CALL fractionalDelayLineProcess(EEL_F **blocks, EEL_F *start, EEL_F *x)
{
	int offs1 = (int)(*start + NSEEL_CLOSEFACTOR);
	EEL_F *fdl = __NSEEL_RAMAlloc(blocks, offs1);
	if (!fdl || fdl == &nseel_ramalloc_onfail)
		return 0;
	int inPtr = (int)(fdl[0] + NSEEL_CLOSEFACTOR);
	fdl[5 + inPtr++] = *x;
	int len = (int)(fdl[2] + NSEEL_CLOSEFACTOR);
	if (inPtr == len) // Check for end condition
		inPtr -= len;
	fdl[0] = (EEL_F)inPtr;
	int outPtr = (int)(fdl[1] + NSEEL_CLOSEFACTOR);
	EEL_F lastOutput = fdl[5 + outPtr++] * fdl[4]; // first 1/2 of interpolation
	if (outPtr < len) // Check for end condition
		lastOutput += fdl[5 + outPtr] * fdl[3]; // second 1/2 of interpolation
	else
	{
		lastOutput += fdl[5] * fdl[3]; // second 1/2 of interpolation
		outPtr -= len;
	}
	fdl[1] = (EEL_F)outPtr;
	return lastOutput;
}
static EEL_F NSEEL_CGEN_CALL FIRInit(EEL_F **blocks, EEL_F *start, EEL_F *length)
{
	int offs1 = (int)(*start + NSEEL_CLOSEFACTOR);
	int hlen = (int)(*length + NSEEL_CLOSEFACTOR);
	EEL_F *fir = __NSEEL_RAMAlloc(blocks, offs1);
	if (!fir || fir == &nseel_ramalloc_onfail)
		return 0;
	memset(fir, 0, (hlen + 2) * sizeof(EEL_F));
	fir[1] = hlen;
	return (EEL_F)(hlen + 2);
}
static EEL_F NSEEL_CGEN_CALL FIRProcess(EEL_F **blocks, EEL_F *start, EEL_F *x, EEL_F *coe)
{
	int offs1 = (int)(*start + NSEEL_CLOSEFACTOR);
	EEL_F *fir = __NSEEL_RAMAlloc(blocks, offs1);
	if (!fir || fir == &nseel_ramalloc_onfail)
		return 0;
	int offs2 = (int)(*coe + NSEEL_CLOSEFACTOR);
	EEL_F *coeffs = __NSEEL_RAMAlloc(blocks, offs2);
	if (!coeffs || coeffs == &nseel_ramalloc_onfail)
		return 0;
	EEL_F *coeff = coeffs;
	int coeffslength = (int)(fir[1] + NSEEL_CLOSEFACTOR);
	EEL_F *coeff_end = coeffs + coeffslength;
	int pos = (int)(fir[0] + NSEEL_CLOSEFACTOR);
	EEL_F *dline = &fir[2];
	EEL_F *buf_val = dline + pos;
	*buf_val = *x;
	EEL_F y = 0.0f;
	while (buf_val >= dline)
		y += *buf_val-- * *coeff++;
	buf_val = dline + coeffslength - 1;
	while (coeff < coeff_end)
		y += *buf_val-- * *coeff++;
	if (++pos >= coeffslength)
		fir[0] = (EEL_F)0;
	else
		fir[0] = (EEL_F)pos;
	return y;
}
#define NRAND 624
#define MRAND 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
EEL_F NSEEL_CGEN_CALL nseel_int_rand(EEL_F amplitude)
{
	unsigned int y;
	static unsigned int mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */
	static unsigned int mt[NRAND]; /* the array for the state vector  */
	static unsigned int __idx;
	unsigned int mti = __idx;
	if (!mti)
	{
		unsigned int s = 0x4141f00d;
		mt[0] = s & 0xffffffffUL;
		for (mti = 1; mti < NRAND; mti++)
		{
			mt[mti] =
				(1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
			/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
			/* In the previous versions, MSBs of the seed affect   */
			/* only MSBs of the array mt[].                        */
			/* 2002/01/09 modified by Makoto Matsumoto             */
			mt[mti] &= 0xffffffffUL;
			/* for >32 bit machines */
		}
		__idx = NRAND; // mti = N (from loop)
	}
	if (mti >= NRAND) { /* generate N words at one time */
		int kk;
		__idx = 1;

		for (kk = 0; kk < NRAND - MRAND; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MRAND] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk < NRAND - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MRAND - NRAND)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[NRAND - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[NRAND - 1] = mt[MRAND - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}
	else
		__idx++;
	y = mt[mti];
	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);
	return (EEL_F)(y*(1.0 / (double)0xFFFFFFFF) * fabs(amplitude));
}
static functionType fnTable1[] = {
   { "sin",    nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_WONTMAKEDENORMAL, {&sin} },
   { "cos",    nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&cos} },
   { "tan",    nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&tan}  },
   { "sqrt",   nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_WONTMAKEDENORMAL, {&sqrt}, },
   { "log",    nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&log} },
   { "log10",  nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&log10} },{ "asin",   nseel_asm_1pdd,nseel_asm_1pdd_end,  1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&asin}, },
   { "acos",   nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&acos}, },
   { "atan",   nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&atan}, },
   { "atan2",  nseel_asm_2pdd,nseel_asm_2pdd_end, 2 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_TWOPARMSONFPSTACK, {&atan2}, },
   { "hypotFast",  nseel_asm_2pdd,nseel_asm_2pdd_end, 2 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_TWOPARMSONFPSTACK, {&hypotFast}, },
   { "hypot",  nseel_asm_2pdd,nseel_asm_2pdd_end, 2 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_TWOPARMSONFPSTACK, {&hypot}, },
   { "pow",    nseel_asm_2pdd,nseel_asm_2pdd_end, 2 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_TWOPARMSONFPSTACK, {&pow}, },
   { "exp",    nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK, {&exp}, },
   { "abs",    nseel_asm_abs,nseel_asm_abs_end,   1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(0) | BIF_WONTMAKEDENORMAL },
   { "sqr",    nseel_asm_sqr,nseel_asm_sqr_end,   1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(1) },
   { "min",    nseel_asm_min,nseel_asm_min_end,   2 | NSEEL_NPARAMS_FLAG_CONST | BIF_FPSTACKUSE(3) | BIF_WONTMAKEDENORMAL },
   { "max",    nseel_asm_max,nseel_asm_max_end,   2 | NSEEL_NPARAMS_FLAG_CONST | BIF_FPSTACKUSE(3) | BIF_WONTMAKEDENORMAL },
   { "sign",   nseel_asm_sign,nseel_asm_sign_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(2) | BIF_CLEARDENORMAL, },
   { "rand",   nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&nseel_int_rand}, },
   { "round",  nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&round} },
   { "floor",  nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&floor} },
   { "ceil",   nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&ceil} },
   { "expint", nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&expint} },
   { "expintFast",nseel_asm_1pdd,nseel_asm_1pdd_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_CLEARDENORMAL, {&expint_interpolation} },
   { "invsqrt",nseel_asm_1pdd,nseel_asm_1pdd_end,1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(3), {&invsqrt} },
   { "invsqrtFast",nseel_asm_invsqrt,nseel_asm_invsqrt_end,1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(3), },
   { "circshift",_asm_generic3parm,_asm_generic3parm_end,3,{&__NSEEL_circshift},NSEEL_PProc_RAM},
   { "convolve_c",_asm_generic3parm,_asm_generic3parm_end,3,{&eel_convolve_c},NSEEL_PProc_RAM},
   { "fft",_asm_generic2parm,_asm_generic2parm_end,2,{&eel_fft},NSEEL_PProc_RAM},
   { "ifft",_asm_generic2parm,_asm_generic2parm_end,2,{&eel_ifft},NSEEL_PProc_RAM},
   { "fft_real",_asm_generic2parm,_asm_generic2parm_end,2,{&eel_fft_real},NSEEL_PProc_RAM},
   { "ifft_real",_asm_generic2parm,_asm_generic2parm_end,2,{&eel_ifft_real},NSEEL_PProc_RAM},
   { "fft_permute",_asm_generic2parm,_asm_generic2parm_end,2,{&eel_fft_permute},NSEEL_PProc_RAM},
   { "fft_ipermute",_asm_generic2parm,_asm_generic2parm_end,2,{&eel_ifft_permute},NSEEL_PProc_RAM},
   { "__dbg_getstackptr",nseel_asm_dbg_getstackptr,nseel_asm_dbg_getstackptr_end, 1 | NSEEL_NPARAMS_FLAG_CONST | BIF_RETURNSONSTACK | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(1),  },
  {"freembuf",_asm_generic1parm,_asm_generic1parm_end,1,{&__NSEEL_RAM_MemFree},NSEEL_PProc_RAM},
  {"memcpy",_asm_generic3parm,_asm_generic3parm_end,3,{&__NSEEL_RAM_MemCpy},NSEEL_PProc_RAM},
  {"memset",_asm_generic3parm,_asm_generic3parm_end,3,{&__NSEEL_RAM_MemSet},NSEEL_PProc_RAM},
  {"__memtop",_asm_generic1parm,_asm_generic1parm_end,1,{&__NSEEL_RAM_MemTop},NSEEL_PProc_RAM},
  {"mem_set_values",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&__NSEEL_RAM_Mem_SetValues},NSEEL_PProc_RAM},
  {"mem_get_values",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&__NSEEL_RAM_Mem_GetValues},NSEEL_PProc_RAM},
  {"stack_push",nseel_asm_stack_push,nseel_asm_stack_push_end,1 | BIF_FPSTACKUSE(0),{0,},NSEEL_PProc_Stack},
  {"stack_pop",nseel_asm_stack_pop,nseel_asm_stack_pop_end,1 | BIF_FPSTACKUSE(1),{0,},NSEEL_PProc_Stack},
  {"stack_peek",nseel_asm_stack_peek,nseel_asm_stack_peek_end,1 | NSEEL_NPARAMS_FLAG_CONST | BIF_LASTPARMONSTACK | BIF_FPSTACKUSE(0),{0,},NSEEL_PProc_Stack},
  {"stack_exch",nseel_asm_stack_exch,nseel_asm_stack_exch_end,1 | BIF_FPSTACKUSE(1), {0,},NSEEL_PProc_Stack_PeekTop},
  {"sleep",_asm_generic1parm_retd,_asm_generic1parm_retd_end,1 | BIF_RETURNSONSTACK,{&_eel_sleep},NSEEL_PProc_THIS},
  {"time",_asm_generic1parm,_asm_generic1parm_end,1,{&_eel_time},NSEEL_PProc_THIS},
  {"time_precise",_asm_generic1parm,_asm_generic1parm_end,1,{&_eel_time_precise},NSEEL_PProc_THIS},
  {"strlen",_asm_generic1parm_retd,_asm_generic1parm_retd_end,1 | BIF_RETURNSONSTACK,{&_eel_strlen},NSEEL_PProc_THIS},
  {"strcmp",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&_eel_strcmp},NSEEL_PProc_THIS},
  {"match",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&_eel_match},NSEEL_PProc_THIS},
  {"matchi",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&_eel_matchi},NSEEL_PProc_THIS},
  {"stricmp",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&_eel_stricmp},NSEEL_PProc_THIS},
  {"strncmp",_asm_generic3parm_retd,_asm_generic3parm_retd_end,3 | BIF_RETURNSONSTACK,{&_eel_strncmp},NSEEL_PProc_THIS},
  {"strnicmp",_asm_generic3parm_retd,_asm_generic3parm_retd_end,3 | BIF_RETURNSONSTACK,{&_eel_strnicmp},NSEEL_PProc_THIS},
  {"printf",_asm_generic2parm_retd,_asm_generic2parm_retd_end,1 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&_eel_printf},NSEEL_PProc_THIS},
  {"sprintf",_asm_generic2parm_retd,_asm_generic2parm_retd_end,1 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&_eel_sprintf},NSEEL_PProc_THIS},
  {"resetStringContainers",_asm_generic1parm_retd,_asm_generic1parm_retd_end,1 | BIF_TAKES_VARPARM_EX | BIF_RETURNSONSTACK,{&_eel_delete_all_strings},NSEEL_PProc_THIS},
  {"resetAllSinks",_asm_generic1parm_retd,_asm_generic1parm_retd_end,1 | BIF_TAKES_VARPARM_EX | BIF_RETURNSONSTACK,{&_eel_delete_all_sinks},NSEEL_PProc_THIS},
  {"loadFileToSink",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&_eel_openText2Sink},NSEEL_PProc_THIS},
  {"getValuesFromSink",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&_eel_sink2Val},NSEEL_PProc_THIS},
  {"importFLTFromStr",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&_eel_importFloatArrayFromString},NSEEL_PProc_THIS},
  {"stftCheckMemoryRequirement",_asm_generic2parm_retd,_asm_generic2parm_retd_end,1 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&stftCheckMemoryRequirement},NSEEL_PProc_THIS},
  {"stftInit",_asm_generic2parm_retd,_asm_generic2parm_retd_end,1 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&stftInit},NSEEL_PProc_THIS},
  {"stftForward",_asm_generic2parm_retd,_asm_generic2parm_retd_end,1 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&stftForward},NSEEL_PProc_THIS},
  {"stftBackward",_asm_generic2parm_retd,_asm_generic2parm_retd_end,1 | BIF_TAKES_VARPARM | BIF_RETURNSONSTACK,{&stftBackward},NSEEL_PProc_THIS},
  {"FIRInit",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&FIRInit},NSEEL_PProc_RAM},
  {"FIRProcess",_asm_generic3parm_retd,_asm_generic3parm_retd_end,3 | BIF_RETURNSONSTACK,{&FIRProcess},NSEEL_PProc_RAM},
  {"fractionalDelayLineInit",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&fractionalDelayLineInit},NSEEL_PProc_RAM},
  {"fractionalDelayLineClear",_asm_generic1parm_retd,_asm_generic1parm_retd_end,1 | BIF_RETURNSONSTACK,{&fractionalDelayLineClear},NSEEL_PProc_RAM},
  {"fractionalDelayLineSetDelay",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&fractionalDelayLineSetDelay},NSEEL_PProc_RAM},
  {"fractionalDelayLineProcess",_asm_generic2parm_retd,_asm_generic2parm_retd_end,2 | BIF_RETURNSONSTACK,{&fractionalDelayLineProcess},NSEEL_PProc_RAM},
};
void NSEEL_quit(NSEEL_VMCTX _ctx)
{
	compileContext *ctx = (compileContext*)_ctx;
}
//---------------------------------------------------------------------------------------------------------------
static void freeBlocks(llBlock **start)
{
  llBlock *s=*start;
  *start=0;
  while (s)
  {
    llBlock *llB = s->next;
    free(s);
    s=llB;
  }
}
//---------------------------------------------------------------------------------------------------------------
static void *__newBlock(llBlock **start, int size)
{
  llBlock *llb;
  int alloc_size;
  if (*start && (LLB_DSIZE - (*start)->sizeused) >= size)
  {
    void *t=(*start)->block+(*start)->sizeused;
    (*start)->sizeused+=(size+7)&~7;
    return t;
  }
  alloc_size=sizeof(llBlock);
  if ((int)size > LLB_DSIZE) alloc_size += size - LLB_DSIZE;
  llb = (llBlock *)malloc(alloc_size); // grab bigger block if absolutely necessary (heh)
  if (!llb) return NULL;
  llb->sizeused=(size+7)&~7;
  llb->next = *start;  
  *start = llb;
  return llb->block;
}
//---------------------------------------------------------------------------------------------------------------
opcodeRec *nseel_createCompiledValue(compileContext *ctx, EEL_F value)
{
  opcodeRec *r=newOpCode(ctx,NULL,OPCODETYPE_DIRECTVALUE);
  if (r)
  {
    r->parms.dv.directValue = value; 
  }
  return r;
}
opcodeRec *nseel_createCompiledValuePtr(compileContext *ctx, EEL_F *addrValue, const char *namestr)
{
  opcodeRec *r=newOpCode(ctx,namestr,OPCODETYPE_VARPTR);
  if (!r) return 0;
  r->parms.dv.valuePtr=addrValue;
  return r;
}
static int validate_varname_for_function(compileContext *ctx, const char *name)
{
  if (!ctx->function_curName || !ctx->function_globalFlag) return 1;
  if (ctx->function_localTable_Size[2] > 0 && ctx->function_localTable_Names[2])
  {
    char * const * const namelist = ctx->function_localTable_Names[2];
    const int namelist_sz = ctx->function_localTable_Size[2];
    int i;
    const size_t name_len = strlen(name);
    for (i=0;i<namelist_sz;i++) 
    {
      const char *nmchk=namelist[i];
      const size_t l = strlen(nmchk);
      if (l > 1 && nmchk[l-1] == '*')
      {
        if (name_len >= l && !strncmp(nmchk,name,l-1) && name[l-1]=='.')  return 1;
      }
      else
      {
        if (name_len == l && !strcmp(nmchk,name)) return 1;
      }
    }
  }
  return 0;
}
opcodeRec *nseel_resolve_named_symbol(compileContext *ctx, opcodeRec *rec, int parmcnt, int *errOut)
{
  const int isFunctionMode = parmcnt >= 0;
  int rel_prefix_len=0;
  int rel_prefix_idx=-2;
  int i;    
  char match_parmcnt[4]={-1,-1,-1,-1}; // [3] is guess
  unsigned char match_parmcnt_pos=0;
  char *sname = (char *)rec->relname;
  const char *prevent_function_calls = NULL;
  if (errOut) *errOut = 0;
  if (rec->opcodeType != OPCODETYPE_VARPTR || !sname || !sname[0]) return NULL;
  if (ctx->function_curName)
  {
    if (!strncmp(sname,"this.",5))
    {
      rel_prefix_len=5;
      rel_prefix_idx=-1;
    } 
    else if (!strcmp(sname,"this"))
    {
      rel_prefix_len=4;
      rel_prefix_idx=-1;
    } 
    // scan for parameters/local variables before user functions   
    if (rel_prefix_idx < -1 && ctx->function_localTable_Size[0] > 0 && ctx->function_localTable_Names[0] && ctx->function_localTable_ValuePtrs)
    {
      char * const * const namelist = ctx->function_localTable_Names[0];
      const int namelist_sz = ctx->function_localTable_Size[0];
      for (i=0; i < namelist_sz; i++)
      {
        const char *p = namelist[i];
        if (p)
        {
          if (!isFunctionMode && !strncmp(p,sname,NSEEL_MAX_VARIABLE_NAMELEN))
          {
            rec->opcodeType = OPCODETYPE_VARPTRPTR;
            rec->parms.dv.valuePtr=(EEL_F *)(ctx->function_localTable_ValuePtrs+i);
            rec->parms.dv.directValue=0.0;
            return rec;
          }
          else 
          {
            const size_t plen = strlen(p);
            if (plen > 1 && p[plen-1] == '*' && !strncmp(p,sname,plen-1) && ((sname[plen-1] == '.'&&sname[plen]) || !sname[plen-1]))
            {
              rel_prefix_len=(int) (sname[plen-1] ? plen : plen-1);
              rel_prefix_idx=i;
              break;
            }
          }
        }
      }
    }
    // if instance name set, translate sname or sname.* into "this.sname.*"
    if (rel_prefix_idx < -1 && ctx->function_localTable_Size[1] > 0 && ctx->function_localTable_Names[1])
    {
      char * const * const namelist = ctx->function_localTable_Names[1];
      const int namelist_sz = ctx->function_localTable_Size[1];
      const char *full_sname = rec->relname; // include # in checks
      for (i=0; i < namelist_sz; i++)
      {
        const char *p = namelist[i];
        if (p && *p)
        {
          const size_t tl = strlen(p);     
          if (!strncmp(p,full_sname,tl) && (full_sname[tl] == 0 || full_sname[tl] == '.'))
          {
            rel_prefix_len=0; // treat as though this. prefixes is present
            rel_prefix_idx=-1;
            break;
          }
        }
      }
    }
    if (rel_prefix_idx >= -1) 
    {
      ctx->function_usesNamespaces=1;
    }
  } // ctx->function_curName
  if (!isFunctionMode)
  {
    // instance variables
    if (rel_prefix_idx >= -1) 
    {
      rec->opcodeType = OPCODETYPE_VALUE_FROM_NAMESPACENAME;
      rec->namespaceidx = rel_prefix_idx;
      if (rel_prefix_len > 0) 
      {
        memmove(sname, sname+rel_prefix_len, strlen(sname + rel_prefix_len) + 1);
      }
    }
    else 
    {
      // no namespace index, so it must be a global
      if (!validate_varname_for_function(ctx,rec->relname)) 
      {
        if (errOut) *errOut = 1;
        if (ctx->last_error_string[0]) lstrcatn(ctx->last_error_string, ", ", sizeof(ctx->last_error_string));
        snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"global '%s' inaccessible",rec->relname);
        return NULL;
      }
    }
    return rec;
  }
  if (ctx->func_check)
    prevent_function_calls = ctx->func_check(sname,ctx->func_check_user);
  ////////// function mode
  // first off, while() and loop() are special and can't be overridden
  //
  if (parmcnt == 1 && !strcmp("while",sname) && !prevent_function_calls)
  {
    rec->opcodeType = OPCODETYPE_FUNC1;
    rec->fntype = FN_WHILE;
    return rec;
  }
  if (parmcnt == 2 && !strcmp("loop",sname) && !prevent_function_calls)
  {
    rec->opcodeType = OPCODETYPE_FUNC2;
    rec->fntype = FN_LOOP;
    return rec;
  }
  //
  // resolve user function names before builtin functions -- this allows the user to override default functions
  {
    _codeHandleFunctionRec *best=NULL;
    size_t bestlen=0;
    const char * const ourcall = sname+rel_prefix_len;
    const size_t ourcall_len = strlen(ourcall);
    int pass;
    for (pass=0;pass<2;pass++)
    {
      _codeHandleFunctionRec *fr = pass ? ctx->functions_common : ctx->functions_local;
      // sname is [namespace.[ns.]]function, find best match of function that matches the right end   
      while (fr)
      {
        int this_np = fr->num_params;
        const char *thisfunc = fr->fname;
        const size_t thisfunc_len = strlen(thisfunc);
        if (this_np < 1) this_np=1;
        if (thisfunc_len == ourcall_len && !strcmp(thisfunc,ourcall))
        {
          if (this_np == parmcnt)
          {
            bestlen = thisfunc_len;
            best = fr;
            break; // found exact match, finished
          }
          else
          {
            if (match_parmcnt_pos < 3) match_parmcnt[match_parmcnt_pos++] = fr->num_params;
          }
        }
        if (thisfunc_len > bestlen && thisfunc_len < ourcall_len && ourcall[ourcall_len - thisfunc_len - 1] == '.' && !strcmp(thisfunc,ourcall + ourcall_len - thisfunc_len))
        {
          if (this_np == parmcnt) 
          {
            bestlen = thisfunc_len;
            best = fr;
          }
          else
            if (match_parmcnt[3]<0) match_parmcnt[3]=fr->num_params;
        }
        fr=fr->next;
      }
      if (fr) break; // found exact match, finished
    }
    if (best)
    {
      switch (parmcnt)
      {
        case 0:
        case 1: rec->opcodeType = OPCODETYPE_FUNC1; break;
        case 2: rec->opcodeType = OPCODETYPE_FUNC2; break;
        case 3: rec->opcodeType = OPCODETYPE_FUNC3; break;
        default: rec->opcodeType = OPCODETYPE_FUNCX; break;
      }
      if (ourcall != rec->relname) memmove((char *)rec->relname, ourcall, strlen(ourcall)+1);
      if (ctx->function_curName && rel_prefix_idx<0)
      {
        // if no namespace specified, and this.commonprefix.func() called, remove common prefixes and set prefixidx to be this
        const char *p=ctx->function_curName;
        if (*p) p++;
        while (*p && *p != '.')  p++;
        if (*p && p[1]) // we have a dot!
        {
          while (p[1]) p++; // go to last char of string, which doesn't allow possible trailing dot to be checked
          while (--p > ctx->function_curName) // do not check possible leading dot
          {            
            if (*p == '.')
            {
              const size_t cmplen = p+1-ctx->function_curName;
              if (!strncmp(rec->relname,ctx->function_curName,cmplen) && rec->relname[cmplen])
              {
                const char *src=rec->relname + cmplen;
                memmove((char *)rec->relname, src, strlen(src)+1);
                rel_prefix_idx=-1; 
                ctx->function_usesNamespaces=1;
                break;
              }
            }
          }
        }
      }
      if (ctx->function_curName && rel_prefix_idx < -1 && 
          strchr(rec->relname,'.') && !validate_varname_for_function(ctx,rec->relname))
      {
        if (errOut) *errOut = 1;
        if (ctx->last_error_string[0]) lstrcatn(ctx->last_error_string, ", ", sizeof(ctx->last_error_string));
        snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"namespaced function '%s' inaccessible",rec->relname);
        return NULL;
      }
      rec->namespaceidx = rel_prefix_idx;
      rec->fntype = FUNCTYPE_EELFUNC;
      rec->fn = best;
      return rec;
    }    
  }
  if (prevent_function_calls)
  {
    if (ctx->last_error_string[0]) lstrcatn(ctx->last_error_string, ", ", sizeof(ctx->last_error_string));
    snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"'%.30s': %s",sname, prevent_function_calls);
    if (errOut) *errOut = 0;
    return NULL;
  }
  // convert legacy pow() to FN_POW
  if (!strcmp("pow",sname))
  {
    if (parmcnt == 2)
    {
      rec->opcodeType = OPCODETYPE_FUNC2;
      rec->fntype = FN_POW;
      return rec;
    }
    if (match_parmcnt_pos < 3) match_parmcnt[match_parmcnt_pos++] = 2;
  }
  else if (!strcmp("__denormal_likely",sname) || !strcmp("__denormal_unlikely",sname))
  {
    if (parmcnt == 1)
    {
      rec->opcodeType = OPCODETYPE_FUNC1;
      rec->fntype = !strcmp("__denormal_likely",sname) ? FN_DENORMAL_LIKELY : FN_DENORMAL_UNLIKELY;
      return rec;
    }
  }
  for (i = 0; fnTable1 + i; i++)
  {
	  if (i >= (int)(sizeof(fnTable1) / sizeof(fnTable1[0])))
		  break;
    functionType *f = fnTable1 + i;
    if (!strcmp(f->name, sname))
    {
      const int pc_needed=(f->nParams&FUNCTIONTYPE_PARAMETERCOUNTMASK);
      if ((f->nParams&BIF_TAKES_VARPARM_EX)==BIF_TAKES_VARPARM ? (parmcnt >= pc_needed) : (parmcnt == pc_needed))
      {
        rec->fntype = FUNCTYPE_FUNCTIONTYPEREC;
        rec->fn = (void *)f;
        switch (parmcnt)
        {
          case 0:
          case 1: rec->opcodeType = OPCODETYPE_FUNC1; break;
          case 2: rec->opcodeType = OPCODETYPE_FUNC2; break;
          case 3: rec->opcodeType = OPCODETYPE_FUNC3; break;
          default: rec->opcodeType = OPCODETYPE_FUNCX; break;
        }
        return rec;
      }
      if (match_parmcnt_pos < 3) match_parmcnt[match_parmcnt_pos++] = (f->nParams&FUNCTIONTYPE_PARAMETERCOUNTMASK);
    }
  }
  if (ctx->last_error_string[0]) lstrcatn(ctx->last_error_string, ", ", sizeof(ctx->last_error_string));
  if (match_parmcnt[3] >= 0)
  {
    if (match_parmcnt_pos<3) match_parmcnt[match_parmcnt_pos] = match_parmcnt[3];
    match_parmcnt_pos++;
  }
  if (!match_parmcnt_pos)
    snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"'%.30s' undefined",sname);
  else
  {
    int x;
    snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"'%.30s' needs ",sname);
    for (x = 0; x < match_parmcnt_pos; x++)
      snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"%s%d",x==0?"" : x == match_parmcnt_pos-1?" or ":",",match_parmcnt[x]);
    lstrcatn(ctx->last_error_string," parms",sizeof(ctx->last_error_string));
  }
  if (errOut) *errOut = match_parmcnt_pos > 0 ? parmcnt<match_parmcnt[0]?2:(match_parmcnt[0] < 2 ? 4:1) : 0;
  return NULL;
}
opcodeRec *nseel_setCompiledFunctionCallParameters(compileContext *ctx, opcodeRec *fn, opcodeRec *code1, opcodeRec *code2, opcodeRec *code3, opcodeRec *postCode, int *errOut)
{
  opcodeRec *r;
  int np=0,x;
  if (!fn || fn->opcodeType != OPCODETYPE_VARPTR || !fn->relname || !fn->relname[0]) 
  {
    return NULL;
  }
  fn->parms.parms[0] = code1;
  fn->parms.parms[1] = code2;
  fn->parms.parms[2] = code3;
  for (x=0;x<3;x++)
  {
    opcodeRec *prni=fn->parms.parms[x];
    while (prni && np < NSEEL_MAX_EELFUNC_PARAMETERS)
    {
      const int isMP = prni->opcodeType == OPCODETYPE_MOREPARAMS;
      np++;
      if (!isMP) break;
      prni = prni->parms.parms[1];
    }
  }
  r = nseel_resolve_named_symbol(ctx, fn, np<1 ? 1 : np ,errOut);
  if (postCode && r)
  {
    if (code1 && r->opcodeType == OPCODETYPE_FUNC1 && r->fntype == FN_WHILE)
    {
      // change while(x) (postcode) to be 
      // while ((x) ? (postcode;1) : 0);
      r->parms.parms[0] = 
        nseel_createIfElse(ctx,r->parms.parms[0],
                               nseel_createSimpleCompiledFunction(ctx,FN_JOIN_STATEMENTS,2,postCode,nseel_createCompiledValue(ctx,1.0f)),
                               NULL); // NULL defaults to 0.0
    }
    else
    {
      snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"syntax error following function");
      *errOut = -1;
      return NULL;
    }
  }
  return r;
}
eelStringSegmentRec *nseel_createStringSegmentRec(compileContext *ctx, const char *str, int len)
{
  eelStringSegmentRec *r = newTmpBlock(ctx,sizeof(eelStringSegmentRec));
  if (r)
  {
    r->_next=0;
    r->str_start=str;
    r->str_len = len;
  }
  return r;
}
opcodeRec *nseel_eelMakeOpcodeFromStringSegments(compileContext *ctx, eelStringSegmentRec *rec)
{
  if (ctx && ctx->onString)
  {
    return nseel_createCompiledValue(ctx, ctx->onString(ctx->caller_this,rec));
  }
  return NULL;
}
opcodeRec *nseel_createMoreParametersOpcode(compileContext *ctx, opcodeRec *code1, opcodeRec *code2)
{
  opcodeRec *r=code1 && code2 ? newOpCode(ctx,NULL,OPCODETYPE_MOREPARAMS) : NULL;
  if (r)
  {
    r->parms.parms[0] = code1;
    r->parms.parms[1] = code2;
  }
  return r;
}
opcodeRec *nseel_createIfElse(compileContext *ctx, opcodeRec *code1, opcodeRec *code2, opcodeRec *code3)
{
  opcodeRec *r=code1 ? newOpCode(ctx,NULL,OPCODETYPE_FUNC3) : NULL;
  if (r)
  {
    if (!code2) code2 = nseel_createCompiledValue(ctx,0.0);
    if (!code3) code3 = nseel_createCompiledValue(ctx,0.0);
    if (!code2||!code3) return NULL;
    r->fntype = FN_IF_ELSE;
    r->parms.parms[0] = code1;
    r->parms.parms[1] = code2;
    r->parms.parms[2] = code3;
  }
  return r;
}
opcodeRec *nseel_createMemoryAccess(compileContext *ctx, opcodeRec *code1, opcodeRec *code2)
{
  if (code2 && (code2->opcodeType != OPCODETYPE_DIRECTVALUE || code2->parms.dv.directValue != 0.0))
  {
    code1 = nseel_createSimpleCompiledFunction(ctx,FN_ADD,2,code1,code2);
  }
  return nseel_createSimpleCompiledFunction(ctx, FN_MEMORY,1,code1,0);
}
opcodeRec *nseel_createSimpleCompiledFunction(compileContext *ctx, int fn, int np, opcodeRec *code1, opcodeRec *code2)
{
  opcodeRec *r=code1 && (np<2 || code2) ? newOpCode(ctx,NULL,np>=2 ? OPCODETYPE_FUNC2:OPCODETYPE_FUNC1) : NULL;
  if (r)
  {
    r->fntype = fn;
    r->parms.parms[0] = code1;
    r->parms.parms[1] = code2;
    if (fn == FN_JOIN_STATEMENTS)
    {
      r->fn = r; // for joins, fn is temporarily used for tail pointers
      if (code1 && code1->opcodeType == OPCODETYPE_FUNC2 && code1->fntype == fn)
      {
        opcodeRec *t = (opcodeRec *)code1->fn;
        // keep joins in the form of dosomething->morestuff. 
        // in this instance, code1 is previous stuff to do, code2 is new stuff to do
        r->parms.parms[0] = t->parms.parms[1];
        code1->fn = (t->parms.parms[1] = r);
        return code1;
      }
    }
  }
  return r;  
}
// these are bitmasks; on request you can tell what is supported, and compileOpcodes will return one of them
#define RETURNVALUE_IGNORE 0 // ignore return value
#define RETURNVALUE_NORMAL 1 // pointer
#define RETURNVALUE_FPSTACK 2
#define RETURNVALUE_BOOL 4 // P1 is nonzero if true
#define RETURNVALUE_BOOL_REVERSED 8 // P1 is zero if true
static int compileOpcodes(compileContext *ctx, opcodeRec *op, unsigned char *bufOut, int bufOut_len, int *computTable, const namespaceInformation *namespacePathToThis, 
                          int supportedReturnValues, int *rvType, int *fpStackUsage, int *canHaveDenormalOutput);
static unsigned char *compileCodeBlockWithRet(compileContext *ctx, opcodeRec *rec, int *computTableSize, const namespaceInformation *namespacePathToThis, 
                                              int supportedReturnValues, int *rvType, int *fpStackUse, int *canHaveDenormalOutput);
_codeHandleFunctionRec *eel_createFunctionNamespacedInstance(compileContext *ctx, _codeHandleFunctionRec *fr, const char *nameptr)
{
  size_t n;
  _codeHandleFunctionRec *subfr = 
    fr->isCommonFunction ? 
      ctx->isSharedFunctions ? newDataBlock(sizeof(_codeHandleFunctionRec),8) : 
      newCtxDataBlock(sizeof(_codeHandleFunctionRec),8) :  // if common function, but derived version is in non-common context, set ownership to VM rather than us
    newTmpBlock(ctx,sizeof(_codeHandleFunctionRec));
  if (!subfr) return 0;
  // fr points to functionname()'s rec, nameptr to blah.functionname()
  *subfr = *fr;
  n = strlen(nameptr);
  if (n > sizeof(subfr->fname)-1) n=sizeof(subfr->fname)-1;
  memcpy(subfr->fname,nameptr,n);
  subfr->fname[n]=0;
  subfr->next = NULL;
  subfr->startptr=0; // make sure this code gets recompiled (with correct member ptrs) for this instance!
  // subfr->derivedCopies already points to the right place
  fr->derivedCopies = subfr; 
  return subfr;
}
static void combineNamespaceFields(char *nm, const namespaceInformation *namespaceInfo, const char *relname, int thisctx) // nm must be NSEEL_MAX_VARIABLE_NAMELEN+1 bytes
{
  const char *prefix = namespaceInfo ? 
                          thisctx<0 ? (thisctx == -1 ? namespaceInfo->namespacePathToThis : NULL) :  (thisctx < MAX_SUB_NAMESPACES ? namespaceInfo->subParmInfo[thisctx] : NULL)
                        : NULL;
  int lfp = 0, lrn=relname ? (int)strlen(relname) : 0;
  if (prefix) while (prefix[lfp] && prefix[lfp] != ':' && lfp < NSEEL_MAX_VARIABLE_NAMELEN) lfp++;
  if (!relname) relname = "";
  while (*relname == '.') // if relname begins with ., then remove a chunk of context from prefix
  {
    relname++;
    while (lfp>0 && prefix[lfp-1] != '.') lfp--;
    if (lfp>0) lfp--;       
  }
  if (lfp > NSEEL_MAX_VARIABLE_NAMELEN-3) lfp=NSEEL_MAX_VARIABLE_NAMELEN-3;
  if (lfp>0) memcpy(nm,prefix,lfp);
  if (lrn > NSEEL_MAX_VARIABLE_NAMELEN - lfp - (lfp>0)) lrn=NSEEL_MAX_VARIABLE_NAMELEN - lfp - (lfp>0);
  if (lrn > 0)
  {
    if (lfp>0) nm[lfp++] = '.';
    memcpy(nm+lfp,relname,lrn);
    lfp+=lrn;
  }
  nm[lfp++]=0;
}
//---------------------------------------------------------------------------------------------------------------
static void *nseel_getBuiltinFunctionAddress(compileContext *ctx, 
      int fntype, void *fn, 
      NSEEL_PPPROC *pProc, void ***replList, 
      void **endP, int *abiInfo, int preferredReturnValues, const EEL_F *hasConstParm1, const EEL_F *hasConstParm2)
{
  const EEL_F *firstConstParm = hasConstParm1 ? hasConstParm1 : hasConstParm2;
  static void *pow_replptrs[4]={&pow,};      
  switch (fntype)
  {
#define RF(x) *endP = nseel_asm_##x##_end; return (void*)nseel_asm_##x
    case FN_MUL_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(mul_op);
    case FN_DIV_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(div_op);
    case FN_OR_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(or_op);
    case FN_XOR_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(xor_op);
    case FN_AND_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(and_op);
    case FN_MOD_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(mod_op);
    case FN_ADD_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(add_op);
    case FN_SUB_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(sub_op);
    case FN_POW_OP:
      *abiInfo=BIF_LASTPARMONSTACK|BIF_CLEARDENORMAL;
      *replList = pow_replptrs;
    RF(2pdds);
    case FN_POW: 
      *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK;//BIF_FPSTACKUSE(2) might be safe, need to look at pow()'s implementation, but safer bet is to disallow fp stack caching for this expression
      *replList = pow_replptrs;
    RF(2pdd);
    case FN_ADD: 
       *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2)|BIF_WONTMAKEDENORMAL;
        // for x +- non-denormal-constant,  we can set BIF_CLEARDENORMAL
       if (firstConstParm && fabs(*firstConstParm) > 1.0e-10) *abiInfo |= BIF_CLEARDENORMAL;
    RF(add);
    case FN_SUB: 
       *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK|BIF_FPSTACKUSE(2)|BIF_WONTMAKEDENORMAL; 
        // for x +- non-denormal-constant,  we can set BIF_CLEARDENORMAL
       if (firstConstParm && fabs(*firstConstParm) > 1.0e-10) *abiInfo |= BIF_CLEARDENORMAL;
    RF(sub);
    case FN_MULTIPLY: 
        *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2); 
         // for x*constant-greater-than-eq-1, we can set BIF_WONTMAKEDENORMAL
        if (firstConstParm && fabs(*firstConstParm) >= 1.0) *abiInfo |= BIF_WONTMAKEDENORMAL;
    RF(mul);
    case FN_DIVIDE: 
        *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK|BIF_FPSTACKUSE(2); 
        // for x/constant-less-than-eq-1, we can set BIF_WONTMAKEDENORMAL
        if (firstConstParm && fabs(*firstConstParm) <= 1.0) *abiInfo |= BIF_WONTMAKEDENORMAL;
    RF(div);
    case FN_MOD:
      *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK|BIF_FPSTACKUSE(1)|BIF_CLEARDENORMAL;
    RF(mod);
    case FN_ASSIGN:
      *abiInfo = BIF_FPSTACKUSE(1)|BIF_CLEARDENORMAL;
    RF(assign);
    case FN_AND: *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL; RF(and);
    case FN_OR: *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL; RF(or);
    case FN_XOR:
      *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(xor);
    case FN_SHR:
      *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(shr);
    case FN_SHL:
      *abiInfo = BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK|BIF_FPSTACKUSE(2)|BIF_CLEARDENORMAL;
    RF(shl);
	case FN_NOTNOT: *abiInfo = BIF_LASTPARM_ASBOOL | BIF_RETURNSBOOL | BIF_FPSTACKUSE(1); RF(bnotnot);
    case FN_UMINUS: *abiInfo = BIF_RETURNSONSTACK|BIF_LASTPARMONSTACK|BIF_WONTMAKEDENORMAL; RF(uminus);
    case FN_NOT: *abiInfo = BIF_LASTPARM_ASBOOL|BIF_RETURNSBOOL|BIF_FPSTACKUSE(1); RF(bnot);
    case FN_EQ:
      *abiInfo = BIF_TWOPARMSONFPSTACK_LAZY|BIF_RETURNSBOOL|BIF_FPSTACKUSE(2);
    RF(equal);
    case FN_EQ_EXACT:
      *abiInfo=BIF_TWOPARMSONFPSTACK_LAZY|BIF_RETURNSBOOL|BIF_FPSTACKUSE(2);
    RF(equal_exact);
    case FN_NE:
      *abiInfo=BIF_TWOPARMSONFPSTACK_LAZY|BIF_RETURNSBOOL|BIF_FPSTACKUSE(2);
    RF(notequal);
    case FN_NE_EXACT:
      *abiInfo=BIF_TWOPARMSONFPSTACK_LAZY|BIF_RETURNSBOOL|BIF_FPSTACKUSE(2);
    RF(notequal_exact);
    case FN_LOGICAL_AND:
      *abiInfo = BIF_RETURNSBOOL;
    RF(band);
    case FN_LOGICAL_OR:
      *abiInfo = BIF_RETURNSBOOL;
    RF(bor);
#ifdef GLUE_HAS_FXCH
    case FN_GT:
      *abiInfo = BIF_TWOPARMSONFPSTACK|BIF_RETURNSBOOL|BIF_FPSTACKUSE(2);
    RF(above);
    case FN_GTE:
      *abiInfo = BIF_TWOPARMSONFPSTACK|BIF_RETURNSBOOL|BIF_REVERSEFPORDER|BIF_FPSTACKUSE(2);
    RF(beloweq);
    case FN_LT:
      *abiInfo = BIF_TWOPARMSONFPSTACK|BIF_RETURNSBOOL|BIF_REVERSEFPORDER|BIF_FPSTACKUSE(2);
    RF(above);
    case FN_LTE:
      *abiInfo = BIF_TWOPARMSONFPSTACK|BIF_RETURNSBOOL|BIF_FPSTACKUSE(2);
    RF(beloweq);
#else
    case FN_GT:
      *abiInfo = BIF_RETURNSBOOL|BIF_LASTPARMONSTACK;
    RF(above);
    case FN_GTE:
      *abiInfo = BIF_RETURNSBOOL|BIF_LASTPARMONSTACK;
    RF(aboveeq);
    case FN_LT:
      *abiInfo = BIF_RETURNSBOOL|BIF_LASTPARMONSTACK;
    RF(below);
    case FN_LTE:
      *abiInfo = BIF_RETURNSBOOL|BIF_LASTPARMONSTACK;
    RF(beloweq);
#endif
#undef RF
#define RF(x) *endP = _asm_##x##_end; return (void*)_asm_##x
    case FN_MEMORY:
      {
        static void *replptrs[4]={&__NSEEL_RAMAlloc,};      
        *replList = replptrs;
        *abiInfo = BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(1)|BIF_CLEARDENORMAL;
        #ifdef GLUE_MEM_NEEDS_PPROC
          *pProc = NSEEL_PProc_RAM;
        #endif
        RF(megabuf);
      }
    break;
#undef RF
    case FUNCTYPE_FUNCTIONTYPEREC:
      if (fn)
      {
        functionType *p=(functionType *)fn;
        // if prefers fpstack or int, or ignoring value, then use fp-stack versions
        if ((preferredReturnValues&(RETURNVALUE_BOOL|RETURNVALUE_FPSTACK)) || !preferredReturnValues)
        {
          static functionType min2={ "min",    nseel_asm_min_fp,nseel_asm_min_fp_end,   2|NSEEL_NPARAMS_FLAG_CONST|BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2)|BIF_WONTMAKEDENORMAL };
          static functionType max2={ "max",    nseel_asm_max_fp,nseel_asm_max_fp_end,   2|NSEEL_NPARAMS_FLAG_CONST|BIF_RETURNSONSTACK|BIF_TWOPARMSONFPSTACK_LAZY|BIF_FPSTACKUSE(2)|BIF_WONTMAKEDENORMAL };
          if (p->afunc == (void*)nseel_asm_min) p = &min2;
          else if (p->afunc == (void*)nseel_asm_max) p = &max2;
        }
        *replList=p->replptrs;
        *pProc=p->pProc;
        *endP = p->func_e;
        *abiInfo = p->nParams & BIF_NPARAMS_MASK;
        if (firstConstParm)
        {
          const char *name=p->name;
          if (!strcmp(name,"min") && *firstConstParm < -1.0e-10) *abiInfo |= BIF_CLEARDENORMAL;
          else if (!strcmp(name,"max") && *firstConstParm > 1.0e-10) *abiInfo |= BIF_CLEARDENORMAL;
        }
        return p->afunc; 
      }
    break;
  }
  return 0;
}
static void *nseel_getEELFunctionAddress(compileContext *ctx, 
      opcodeRec *op,
      int *customFuncParmSize, int *customFuncLocalStorageSize,
      EEL_F ***customFuncLocalStorage, int *computTableTop, 
      void **endP, int *isRaw, int wantCodeGenerated,
      const namespaceInformation *namespacePathToThis, int *rvMode, int *fpStackUse, int *canHaveDenormalOutput,
      opcodeRec **ordered_parmptrs, int num_ordered_parmptrs      
      ) // if wantCodeGenerated is false, can return bogus pointers in raw mode
{
  _codeHandleFunctionRec *fn = (_codeHandleFunctionRec*)op->fn;
  namespaceInformation local_namespace={NULL};
  char prefix_buf[NSEEL_MAX_VARIABLE_NAMELEN+1], nm[NSEEL_MAX_FUNCSIG_NAME+1];
  if (!fn) return NULL;
  // op->relname ptr is [whatever.]funcname
  if (fn->parameterAsNamespaceMask || fn->usesNamespaces)
  {
    if (wantCodeGenerated)
    {
      char *p = prefix_buf;
      combineNamespaceFields(nm,namespacePathToThis,op->relname,op->namespaceidx);
      lstrcpyn_safe(prefix_buf,nm,sizeof(prefix_buf));
      local_namespace.namespacePathToThis = prefix_buf;
      // nm is full path of function, prefix_buf will be the path not including function name (unless function name only)
      while (*p) p++;
      while (p >= prefix_buf && *p != '.') p--;
      if (p > prefix_buf) *p=0;
    }
    if (fn->parameterAsNamespaceMask)
    {
      int x;
      for(x=0;x<MAX_SUB_NAMESPACES && x < fn->num_params;x++)
      {
        if (fn->parameterAsNamespaceMask & (((unsigned int)1)<<x))
        {
          if (wantCodeGenerated)
          {
            const char *rn=NULL;
            char tmp[NSEEL_MAX_VARIABLE_NAMELEN+1];
            if (x < num_ordered_parmptrs && ordered_parmptrs[x]) 
            {
              if (ordered_parmptrs[x]->opcodeType == OPCODETYPE_VARPTR) 
              {
                rn=ordered_parmptrs[x]->relname;
              }
              else if (ordered_parmptrs[x]->opcodeType == OPCODETYPE_VALUE_FROM_NAMESPACENAME)
              {
                const char *p=ordered_parmptrs[x]->relname;
                if (*p == '#') p++;
                combineNamespaceFields(tmp,namespacePathToThis,p,ordered_parmptrs[x]->namespaceidx);
                rn = tmp;
              }
            }
            if (!rn) 
            {
              // todo: figure out how to give correct line number/offset (ugh)
              WDL_snprintf(ctx->last_error_string,sizeof(ctx->last_error_string),"parameter %d to %s() must be namespace",x+1,fn->fname);
              return NULL;
            }
            lstrcatn(nm,":",sizeof(nm));
            local_namespace.subParmInfo[x] = nm+strlen(nm);
            lstrcatn(nm,rn,sizeof(nm));
          }         
          ordered_parmptrs[x] = NULL; // prevent caller from bothering generating parameters
        }
      }
    }
    if (wantCodeGenerated)
    {
      _codeHandleFunctionRec *fr = fn;
      // find namespace-adjusted function (if generating code, otherwise assume size is the same)
      fn = 0; // if this gets re-set, it will be the new function
      while (fr && !fn)
      {
        if (!strcmp(fr->fname,nm)) fn = fr;
        fr=fr->derivedCopies;
      }
      if (!fn) // generate copy of function
      {
        fn = eel_createFunctionNamespacedInstance(ctx,(_codeHandleFunctionRec*)op->fn,nm);
      }
    }
  }
  if (!fn) return NULL;
  if (!fn->startptr && fn->opcodes && fn->startptr_size > 0)
  {
    int sz;
    fn->tmpspace_req=0;
    fn->rvMode = RETURNVALUE_IGNORE;
    fn->canHaveDenormalOutput=0;
    sz=compileOpcodes(ctx,fn->opcodes,NULL,128*1024*1024,&fn->tmpspace_req,wantCodeGenerated ? &local_namespace : NULL,RETURNVALUE_NORMAL|RETURNVALUE_FPSTACK,&fn->rvMode,&fn->fpStackUsage,&fn->canHaveDenormalOutput);
    if (!wantCodeGenerated)
    {
      // don't compile anything for now, just give stats
      if (computTableTop) *computTableTop += fn->tmpspace_req;
      *customFuncParmSize = fn->num_params;
      *customFuncLocalStorage = fn->localstorage;
      *customFuncLocalStorageSize = fn->localstorage_size;
      *rvMode = fn->rvMode;
      *fpStackUse = fn->fpStackUsage;
      if (canHaveDenormalOutput) *canHaveDenormalOutput=fn->canHaveDenormalOutput;
      if (sz <= NSEEL_MAX_FUNCTION_SIZE_FOR_INLINE && !(0&OPTFLAG_NO_INLINEFUNC))
      {
        *isRaw = 1;
        *endP = ((char *)1) + sz;
        return (char *)1;
      }
      *endP = (void*)nseel_asm_fcall_end;
      return (void*)nseel_asm_fcall;
    }
    if (sz <= NSEEL_MAX_FUNCTION_SIZE_FOR_INLINE && !(0&OPTFLAG_NO_INLINEFUNC))
    {
      void *p=newTmpBlock(ctx,sz);
      fn->tmpspace_req=0;
      if (p)
      {
        fn->canHaveDenormalOutput=0;
        if (fn->isCommonFunction) ctx->isGeneratingCommonFunction++;
        sz=compileOpcodes(ctx,fn->opcodes,(unsigned char*)p,sz,&fn->tmpspace_req,&local_namespace,RETURNVALUE_NORMAL|RETURNVALUE_FPSTACK,&fn->rvMode,&fn->fpStackUsage,&fn->canHaveDenormalOutput);
        if (fn->isCommonFunction) ctx->isGeneratingCommonFunction--;
        // recompile function with native context pointers
        if (sz>0)
        {
          fn->startptr_size=sz;
          fn->startptr=p;
        }
      }
    }
    else
    {
      unsigned char *codeCall;
      fn->tmpspace_req=0;
      fn->fpStackUsage=0;
      fn->canHaveDenormalOutput=0;
      if (fn->isCommonFunction) ctx->isGeneratingCommonFunction++;
      codeCall=compileCodeBlockWithRet(ctx,fn->opcodes,&fn->tmpspace_req,&local_namespace,RETURNVALUE_NORMAL|RETURNVALUE_FPSTACK,&fn->rvMode,&fn->fpStackUsage,&fn->canHaveDenormalOutput);
      if (fn->isCommonFunction) ctx->isGeneratingCommonFunction--;
      if (codeCall)
      {
        void *f=GLUE_realAddress(nseel_asm_fcall,nseel_asm_fcall_end,&sz);
        fn->startptr = newTmpBlock(ctx,sz);
        if (fn->startptr)
        {
          memcpy(fn->startptr,f,sz);
          EEL_GLUE_set_immediate(fn->startptr,(INT_PTR)codeCall);
          fn->startptr_size = sz;
        }
      }
    }
  }
  if (fn->startptr)
  {
    if (computTableTop) *computTableTop += fn->tmpspace_req;
    *customFuncParmSize = fn->num_params;
    *customFuncLocalStorage = fn->localstorage;
    *customFuncLocalStorageSize = fn->localstorage_size;
    *rvMode = fn->rvMode;
    *fpStackUse = fn->fpStackUsage;
    if (canHaveDenormalOutput) *canHaveDenormalOutput= fn->canHaveDenormalOutput;
    *endP = (char*)fn->startptr + fn->startptr_size;
    *isRaw=1;
    return fn->startptr;
  }
  return 0;
}
// returns true if does something (other than calculating and throwing away a value)
static char optimizeOpcodes(compileContext *ctx, opcodeRec *op, int needsResult)
{
  opcodeRec *lastJoinOp=NULL;
  char retv, retv_parm[3], joined_retv=0;
  while (op && op->opcodeType == OPCODETYPE_FUNC2 && op->fntype == FN_JOIN_STATEMENTS)
  {
    if (!optimizeOpcodes(ctx,op->parms.parms[0], 0) || OPCODE_IS_TRIVIAL(op->parms.parms[0]))
    {
      // direct value, can skip ourselves
      memcpy(op,op->parms.parms[1],sizeof(*op));
    }
    else
    {
      joined_retv |= 1;
      lastJoinOp = op;
      op = op->parms.parms[1];
    }
  }
goto start_over;
#define RESTART_DIRECTVALUE(X) { op->parms.dv.directValue = (X); goto start_over_directvalue; }
start_over_directvalue:
  op->opcodeType = OPCODETYPE_DIRECTVALUE;
  op->parms.dv.valuePtr=NULL;
start_over: // when an opcode changed substantially in optimization, goto here to reprocess it
  retv = retv_parm[0]=retv_parm[1]=retv_parm[2]=0;
  if (!op || // should never really happen
      OPCODE_IS_TRIVIAL(op) || // should happen often (vars)
      op->opcodeType < 0 || op->opcodeType >= OPCODETYPE_INVALID // should never happen (assert would be appropriate heh)
      ) return joined_retv;
  if (!needsResult)
  {
    if (op->fntype == FUNCTYPE_EELFUNC) 
    {
      needsResult=1; // assume eel functions are non-const for now
    }
    else if (op->fntype == FUNCTYPE_FUNCTIONTYPEREC)
    {
      functionType  *pfn = (functionType *)op->fn;
      if (!pfn || !(pfn->nParams&NSEEL_NPARAMS_FLAG_CONST)) needsResult=1;
    }
    else if (op->fntype >= FN_NONCONST_BEGIN && op->fntype < FUNCTYPE_SIMPLEMAX)
    {
      needsResult=1;
    }
  }
  if (op->opcodeType>=OPCODETYPE_FUNC2) retv_parm[1] = optimizeOpcodes(ctx,op->parms.parms[1], needsResult);
  if (op->opcodeType>=OPCODETYPE_FUNC3) retv_parm[2] = optimizeOpcodes(ctx,op->parms.parms[2], needsResult);
  retv_parm[0] = optimizeOpcodes(ctx,op->parms.parms[0], needsResult || 
      (FNPTR_HAS_CONDITIONAL_EXEC(op) && (retv_parm[1] || retv_parm[2] || op->opcodeType <= OPCODETYPE_FUNC1)) );
  if (op->opcodeType != OPCODETYPE_MOREPARAMS)
  {
    if (op->fntype >= 0 && op->fntype < FUNCTYPE_SIMPLEMAX)
    {
      if (op->opcodeType == OPCODETYPE_FUNC1) // within FUNCTYPE_SIMPLE
      {
        if (op->parms.parms[0]->opcodeType == OPCODETYPE_DIRECTVALUE)
        {
          switch (op->fntype)
          {
            case FN_NOTNOT: RESTART_DIRECTVALUE(fabs(op->parms.parms[0]->parms.dv.directValue)>=NSEEL_CLOSEFACTOR ? 1.0 : 0.0);
            case FN_NOT:    RESTART_DIRECTVALUE(fabs(op->parms.parms[0]->parms.dv.directValue)>=NSEEL_CLOSEFACTOR ? 0.0 : 1.0);
            case FN_UMINUS: RESTART_DIRECTVALUE(- op->parms.parms[0]->parms.dv.directValue);
          }
        }
        else if (op->fntype == FN_NOT || op->fntype == FN_NOTNOT)
        {
          if (op->parms.parms[0]->opcodeType == OPCODETYPE_FUNC1)
          {
            switch (op->parms.parms[0]->fntype)
            {
              case FN_UMINUS:
              case FN_NOTNOT: // ignore any NOTNOTs UMINUS or UPLUS, they would have no effect anyway
                op->parms.parms[0] = op->parms.parms[0]->parms.parms[0];
              goto start_over;
              case FN_NOT:
                op->fntype = op->fntype==FN_NOT ? FN_NOTNOT : FN_NOT; // switch between FN_NOT and FN_NOTNOT
                op->parms.parms[0] = op->parms.parms[0]->parms.parms[0];
              goto start_over;
            }
          }
          else if (op->parms.parms[0]->opcodeType == OPCODETYPE_FUNC2)
          {
            int repl_type = -1;
            switch (op->parms.parms[0]->fntype)
            {
              case FN_EQ: repl_type = FN_NE; break;
              case FN_NE: repl_type = FN_EQ; break;
              case FN_EQ_EXACT: repl_type = FN_NE_EXACT; break;
              case FN_NE_EXACT: repl_type = FN_EQ_EXACT; break;
              case FN_LT:  repl_type = FN_GTE; break;
              case FN_LTE: repl_type = FN_GT; break;
              case FN_GT:  repl_type = FN_LTE; break;
              case FN_GTE: repl_type = FN_LT; break;
            }
            if (repl_type != -1)
            {
              const int oldtype = op->fntype;
              memcpy(op,op->parms.parms[0],sizeof(*op));
              if (oldtype == FN_NOT) op->fntype = repl_type;
              goto start_over;
            }
          }
        }
      }
      else if (op->opcodeType == OPCODETYPE_FUNC2)  // within FUNCTYPE_SIMPLE
      {
        const int dv0 = op->parms.parms[0]->opcodeType == OPCODETYPE_DIRECTVALUE;
        const int dv1 = op->parms.parms[1]->opcodeType == OPCODETYPE_DIRECTVALUE;
        if (dv0 && dv1)
        {
          int reval = -1;
          switch (op->fntype)
          {
            case FN_MOD:
              {
                int a = (int) op->parms.parms[1]->parms.dv.directValue;
                if (a) 
                {
                  a = (int) op->parms.parms[0]->parms.dv.directValue % a;
                  if (a<0) a=-a;
                }
                RESTART_DIRECTVALUE((EEL_F)a);
              }
            break;
            case FN_SHL:      RESTART_DIRECTVALUE(((int)op->parms.parms[0]->parms.dv.directValue) << ((int)op->parms.parms[1]->parms.dv.directValue));
            case FN_SHR:      RESTART_DIRECTVALUE(((int)op->parms.parms[0]->parms.dv.directValue) >> ((int)op->parms.parms[1]->parms.dv.directValue));
            case FN_POW:      RESTART_DIRECTVALUE(pow(op->parms.parms[0]->parms.dv.directValue, op->parms.parms[1]->parms.dv.directValue));
            case FN_DIVIDE:   RESTART_DIRECTVALUE(op->parms.parms[0]->parms.dv.directValue / op->parms.parms[1]->parms.dv.directValue);
            case FN_MULTIPLY: RESTART_DIRECTVALUE(op->parms.parms[0]->parms.dv.directValue * op->parms.parms[1]->parms.dv.directValue);
            case FN_ADD:      RESTART_DIRECTVALUE(op->parms.parms[0]->parms.dv.directValue + op->parms.parms[1]->parms.dv.directValue);
            case FN_SUB:      RESTART_DIRECTVALUE(op->parms.parms[0]->parms.dv.directValue - op->parms.parms[1]->parms.dv.directValue);
            case FN_AND:      RESTART_DIRECTVALUE((double) (((WDL_INT64)op->parms.parms[0]->parms.dv.directValue) & ((WDL_INT64)op->parms.parms[1]->parms.dv.directValue)));
            case FN_OR:       RESTART_DIRECTVALUE((double) (((WDL_INT64)op->parms.parms[0]->parms.dv.directValue) | ((WDL_INT64)op->parms.parms[1]->parms.dv.directValue)));
            case FN_XOR:      RESTART_DIRECTVALUE((double) (((WDL_INT64)op->parms.parms[0]->parms.dv.directValue) ^ ((WDL_INT64)op->parms.parms[1]->parms.dv.directValue)));
            case FN_EQ:       reval = fabs(op->parms.parms[0]->parms.dv.directValue - op->parms.parms[1]->parms.dv.directValue) < NSEEL_CLOSEFACTOR; break;
            case FN_NE:       reval = fabs(op->parms.parms[0]->parms.dv.directValue - op->parms.parms[1]->parms.dv.directValue) >= NSEEL_CLOSEFACTOR; break;
            case FN_EQ_EXACT: reval = op->parms.parms[0]->parms.dv.directValue == op->parms.parms[1]->parms.dv.directValue; break;
            case FN_NE_EXACT: reval = op->parms.parms[0]->parms.dv.directValue != op->parms.parms[1]->parms.dv.directValue; break;
            case FN_LT:       reval = op->parms.parms[0]->parms.dv.directValue < op->parms.parms[1]->parms.dv.directValue; break;
            case FN_LTE:      reval = op->parms.parms[0]->parms.dv.directValue <= op->parms.parms[1]->parms.dv.directValue; break;
            case FN_GT:       reval = op->parms.parms[0]->parms.dv.directValue > op->parms.parms[1]->parms.dv.directValue; break;
            case FN_GTE:      reval = op->parms.parms[0]->parms.dv.directValue >= op->parms.parms[1]->parms.dv.directValue; break;
            case FN_LOGICAL_AND: reval = fabs(op->parms.parms[0]->parms.dv.directValue) >= NSEEL_CLOSEFACTOR && fabs(op->parms.parms[1]->parms.dv.directValue) >= NSEEL_CLOSEFACTOR; break;
            case FN_LOGICAL_OR:  reval = fabs(op->parms.parms[0]->parms.dv.directValue) >= NSEEL_CLOSEFACTOR || fabs(op->parms.parms[1]->parms.dv.directValue) >= NSEEL_CLOSEFACTOR; break;
          }
          if (reval >= 0) RESTART_DIRECTVALUE((EEL_F) reval);
        }
        else if (dv0 || dv1)
        {
          double dvalue = op->parms.parms[!dv0]->parms.dv.directValue;
          switch (op->fntype)
          {
            case FN_OR:
            case FN_XOR:
              if (!(WDL_INT64)dvalue)
              {
                // replace with or0
                static functionType fr={"or0",nseel_asm_or0, nseel_asm_or0_end, 1|NSEEL_NPARAMS_FLAG_CONST|BIF_LASTPARMONSTACK|BIF_RETURNSONSTACK|BIF_CLEARDENORMAL, {0}, NULL};
                op->opcodeType = OPCODETYPE_FUNC1;
                op->fntype = FUNCTYPE_FUNCTIONTYPEREC;
                op->fn = &fr;
                if (dv0) op->parms.parms[0] = op->parms.parms[1];
                goto start_over;
              }
            break;
            case FN_SUB:
              if (dv0) 
              {
                if (dvalue == 0.0)
                {
                  op->opcodeType = OPCODETYPE_FUNC1;
                  op->fntype = FN_UMINUS;
                  op->parms.parms[0] = op->parms.parms[1];
                  goto start_over;
                }
                break;
              }
              // fall through, if dv1 we can remove +0.0
            case FN_ADD:
              if (dvalue == 0.0) 
              {
                memcpy(op,op->parms.parms[!!dv0],sizeof(*op));
                goto start_over;
              }
            break;
            case FN_AND:
              if ((WDL_INT64)dvalue) break;
              dvalue = 0.0; // treat x&0 as x*0, which optimizes to 0
              // fall through
            case FN_MULTIPLY:
              if (dvalue == 0.0) // remove multiply by 0.0 (using 0.0 direct value as replacement), unless the nonzero side did something
              {
                if (!retv_parm[!!dv0]) 
                {
                  memcpy(op,op->parms.parms[!dv0],sizeof(*op)); // set to 0 if other action wouldn't do anything
                  goto start_over;
                }
                else
                {
                  // this is 0.0 * oldexpressionthatmustbeprocessed or oldexpressionthatmustbeprocessed*0.0
                  op->fntype = FN_JOIN_STATEMENTS;
                  if (dv0) // 0.0*oldexpression, reverse the order so that 0 is returned
                  {
                    // set to (oldexpression;0)
                    opcodeRec *tmp = op->parms.parms[1];
                    op->parms.parms[1] = op->parms.parms[0];
                    op->parms.parms[0] = tmp;
                  }
                  goto start_over;
                }
              }
              else if (dvalue == 1.0) // remove multiply by 1.0 (using non-1.0 value as replacement)
              {
                memcpy(op,op->parms.parms[!!dv0],sizeof(*op));
                goto start_over;
              }
            break;
            case FN_POW:
              if (dv1)
              {
                // x^0 = 1
                if (fabs(dvalue) < 1e-30)
                {
                  RESTART_DIRECTVALUE(1.0);
                }
                // x^1 = x
                if (fabs(dvalue-1.0) < 1e-30)
                {
                  memcpy(op,op->parms.parms[0],sizeof(*op));
                  goto start_over;
                }
              }
              else if (dv0)
              {
                // pow(constant, x) = exp((x) * ln(constant)), if constant>0
                // opcodeRec *parm0 = op->parms.parms[0];
                if (dvalue > 0.0)
                {
                  static functionType expcpy={ "exp",    nseel_asm_1pdd,nseel_asm_1pdd_end,   1|NSEEL_NPARAMS_FLAG_CONST|BIF_RETURNSONSTACK|BIF_LASTPARMONSTACK, {&exp}, };
                  // 1^x = 1
                  if (fabs(dvalue-1.0) < 1e-30)
                  {
                    RESTART_DIRECTVALUE(1.0);
                  }
                  dvalue=log(dvalue);
                  if (fabs(dvalue-1.0) < 1e-9)
                  {
                    // caller wanted e^x
                    op->parms.parms[0]=op->parms.parms[1];
                  }
                  else
                  {
                    // it would be nice to replace 10^x with exp(log(10)*x) or 2^x with exp(log(2),x), but 
                    // doing so breaks rounding. we could maybe only allow 10^x, which is used for dB conversion,
                    // but for now we should just force the programmer do it exp(log(10)*x) themselves.
                    break;
                    /* 
                    parm0->opcodeType = OPCODETYPE_FUNC2;
                    parm0->fntype = FN_MULTIPLY;
                    parm0->parms.parms[0] = nseel_createCompiledValue(ctx,dvalue);
                    parm0->parms.parms[1] = op->parms.parms[1];
                    */
                  }
                  op->opcodeType = OPCODETYPE_FUNC1;
                  op->fntype = FUNCTYPE_FUNCTIONTYPEREC;
                  op->fn = &expcpy;
                  goto start_over;
                }
              }
            break;
            case FN_MOD:
              if (dv1)
              {
                const int a = (int) dvalue;
                if (!a) 
                {
                  RESTART_DIRECTVALUE(0.0);
                }
              }
            break;
            case FN_DIVIDE:
              if (dv1)
              {
                if (dvalue == 1.0)  // remove divide by 1.0  (using non-1.0 value as replacement)
                {
                  memcpy(op,op->parms.parms[!!dv0],sizeof(*op));
                  goto start_over;
                }
                else
                {
                  // change to a multiply
                  if (dvalue == 0.0)
                  {
                    op->fntype = FN_MULTIPLY;
                    goto start_over;
                  }
                  else
                  {
					  op->fntype = FN_MULTIPLY;
					  op->parms.parms[1]->parms.dv.directValue = 1.0 / dvalue;
					  op->parms.parms[1]->parms.dv.valuePtr = NULL;
					  goto start_over;
                  }
                }
              }
              else if (dvalue == 0.0)
              {
                if (!retv_parm[!!dv0])
                {
                  // if 0/x set to always 0.
                  // this is 0.0 / (oldexpression that can be eliminated)
                  memcpy(op,op->parms.parms[!dv0],sizeof(*op)); // set to 0 if other action wouldn't do anything
                }
                else
                {
                  opcodeRec *tmp;
                  // this is 0.0 / oldexpressionthatmustbeprocessed
                  op->fntype = FN_JOIN_STATEMENTS;
                  tmp = op->parms.parms[1];
                  op->parms.parms[1] = op->parms.parms[0];
                  op->parms.parms[0] = tmp;
                  // set to (oldexpression;0)
                }
                goto start_over;
              }
            break;
            case FN_EQ:
              if (dvalue == 0.0)
              {
                // convert x == 0.0 to !x
                op->opcodeType=OPCODETYPE_FUNC1;
                op->fntype = FN_NOT;
                if (dv0) op->parms.parms[0]=op->parms.parms[1];
                goto start_over;
              }
            break;
            case FN_NE:
              if (dvalue == 0.0)
              {
                // convert x != 0.0 to !!
                op->opcodeType=OPCODETYPE_FUNC1;
                op->fntype = FN_NOTNOT;
                if (dv0) op->parms.parms[0]=op->parms.parms[1];
                goto start_over;
              }
            break;
            case FN_LOGICAL_AND:
              if (dv0)
              {
                // dvalue && expr
                if (fabs(dvalue) < NSEEL_CLOSEFACTOR)
                {
                  // 0 && expr, replace with 0
                  RESTART_DIRECTVALUE(0.0);
                }
                else
                {
                  // 1 && expr, replace with 0 != expr
                  op->fntype = FN_NE;
                  op->parms.parms[0]->parms.dv.valuePtr=NULL;
                  op->parms.parms[0]->parms.dv.directValue = 0.0;
                }
              }
              else
              {
                // expr && dvalue
                if (fabs(dvalue) < NSEEL_CLOSEFACTOR)
                {
                  // expr && 0
                  if (!retv_parm[0]) 
                  {
                    // expr has no consequence, drop it
                    RESTART_DIRECTVALUE(0.0);
                  }
                  else
                  {
                    // replace with (expr; 0)
                    op->fntype = FN_JOIN_STATEMENTS;
                    op->parms.parms[1]->parms.dv.valuePtr=NULL;
                    op->parms.parms[1]->parms.dv.directValue = 0.0;
                  }
                }
                else
                {
                  // expr && 1, replace with expr != 0
                  op->fntype = FN_NE;
                  op->parms.parms[1]->parms.dv.valuePtr=NULL;
                  op->parms.parms[1]->parms.dv.directValue = 0.0;
                }
              }
            goto start_over;
            case FN_LOGICAL_OR:
              if (dv0)
              {
                // dvalue || expr
                if (fabs(dvalue) >= NSEEL_CLOSEFACTOR)
                {
                  // 1 || expr, replace with 1
                  RESTART_DIRECTVALUE(1.0);
                }
                else
                {
                  // 0 || expr, replace with 0 != expr
                  op->fntype = FN_NE;
                  op->parms.parms[0]->parms.dv.valuePtr=NULL;
                  op->parms.parms[0]->parms.dv.directValue = 0.0;
                }
              }
              else
              {
                // expr || dvalue
                if (fabs(dvalue) >= NSEEL_CLOSEFACTOR)
                {
                  // expr || 1
                  if (!retv_parm[0]) 
                  {
                    // expr has no consequence, drop it and return 1
                    RESTART_DIRECTVALUE(1.0);
                  }
                  else
                  {
                    // replace with (expr; 1)
                    op->fntype = FN_JOIN_STATEMENTS;
                    op->parms.parms[1]->parms.dv.valuePtr=NULL;
                    op->parms.parms[1]->parms.dv.directValue = 1.0;
                  }
                }
                else
                {
                  // expr || 0, replace with expr != 0
                  op->fntype = FN_NE;
                  op->parms.parms[1]->parms.dv.valuePtr=NULL;
                  op->parms.parms[1]->parms.dv.directValue = 0.0;
                }
              }
            goto start_over;
          }
        } // dv0 || dv1
        // general optimization of two parameters
        switch (op->fntype)
        {
          case FN_MULTIPLY:
          {
            opcodeRec *first_parm = op->parms.parms[0],*second_parm = op->parms.parms[1];
            if (second_parm->opcodeType == first_parm->opcodeType) 
            {
              switch(first_parm->opcodeType)
              {
                case OPCODETYPE_VALUE_FROM_NAMESPACENAME:
                  if (first_parm->namespaceidx != second_parm->namespaceidx) break;
                  // fall through
                case OPCODETYPE_VARPTR:
                  if (first_parm->relname && second_parm->relname && !strcmp(second_parm->relname,first_parm->relname)) second_parm=NULL;
                break;
                case OPCODETYPE_VARPTRPTR:
                  if (first_parm->parms.dv.valuePtr && first_parm->parms.dv.valuePtr==second_parm->parms.dv.valuePtr) second_parm=NULL;
                break;
              }
              if (!second_parm) // switch from x*x to sqr(x)
              {
                static functionType sqrcpy={ "sqr",    nseel_asm_sqr,nseel_asm_sqr_end,   1|NSEEL_NPARAMS_FLAG_CONST|BIF_RETURNSONSTACK|BIF_LASTPARMONSTACK|BIF_FPSTACKUSE(1) };
                op->opcodeType = OPCODETYPE_FUNC1;
                op->fntype = FUNCTYPE_FUNCTIONTYPEREC;
                op->fn = &sqrcpy;
                goto start_over;
              }
            }
          }
          break;
          case FN_POW:
            {
              opcodeRec *first_parm = op->parms.parms[0];
              if (first_parm->opcodeType == op->opcodeType && first_parm->fntype == FN_POW)
              {
                // since first_parm is a pow too, we can multiply the exponents.
                // set our base to be the base of the inner pow
                op->parms.parms[0] = first_parm->parms.parms[0];
                // make the old extra pow be a multiply of the exponents
                first_parm->fntype = FN_MULTIPLY;
                first_parm->parms.parms[0] = op->parms.parms[1];
                // put that as the exponent
                op->parms.parms[1] = first_parm;
                goto start_over;
              }
            }
          break;
          case FN_LOGICAL_AND:
          case FN_LOGICAL_OR:
            if (op->parms.parms[0]->fntype == FN_NOTNOT)
            {
              // remove notnot, unnecessary for input to &&/|| operators
              op->parms.parms[0] = op->parms.parms[0]->parms.parms[0];
              goto start_over;
            }
            if (op->parms.parms[1]->fntype == FN_NOTNOT)
            {
              // remove notnot, unnecessary for input to &&/|| operators
              op->parms.parms[1] = op->parms.parms[1]->parms.parms[0];
              goto start_over;
            }        
          break;
        }
      }
      else if (op->opcodeType==OPCODETYPE_FUNC3)  // within FUNCTYPE_SIMPLE
      {
        if (op->fntype == FN_IF_ELSE)
        {
          if (op->parms.parms[0]->opcodeType == OPCODETYPE_DIRECTVALUE)
          {
            int s = fabs(op->parms.parms[0]->parms.dv.directValue) >= NSEEL_CLOSEFACTOR;
            memcpy(op,op->parms.parms[s ? 1 : 2],sizeof(opcodeRec));
            goto start_over;
          }
          if (op->parms.parms[0]->opcodeType == OPCODETYPE_FUNC1)
          {
            if (op->parms.parms[0]->fntype == FN_NOTNOT)
            {
              // remove notnot, unnecessary for input to ? operator
              op->parms.parms[0] = op->parms.parms[0]->parms.parms[0];
              goto start_over;
            }
          }
        }
      }
      if (op->fntype >= FN_NONCONST_BEGIN && op->fntype < FUNCTYPE_SIMPLEMAX) retv|=1;
      // FUNCTYPE_SIMPLE
    }   
    else if (op->fntype == FUNCTYPE_FUNCTIONTYPEREC && op->fn)
    {
      /*
      probably worth doing reduction on:
      _divop (constant change to multiply)
      _and
      _or
      abs
      maybe:
      min
      max
      also, optimize should (recursively or maybe iteratively?) search transitive functions (mul/div) for more constant reduction possibilities
      */
      functionType  *pfn = (functionType *)op->fn;
      if (!(pfn->nParams&NSEEL_NPARAMS_FLAG_CONST)) retv|=1;
      if (op->opcodeType==OPCODETYPE_FUNC1) // within FUNCTYPE_FUNCTIONTYPEREC
      {
        if (op->parms.parms[0]->opcodeType == OPCODETYPE_DIRECTVALUE)
        {
          int suc=1;
          EEL_F v = op->parms.parms[0]->parms.dv.directValue;
  #define DOF(x) if (!strcmp(pfn->name,#x)) v = x(v); else
  #define DOF2(x,y) if (!strcmp(pfn->name,#x)) v = x(y); else
          DOF(sin)
          DOF(cos)
          DOF(tan)
          DOF(asin)
          DOF(acos)
          DOF(atan)
          DOF2(sqrt, fabs(v))
          DOF(exp)
          DOF(log)
          DOF(log10)
          /* else */ suc=0;
  #undef DOF
  #undef DOF2
          if (suc)
          {
            RESTART_DIRECTVALUE(v);
          }
        }
      }
      else if (op->opcodeType==OPCODETYPE_FUNC2)  // within FUNCTYPE_FUNCTIONTYPEREC
      {
        const int dv0=op->parms.parms[0]->opcodeType == OPCODETYPE_DIRECTVALUE;
        const int dv1=op->parms.parms[1]->opcodeType == OPCODETYPE_DIRECTVALUE;
        if (dv0 && dv1)
        {
          if (!strcmp(pfn->name,"atan2")) 
          {
            RESTART_DIRECTVALUE(atan2(op->parms.parms[0]->parms.dv.directValue, op->parms.parms[1]->parms.dv.directValue));
          }
        }
      }
      // FUNCTYPE_FUNCTIONTYPEREC
    }
    else
    {
      // unknown or eel func, assume non-const
      retv |= 1;
    }
  }
  // if we need results, or our function has effects itself, then finish
  if (retv || needsResult)
  {
    return retv || joined_retv || retv_parm[0] || retv_parm[1] || retv_parm[2];
  }
  // we don't need results here, and our function is const, which means we can remove it
  {
    int cnt=0, idx1=0, idx2=0, x;
    for (x=0;x<3;x++) if (retv_parm[x]) { if (!cnt++) idx1=x; else idx2=x; }
    if (!cnt) // none of the parameters do anything, remove this opcode
    {
      if (lastJoinOp)
      {
        // replace previous join with its first linked opcode, removing this opcode completely
        memcpy(lastJoinOp,lastJoinOp->parms.parms[0],sizeof(*lastJoinOp));
      }
      else if (op->opcodeType != OPCODETYPE_DIRECTVALUE)
      {
        // allow caller to easily detect this as trivial and remove it
        op->opcodeType = OPCODETYPE_DIRECTVALUE;
        op->parms.dv.valuePtr=NULL;
        op->parms.dv.directValue=0.0;
      }
      // return joined_retv below
    }
    else
    {
      // if parameters are non-const, and we're a conditional, preserve our function
      if (FNPTR_HAS_CONDITIONAL_EXEC(op)) return 1;
      // otherwise, condense into either the non-const statement, or a join
      if (cnt==1)
      {
        memcpy(op,op->parms.parms[idx1],sizeof(*op));
      }
      else if (cnt == 2)
      {
        op->opcodeType = OPCODETYPE_FUNC2;
        op->fntype = FN_JOIN_STATEMENTS;
        op->fn = op;
        op->parms.parms[0] = op->parms.parms[idx1];
        op->parms.parms[1] = op->parms.parms[idx2];
        op->parms.parms[2] = NULL;
      }
      else
      {
        // todo need to create a new opcodeRec here, for now just leave as is 
        // (non-conditional const 3 parameter functions are rare anyway)
      }
      return 1;
    }
  }
  return joined_retv;
}
static int generateValueToReg(compileContext *ctx, opcodeRec *op, unsigned char *bufOut, int whichReg, const namespaceInformation *functionPrefix, int allowCache)
{
  EEL_F *b=NULL;
  if (op->opcodeType==OPCODETYPE_VALUE_FROM_NAMESPACENAME)
  {
    char nm[NSEEL_MAX_VARIABLE_NAMELEN+1];
    const char *p = op->relname;
    combineNamespaceFields(nm,functionPrefix,p+(*p == '#'),op->namespaceidx);
    if (!nm[0]) return -1;
	b = nseel_int_register_var(ctx, nm, 0, NULL);
	if (!b) RET_MINUS1_FAIL("error registering var")
  }
  else
  {
    if (op->opcodeType != OPCODETYPE_DIRECTVALUE) allowCache=0;
    b=op->parms.dv.valuePtr;
    if (!b && op->opcodeType == OPCODETYPE_VARPTR && op->relname && op->relname[0]) 
    {
      op->parms.dv.valuePtr = b = nseel_int_register_var(ctx,op->relname,0,NULL);
    }
    if (b && op->opcodeType == OPCODETYPE_VARPTRPTR) b = *(EEL_F **)b;
    if (!b && allowCache)
    {
      int n=50; // only scan last X items
      opcodeRec *r = ctx->directValueCache;
      while (r && n--)
      {
        if (r->parms.dv.directValue == op->parms.dv.directValue && (b=r->parms.dv.valuePtr)) break;
        r=(opcodeRec*)r->fn;
      }
    }
    if (!b)
    {
      ctx->l_stats[3]++;
      if (ctx->isGeneratingCommonFunction)
        b = newCtxDataBlock(sizeof(EEL_F),sizeof(EEL_F));
      else
        b = newDataBlock(sizeof(EEL_F),sizeof(EEL_F));
      if (!b) RET_MINUS1_FAIL("error allocating data block")
      if (op->opcodeType != OPCODETYPE_VARPTRPTR) op->parms.dv.valuePtr = b;
      #if EEL_F_SIZE == 8
        *b = op->parms.dv.directValue;
      #else
        *b = op->parms.dv.directValue;
      #endif
      if (allowCache)
      {
        op->fn = ctx->directValueCache;
        ctx->directValueCache = op;
      }
    }
  }
  GLUE_MOV_PX_DIRECTVALUE_GEN(bufOut,(INT_PTR)b,whichReg);
  return GLUE_MOV_PX_DIRECTVALUE_SIZE;
}
unsigned char *compileCodeBlockWithRet(compileContext *ctx, opcodeRec *rec, int *computTableSize, const namespaceInformation *namespacePathToThis, 
                                       int supportedReturnValues, int *rvType, int *fpStackUsage, int *canHaveDenormalOutput)
{
  unsigned char *p, *newblock2;
  // generate code call
  int funcsz=compileOpcodes(ctx,rec,NULL,1024*1024*128,NULL,namespacePathToThis,supportedReturnValues, rvType,fpStackUsage, NULL);
  if (funcsz<0) return NULL;
  p = newblock2 = newCodeBlock(funcsz+ sizeof(GLUE_RET),32);
  if (!newblock2) return NULL;
  *fpStackUsage=0;
  funcsz=compileOpcodes(ctx,rec,p, funcsz, computTableSize,namespacePathToThis,supportedReturnValues, rvType,fpStackUsage, canHaveDenormalOutput);         
  if (funcsz<0) return NULL;
  p+=funcsz;
  memcpy(p,&GLUE_RET,sizeof(GLUE_RET)); p+=sizeof(GLUE_RET);
#ifdef __arm__
  __clear_cache(newblock2,p);
#endif
  ctx->l_stats[2]+=funcsz+2;
  return newblock2;
}      
static int compileNativeFunctionCall(compileContext *ctx, opcodeRec *op, unsigned char *bufOut, int bufOut_len, int *computTableSize, const namespaceInformation *namespacePathToThis, 
                                     int *rvMode, int *fpStackUsage, int preferredReturnValues, int *canHaveDenormalOutput)
{
  // builtin function generation
  int func_size=0;
  int cfunc_abiinfo=0;
  int local_fpstack_use=0; // how many items we have pushed onto the fp stack
  int parm_size=0;
  int restore_stack_amt=0;
  void *func_e=NULL;
  NSEEL_PPPROC preProc=0;
  void **repl=NULL;
  int n_params= 1 + op->opcodeType - OPCODETYPE_FUNC1;
  const int parm0_dv = op->parms.parms[0]->opcodeType == OPCODETYPE_DIRECTVALUE;
  const int parm1_dv = n_params > 1 && op->parms.parms[1]->opcodeType == OPCODETYPE_DIRECTVALUE;
  void *func = nseel_getBuiltinFunctionAddress(ctx, op->fntype, op->fn, &preProc,&repl,&func_e,&cfunc_abiinfo,preferredReturnValues, 
       parm0_dv ? &op->parms.parms[0]->parms.dv.directValue : NULL,
       parm1_dv ? &op->parms.parms[1]->parms.dv.directValue : NULL
       );
  if (!func) RET_MINUS1_FAIL("error getting funcaddr")
  *fpStackUsage=BIF_GETFPSTACKUSE(cfunc_abiinfo);
  *rvMode = RETURNVALUE_NORMAL;
  if (cfunc_abiinfo & BIF_TAKES_VARPARM)
  {
#if defined(__arm__) || defined(__ppc__) || (defined (_M_ARM) && _M_ARM == 7)
    const int max_params=4096; // 32kb max offset addressing for stack, so 4096*4 = 16384, should be safe
#else
    const int max_params=32768; // sanity check, the stack is free to grow on x86/x86-64
#endif
    int x;
    // this mode is less efficient in that it creates a list of pointers on the stack to pass to the function
    // but it is more flexible and works for >3 parameters.
    if (op->opcodeType == OPCODETYPE_FUNCX)
    {
      n_params=0;
      for (x=0;x<3;x++)
      {
        opcodeRec *prni=op->parms.parms[x];
        while (prni)
        {
          const int isMP = prni->opcodeType == OPCODETYPE_MOREPARAMS;
          n_params++;
          if (!isMP||n_params>=max_params) break;
          prni = prni->parms.parms[1];
        }
      }
    }
    restore_stack_amt = (sizeof(void *) * n_params + 15)&~15;
    if (restore_stack_amt)
    {
      int offs = restore_stack_amt;
      while (offs > 0)
      {
        int amt = offs;
        if (amt > 4096) amt=4096;
        if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_MOVE_STACK_SIZE)) RET_MINUS1_FAIL("insufficient size for varparm")
        if (bufOut) GLUE_MOVE_STACK(bufOut+parm_size, - amt);
        parm_size += GLUE_MOVE_STACK_SIZE;
        offs -= amt;
        if (offs>0) // make sure this page is in memory
        {
          if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_STORE_P1_TO_STACK_AT_OFFS_SIZE))
            RET_MINUS1_FAIL("insufficient size for varparm stackchk")
          if (bufOut) GLUE_STORE_P1_TO_STACK_AT_OFFS(bufOut+parm_size,0);
          parm_size += GLUE_STORE_P1_TO_STACK_AT_OFFS_SIZE;
        }
      }
    }
    if (op->opcodeType == OPCODETYPE_FUNCX)
    {     
      n_params=0;
      for (x=0;x<3;x++)
      {
        opcodeRec *prni=op->parms.parms[x];
        while (prni)
        {
          const int isMP = prni->opcodeType == OPCODETYPE_MOREPARAMS;
          opcodeRec *r = isMP ? prni->parms.parms[0] : prni;
          if (r)
          {
            int canHaveDenorm=0;
            int rvt=RETURNVALUE_NORMAL;
            int subfpstackuse=0;
            int lsz = compileOpcodes(ctx,r,bufOut ? bufOut + parm_size : NULL,bufOut_len - parm_size, computTableSize, namespacePathToThis, rvt,&rvt, &subfpstackuse, &canHaveDenorm);
            if (canHaveDenorm && canHaveDenormalOutput) *canHaveDenormalOutput = 1;
            if (lsz<0) RET_MINUS1_FAIL("call coc for varparmX failed")
            if (rvt != RETURNVALUE_NORMAL) RET_MINUS1_FAIL("call coc for varparmX gave bad type back");
            parm_size += lsz;            
            if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_STORE_P1_TO_STACK_AT_OFFS_SIZE)) RET_MINUS1_FAIL("call coc for varparmX size");
            if (bufOut) GLUE_STORE_P1_TO_STACK_AT_OFFS(bufOut + parm_size, n_params*sizeof(void *));
            parm_size+=GLUE_STORE_P1_TO_STACK_AT_OFFS_SIZE;
            if (subfpstackuse+local_fpstack_use > *fpStackUsage) *fpStackUsage = subfpstackuse+local_fpstack_use;
          }
          else RET_MINUS1_FAIL("zero parameter varparmX")
          n_params++;
          if (!isMP||n_params>=max_params) break;
          prni = prni->parms.parms[1];
        }
      }
    }
    else for (x=0;x<n_params;x++)
    {
      opcodeRec *r = op->parms.parms[x];
      if (r)
      {
        int canHaveDenorm=0;
        int subfpstackuse=0;
        int rvt=RETURNVALUE_NORMAL;
        int lsz = compileOpcodes(ctx,r,bufOut ? bufOut + parm_size : NULL,bufOut_len - parm_size, computTableSize, namespacePathToThis, rvt,&rvt, &subfpstackuse, &canHaveDenorm);
        if (canHaveDenorm && canHaveDenormalOutput) *canHaveDenormalOutput = 1;
        if (lsz<0) RET_MINUS1_FAIL("call coc for varparm123 failed")
        if (rvt != RETURNVALUE_NORMAL) RET_MINUS1_FAIL("call coc for varparm123 gave bad type back");
        parm_size += lsz;
        if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_STORE_P1_TO_STACK_AT_OFFS_SIZE)) RET_MINUS1_FAIL("call coc for varparm123 size");
        if (bufOut) GLUE_STORE_P1_TO_STACK_AT_OFFS(bufOut + parm_size, x*sizeof(void *));
        parm_size+=GLUE_STORE_P1_TO_STACK_AT_OFFS_SIZE;
        if (subfpstackuse+local_fpstack_use > *fpStackUsage) *fpStackUsage = subfpstackuse+local_fpstack_use;
      }
      else RET_MINUS1_FAIL("zero parameter for varparm123");
    }
    if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_MOV_PX_DIRECTVALUE_SIZE+GLUE_MOVE_PX_STACKPTR_SIZE)) RET_MINUS1_FAIL("insufficient size for varparm p1")
    if (bufOut) GLUE_MOV_PX_DIRECTVALUE_GEN(bufOut+parm_size, (INT_PTR)n_params,1);
    parm_size+=GLUE_MOV_PX_DIRECTVALUE_SIZE;
    if (bufOut) GLUE_MOVE_PX_STACKPTR_GEN(bufOut+parm_size, 0);
    parm_size+=GLUE_MOVE_PX_STACKPTR_SIZE;
  }
  else // not varparm
  {
    int pn;
  #ifdef GLUE_HAS_FXCH
    int need_fxch=0;
  #endif
    int last_nt_parm=-1, last_nt_parm_type;
    if (op->opcodeType == OPCODETYPE_FUNCX)
    {
      // this is not yet supported (calling conventions will need to be sorted, among other things)
      RET_MINUS1_FAIL("funcx for native functions requires BIF_TAKES_VARPARM or BIF_TAKES_VARPARM_EX")
    }
    if (parm0_dv) 
    {
      if (func == nseel_asm_stack_pop)
      {
        func = GLUE_realAddress(nseel_asm_stack_pop_fast,nseel_asm_stack_pop_fast_end,&func_size);
        if (!func || bufOut_len < func_size) RET_MINUS1_FAIL(func?"failed on popfast size":"failed on popfast addr")
        if (bufOut) 
        {
          memcpy(bufOut,func,func_size);
          NSEEL_PProc_Stack(bufOut,func_size,ctx);
        }
        return func_size;            
      }
      else if (func == nseel_asm_stack_peek)
      {
        int f = (int) op->parms.parms[0]->parms.dv.directValue;
        if (!f)
        {
          func = GLUE_realAddress(nseel_asm_stack_peek_top,nseel_asm_stack_peek_top_end,&func_size);
          if (!func || bufOut_len < func_size) RET_MINUS1_FAIL(func?"failed on peek size":"failed on peek addr")
          if (bufOut) 
          {
            memcpy(bufOut,func,func_size);
            NSEEL_PProc_Stack_PeekTop(bufOut,func_size,ctx);
          }
          return func_size;
        }
        else
        {
          func = GLUE_realAddress(nseel_asm_stack_peek_int,nseel_asm_stack_peek_int_end,&func_size);
          if (!func || bufOut_len < func_size) RET_MINUS1_FAIL(func?"failed on peekint size":"failed on peekint addr")
          if (bufOut)
          {
            memcpy(bufOut,func,func_size);
            NSEEL_PProc_Stack_PeekInt(bufOut,func_size,ctx,f*sizeof(EEL_F));
          }
          return func_size;
        }
      }
    }
    // end of built-in function specific special casing
    // first pass, calculate any non-trivial parameters
    for (pn=0; pn < n_params; pn++)
    { 
      if (!OPCODE_IS_TRIVIAL(op->parms.parms[pn]))
      {
        int canHaveDenorm=0;
        int subfpstackuse=0;
        int lsz=0; 
        int rvt=RETURNVALUE_NORMAL;
        int may_need_fppush=-1;
        if (last_nt_parm>=0)
        {
          if (last_nt_parm_type==RETURNVALUE_FPSTACK)
          {          
            may_need_fppush= parm_size;
          }
          else
          {
            // push last result
            if (bufOut_len < parm_size + (int)sizeof(GLUE_PUSH_P1)) RET_MINUS1_FAIL("failed on size, pushp1")
            if (bufOut) memcpy(bufOut + parm_size, &GLUE_PUSH_P1, sizeof(GLUE_PUSH_P1));
            parm_size += sizeof(GLUE_PUSH_P1);
          }
        }         
        if (func == nseel_asm_bnot) rvt=RETURNVALUE_BOOL_REVERSED|RETURNVALUE_BOOL;
        else if (pn == n_params - 1)
        {
          if (cfunc_abiinfo&BIF_LASTPARMONSTACK) rvt=RETURNVALUE_FPSTACK;
          else if (cfunc_abiinfo&BIF_LASTPARM_ASBOOL) rvt=RETURNVALUE_BOOL;
          else if (func == nseel_asm_assign) rvt=RETURNVALUE_FPSTACK|RETURNVALUE_NORMAL;
        }
        else if (pn == n_params -2 && (cfunc_abiinfo&BIF_SECONDLASTPARMST))
        {
          rvt=RETURNVALUE_FPSTACK;
        }
        lsz = compileOpcodes(ctx,op->parms.parms[pn],bufOut ? bufOut + parm_size : NULL,bufOut_len - parm_size, computTableSize, namespacePathToThis, rvt,&rvt, &subfpstackuse, &canHaveDenorm);
        if (lsz<0) RET_MINUS1_FAIL("call coc failed")
        if (func == nseel_asm_bnot && rvt==RETURNVALUE_BOOL_REVERSED)
        {
          // remove bnot, compileOpcodes() used fptobool_rev
			func = nseel_asm_bnotnot;
			func_e = nseel_asm_bnotnot_end;
          rvt = RETURNVALUE_BOOL;
        }
        if (canHaveDenorm && canHaveDenormalOutput) *canHaveDenormalOutput = 1;
        parm_size += lsz;            
        if (may_need_fppush>=0)
        {
          if (local_fpstack_use+subfpstackuse >= (GLUE_MAX_FPSTACK_SIZE-1) || (0&OPTFLAG_NO_FPSTACK))
          {
            if (bufOut_len < parm_size + (int)sizeof(GLUE_POP_FPSTACK_TOSTACK)) 
              RET_MINUS1_FAIL("failed on size, popfpstacktostack")
            if (bufOut) 
            {
              memmove(bufOut + may_need_fppush + sizeof(GLUE_POP_FPSTACK_TOSTACK), bufOut + may_need_fppush, parm_size - may_need_fppush);
              memcpy(bufOut + may_need_fppush, &GLUE_POP_FPSTACK_TOSTACK, sizeof(GLUE_POP_FPSTACK_TOSTACK));
            }
            parm_size += sizeof(GLUE_POP_FPSTACK_TOSTACK);
          }
          else
          {
            local_fpstack_use++;
          }
        }
        if (subfpstackuse+local_fpstack_use > *fpStackUsage)
			*fpStackUsage = subfpstackuse+local_fpstack_use;
        last_nt_parm = pn;
        last_nt_parm_type = rvt;
        if (pn == n_params - 1 && func == nseel_asm_assign)
        {
          if (!canHaveDenorm)
          {
            if (rvt == RETURNVALUE_FPSTACK)
            {
              cfunc_abiinfo |= BIF_LASTPARMONSTACK;
              func = nseel_asm_assign_fast_fromfp;
              func_e = nseel_asm_assign_fast_fromfp_end;
            }
            else
            {
              func = nseel_asm_assign_fast;
              func_e = nseel_asm_assign_fast_end;
            }
          }
          else
          {
            if (rvt == RETURNVALUE_FPSTACK)
            {
              cfunc_abiinfo |= BIF_LASTPARMONSTACK;
              func = nseel_asm_assign_fromfp;
              func_e = nseel_asm_assign_fromfp_end;
            }
          }
        }
      }
    }
    pn = last_nt_parm;
    if (pn >= 0) // if the last thing executed doesn't go to the last parameter, move it there
    {
      if ((cfunc_abiinfo&BIF_SECONDLASTPARMST) && pn == n_params-2)
      {
        // do nothing, things are in the right place
      }
      else if (pn != n_params-1)
      {
        // generate mov p1->pX
        if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_SET_PX_FROM_P1_SIZE)) RET_MINUS1_FAIL("size, pxfromp1")
        if (bufOut) GLUE_SET_PX_FROM_P1(bufOut + parm_size,n_params - 1 - pn);
        parm_size += GLUE_SET_PX_FROM_P1_SIZE;
      }
    }
    // pop any pushed parameters
    while (--pn >= 0)
    { 
      if (!OPCODE_IS_TRIVIAL(op->parms.parms[pn]))
      {
        if ((cfunc_abiinfo&BIF_SECONDLASTPARMST) && pn == n_params-2)
        {
          if (!local_fpstack_use)
          {
            if (bufOut_len < parm_size + (int)sizeof(GLUE_POP_STACK_TO_FPSTACK)) RET_MINUS1_FAIL("size, popstacktofpstack 2")
            if (bufOut) memcpy(bufOut+parm_size,GLUE_POP_STACK_TO_FPSTACK,sizeof(GLUE_POP_STACK_TO_FPSTACK));
            parm_size += sizeof(GLUE_POP_STACK_TO_FPSTACK);
            #ifdef GLUE_HAS_FXCH
              need_fxch = 1;
            #endif
          }
          else
          {
            local_fpstack_use--;
          }
        }
        else
        {
          if ((size_t)bufOut_len < (size_t)(parm_size) + GLUE_POP_PX_SIZE) RET_MINUS1_FAIL("size, poppx")
          if (bufOut) GLUE_POP_PX(bufOut + parm_size,n_params - 1 - pn);
          parm_size += GLUE_POP_PX_SIZE;
        }
      }
    }
    // finally, set trivial pointers
    for (pn=0; pn < n_params; pn++)
    { 
      if (OPCODE_IS_TRIVIAL(op->parms.parms[pn]))
      {
        if (pn == n_params-2 && (cfunc_abiinfo&(BIF_SECONDLASTPARMST)))  // second to last parameter
        {
          int a = compileOpcodes(ctx,op->parms.parms[pn],bufOut ? bufOut+parm_size : NULL,bufOut_len - parm_size,computTableSize,namespacePathToThis,
                                  RETURNVALUE_FPSTACK,NULL,NULL,canHaveDenormalOutput);
          if (a<0) RET_MINUS1_FAIL("coc call here 2")
          parm_size+=a;
          #ifdef GLUE_HAS_FXCH
            need_fxch = 1;
          #endif
        }
        else if (pn == n_params-1)  // last parameter, but we should call compileOpcodes to get it in the right format (compileOpcodes can optimize that process if it needs to)
        {
          int rvt=0, a;
          int wantFpStack = func == nseel_asm_assign;
  #ifdef GLUE_PREFER_NONFP_DV_ASSIGNS // x86-64, and maybe others, prefer to avoid the fp stack for a simple copy
          if (wantFpStack &&
              (op->parms.parms[pn]->opcodeType != OPCODETYPE_DIRECTVALUE ||
              (op->parms.parms[pn]->parms.dv.directValue != 1.0 && op->parms.parms[pn]->parms.dv.directValue != 0.0)))
          {
            wantFpStack=0;
          }
  #endif
          a = compileOpcodes(ctx,op->parms.parms[pn],bufOut ? bufOut+parm_size : NULL,bufOut_len - parm_size,computTableSize,namespacePathToThis,
            func == nseel_asm_bnot ? (RETURNVALUE_BOOL_REVERSED|RETURNVALUE_BOOL) :
              (cfunc_abiinfo & BIF_LASTPARMONSTACK) ? RETURNVALUE_FPSTACK : 
              (cfunc_abiinfo & BIF_LASTPARM_ASBOOL) ? RETURNVALUE_BOOL : 
              wantFpStack ? (RETURNVALUE_FPSTACK|RETURNVALUE_NORMAL) : 
              RETURNVALUE_NORMAL,       
            &rvt, NULL,canHaveDenormalOutput);
          if (a<0) RET_MINUS1_FAIL("coc call here 3")
          if (func == nseel_asm_bnot && rvt == RETURNVALUE_BOOL_REVERSED)
          {
            // remove bnot, compileOpcodes() used fptobool_rev
			  func = nseel_asm_bnotnot;
			  func_e = nseel_asm_bnotnot_end;
            rvt = RETURNVALUE_BOOL;
          }
          parm_size+=a;
          #ifdef GLUE_HAS_FXCH
            need_fxch = 0;
          #endif
          if (func == nseel_asm_assign)
          {
            if (rvt == RETURNVALUE_FPSTACK)
            {
				func = nseel_asm_assign_fast_fromfp;
				func_e = nseel_asm_assign_fast_fromfp_end;
            }
            else
            {
               // assigning a value (from a variable or other non-computer), can use a fast assign (no denormal/result checking)
              func = nseel_asm_assign_fast;
              func_e = nseel_asm_assign_fast_end;
            }
          }
        }
        else
        {
          if ((size_t)bufOut_len < ((size_t)(parm_size) + GLUE_MOV_PX_DIRECTVALUE_SIZE)) RET_MINUS1_FAIL("size, pxdvsz")
          if (bufOut) 
          {
            if (generateValueToReg(ctx,op->parms.parms[pn],bufOut + parm_size,n_params - 1 - pn,namespacePathToThis, 0/*nocaching, function gets pointer*/)<0) RET_MINUS1_FAIL("gvtr")
          }
          parm_size += GLUE_MOV_PX_DIRECTVALUE_SIZE;
        }
      }
    }
  #ifdef GLUE_HAS_FXCH
    if ((cfunc_abiinfo&(BIF_SECONDLASTPARMST)) && !(cfunc_abiinfo&(BIF_LAZYPARMORDERING))&&
        ((!!need_fxch)^!!(cfunc_abiinfo&BIF_REVERSEFPORDER)) 
        )
    {
      // emit fxch
      if ((size_t)bufOut_len < sizeof(GLUE_FXCH)) RET_MINUS1_FAIL("len,fxch")
      if (bufOut) 
      { 
        memcpy(bufOut+parm_size,GLUE_FXCH,sizeof(GLUE_FXCH));
      }
      parm_size+=sizeof(GLUE_FXCH);
    }
  #endif
    if (!*canHaveDenormalOutput)
    {
      // if add_op or sub_op, and non-denormal input, safe to omit denormal checks
      if (func == (void*)nseel_asm_add_op)
      {
        func = nseel_asm_add_op_fast;
        func_e = nseel_asm_add_op_fast_end;
      }
      else if (func == (void*)nseel_asm_sub_op)
      {
        func = nseel_asm_sub_op_fast;
        func_e = nseel_asm_sub_op_fast_end;
      }
      // or if mul/div by a fixed value of >= or <= 1.0
      else if (func == (void *)nseel_asm_mul_op && parm1_dv && fabs(op->parms.parms[1]->parms.dv.directValue) >= 1.0)
      {
        func = nseel_asm_mul_op_fast;
        func_e = nseel_asm_mul_op_fast_end;
      }
      else if (func == (void *)nseel_asm_div_op && parm1_dv && fabs(op->parms.parms[1]->parms.dv.directValue) <= 1.0)
      {
        func = nseel_asm_div_op_fast;
        func_e = nseel_asm_div_op_fast_end;
      }
    }
  } // not varparm
  if (cfunc_abiinfo & (BIF_CLEARDENORMAL | BIF_RETURNSBOOL) ) *canHaveDenormalOutput=0;
  else if (!(cfunc_abiinfo & BIF_WONTMAKEDENORMAL)) *canHaveDenormalOutput=1;
  func = GLUE_realAddress(func,func_e,&func_size);
  if (!func) RET_MINUS1_FAIL("failrealladdrfunc")
  if (bufOut_len < parm_size + func_size) RET_MINUS1_FAIL("funcsz")
  if (bufOut)
  {
    unsigned char *p=bufOut + parm_size;
    memcpy(p, func, func_size);
    if (preProc) p=preProc(p,func_size,ctx);
    if (repl)
    {
      if (repl[0]) p=EEL_GLUE_set_immediate(p,(INT_PTR)repl[0]);
      if (repl[1]) p=EEL_GLUE_set_immediate(p,(INT_PTR)repl[1]);
      if (repl[2]) p=EEL_GLUE_set_immediate(p,(INT_PTR)repl[2]);
      if (repl[3]) p=EEL_GLUE_set_immediate(p,(INT_PTR)repl[3]);
    }
  }
  if (restore_stack_amt)
  {
    if ((size_t)bufOut_len < ((size_t)(parm_size + func_size) + GLUE_MOVE_STACK_SIZE)) RET_MINUS1_FAIL("insufficient size for varparm")
    if (bufOut) GLUE_MOVE_STACK(bufOut + parm_size + func_size, restore_stack_amt); 
    parm_size += GLUE_MOVE_STACK_SIZE;
  }
  if (cfunc_abiinfo&BIF_RETURNSONSTACK) *rvMode = RETURNVALUE_FPSTACK;
  else if (cfunc_abiinfo&BIF_RETURNSBOOL) *rvMode=RETURNVALUE_BOOL;
  return parm_size + func_size;
}
static int compileEelFunctionCall(compileContext *ctx, opcodeRec *op, unsigned char *bufOut, int bufOut_len, int *computTableSize, const namespaceInformation *namespacePathToThis, 
                                  int *rvMode, int *fpStackUse, int *canHaveDenormalOutput)
{
  int func_size=0, parm_size=0;
  int pn;
  int last_nt_parm=-1,last_nt_parm_mode=0;
  void *func_e=NULL;
  int n_params;
  opcodeRec *parmptrs[NSEEL_MAX_EELFUNC_PARAMETERS];
  int cfp_numparams=-1;
  int cfp_statesize=0;
  EEL_F **cfp_ptrs=NULL;
  int func_raw=0;
  int do_parms;
  int x;
  void *func;
  for (x=0; x < 3; x ++) parmptrs[x] = op->parms.parms[x];
  if (op->opcodeType == OPCODETYPE_FUNCX)
  {
    n_params=0;
    for (x=0;x<3;x++)
    {
      opcodeRec *prni=op->parms.parms[x];
      while (prni && n_params < NSEEL_MAX_EELFUNC_PARAMETERS)
      {
        const int isMP = prni->opcodeType == OPCODETYPE_MOREPARAMS;
        parmptrs[n_params++] = isMP ? prni->parms.parms[0] : prni;
        if (!isMP) break;
        prni = prni->parms.parms[1];
      }
    }
  }
  else 
  {
    n_params = 1 + op->opcodeType - OPCODETYPE_FUNC1;
  }
  *fpStackUse = 0;
  func = nseel_getEELFunctionAddress(ctx, op,
                                      &cfp_numparams,&cfp_statesize,&cfp_ptrs, 
                                      computTableSize, 
                                      &func_e, &func_raw,                                              
                                      !!bufOut,namespacePathToThis,rvMode,fpStackUse,canHaveDenormalOutput, parmptrs, n_params);
  if (func_raw) func_size = (int) ((char*)func_e  - (char*)func);
  else if (func) func = GLUE_realAddress(func,func_e,&func_size);
  if (!func) RET_MINUS1_FAIL("eelfuncaddr")
  *fpStackUse += 1;
  if (cfp_numparams>0 && n_params != cfp_numparams)
  {
    RET_MINUS1_FAIL("eelfuncnp")
  }
  // user defined function
  do_parms = cfp_numparams>0 && cfp_ptrs && cfp_statesize>0;
  // if function local/parameter state is zero, we need to allocate storage for it
  if (cfp_statesize>0 && cfp_ptrs && !cfp_ptrs[0])
  {
    EEL_F *pstate = newDataBlock(sizeof(EEL_F)*cfp_statesize,8);
    if (!pstate) RET_MINUS1_FAIL("eelfuncdb")
    for (pn=0;pn<cfp_statesize;pn++)
    {
      pstate[pn]=0;
      cfp_ptrs[pn] = pstate + pn;
    }
  }
  // first process parameters that are non-trivial
  for (pn=0; pn < n_params; pn++)
  { 
    int needDenorm=0;
    int lsz,sUse=0;                      
    if (!parmptrs[pn] || OPCODE_IS_TRIVIAL(parmptrs[pn])) continue; // skip and process after
    if (last_nt_parm >= 0 && do_parms)
    {
      if (last_nt_parm_mode == RETURNVALUE_FPSTACK)
      {
        if (bufOut_len < parm_size + (int)sizeof(GLUE_POP_FPSTACK_TOSTACK)) RET_MINUS1_FAIL("eelfunc_size popfpstacktostack")
        if (bufOut) memcpy(bufOut + parm_size,GLUE_POP_FPSTACK_TOSTACK,sizeof(GLUE_POP_FPSTACK_TOSTACK));
        parm_size+=sizeof(GLUE_POP_FPSTACK_TOSTACK);
      }
      else
      {
        if (bufOut_len < parm_size + (int)sizeof(GLUE_PUSH_P1PTR_AS_VALUE)) RET_MINUS1_FAIL("eelfunc_size pushp1ptrasval")
        // push
        if (bufOut) memcpy(bufOut + parm_size,&GLUE_PUSH_P1PTR_AS_VALUE,sizeof(GLUE_PUSH_P1PTR_AS_VALUE));
        parm_size+=sizeof(GLUE_PUSH_P1PTR_AS_VALUE);
      }
    }
    last_nt_parm_mode=0;
    lsz = compileOpcodes(ctx,parmptrs[pn],bufOut ? bufOut + parm_size : NULL,bufOut_len - parm_size, computTableSize, namespacePathToThis,
      do_parms ? (RETURNVALUE_FPSTACK|RETURNVALUE_NORMAL) : RETURNVALUE_IGNORE,&last_nt_parm_mode,&sUse, &needDenorm);
    // todo: if needDenorm, denorm convert when copying parameter
    if (lsz<0) RET_MINUS1_FAIL("eelfunc, coc fail")
    if (last_nt_parm_mode == RETURNVALUE_FPSTACK) sUse++;
    if (sUse > *fpStackUse) *fpStackUse=sUse;
    parm_size += lsz;
    last_nt_parm = pn;
  }
  // pop non-trivial results into place
  if (last_nt_parm >=0 && do_parms)
  {
    while (--pn >= 0)
    { 
      if (!parmptrs[pn] || OPCODE_IS_TRIVIAL(parmptrs[pn])) continue; // skip and process after
      if (pn == last_nt_parm)
      {
        if (last_nt_parm_mode == RETURNVALUE_FPSTACK)
        {
          // pop to memory directly
          const int cpsize = GLUE_POP_FPSTACK_TO_PTR(NULL,NULL);
          if (bufOut_len < parm_size + cpsize) RET_MINUS1_FAIL("eelfunc size popfpstacktoptr")
          if (bufOut) GLUE_POP_FPSTACK_TO_PTR((unsigned char *)bufOut + parm_size,cfp_ptrs[pn]);
          parm_size += cpsize;
        }
        else
        {
          // copy direct p1ptr to mem
          const int cpsize = GLUE_COPY_VALUE_AT_P1_TO_PTR(NULL,NULL);
          if (bufOut_len < parm_size + cpsize) RET_MINUS1_FAIL("eelfunc size copyvalueatp1toptr")
          if (bufOut) GLUE_COPY_VALUE_AT_P1_TO_PTR((unsigned char *)bufOut + parm_size,cfp_ptrs[pn]);
          parm_size += cpsize;
        }
      }
      else
      {
        const int popsize =  GLUE_POP_VALUE_TO_ADDR(NULL,NULL);
        if (bufOut_len < parm_size + popsize) RET_MINUS1_FAIL("eelfunc size pop value to addr")
        if (bufOut) GLUE_POP_VALUE_TO_ADDR((unsigned char *)bufOut + parm_size,cfp_ptrs[pn]);
        parm_size+=popsize;
      }
    }
  }
  // finally, set any trivial parameters
  if (do_parms)
  {
    const int cpsize = GLUE_MOV_PX_DIRECTVALUE_SIZE + GLUE_COPY_VALUE_AT_P1_TO_PTR(NULL,NULL);
    for (pn=0; pn < n_params; pn++)
    { 
      if (!parmptrs[pn] || !OPCODE_IS_TRIVIAL(parmptrs[pn])) continue; // set trivial values, we already set nontrivials
      if (bufOut_len < parm_size + cpsize) RET_MINUS1_FAIL("eelfunc size trivial set")
      if (bufOut) 
      {
        if (generateValueToReg(ctx,parmptrs[pn],bufOut + parm_size,0,namespacePathToThis, 1)<0) RET_MINUS1_FAIL("eelfunc gvr fail")
        GLUE_COPY_VALUE_AT_P1_TO_PTR(bufOut + parm_size + GLUE_MOV_PX_DIRECTVALUE_SIZE,cfp_ptrs[pn]);
      }
      parm_size += cpsize;
    }
  }
  if (bufOut_len < parm_size + func_size) RET_MINUS1_FAIL("eelfunc size combined")
  if (bufOut) memcpy(bufOut + parm_size, func, func_size);
  return parm_size + func_size;
  // end of EEL function generation
}
#define CHECK_SIZE_FORJMP(x,y)
#define RET_MINUS1_FAIL_FALLBACK(err,j) RET_MINUS1_FAIL(err)
static int compileOpcodesInternal(compileContext *ctx, opcodeRec *op, unsigned char *bufOut, int bufOut_len, int *computTableSize, const namespaceInformation *namespacePathToThis, int *calledRvType, int preferredReturnValues, int *fpStackUse, int *canHaveDenormalOutput)
{
  int rv_offset=0, denormal_force=-1;
  if (!op) RET_MINUS1_FAIL("coi !op")
  *fpStackUse=0;
  for (;;)
  {
    // special case: statement delimiting means we can process the left side into place, and iteratively do the second parameter without recursing
    // also we don't need to save/restore anything to the stack (which the normal 2 parameter function processing does)
    if (op->opcodeType == OPCODETYPE_FUNC2 && op->fntype == FN_JOIN_STATEMENTS)
    {
      int fUse1;
      int parm_size = compileOpcodes(ctx,op->parms.parms[0],bufOut,bufOut_len, computTableSize, namespacePathToThis, RETURNVALUE_IGNORE, NULL,&fUse1,NULL);
      if (parm_size < 0) RET_MINUS1_FAIL("coc join fail")
      op = op->parms.parms[1];
      if (!op) RET_MINUS1_FAIL("join got to null")
      if (fUse1>*fpStackUse) *fpStackUse=fUse1;
      if (bufOut) bufOut += parm_size;
      bufOut_len -= parm_size;
      rv_offset += parm_size;
      denormal_force=-1;
    }
    // special case: __denormal_likely(), __denormal_unlikely()
    else if (op->opcodeType == OPCODETYPE_FUNC1 && (op->fntype == FN_DENORMAL_LIKELY || op->fntype == FN_DENORMAL_UNLIKELY))
    {
      denormal_force = op->fntype == FN_DENORMAL_LIKELY;
      op = op->parms.parms[0];
    }
    else 
    {
      break;
    }
  }
  if (denormal_force >= 0 && canHaveDenormalOutput)
  {  
    *canHaveDenormalOutput = denormal_force;
    canHaveDenormalOutput = &denormal_force; // prevent it from being changed by functions below
  }
  // special case: BAND/BOR
  if (op->opcodeType == OPCODETYPE_FUNC2 && (op->fntype == FN_LOGICAL_AND || op->fntype == FN_LOGICAL_OR))
  {
    int fUse1=0;
    int parm_size;
    int retType=RETURNVALUE_IGNORE;
    if (preferredReturnValues != RETURNVALUE_IGNORE) retType = RETURNVALUE_BOOL;
    *calledRvType = retType;
    parm_size = compileOpcodes(ctx,op->parms.parms[0],bufOut,bufOut_len, computTableSize, namespacePathToThis, RETURNVALUE_BOOL, NULL, &fUse1, NULL);
    if (parm_size < 0) RET_MINUS1_FAIL("loop band/bor coc fail")
    if (fUse1 > *fpStackUse) *fpStackUse=fUse1;
    {
      int sz2, fUse2=0;
      unsigned char *destbuf;
      const int testsz=op->fntype == FN_LOGICAL_OR ? sizeof(GLUE_JMP_IF_P1_NZ) : sizeof(GLUE_JMP_IF_P1_Z);
      if (bufOut_len < parm_size+testsz) RET_MINUS1_FAIL_FALLBACK("band/bor size fail",doNonInlinedAndOr_)
      if (bufOut)  memcpy(bufOut+parm_size,op->fntype == FN_LOGICAL_OR ? GLUE_JMP_IF_P1_NZ : GLUE_JMP_IF_P1_Z,testsz); 
      parm_size += testsz;
      destbuf = bufOut + parm_size;
      sz2= compileOpcodes(ctx,op->parms.parms[1],bufOut?bufOut+parm_size:NULL,bufOut_len-parm_size, computTableSize, namespacePathToThis, retType, NULL,&fUse2, NULL);
      CHECK_SIZE_FORJMP(sz2,doNonInlinedAndOr_)
      if (sz2<0) RET_MINUS1_FAIL("band/bor coc fail")
      parm_size+=sz2;
      if (bufOut) GLUE_JMP_SET_OFFSET(destbuf, (bufOut + parm_size) - destbuf);
      if (fUse2 > *fpStackUse) *fpStackUse=fUse2;
      return rv_offset + parm_size;
    }
  }  
  if (op->opcodeType == OPCODETYPE_FUNC3 && op->fntype == FN_IF_ELSE) // special case: IF
  {
    int fUse1=0;
    int use_rv = RETURNVALUE_IGNORE;
    int rvMode=0;
    int parm_size = compileOpcodes(ctx,op->parms.parms[0],bufOut,bufOut_len, computTableSize, namespacePathToThis, RETURNVALUE_BOOL|RETURNVALUE_BOOL_REVERSED, &rvMode,&fUse1, NULL);
    if (parm_size < 0) RET_MINUS1_FAIL("if coc fail")
    if (fUse1 > *fpStackUse) *fpStackUse=fUse1;
    if (preferredReturnValues & RETURNVALUE_NORMAL) use_rv=RETURNVALUE_NORMAL;
    else if (preferredReturnValues & RETURNVALUE_FPSTACK) use_rv=RETURNVALUE_FPSTACK;
    else if (preferredReturnValues & RETURNVALUE_BOOL) use_rv=RETURNVALUE_BOOL;
    *calledRvType = use_rv;
    {
      int csz,hasSecondHalf;
      if (rvMode & RETURNVALUE_BOOL_REVERSED)
      {
        if (bufOut_len < parm_size + (int)sizeof(GLUE_JMP_IF_P1_NZ)) RET_MINUS1_FAIL_FALLBACK("if size fail",doNonInlineIf_)
        if (bufOut) memcpy(bufOut+parm_size,GLUE_JMP_IF_P1_NZ,sizeof(GLUE_JMP_IF_P1_NZ));
        parm_size += sizeof(GLUE_JMP_IF_P1_NZ);
      }
      else
      {
        if (bufOut_len < parm_size + (int)sizeof(GLUE_JMP_IF_P1_Z)) RET_MINUS1_FAIL_FALLBACK("if size fail",doNonInlineIf_)
        if (bufOut) memcpy(bufOut+parm_size,GLUE_JMP_IF_P1_Z,sizeof(GLUE_JMP_IF_P1_Z));
        parm_size += sizeof(GLUE_JMP_IF_P1_Z);
      }
      csz=compileOpcodes(ctx,op->parms.parms[1],bufOut ? bufOut+parm_size : NULL,bufOut_len - parm_size, computTableSize, namespacePathToThis, use_rv, NULL,&fUse1, canHaveDenormalOutput);
      if (fUse1 > *fpStackUse) *fpStackUse=fUse1;
      hasSecondHalf = preferredReturnValues || !OPCODE_IS_TRIVIAL(op->parms.parms[2]);
      CHECK_SIZE_FORJMP(csz,doNonInlineIf_)
      if (csz<0) RET_MINUS1_FAIL("if coc fial")
      if (bufOut) GLUE_JMP_SET_OFFSET(bufOut + parm_size, csz + (hasSecondHalf?sizeof(GLUE_JMP_NC):0));
      parm_size+=csz;
      if (hasSecondHalf)
      {
        if (bufOut_len < parm_size + (int)sizeof(GLUE_JMP_NC)) RET_MINUS1_FAIL_FALLBACK("if len fail",doNonInlineIf_)
        if (bufOut) memcpy(bufOut+parm_size,GLUE_JMP_NC,sizeof(GLUE_JMP_NC));
        parm_size+=sizeof(GLUE_JMP_NC);
        csz=compileOpcodes(ctx,op->parms.parms[2],bufOut ? bufOut+parm_size : NULL,bufOut_len - parm_size, computTableSize, namespacePathToThis, use_rv, NULL, &fUse1, canHaveDenormalOutput);
        CHECK_SIZE_FORJMP(csz,doNonInlineIf_)
        if (csz<0) RET_MINUS1_FAIL("if coc 2 fail")
        // update jump address
        if (bufOut) GLUE_JMP_SET_OFFSET(bufOut + parm_size,csz); 
        parm_size+=csz;       
        if (fUse1 > *fpStackUse) *fpStackUse=fUse1;
      }
      return rv_offset + parm_size;
    }
  }
  {
    // special case: while
    if (op->opcodeType == OPCODETYPE_FUNC1 && op->fntype == FN_WHILE)
    {
      *calledRvType = RETURNVALUE_BOOL;
	  {
		  unsigned char *jzoutpt;
		  unsigned char *looppt;
		  int parm_size = 0, subsz;
		  if (bufOut_len < parm_size + (int)(GLUE_WHILE_SETUP_SIZE + sizeof(GLUE_WHILE_BEGIN))) RET_MINUS1_FAIL("while size fail 1")
			  if (bufOut) memcpy(bufOut + parm_size, GLUE_WHILE_SETUP, GLUE_WHILE_SETUP_SIZE);
		  parm_size += GLUE_WHILE_SETUP_SIZE;
		  looppt = bufOut + parm_size;
		  if (bufOut) memcpy(bufOut + parm_size, GLUE_WHILE_BEGIN, sizeof(GLUE_WHILE_BEGIN));
		  parm_size += sizeof(GLUE_WHILE_BEGIN);
		  subsz = compileOpcodes(ctx, op->parms.parms[0], bufOut ? (bufOut + parm_size) : NULL, bufOut_len - parm_size, computTableSize, namespacePathToThis, RETURNVALUE_BOOL, NULL, fpStackUse, NULL);
		  if (subsz < 0) RET_MINUS1_FAIL("while coc fail")
			  if (bufOut_len < parm_size + (int)(sizeof(GLUE_WHILE_END) + sizeof(GLUE_WHILE_CHECK_RV))) RET_MINUS1_FAIL("which size fial 2")
				  parm_size += subsz;
		  if (bufOut) memcpy(bufOut + parm_size, GLUE_WHILE_END, sizeof(GLUE_WHILE_END));
		  parm_size += sizeof(GLUE_WHILE_END);
		  jzoutpt = bufOut + parm_size;
		  if (bufOut) memcpy(bufOut + parm_size, GLUE_WHILE_CHECK_RV, sizeof(GLUE_WHILE_CHECK_RV));
		  parm_size += sizeof(GLUE_WHILE_CHECK_RV);
		  if (bufOut)
		  {
			  GLUE_JMP_SET_OFFSET(bufOut + parm_size, (looppt - (bufOut + parm_size)));
			  GLUE_JMP_SET_OFFSET(jzoutpt, (bufOut + parm_size) - jzoutpt);
		  }
		  return rv_offset + parm_size;
	  }
    }
    // special case: loop
    if (op->opcodeType == OPCODETYPE_FUNC2 && op->fntype == FN_LOOP)
    {
      int fUse1;
      int parm_size = compileOpcodes(ctx,op->parms.parms[0],bufOut,bufOut_len, computTableSize, namespacePathToThis, RETURNVALUE_FPSTACK, NULL,&fUse1, NULL);
      if (parm_size < 0) RET_MINUS1_FAIL("loop coc fail")
      *calledRvType = RETURNVALUE_BOOL;
      if (fUse1 > *fpStackUse) *fpStackUse=fUse1;
	  {
		  int subsz;
		  int fUse2 = 0;
		  unsigned char *skipptr1, *loopdest;
		  if (bufOut_len < parm_size + (int)(sizeof(GLUE_LOOP_LOADCNT) + GLUE_LOOP_CLAMPCNT_SIZE + GLUE_LOOP_BEGIN_SIZE)) RET_MINUS1_FAIL("loop size fail")
			  // store, convert to int, compare against 1, if less than, skip to end
			  if (bufOut) memcpy(bufOut + parm_size, GLUE_LOOP_LOADCNT, sizeof(GLUE_LOOP_LOADCNT));
		  parm_size += sizeof(GLUE_LOOP_LOADCNT);
		  skipptr1 = bufOut + parm_size;
		  // compare aginst max loop length, jump to loop start if not above it
		  if (bufOut) memcpy(bufOut + parm_size, GLUE_LOOP_CLAMPCNT, GLUE_LOOP_CLAMPCNT_SIZE);
		  parm_size += GLUE_LOOP_CLAMPCNT_SIZE;
		  // loop code:
		  loopdest = bufOut + parm_size;
		  if (bufOut) memcpy(bufOut + parm_size, GLUE_LOOP_BEGIN, GLUE_LOOP_BEGIN_SIZE);
		  parm_size += GLUE_LOOP_BEGIN_SIZE;
		  subsz = compileOpcodes(ctx, op->parms.parms[1], bufOut ? (bufOut + parm_size) : NULL, bufOut_len - parm_size, computTableSize, namespacePathToThis, RETURNVALUE_IGNORE, NULL, &fUse2, NULL);
		  if (subsz < 0) RET_MINUS1_FAIL("loop coc fail")
			  if (fUse2 > *fpStackUse) *fpStackUse = fUse2;
		  parm_size += subsz;
		  if (bufOut_len < parm_size + (int)sizeof(GLUE_LOOP_END)) RET_MINUS1_FAIL("loop size fail 2")
			  if (bufOut) memcpy(bufOut + parm_size, GLUE_LOOP_END, sizeof(GLUE_LOOP_END));
		  parm_size += sizeof(GLUE_LOOP_END);
		  if (bufOut)
		  {
			  GLUE_JMP_SET_OFFSET(bufOut + parm_size, loopdest - (bufOut + parm_size));
			  GLUE_JMP_SET_OFFSET(skipptr1, (bufOut + parm_size) - skipptr1);
		  }
		  return rv_offset + parm_size;
	}
    }    
  }
  switch (op->opcodeType)
  {
    case OPCODETYPE_DIRECTVALUE:
        if (preferredReturnValues == RETURNVALUE_BOOL)
        {
          int w = fabs(op->parms.dv.directValue) >= NSEEL_CLOSEFACTOR;
          int wsz=(w?sizeof(GLUE_SET_P1_NZ):sizeof(GLUE_SET_P1_Z));
          *calledRvType = RETURNVALUE_BOOL;
          if (bufOut_len < wsz) RET_MINUS1_FAIL("direct int size fail3")
          if (bufOut) memcpy(bufOut,w?GLUE_SET_P1_NZ:GLUE_SET_P1_Z,wsz);
          return rv_offset+wsz;
        }
        else if (preferredReturnValues & RETURNVALUE_FPSTACK)
        {
#ifdef GLUE_HAS_FLDZ
          if (op->parms.dv.directValue == 0.0)
          {
            *fpStackUse = 1;
            *calledRvType = RETURNVALUE_FPSTACK;
            if (bufOut_len < sizeof(GLUE_FLDZ)) RET_MINUS1_FAIL("direct fp fail 1")
            if (bufOut) memcpy(bufOut,GLUE_FLDZ,sizeof(GLUE_FLDZ));
            return rv_offset+sizeof(GLUE_FLDZ);
          }
#endif
#ifdef GLUE_HAS_FLD1
          if (op->parms.dv.directValue == 1.0)
          {
            *fpStackUse = 1;
            *calledRvType = RETURNVALUE_FPSTACK;
            if (bufOut_len < sizeof(GLUE_FLD1)) RET_MINUS1_FAIL("direct fp fail 1")
            if (bufOut) memcpy(bufOut,GLUE_FLD1,sizeof(GLUE_FLD1));
            return rv_offset+sizeof(GLUE_FLD1);
          }
#endif
        }
        // fall through
    case OPCODETYPE_DIRECTVALUE_TEMPSTRING:
    case OPCODETYPE_VALUE_FROM_NAMESPACENAME:
    case OPCODETYPE_VARPTR:
    case OPCODETYPE_VARPTRPTR:
      #ifdef GLUE_MOV_PX_DIRECTVALUE_TOSTACK_SIZE
        if (OPCODE_IS_TRIVIAL(op))
        {
          if (preferredReturnValues & RETURNVALUE_FPSTACK)
          {
            *fpStackUse = 1;
            if (bufOut_len < GLUE_MOV_PX_DIRECTVALUE_TOSTACK_SIZE) RET_MINUS1_FAIL("direct fp fail 2")
            if (bufOut)
            {
              if (generateValueToReg(ctx,op,bufOut,-1,namespacePathToThis, 1 /*allow caching*/)<0) RET_MINUS1_FAIL("direct fp fail gvr")
            }
            *calledRvType = RETURNVALUE_FPSTACK;
            return rv_offset+GLUE_MOV_PX_DIRECTVALUE_TOSTACK_SIZE;
          }
        }
      #endif
      if (bufOut_len < GLUE_MOV_PX_DIRECTVALUE_SIZE) 
      {
        RET_MINUS1_FAIL("direct value fail 1")
      }
      if (bufOut) 
      {
        if (generateValueToReg(ctx,op,bufOut,0,namespacePathToThis, !!(preferredReturnValues&RETURNVALUE_FPSTACK)/*cache if going to the fp stack*/)<0) RET_MINUS1_FAIL("direct value gvr fail3")
      }
    return rv_offset + GLUE_MOV_PX_DIRECTVALUE_SIZE;
    case OPCODETYPE_FUNCX:
    case OPCODETYPE_FUNC1:
    case OPCODETYPE_FUNC2:
    case OPCODETYPE_FUNC3:
      if (op->fntype == FUNCTYPE_EELFUNC)
      {
        int a;
        a = compileEelFunctionCall(ctx,op,bufOut,bufOut_len,computTableSize,namespacePathToThis, calledRvType,fpStackUse,canHaveDenormalOutput);
        if (a<0) return a;
        rv_offset += a;
      }
      else
      {
        int a;
        a = compileNativeFunctionCall(ctx,op,bufOut,bufOut_len,computTableSize,namespacePathToThis, calledRvType,fpStackUse,preferredReturnValues,canHaveDenormalOutput);
        if (a<0)return a;
        rv_offset += a;
      }        
    return rv_offset;
  }
  RET_MINUS1_FAIL("default opcode fail")
}
int compileOpcodes(compileContext *ctx, opcodeRec *op, unsigned char *bufOut, int bufOut_len, int *computTableSize, const namespaceInformation *namespacePathToThis, 
                   int supportedReturnValues, int *rvType, int *fpStackUse, int *canHaveDenormalOutput)
{
  int code_returns=RETURNVALUE_NORMAL;
  int fpsu=0;
  int codesz;
  int denorm=0;
  codesz = compileOpcodesInternal(ctx,op,bufOut,bufOut_len,computTableSize,namespacePathToThis,&code_returns, supportedReturnValues,&fpsu,&denorm);
  if (denorm && canHaveDenormalOutput) *canHaveDenormalOutput=1;
  if (codesz < 0) return codesz;
  if (fpStackUse) *fpStackUse=fpsu;
  if (bufOut) bufOut += codesz;
  bufOut_len -= codesz;
  if (code_returns == RETURNVALUE_BOOL && !(supportedReturnValues & RETURNVALUE_BOOL) && supportedReturnValues)
  {
    int stubsize;
    void *stub = GLUE_realAddress(nseel_asm_booltofp,nseel_asm_booltofp_end,&stubsize);
    if (!stub || bufOut_len < stubsize) RET_MINUS1_FAIL(stub?"booltofp size":"booltfp addr")
    if (bufOut) 
    {
      memcpy(bufOut,stub,stubsize);
      bufOut += stubsize;
    }
    codesz+=stubsize;
    bufOut_len -= stubsize;
    code_returns = RETURNVALUE_FPSTACK;
  }
  // default processing of code_returns to meet return value requirements
  if (supportedReturnValues & code_returns) 
  {
    if (rvType) *rvType = code_returns;
    return codesz;
  }
  if (rvType) *rvType = RETURNVALUE_IGNORE;
  if (code_returns == RETURNVALUE_NORMAL)
  {
    if (supportedReturnValues & (RETURNVALUE_FPSTACK|RETURNVALUE_BOOL))
    {
      if (bufOut_len < GLUE_PUSH_VAL_AT_PX_TO_FPSTACK_SIZE) RET_MINUS1_FAIL("pushvalatpxtofpstack,size")
      if (bufOut) 
      {
        GLUE_PUSH_VAL_AT_PX_TO_FPSTACK(bufOut,0); // always fld qword [eax] but we might change that later
        bufOut += GLUE_PUSH_VAL_AT_PX_TO_FPSTACK_SIZE;
      }
      codesz += GLUE_PUSH_VAL_AT_PX_TO_FPSTACK_SIZE;  
      bufOut_len -= GLUE_PUSH_VAL_AT_PX_TO_FPSTACK_SIZE;
      if (supportedReturnValues & RETURNVALUE_BOOL) 
      {
        code_returns = RETURNVALUE_FPSTACK;
      }
      else
      {
        if (rvType) *rvType = RETURNVALUE_FPSTACK;
      }
    }
  }
  if (code_returns == RETURNVALUE_FPSTACK)
  {
    if (supportedReturnValues & (RETURNVALUE_BOOL|RETURNVALUE_BOOL_REVERSED))
    {
      int stubsize;
      void *stub;
      if (supportedReturnValues & RETURNVALUE_BOOL_REVERSED)
      {
        if (rvType) *rvType = RETURNVALUE_BOOL_REVERSED;
        stub = GLUE_realAddress(nseel_asm_fptobool_rev,nseel_asm_fptobool_rev_end,&stubsize);
      }
      else
      {
        if (rvType) *rvType = RETURNVALUE_BOOL;
        stub = GLUE_realAddress(nseel_asm_fptobool,nseel_asm_fptobool_end,&stubsize);
      }
      if (!stub || bufOut_len < stubsize) RET_MINUS1_FAIL(stub?"fptobool size":"fptobool addr")
      if (bufOut) 
      {
        memcpy(bufOut,stub,stubsize);
        bufOut += stubsize;
      }
      codesz+=stubsize;
      bufOut_len -= stubsize;
    }
    else if (supportedReturnValues & RETURNVALUE_NORMAL)
    {
      if (computTableSize) (*computTableSize) ++;
      if (bufOut_len < GLUE_POP_FPSTACK_TO_WTP_TO_PX_SIZE) RET_MINUS1_FAIL("popfpstacktowtptopxsize")
      // generate fp-pop to temp space
      if (bufOut) GLUE_POP_FPSTACK_TO_WTP_TO_PX(bufOut,0);
      codesz+=GLUE_POP_FPSTACK_TO_WTP_TO_PX_SIZE;
      if (rvType) *rvType = RETURNVALUE_NORMAL;
    }
    else
    {
      // toss return value that will be ignored
      if (bufOut_len < GLUE_POP_FPSTACK_SIZE) RET_MINUS1_FAIL("popfpstack size")
      if (bufOut) memcpy(bufOut,GLUE_POP_FPSTACK,GLUE_POP_FPSTACK_SIZE);   
      codesz+=GLUE_POP_FPSTACK_SIZE;
    }
  }
  return codesz;
}
//------------------------------------------------------------------------------
NSEEL_CODEHANDLE NSEEL_code_compile(NSEEL_VMCTX _ctx, const char *_expression, int lineoffs)
{
  return NSEEL_code_compile_ex(_ctx,_expression,lineoffs,0);
}
typedef struct topLevelCodeSegmentRec {
  struct topLevelCodeSegmentRec *_next;
  void *code;
  int codesz;
  int tmptable_use;
} topLevelCodeSegmentRec;
NSEEL_CODEHANDLE NSEEL_code_compile_ex(NSEEL_VMCTX _ctx, const char *_expression, int lineoffs, int compile_flags)
{
  compileContext *ctx = (compileContext *)_ctx;
  const char *endptr;
  const char *_expression_end;
  codeHandleType *handle;
  topLevelCodeSegmentRec *startpts_tail=NULL;
  topLevelCodeSegmentRec *startpts=NULL;
  _codeHandleFunctionRec *oldCommonFunctionList;
  int curtabptr_sz=0;
  void *curtabptr=NULL;
  int had_err=0;
  if (!ctx) return 0;
  ctx->directValueCache=0;
  ctx->gotEndOfInput=0;
  if (compile_flags & NSEEL_CODE_COMPILE_FLAG_COMMONFUNCS_RESET)
  {
    ctx->functions_common=NULL; // reset common function list
  }
  else
  {
    // reset common compiled function code, forcing a recompile if shared
    _codeHandleFunctionRec *a = ctx->functions_common;
    while (a)
    {
      _codeHandleFunctionRec *b = a->derivedCopies;
      if (a->localstorage) 
      {
        // force local storage actual values to be reallocated if used again
        memset(a->localstorage,0,sizeof(EEL_F *) * a->localstorage_size);
      }
      a->startptr = NULL; // force this copy to be recompiled
      while (b)
      {
        b->startptr = NULL; // force derived copies to get recompiled
        // no need to reset b->localstorage, since it points to a->localstorage
        b=b->derivedCopies;
      }
      a=a->next;
    }
  }
  ctx->last_error_string[0]=0;
  if (!_expression || !*_expression) return 0;
  _expression_end = _expression + strlen(_expression);
  oldCommonFunctionList = ctx->functions_common;
  ctx->isGeneratingCommonFunction=0;
  ctx->isSharedFunctions = !!(compile_flags & NSEEL_CODE_COMPILE_FLAG_COMMONFUNCS);
  ctx->functions_local = NULL;
  freeBlocks(&ctx->tmpblocks_head);  // free blocks
  freeBlocks(&ctx->blocks_head);  // free blocks
  freeBlocks(&ctx->blocks_head_data);  // free blocks
  memset(ctx->l_stats,0,sizeof(ctx->l_stats));
  handle = (codeHandleType*)newDataBlock(sizeof(codeHandleType),8);
  if (!handle) 
  {
    return 0;
  }
  memset(handle,0,sizeof(codeHandleType));
  ctx->l_stats[0] += (int)(_expression_end - _expression);
  ctx->tmpCodeHandle = handle;
  endptr=_expression;
  while (*endptr)
  {
    int computTableTop = 0;
    int startptr_size=0;
    void *startptr=NULL;
    opcodeRec *start_opcode=NULL;
    const char *expr=endptr;
    int function_numparms=0;
    char is_fname[NSEEL_MAX_VARIABLE_NAMELEN+1];
    is_fname[0]=0;
    memset(ctx->function_localTable_Size,0,sizeof(ctx->function_localTable_Size));
    memset(ctx->function_localTable_Names,0,sizeof(ctx->function_localTable_Names));
    ctx->function_localTable_ValuePtrs=0;
    ctx->function_usesNamespaces=0;
    ctx->function_curName=NULL;
    ctx->function_globalFlag=0;
    ctx->errVar=0;
    // single out top level segment
    {
      int had_something = 0, pcnt=0, pcnt2=0;
      int state=0;
      for (;;)
      {
        int l;
        const char *p=nseel_simple_tokenizer(&endptr,_expression_end,&l,&state);
        if (!p) 
        {
          if (pcnt || pcnt2) ctx->gotEndOfInput|=4;
          break;
        }
        if (*p == ';') 
        {
          if (had_something && !pcnt && !pcnt2) break;
        }
        else if (*p == '/' && l > 1 && (p[1] == '/' || p[1] == '*')) 
        {
			if (l > 19 && !strncmp(p, "//#eel-no-optimize:", 19))
			{
				//0 = atoi(p + 19);
			}
        }
        else
        {
          if (!had_something) 
          {
            expr = p;
            had_something = 1;
          }
          if (*p == '(') pcnt++;
          else if (*p == ')')  {  if (--pcnt<0) pcnt=0; }
          else if (*p == '[') pcnt2++;
          else if (*p == ']')  {  if (--pcnt2<0) pcnt2=0; }
        }
      }
      if (!*expr || !had_something) break;
    }
    // parse   
    {
      int tmplen,funcname_len;
      const char *p = expr;
      const char *tok1 = nseel_simple_tokenizer(&p,endptr,&tmplen,NULL);
      const char *funcname = nseel_simple_tokenizer(&p,endptr,&funcname_len,NULL);
      if (tok1 && funcname && tmplen == 8 && !strncmp(tok1,"function",8) && (isalpha(funcname[0]) || funcname[0] == '_'))
      {
        int had_parms_locals=0;
        if (funcname_len > (int)(sizeof(is_fname)-1))
        	funcname_len = (int)(sizeof(is_fname)-1);
        memcpy(is_fname, funcname, funcname_len);
        is_fname[funcname_len]=0;
        ctx->function_curName = is_fname; // only assigned for the duration of the loop, cleared later //-V507
        while (NULL != (tok1 = nseel_simple_tokenizer(&p,endptr,&tmplen,NULL)))
        {
          int is_parms = 0, localTableContext = 0;
          int maxcnt=0;
          const char *sp_save;
          if (tok1[0] == '(')
          {
            if (had_parms_locals) 
            {
              expr = p-1; // begin compilation at this code!
              break;
            }
            is_parms = 1;
          }
          else
          {
            if (tmplen == 5 && !strncmp(tok1,"local",tmplen)) localTableContext=0;
            else if (tmplen == 6 && !strncmp(tok1,"static",tmplen)) localTableContext=0;
            else if (tmplen == 8 && !strncmp(tok1,"instance",tmplen)) localTableContext=1;
            else if ((tmplen == 7 && !strncmp(tok1,"globals",tmplen))  ||
                     (tmplen == 6 && !strncmp(tok1,"global",tmplen)))
            {
              ctx->function_globalFlag = 1;
              localTableContext=2;
            }
            else break; // unknown token!
            tok1 = nseel_simple_tokenizer(&p,endptr,&tmplen,NULL);
            if (!tok1 || tok1[0] != '(') break;
          }
          had_parms_locals = 1;
          sp_save=p;
          while (NULL != (tok1 = nseel_simple_tokenizer(&p,endptr,&tmplen,NULL)))
          {
            if (tok1[0] == ')') break;
            if (*tok1 == '#' && localTableContext!=1 && localTableContext!=2) 
            {
              ctx->errVar = (int) (tok1 - _expression);
              lstrcpyn_safe(ctx->last_error_string,"#string can only be in instance() or globals()",sizeof(ctx->last_error_string));
              goto had_error;
            }
            if (isalpha(*tok1) || *tok1 == '_' || *tok1 == '#') 
            {
              maxcnt++;
              if (p < endptr && *p == '*')
              {
                if (!is_parms && localTableContext!=2)
                {
                  ctx->errVar = (int) (p - _expression);
                  lstrcpyn_safe(ctx->last_error_string,"namespace* can only be used in parameters or globals()",sizeof(ctx->last_error_string));
                  goto had_error;
                }
                p++;
              }
            }
            else if (*tok1 != ',')
            {
              ctx->errVar = (int)(tok1 - _expression);
              lstrcpyn_safe(ctx->last_error_string,"unknown character in function parameters",sizeof(ctx->last_error_string));
              goto had_error;
            }
          }
          if (tok1 && maxcnt > 0)
          {
            char **ot = ctx->function_localTable_Names[localTableContext];
            const int osz = ctx->function_localTable_Size[localTableContext];            
            maxcnt += osz;
            ctx->function_localTable_Names[localTableContext] = (char **)newTmpBlock(ctx,sizeof(char *) * maxcnt);
            if (ctx->function_localTable_Names[localTableContext])
            {
              int i=osz;
              if (osz && ot) memcpy(ctx->function_localTable_Names[localTableContext],ot,sizeof(char *) * osz);
              p=sp_save;
              while (NULL != (tok1 = nseel_simple_tokenizer(&p,endptr,&tmplen,NULL)))
              {
                if (tok1[0] == ')') break;
                if (isalpha(*tok1) || *tok1 == '_' || *tok1 == '#') 
                {
                  char *newstr;
                  int l = tmplen;
                  if (*p == '*')  // xyz* for namespace
                  {
                    p++;
                    l++;
                  }
                  if (l > NSEEL_MAX_VARIABLE_NAMELEN) l = NSEEL_MAX_VARIABLE_NAMELEN;
                  newstr = newTmpBlock(ctx,l+1);
                  if (newstr)
                  {
                    memcpy(newstr,tok1,l);
                    newstr[l]=0;
                    ctx->function_localTable_Names[localTableContext][i++] = newstr;
                  }
                }
              }
              ctx->function_localTable_Size[localTableContext]=i;
              if (is_parms) function_numparms = i;
            }         
          }
        }
      }
    }
    if (ctx->function_localTable_Size[0]>0)
    {
      ctx->function_localTable_ValuePtrs = 
          ctx->isSharedFunctions ? newDataBlock(ctx->function_localTable_Size[0] * sizeof(EEL_F *),8) : 
                                   newTmpBlock(ctx,ctx->function_localTable_Size[0] * sizeof(EEL_F *)); 
      if (!ctx->function_localTable_ValuePtrs)
      {
        ctx->function_localTable_Size[0]=0;
        function_numparms=0;
      }
      else
      {
        memset(ctx->function_localTable_ValuePtrs,0,sizeof(EEL_F *) * ctx->function_localTable_Size[0]); // force values to be allocated
      }
    }
   {
     int nseelparse(compileContext* context);
     void nseelrestart (void *input_file ,void *yyscanner );
     ctx->rdbuf_start = _expression;
	 ctx->rdbuf = expr;
	 ctx->rdbuf_end = endptr;
	 if (!nseelparse(ctx) && !ctx->errVar)
		 start_opcode = ctx->result;
     ctx->rdbuf = NULL;
   }
    if (start_opcode)
    {
      int rvMode=0, fUse=0;
	  optimizeOpcodes(ctx, start_opcode, is_fname[0] ? 1 : 0);
      startptr_size = compileOpcodes(ctx,start_opcode,NULL,1024*1024*256,NULL, NULL, 
        is_fname[0] ? (RETURNVALUE_NORMAL|RETURNVALUE_FPSTACK) : RETURNVALUE_IGNORE, &rvMode, &fUse, NULL); // if not a function, force return value as address (avoid having to pop it ourselves
      // if a function, allow the code to decide how return values are generated
      if (is_fname[0])
      {
        _codeHandleFunctionRec *fr = ctx->isSharedFunctions ? newDataBlock(sizeof(_codeHandleFunctionRec),8) : 
                                        newTmpBlock(ctx,sizeof(_codeHandleFunctionRec)); 
        if (fr)
        {
          memset(fr,0,sizeof(_codeHandleFunctionRec));
          fr->startptr_size = startptr_size;
          fr->opcodes = start_opcode;
          fr->rvMode = rvMode;
          fr->fpStackUsage=fUse;
          fr->tmpspace_req = computTableTop;
          if (ctx->function_localTable_Size[0] > 0 && ctx->function_localTable_ValuePtrs)
          {
            if (ctx->function_localTable_Names[0])
            {
              int i;
              for(i=0;i<function_numparms;i++)
              {
                const char *nptr = ctx->function_localTable_Names[0][i];
                if (nptr && *nptr && nptr[strlen(nptr)-1] == '*') 
                {
                  fr->parameterAsNamespaceMask |= ((unsigned int)1)<<i;
                }
              }
            }
            fr->num_params=function_numparms;
            fr->localstorage = ctx->function_localTable_ValuePtrs;
            fr->localstorage_size = ctx->function_localTable_Size[0];
          }
          fr->usesNamespaces = ctx->function_usesNamespaces;
          fr->isCommonFunction = ctx->isSharedFunctions;
          lstrcpyn_safe(fr->fname,is_fname,sizeof(fr->fname));
          if (ctx->isSharedFunctions)
          {
            fr->next = ctx->functions_common;
            ctx->functions_common = fr;
          }
          else
          {
            fr->next = ctx->functions_local;
            ctx->functions_local = fr;
          }         
        }
        continue;
      }
      if (!startptr_size) continue; // optimized away
      if (startptr_size>0)
      {
        startptr = newTmpBlock(ctx,startptr_size);
        if (startptr)
        {
          startptr_size=compileOpcodes(ctx,start_opcode,(unsigned char*)startptr,startptr_size,&computTableTop, NULL, RETURNVALUE_IGNORE, NULL,NULL, NULL);
          if (startptr_size<=0) startptr = NULL;
        }
      }
    }
    if (!startptr) 
    {  
had_error:
		//if (!ctx->last_error_string[0])
		{
			int byteoffs = ctx->errVar;
			int linenumber;
			char cur_err[sizeof(ctx->last_error_string)];
			lstrcpyn_safe(cur_err, ctx->last_error_string, sizeof(cur_err));
			if (cur_err[0]) lstrcatn(cur_err, ": ", sizeof(cur_err));
			else lstrcpyn_safe(cur_err, "syntax error: ", sizeof(cur_err));
			if (_expression + byteoffs >= _expression_end)
			{
				if (ctx->gotEndOfInput & 4) byteoffs = (int)(expr - _expression);
				else byteoffs = (int)(_expression_end - _expression);
			}
			if (byteoffs < 0) byteoffs = 0;
			linenumber = findLineNumber(_expression, byteoffs) + 1;
			if (ctx->gotEndOfInput & 4)
			{
				WDL_snprintf(ctx->last_error_string, sizeof(ctx->last_error_string), "%d: %smissing ) or ]", linenumber + lineoffs, cur_err);
			}
			else
			{
				const char *p = _expression + byteoffs;
				int x = 0, right_amt_nospace = 0, left_amt_nospace = 0;
				while (x < 32 && p - x > _expression && p[-x] != '\r' && p[-x] != '\n')
				{
					if (!isspace(p[-x])) left_amt_nospace = x;
					x++;
				}
				x = 0;
				while (x < 60 && p[x] && p[x] != '\r' && p[x] != '\n')
				{
					if (!isspace(p[x])) right_amt_nospace = x;
					x++;
				}
				if (right_amt_nospace < 1) right_amt_nospace = 1;
				// display left_amt >>>> right_amt_nospace
				if (left_amt_nospace > 0)
					WDL_snprintf(ctx->last_error_string, sizeof(ctx->last_error_string), "%d: %s'%.*s <!> %.*s'", linenumber + lineoffs, cur_err,
						left_amt_nospace, p - left_amt_nospace,
						right_amt_nospace, p);
				else
					WDL_snprintf(ctx->last_error_string, sizeof(ctx->last_error_string), "%d: %s'%.*s'", linenumber + lineoffs, cur_err, right_amt_nospace, p);
			}
	}
		startpts = NULL;
		startpts_tail = NULL;
		had_err = 1;
		break;
    }
    if (!is_fname[0]) // redundant check (if is_fname[0] is set and we succeeded, it should continue)
                      // but we'll be on the safe side
    {
      topLevelCodeSegmentRec *p = newTmpBlock(ctx,sizeof(topLevelCodeSegmentRec));
      p->_next=0;
      p->code = startptr;
      p->codesz = startptr_size;
      p->tmptable_use = computTableTop;
      if (!startpts_tail) startpts_tail=startpts=p;
      else
      {
        startpts_tail->_next=p;
        startpts_tail=p;
      }
      if (curtabptr_sz < computTableTop)
      {
        curtabptr_sz=computTableTop;
      }
    }
  }
  memset(ctx->function_localTable_Size,0,sizeof(ctx->function_localTable_Size));
  memset(ctx->function_localTable_Names,0,sizeof(ctx->function_localTable_Names));
  ctx->function_localTable_ValuePtrs=0;
  ctx->function_usesNamespaces=0;
  ctx->function_curName=NULL;
  ctx->function_globalFlag=0;
  ctx->tmpCodeHandle = NULL;
  if (handle->want_stack)
  {
    if (!handle->stack) startpts=NULL;
  }
  if (startpts) 
  {
    curtabptr_sz += 2; // many functions use the worktable for temporary storage of up to 2 EEL_F's
    handle->workTable_size = curtabptr_sz;
    handle->workTable = curtabptr = newDataBlock((curtabptr_sz+MIN_COMPUTABLE_SIZE + COMPUTABLE_EXTRA_SPACE) * sizeof(EEL_F),32);
    if (!curtabptr) startpts=NULL;
  }
  if (startpts || (!had_err && (compile_flags & NSEEL_CODE_COMPILE_FLAG_COMMONFUNCS)))
  {
    unsigned char *writeptr;
    topLevelCodeSegmentRec *p=startpts;
    int size=sizeof(GLUE_RET); // for ret at end :)
    int wtpos=0;
    // now we build one big code segment out of our list of them, inserting a mov esi, computable before each item as necessary
    while (p)
    {
      if (wtpos <= 0)
      {
        wtpos=MIN_COMPUTABLE_SIZE;
        size += GLUE_RESET_WTP(NULL,0);
      }
      size+=p->codesz;
      wtpos -= p->tmptable_use;
      p=p->_next;
    }
    handle->code = newCodeBlock(size,32);
    if (handle->code)
    {
      writeptr=(unsigned char *)handle->code;
      p=startpts;
      wtpos=0;
      while (p)
      {
        if (wtpos <= 0)
        {
          wtpos=MIN_COMPUTABLE_SIZE;
          writeptr+=GLUE_RESET_WTP(writeptr,curtabptr);
        }
        memcpy(writeptr,(char*)p->code,p->codesz);
        writeptr += p->codesz;
        wtpos -= p->tmptable_use;
        p=p->_next;
      }
      memcpy(writeptr,&GLUE_RET,sizeof(GLUE_RET)); writeptr += sizeof(GLUE_RET);
      ctx->l_stats[1]=size;
      handle->code_size = (int) (writeptr - (unsigned char *)handle->code);
#ifdef __arm__
      __clear_cache(handle->code,writeptr);
#endif
    }
    handle->blocks = ctx->blocks_head;
    handle->blocks_data = ctx->blocks_head_data;
    ctx->blocks_head=0;
    ctx->blocks_head_data=0;
  }
  else
  {
    // failed compiling, or failed calloc()
    handle=NULL;              // return NULL (after resetting blocks_head)
  }
  ctx->directValueCache=0;
  ctx->functions_local = NULL;
  ctx->isGeneratingCommonFunction=0;
  ctx->isSharedFunctions=0;
  freeBlocks(&ctx->tmpblocks_head);  // free blocks
  freeBlocks(&ctx->blocks_head);  // free blocks of code (will be nonzero only on error)
  freeBlocks(&ctx->blocks_head_data);  // free blocks of data (will be nonzero only on error)
  if (handle)
  {
    handle->ramPtr = ctx->ram_state.blocks;
    memcpy(handle->code_stats,ctx->l_stats,sizeof(ctx->l_stats));
    nseel_evallib_stats[0]+=ctx->l_stats[0];
    nseel_evallib_stats[1]+=ctx->l_stats[1];
    nseel_evallib_stats[2]+=ctx->l_stats[2];
    nseel_evallib_stats[3]+=ctx->l_stats[3];
    nseel_evallib_stats[4]++;
  }
  else
  {
    ctx->functions_common = oldCommonFunctionList; // failed compiling, remove any added common functions from the list
    // remove any derived copies of functions due to error, since we may have added some that have been freed
    while (oldCommonFunctionList)
    {
      oldCommonFunctionList->derivedCopies=NULL;
      oldCommonFunctionList=oldCommonFunctionList->next;
    }
  }
  memset(ctx->l_stats,0,sizeof(ctx->l_stats));
  return (NSEEL_CODEHANDLE)handle;
}
//------------------------------------------------------------------------------
void NSEEL_code_execute(NSEEL_CODEHANDLE code)
{
  codeHandleType *h = (codeHandleType *)code;
  INT_PTR codeptr = (INT_PTR) h->code;
  INT_PTR tabptr = (INT_PTR)h->workTable;
  GLUE_CALL_CODE(tabptr,codeptr,(INT_PTR)h->ramPtr);
}
int NSEEL_code_geterror_flag(NSEEL_VMCTX ctx)
{
  compileContext *c=(compileContext *)ctx;
  if (c) return (c->gotEndOfInput ? 1 : 0);
  return 0;
}
char *NSEEL_code_getcodeerror(NSEEL_VMCTX ctx)
{
  compileContext *c=(compileContext *)ctx;
  if (ctx && c->last_error_string[0]) return c->last_error_string;
  return 0;
}
//------------------------------------------------------------------------------
void NSEEL_code_free(NSEEL_CODEHANDLE code)
{
  codeHandleType *h = (codeHandleType *)code;
  if (h != NULL)
  {
    nseel_evallib_stats[0]-=h->code_stats[0];
    nseel_evallib_stats[1]-=h->code_stats[1];
    nseel_evallib_stats[2]-=h->code_stats[2];
    nseel_evallib_stats[3]-=h->code_stats[3];
    nseel_evallib_stats[4]--;
	freeBlocks(&h->blocks);
    freeBlocks(&h->blocks_data);
  }
}
//------------------------------------------------------------------------------
void NSEEL_VM_freevars(NSEEL_VMCTX _ctx)
{
  if (_ctx)
  {
    compileContext *ctx=(compileContext *)_ctx;
    free(ctx->varTable_Values);
    free(ctx->varTable_Names);
    ctx->varTable_Values=0;
    ctx->varTable_Names=0;
    ctx->varTable_numBlocks=0;
  }
}
NSEEL_VMCTX NSEEL_VM_alloc() // return a handle
{
  compileContext *ctx=calloc(1,sizeof(compileContext));
  if (ctx) 
  {
	  ctx->caller_this = ctx;
	  ctx->scanner = ctx;
    ctx->ram_state.maxblocks = NSEEL_RAM_BLOCKS_DEFAULTMAX;
	ctx->ram_state.closefact = NSEEL_CLOSEFACTOR;
	ctx->onString = addStringCallback;
	ctx->m_string_context = (eel_string_context_state*)malloc(sizeof(eel_string_context_state));
	Initeel_string_context_state(ctx->m_string_context);
  }
  return ctx;
}
int NSEEL_VM_setramsize(NSEEL_VMCTX _ctx, int maxent)
{
  compileContext *ctx = (compileContext *)_ctx;
  if (!ctx) return 0;
  if (maxent > 0)
  {
    maxent = (maxent + NSEEL_RAM_ITEMSPERBLOCK - 1)/NSEEL_RAM_ITEMSPERBLOCK;
    if (maxent > NSEEL_RAM_BLOCKS) maxent = NSEEL_RAM_BLOCKS;
    ctx->ram_state.maxblocks = maxent;
  }
  return ctx->ram_state.maxblocks * NSEEL_RAM_ITEMSPERBLOCK;
}
void NSEEL_VM_SetFunctionValidator(NSEEL_VMCTX _ctx, const char * (*validateFunc)(const char *fn_name, void *user), void *user)
{
  if (_ctx)
  {
    compileContext *ctx = (compileContext *)_ctx;
    ctx->func_check = validateFunc;
    ctx->func_check_user = user;
  }
}
void NSEEL_VM_free(NSEEL_VMCTX _ctx) // free when done with a VM and ALL of its code have been freed, as well
{
  if (_ctx)
  {
    compileContext *ctx=(compileContext *)_ctx;
	if (ctx->numSinks)
	{
		for (int i = 0; i < ctx->numSinks; i++)
			if (ctx->inputSink[i])
				free(ctx->inputSink[i]);
		free(ctx->sinksMap);
		free(ctx->sinksLength);
		free(ctx->inputSink);
	}
    NSEEL_VM_freevars(_ctx);
    NSEEL_VM_freeRAM(_ctx);
    freeBlocks(&ctx->pblocks);
    // these should be 0 normally but just in case
    freeBlocks(&ctx->tmpblocks_head);  // free blocks
    freeBlocks(&ctx->blocks_head);  // free blocks
    freeBlocks(&ctx->blocks_head_data);  // free blocks
    ctx->scanner=0;
	Freeel_string_context_state(ctx->m_string_context);
	free(ctx->m_string_context);
    free(ctx);
  }
}
int *NSEEL_code_getstats(NSEEL_CODEHANDLE code)
{
  codeHandleType *h = (codeHandleType *)code;
  if (h)
  {
    return h->code_stats;
  }
  return 0;
}
void NSEEL_VM_SetStringFunc(NSEEL_VMCTX ctx, EEL_F (*onString)(void *caller_this, eelStringSegmentRec *list))
{
  if (ctx)
  {
    compileContext *c=(compileContext*)ctx;
    c->onString = onString;
  }
}
void *NSEEL_PProc_RAM(void *data, int data_size, compileContext *ctx)
{
  if (data_size>0) data=EEL_GLUE_set_immediate(data, (INT_PTR)ctx->ram_state.blocks); 
  return data;
}
void *NSEEL_PProc_THIS(void *data, int data_size, compileContext *ctx)
{
  if (data_size>0) data=EEL_GLUE_set_immediate(data, (INT_PTR)ctx->caller_this);
  return data;
}
void NSEEL_VM_remove_unused_vars(NSEEL_VMCTX _ctx)
{
  compileContext *ctx = (compileContext *)_ctx;
  int wb;
  if (ctx) for (wb = 0; wb < ctx->varTable_numBlocks; wb ++)
  {
    int ti;
    char **plist=ctx->varTable_Names[wb];
    if (!plist) break;
    for (ti = 0; ti < NSEEL_VARS_PER_BLOCK; ti ++)
    {        
      if (plist[ti])
      {
        varNameHdr *v = ((varNameHdr*)plist[ti])-1;
        if (!v->refcnt && !v->isreg) 
        {
          plist[ti]=NULL;
        }
      }
    }
  }
}
void NSEEL_VM_remove_all_nonreg_vars(NSEEL_VMCTX _ctx)
{
  compileContext *ctx = (compileContext *)_ctx;
  int wb;
  if (ctx) for (wb = 0; wb < ctx->varTable_numBlocks; wb ++)
  {
    int ti;
    char **plist=ctx->varTable_Names[wb];
    if (!plist) break;
    for (ti = 0; ti < NSEEL_VARS_PER_BLOCK; ti ++)
    {        
      if (plist[ti])
      {
        varNameHdr *v = ((varNameHdr*)plist[ti])-1;
        if (!v->isreg) 
        {
          plist[ti]=NULL;
        }
      }
    }
  }
}
void NSEEL_VM_clear_var_refcnts(NSEEL_VMCTX _ctx)
{
  compileContext *ctx = (compileContext *)_ctx;
  int wb;
  if (ctx) for (wb = 0; wb < ctx->varTable_numBlocks; wb ++)
  {
    int ti;
    char **plist=ctx->varTable_Names[wb];
    if (!plist) break;
    for (ti = 0; ti < NSEEL_VARS_PER_BLOCK; ti ++)
    {        
      if (plist[ti])
      {
        varNameHdr *v = ((varNameHdr*)plist[ti])-1;
        v->refcnt=0;
      }
    }
  }
}
EEL_F *nseel_int_register_var(compileContext *ctx, const char *name, int isReg, const char **namePtrOut)
{
  int match_wb = -1, match_ti=-1;
  int wb;
  int ti=0;
  for (wb = 0; wb < ctx->varTable_numBlocks; wb ++)
  {
    char **plist=ctx->varTable_Names[wb];
    if (!plist) return NULL; // error!
    for (ti = 0; ti < NSEEL_VARS_PER_BLOCK; ti ++)
    { 
      if (!plist[ti])
      {
        if (match_wb < 0)
        {
          match_wb=wb;
          match_ti=ti;
        }
      }
      else if (!strncmp(plist[ti],name,NSEEL_MAX_VARIABLE_NAMELEN))
      {
        varNameHdr *v = ((varNameHdr*)plist[ti])-1;
        if (isReg < 0)
        {
          EEL_F *p; 
          return (ctx->varTable_Values && NULL != (p = ctx->varTable_Values[wb])) ? p + ti : NULL;
        }
        v->refcnt++;
        if (isReg) v->isreg=isReg;
        if (namePtrOut) *namePtrOut = plist[ti];
        break;
      }
    }
    if (ti < NSEEL_VARS_PER_BLOCK) break;
  }
  if (isReg < 0) return NULL;
  if (wb == ctx->varTable_numBlocks && match_wb >=0 && match_ti >= 0)
  {
    wb = match_wb;
    ti = match_ti;
  }
  if (wb == ctx->varTable_numBlocks)
  {
    ti=0;
    // add new block
    if (!(ctx->varTable_numBlocks&(NSEEL_VARS_MALLOC_CHUNKSIZE-1)) || !ctx->varTable_Values || !ctx->varTable_Names )
    {
      void *nv = realloc(ctx->varTable_Values,(ctx->varTable_numBlocks+NSEEL_VARS_MALLOC_CHUNKSIZE) * sizeof(EEL_F *));
      if (!nv) return NULL;
      ctx->varTable_Values = (EEL_F **)nv;
      nv = realloc(ctx->varTable_Names,(ctx->varTable_numBlocks+NSEEL_VARS_MALLOC_CHUNKSIZE) * sizeof(char **));
      if (!nv) return NULL;
      ctx->varTable_Names = (char ***)nv;
    }
    ctx->varTable_numBlocks++;
    ctx->varTable_Values[wb] = (EEL_F *)newCtxDataBlock(sizeof(EEL_F)*NSEEL_VARS_PER_BLOCK,8);
    ctx->varTable_Names[wb] = (char **)newCtxDataBlock(sizeof(char *)*NSEEL_VARS_PER_BLOCK,1);
    if (ctx->varTable_Values[wb])
    {
      memset(ctx->varTable_Values[wb],0,sizeof(EEL_F)*NSEEL_VARS_PER_BLOCK);
    }
    if (ctx->varTable_Names[wb])
    {
      memset(ctx->varTable_Names[wb],0,sizeof(char *)*NSEEL_VARS_PER_BLOCK);
    }
  }
  if (!ctx->varTable_Names[wb] || !ctx->varTable_Values[wb]) return NULL;
  if (!ctx->varTable_Names[wb][ti])
  {
    size_t l = strlen(name);
    char *b;
    varNameHdr *vh;
    if (l > NSEEL_MAX_VARIABLE_NAMELEN) l = NSEEL_MAX_VARIABLE_NAMELEN;
    b=newCtxDataBlock( (int) (sizeof(varNameHdr) + l+1),1);
    if (!b) return NULL; // malloc fail
    vh=(varNameHdr *)b;
    vh->refcnt=1;
    vh->isreg=isReg;
    b+=sizeof(varNameHdr);
    memcpy(b,name,l);
    b[l] = 0;
    ctx->varTable_Names[wb][ti] = b;
    ctx->varTable_Values[wb][ti]=0.0;
    if (namePtrOut) *namePtrOut = b;
  }
  return ctx->varTable_Values[wb] + ti;
}
//------------------------------------------------------------------------------
void NSEEL_VM_enumallvars(NSEEL_VMCTX ctx, int (*func)(const char *name, EEL_F *val, void *ctx), void *userctx)
{
  compileContext *tctx = (compileContext *) ctx;
  int wb;
  if (!tctx) return;
  for (wb = 0; wb < tctx->varTable_numBlocks; wb ++)
  {
    int ti;
    char **plist=tctx->varTable_Names[wb];
    if (!plist) break;
    for (ti = 0; ti < NSEEL_VARS_PER_BLOCK; ti ++)
    {              
      if (plist[ti] && !func(plist[ti],tctx->varTable_Values[wb] + ti,userctx)) break;
    }
    if (ti < NSEEL_VARS_PER_BLOCK)
      break;
  }
}
//------------------------------------------------------------------------------
EEL_F *NSEEL_VM_regvar(NSEEL_VMCTX _ctx, const char *var)
{
  compileContext *ctx = (compileContext *)_ctx;
  if (!ctx) return 0;
  return nseel_int_register_var(ctx,var,1,NULL);
}
EEL_F *NSEEL_VM_getvar(NSEEL_VMCTX _ctx, const char *var)
{
  compileContext *ctx = (compileContext *)_ctx;
  if (!ctx) return 0;
  return nseel_int_register_var(ctx,var,-1,NULL);
}
int  NSEEL_VM_get_var_refcnt(NSEEL_VMCTX _ctx, const char *name)
{
  compileContext *ctx = (compileContext *)_ctx;
  int wb;
  if (!ctx) return -1;
  for (wb = 0; wb < ctx->varTable_numBlocks; wb ++)
  {
    int ti;
    if (!ctx->varTable_Values[wb] || !ctx->varTable_Names[wb]) break;
    for (ti = 0; ti < NSEEL_VARS_PER_BLOCK; ti ++)
    {        
      if (ctx->varTable_Names[wb][ti] && !stricmp(ctx->varTable_Names[wb][ti],name)) 
      {
        varNameHdr *h = ((varNameHdr *)ctx->varTable_Names[wb][ti])-1;
        return h->refcnt;
      }
    }
  }
  return -1;
}
int NSEEL_VM_openText2Sink(NSEEL_VMCTX _ctx, char *filename, int slot)
{
	compileContext *c = (compileContext *)_ctx;
	int idx = arySearch(c->sinksMap, c->numSinks, slot);
	if (idx >= 0)
	{
		if (c->sinksMap[idx])
			free(c->inputSink[idx]);
	}
	else
	{
		c->numSinks++;
		c->sinksMap = (int*)realloc(c->sinksMap, c->numSinks * sizeof(int));
		c->sinksLength = (int*)realloc(c->sinksLength, c->numSinks * sizeof(int));
		c->inputSink = (EEL_F**)realloc(c->inputSink, c->numSinks * sizeof(EEL_F*));
		idx = c->numSinks - 1;
	}
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("failed to open file\n");
		return 0;
	}
	c->sinksMap[idx] = slot;
	int counter = 0;
	double fval;
	while (fscanf(fp, "%lf", &fval) != EOF)
		counter++;
	c->inputSink[idx] = (EEL_F*)malloc(counter * sizeof(EEL_F));
	c->sinksLength[idx] = counter;
	fseek(fp, 0, SEEK_SET);
	for (int i = 0; i < counter; i++)
		fscanf(fp, "%lf", &c->inputSink[idx][i]);
	fclose(fp);
	return counter;
}
opcodeRec *nseel_createFunctionByName(compileContext *ctx, const char *name, int np, opcodeRec *code1, opcodeRec *code2, opcodeRec *code3)
{
  int i;
  for (i = 0; fnTable1 + i; i++)
  {
	  if (i >= (int)(sizeof(fnTable1) / sizeof(fnTable1[0])))
		  break;
    functionType *f = fnTable1 + i;
    if ((f->nParams&FUNCTIONTYPE_PARAMETERCOUNTMASK) == np && !strcmp(f->name, name))
    {
      opcodeRec *o=newOpCode(ctx,NULL, np==3?OPCODETYPE_FUNC3:np==2?OPCODETYPE_FUNC2:OPCODETYPE_FUNC1);
      if (o) 
      {
        o->fntype = FUNCTYPE_FUNCTIONTYPEREC;
        o->fn = f;
        o->parms.parms[0]=code1;
        o->parms.parms[1]=code2;
        o->parms.parms[2]=code3;
      }
      return o;
    }
  }
  return NULL;
}
#include <float.h>
//------------------------------------------------------------------------------
opcodeRec *nseel_translate(compileContext *ctx, const char *tmp, size_t tmplen) // tmplen 0 = null term
{
  // this depends on the string being nul terminated eventually, tmplen is used more as a hint than anything else
  if ((tmp[0] == '0' || tmp[0] == '$') && toupper(tmp[1])=='X')
  {
    char *p;
    return nseel_createCompiledValue(ctx,(EEL_F)strtoul(tmp+2,&p,16));
  }
  else if (tmp[0] == '$')
  {
    if (tmp[1] == '~')
    {
      char *p=(char*)tmp+2;
      unsigned int v=strtoul(tmp+2,&p,10);
      if (v>53) v=53;
      return nseel_createCompiledValue(ctx,(EEL_F)((((WDL_INT64)1) << v) - 1));
    }
    else if (!tmplen ? !stricmp(tmp,"$E") : (tmplen == 2 && !strnicmp(tmp,"$E",2)))
      return nseel_createCompiledValue(ctx,(EEL_F)2.71828182845904523536);
    else if (!tmplen ? !stricmp(tmp, "$PI") : (tmplen == 3 && !strnicmp(tmp, "$PI", 3)))
      return nseel_createCompiledValue(ctx,(EEL_F)3.141592653589793238463);
    else if (!tmplen ? !stricmp(tmp, "$PHI") : (tmplen == 4 && !strnicmp(tmp, "$PHI", 4)))
      return nseel_createCompiledValue(ctx,(EEL_F)1.618033988749894848205);
	else if (!tmplen ? !stricmp(tmp, "$EPS") : (tmplen == 4 && !strnicmp(tmp, "$EPS", 4)))
		return nseel_createCompiledValue(ctx, (EEL_F)DBL_EPSILON);
    else if ((!tmplen || tmplen == 4) && tmp[1] == '\'' && tmp[2] && tmp[3] == '\'')
      return nseel_createCompiledValue(ctx,(EEL_F)tmp[2]);      
    else return NULL;
  }
  else if (tmp[0] == '\'')
  {
    char b[64];
    int x,sz;
    unsigned int rv=0;
    if (!tmplen) // nul terminated tmplen, calculate a workable length
    {
      // faster than strlen(tmp) if tmp is large, we'll never need more than ~18 chars anyway
      while (tmplen < 32 && tmp[tmplen]) tmplen++;
    }
    sz = tmplen > 0 ? nseel_filter_escaped_string(b,sizeof(b),tmp+1, tmplen - 1, '\'') : 0;
    if (sz > 4) 
    {
      if (ctx->last_error_string[0]) lstrcatn(ctx->last_error_string, ", ", sizeof(ctx->last_error_string));
      snprintf_append(ctx->last_error_string,sizeof(ctx->last_error_string),"multi-byte character '%.5s...' too long",b);
      return NULL; // do not allow 'xyzxy', limit to 4 bytes
    }
    for (x=0;x<sz;x++) rv = (rv<<8) + ((unsigned char*)b)[x];
    return nseel_createCompiledValue(ctx,(EEL_F)rv);
  }
  return nseel_createCompiledValue(ctx,(EEL_F)atof(tmp));
}
