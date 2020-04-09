/*
  Nullsoft Expression Evaluator Library (NS-EEL)
  Copyright (C) 1999-2003 Nullsoft, Inc.
  ns-eel.h: main application interface header
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
#ifndef __NS_EEL_H__
#define __NS_EEL_H__
// put standard includes here
#include <stdlib.h>
#include <stdio.h>
//#define EEL_F_SIZE 4
#ifndef EEL_F_SIZE
#define EEL_F_SIZE 8
#endif
#ifdef _MSC_VER
typedef __int64 WDL_INT64;
typedef unsigned __int64 WDL_UINT64;
#else
typedef long long WDL_INT64;
typedef unsigned long long WDL_UINT64;
#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <stdint.h>
typedef intptr_t INT_PTR;
typedef uintptr_t UINT_PTR;
#endif
#ifdef __GNUC__
// for structures that contain doubles, or doubles in structures that are after stuff of questionable alignment (for OSX/linux)
#define WDL_FIXALIGN  __attribute__ ((aligned (8)))
// usage: void func(int a, const char *fmt, ...) WDL_VARARG_WARN(printf,2,3); // note: if member function, this pointer is counted as well, so as member function that would be 3,4
#define WDL_VARARG_WARN(x,n,s) __attribute__ ((format (x,n,s)))
#else
#define WDL_FIXALIGN 
#define WDL_VARARG_WARN(x,n,s)
#endif
#if !defined(max) && defined(WDL_DEFINE_MINMAX)
#define max(x,y) ((x)<(y)?(y):(x))
#define min(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef wdl_max
#define wdl_max(x,y) ((x)<(y)?(y):(x))
#define wdl_min(x,y) ((x)<(y)?(x):(y))
#define wdl_abs(x) ((x)<0 ? -(x) : (x))
#endif
#ifndef _WIN32
#ifndef strnicmp 
#define strnicmp(x,y,z) strncasecmp(x,y,z)
#endif
#ifndef stricmp 
#define stricmp(x,y) strcasecmp(x,y)
#endif
#endif
#if EEL_F_SIZE == 4
typedef float EEL_F;
typedef float *EEL_F_PTR;
#else
typedef double EEL_F WDL_FIXALIGN;
typedef double *EEL_F_PTR;
#endif
#ifdef _MSC_VER
#define NSEEL_CGEN_CALL __cdecl
#else
#define NSEEL_CGEN_CALL 
#endif
#ifdef __cplusplus
extern "C" {
#endif
// host should implement these (can be empty stub functions if no VM will execute code in multiple threads at once)
  // implement if you will be running the code in same VM from multiple threads, 
  // or VMs that have the same GRAM pointer from different threads, or multiple
  // VMs that have a NULL GRAM pointer from multiple threads.
  // if you give each VM it's own unique GRAM and only run each VM in one thread, then you can leave it blank.
  // or if you're daring....
void NSEEL_HOSTSTUB_EnterMutex();
void NSEEL_HOSTSTUB_LeaveMutex();
typedef void *NSEEL_VMCTX;
typedef void *NSEEL_CODEHANDLE;
void NSEEL_quit(NSEEL_VMCTX _ctx); // clears any added functions
int NSEEL_VM_openText2Sink(NSEEL_VMCTX _ctx, char *filename, int slot);
int *NSEEL_getstats(); // returns a pointer to 5 ints... source bytes, static code bytes, call code bytes, data bytes, number of code handles
void NSEEL_VM_freevars(NSEEL_VMCTX _ctx);
NSEEL_VMCTX NSEEL_VM_alloc(); // return a handle
void NSEEL_VM_free(NSEEL_VMCTX ctx); // free when done with a VM and ALL of its code have been freed, as well
// validateFunc can return error message if not permitted
void NSEEL_VM_SetFunctionValidator(NSEEL_VMCTX, const char * (*validateFunc)(const char *fn_name, void *user), void *user);
void NSEEL_VM_remove_unused_vars(NSEEL_VMCTX _ctx);
void NSEEL_VM_clear_var_refcnts(NSEEL_VMCTX _ctx);
void NSEEL_VM_remove_all_nonreg_vars(NSEEL_VMCTX _ctx);
void NSEEL_VM_enumallvars(NSEEL_VMCTX ctx, int (*func)(const char *name, EEL_F *val, void *ctx), void *userctx); // return false from func to stop
EEL_F *NSEEL_VM_regvar(NSEEL_VMCTX ctx, const char *name); // register a variable (before compilation)
EEL_F *NSEEL_VM_getvar(NSEEL_VMCTX ctx, const char *name); // get a variable (if registered or created by code)
int  NSEEL_VM_get_var_refcnt(NSEEL_VMCTX _ctx, const char *name); // returns -1 if not registered, or >=0
void NSEEL_VM_freeRAM(NSEEL_VMCTX ctx); // clears and frees all (VM) RAM used
void NSEEL_VM_freeRAMIfCodeRequested(NSEEL_VMCTX); // call after code to free the script-requested memory
int NSEEL_VM_wantfreeRAM(NSEEL_VMCTX ctx); // want NSEEL_VM_freeRAMIfCodeRequested?
EEL_F *NSEEL_VM_getramptr(NSEEL_VMCTX ctx, unsigned int offs, int *validCount);
EEL_F *NSEEL_VM_getramptr_noalloc(NSEEL_VMCTX ctx, unsigned int offs, int *validCount);
// set 0 to query. returns actual value used (limits, granularity apply -- see NSEEL_RAM_BLOCKS)
int NSEEL_VM_setramsize(NSEEL_VMCTX ctx, int maxent);
typedef struct eelstrSegRec
{
  struct eelstrSegRec *_next;
  const char *str_start; // escaped characters, including opening/trailing characters
  int str_len; 
} eelStringSegmentRec;
void NSEEL_VM_SetStringFunc(NSEEL_VMCTX ctx, EEL_F (*onString)(void *caller_this, eelStringSegmentRec *list));
NSEEL_CODEHANDLE NSEEL_code_compile(NSEEL_VMCTX ctx, const char *code, int lineoffs);
#define NSEEL_CODE_COMPILE_FLAG_COMMONFUNCS 1 // allows that code's functions to be used in other code (note you shouldn't destroy that codehandle without destroying others first if used)
#define NSEEL_CODE_COMPILE_FLAG_COMMONFUNCS_RESET 2 // resets common code functions
NSEEL_CODEHANDLE NSEEL_code_compile_ex(NSEEL_VMCTX ctx, const char *code, int lineoffs, int flags);
char *NSEEL_code_getcodeerror(NSEEL_VMCTX ctx);
int NSEEL_code_geterror_flag(NSEEL_VMCTX ctx);
void NSEEL_code_execute(NSEEL_CODEHANDLE code);
void NSEEL_code_free(NSEEL_CODEHANDLE code);
int *NSEEL_code_getstats(NSEEL_CODEHANDLE code); // 4 ints...source bytes, static code bytes, call code bytes, data bytes
// global memory control/view
extern int NSEEL_RAM_memused_errors;
// configuration:
#define NSEEL_MAX_VARIABLE_NAMELEN 128  // define this to override the max variable length
#define NSEEL_MAX_EELFUNC_PARAMETERS 40
#define NSEEL_MAX_FUNCSIG_NAME 2048 // longer than variable maxlen, due to multiple namespaces
#include <limits.h>
#define NSEEL_LOOPFUNC_SUPPORT_MAXLEN INT_MAX
#define NSEEL_MAX_FUNCTION_SIZE_FOR_INLINE 2048
//#define EEL_DUMP_OPS // used for testing frontend parser/logic changes
// 512 * 65536 = 32 million entries maximum (256MB RAM)
// default is limited to 128 * 65536 = 8 million entries (64MB RAM)
// default to 8 million entries, use NSEEL_VM_setramsize() to change at runtime
#define NSEEL_RAM_BLOCKS_DEFAULTMAX 128
// 512 entry block table maximum (2k/4k per VM)
#define NSEEL_RAM_BLOCKS_LOG2 9
 // 262144 items per block (2048KB)
#define NSEEL_RAM_ITEMSPERBLOCK_LOG2 18
#define NSEEL_RAM_BLOCKS (1 << NSEEL_RAM_BLOCKS_LOG2)
#define NSEEL_RAM_ITEMSPERBLOCK (1<<NSEEL_RAM_ITEMSPERBLOCK_LOG2)
#define NSEEL_STACK_SIZE 32768 // about 64k overhead if the stack functions are used in a given code handle
#define EEL_BC_TYPE int
#ifdef __cplusplus
}
#endif
#endif//__NS_EEL_H__
