#ifndef _WDL_FFT_H_
#define _WDL_FFT_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "eelCommon.h"

typedef struct {
  EEL_F re;
  EEL_F im;
} WDL_FFT_COMPLEX;

extern void WDL_fft_init();

extern void WDL_fft_complexmul(WDL_FFT_COMPLEX *dest, WDL_FFT_COMPLEX *src, int len);
extern void WDL_fft_complexmul2(WDL_FFT_COMPLEX *dest, WDL_FFT_COMPLEX *src, WDL_FFT_COMPLEX *src2, int len);
extern void WDL_fft_complexmul3(WDL_FFT_COMPLEX *destAdd, WDL_FFT_COMPLEX *src, WDL_FFT_COMPLEX *src2, int len);

/* Expects WDL_FFT_COMPLEX input[0..len-1] scaled by 1.0/len, returns
WDL_FFT_COMPLEX output[0..len-1] order by WDL_fft_permute(len). */
extern void WDL_fft(WDL_FFT_COMPLEX *, int len, int isInverse);

/* Expects EEL_F input[0..len-1] scaled by 0.5/len, returns
WDL_FFT_COMPLEX output[0..len/2-1], for len >= 4 order by
WDL_fft_permute(len/2). Note that output[len/2].re is stored in
output[0].im. */
extern void WDL_real_fft(EEL_F *, int len, int isInverse);

extern int WDL_fft_permute(int fftsize, int idx);
extern int *WDL_fft_permute_tab(int fftsize);

#ifndef EEL_FFT_MINBITLEN
#define EEL_FFT_MINBITLEN 4
#endif

#ifndef EEL_FFT_MAXBITLEN
#define EEL_FFT_MAXBITLEN 15
#endif

#ifndef EEL_FFT_MINBITLEN_REORDER
#define EEL_FFT_MINBITLEN_REORDER (EEL_FFT_MINBITLEN-1)
#endif
extern int *fft_reorder_table_for_bitsize(int bitsz);
extern void initFFTData();
#ifdef __cplusplus
};
#endif

#endif