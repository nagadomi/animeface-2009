#ifndef __NV_PORTABLE_H
#define __NV_PORTABLE_H

#ifdef __cplusplus
extern "C" {
#endif

// 標準Cヘッダー

//#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>

// OpenMP 
#if (_MSC_VER >= 1400 && defined(_OPENMP))
#ifdef _DEBUG
    #undef _DEBUG
    #include <omp.h>
    #define _DEBUG
#else
    #include <omp.h>
#endif
#endif

#if NV_ENABLE_CLOCK
unsigned long nv_clock(void);
#endif
#if NV_ENABLE_SLEEP
void nv_sleep(unsigned int msec);
#endif

// aligne
#if __GNUC__
#define NV_ALIGNED(typed, variable, align_size) typed variable __attribute__((aligned(align_size)))
#else // MSC
#define NV_ALIGNED(typed, variable, align_size) __declspec(align(align_size)) typed variable 
#endif

// snprintf
#ifdef _MSC_VER
#define snprintf _snprintf
#endif

// roundf
#define NV_ROUND_INT(x) ((int)(0.5f + x))

// min max
#ifndef max
#define max(a, b) ((a) > (b) ? (a):(b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a):(b))
#endif




#ifdef __cplusplus
}
#endif


#endif
