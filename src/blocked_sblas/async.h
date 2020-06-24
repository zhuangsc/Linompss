#ifndef __ASYNC_H__
#define __ASYNC_H__


#ifdef SINGLE_PRECISION

#define ASYNC_BREAK		sasync_break
#define ASYNC_CONV		sasync_conv

#else

#define ASYNC_BREAK		dasync_break
#define ASYNC_CONV		dasync_conv

#endif


#include "async_struct.h"


#define LIBSBBLAS_EXPORT __attribute__((__visibility__("default")))


extern LIBSBBLAS_EXPORT async_stat_t sasync_break(int it, float *p, float *c, volatile int *state, int force);
extern LIBSBBLAS_EXPORT async_stat_t dasync_break(int it, double *p, double *c, volatile int *state, int force);

extern LIBSBBLAS_EXPORT async_stat_t sasync_conv(int it, float res, float *c, float div, async_t *sync, int force);
extern LIBSBBLAS_EXPORT async_stat_t dasync_conv(int it, double res, double *c, double div, async_t *sync, int force);


#endif // __ASYNC_H__
