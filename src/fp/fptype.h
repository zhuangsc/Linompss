#ifndef __FPT_H__
#define __FPT_H__


#include <math.h>
#include <float.h>


#define I_ONE		i_one
#define I_MONE		i_mone
#define FP_INF		INFINITY


#ifdef SINGLE_PRECISION

typedef float fp_t;

#define FP_ONE 		s_one
#define FP_MONE 	s_mone
#define FP_NOUGHT	s_zero
#define FP_MAXVAL	FLT_MAX

#define FP_RAND 	drand48
#define FP_SEED 	srand48
#define FP_EXP 		frexpf
#define FP_LOG10 	log10f
#define FP_SQRT 	sqrtf
#define FP_POW 		powf
#define FP_ABS 		fabsf	

#define FP_SCANSPEC	scan_sconspec
#endif

#ifdef DOUBLE_PRECISION

typedef double fp_t;

#define FP_ONE 		d_one
#define FP_MONE 	d_mone
#define FP_NOUGHT 	d_zero
#define FP_MAXVAL	DBL_MAX

#define FP_SQRT 	sqrt
#define FP_RAND 	drand48
#define FP_SEED 	srand48
#define FP_ABS 		fabs	
#define FP_EXP 		frexp
#define FP_LOG10 	log10
#define FP_POW 		pow

#define FP_SCANSPEC	scan_dconspec

#endif


extern int i_one;
extern int i_mone;

extern double d_one;
extern double d_zero;
extern double d_mone;

extern float s_one;
extern float s_mone;
extern float s_zero;

extern const char *scan_dconspec;
extern const char *scan_sconspec;


#endif // __FPT_H__
