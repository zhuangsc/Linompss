#ifndef __CONVERTMAT_H__
#define __CONVERTMAT_H__


#if DOUBLE_PRECISION

#define MATCONVERT		convert_d2s

#else

#define MATCONVERT		convert_s2d

#endif


float *convert_d2s(int m, int n, double *D, int ldimD);  
double *convert_s2d(int m, int n, float *S, int ldimS);


#endif // __CONVERTMAT_H__
