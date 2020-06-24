#ifndef __DENSEUTIL_H__
#define __DENSEUTIL_H__


#include "fpmatr.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef SINGLE_PRECISION

#define DMAT_CP				dmat_scp
#define DMAT_RELERR			dmat_srelerr
#define VECTOR_2NORM		svector_2norm
#define MAT_NORMDIFF		dmat_snormdiff
#define MAT_1NORM			smat_1norm
#define CMH2CM             	scmh2cm
#define MAT_ENERGY_SYM		smat_energy_sym
#define VEC_SANITY_CHECK	svec_sanity_check

#else

#define DMAT_CP				dmat_dcp
#define DMAT_RELERR			dmat_drelerr
#define VECTOR_2NORM		dvector_2norm
#define MAT_NORMDIFF		dmat_dnormdiff
#define MAT_1NORM			dmat_1norm
#define CMH2CM             	dcmh2cm
#define MAT_ENERGY_SYM		dmat_energy_sym
#define VEC_SANITY_CHECK	dvec_sanity_check

#endif


double* dmat_dcp(int m, int n, double *A, int lda);
float* dmat_scp(int m, int n, float *A, int lda);

float  dmat_srelerr(ompssblas_t trian, int m, int n, float *A, float* B);
double dmat_drelerr(ompssblas_t trian, int m, int n, double *A, double* B);

float  svector_2norm(float *v, int length);
double dvector_2norm(double *v, int length);

float dmat_snormdiff(char norm, int m, int n, float *A, int lda, float *A0, int lda0);
double dmat_dnormdiff(char norm, int m, int n, double *A, int lda, double *A0, int lda0);

float* scmh2cm(int m, int n, int tr, int tc, float *Ah, int ldAh); 
double* dcmh2cm(int m, int n, int tr, int tc, double *Ah, int ldAh); 

float smat_1norm(float *mat, int m, int n);
double dmat_1norm(double *mat, int m, int n);


float smat_energy_sym(int n, int lda, int bs, float *mat, float *dist);
double dmat_energy_sym(int n, int lda, int bs, double *mat, double *dist);


void svec_sanity_check(char *n, float *vec, int len);
void dvec_sanity_check(char *n, double *vec, int len);

#endif // __DENSEUTIL_H__
