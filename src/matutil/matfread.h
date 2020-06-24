#ifndef __MATFREAD_H__
#define __MATFREAD_H__


#include <stdio.h>


#if DOUBLE_PRECISION
#define READ_MM2DENSE		read_mm2ddense
#else
#define READ_MM2DENSE		read_mm2sdense
#endif

int read_mm2sdense(FILE *f, int m, int n, float *A);
int read_mm2ddense(FILE *f, int m, int n, double *A);


#endif // __MATFREAD_H__

