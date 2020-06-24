#ifndef __GEQRF_MAIN_H__
#define __GEQRF_MAIN_H__


void sgeqrf(int diagl, int tr, int tc, int bs, int mt, int nt, float **A, float **T, float **S);
void dgeqrf(int diagl, int tr, int tc, int bs, int mt, int nt, double **A, double **T, double **S);


#endif // __GEQRF_MAIN_H__
