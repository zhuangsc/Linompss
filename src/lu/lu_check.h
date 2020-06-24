#ifndef __LU_CHECK_H__
#define __LU_CHECK_H__


#include <stdio.h>
#include "fptype.h"


#define LU_OCTAVE 1


#if LU_OCTAVE


static inline __attribute__((always_inline)) void print_matrix(int m, int mr, int n, int nr, fp_t *A)
{
  int j;
  for(j=0;j<n;j++)
  {
    int i;
    for(i=0;i<m;i++)
    {
      printf("%.4f ",A[j*mr+i]);
    }
    printf("\n");
  }
}


static inline __attribute__((always_inline)) void octave_matrix_def(char *file, char *name,int m, int mr, int n, int nr, fp_t *A)
{
  FILE *fstr=fopen(file,"a");

  fprintf(fstr,"%s = [];\n",name);

  int j;
  for(j=0;j<n;j++)
  {
    int i;
    for(i=0;i<m;i++)
    {
      fprintf(fstr,"%s(%i,%i) = %.6f;\n",name,i+1,j+1,A[j*mr+i]);
    }
  }

  fclose(fstr);
}


static inline __attribute__((always_inline)) void octave_fin(char *file, char *a,char *q)
{
  FILE *fstr=fopen(file,"a");

  fprintf(fstr,"\n[l,u,p] = lu(%s);\n",a);
  fprintf(fstr,"norm(triu(abs(%s)-abs(u)),1)/norm(u,1)\n",q);

  fclose(fstr);
}



#else


static inline __attribute__((always_inline)) void octave_matrix_def(char *file, char *name, int m, int mr, int n, int nr, fp_t *A) {}
static inline __attribute__((always_inline)) void octave_fin(char *file, char *a,char *q) {}


#endif


fp_t lu_check(int check, int m, int n, fp_t *A, fp_t *Aorig);


#endif // __LU_CHECK_H__
