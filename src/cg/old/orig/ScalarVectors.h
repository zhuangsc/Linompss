#ifndef ScalarVector

#define ScalarVector 1

#include <stdio.h>

extern int CreateInts (int **vint, int num);

extern int InitInts (int *vint, int n, int frst, int incr);

extern int CopyShiftInts (int *src, int *dst, int n, int shft);

extern void GetIntFromString (char *string, int *pnum, int numC, int shft);

extern void GetIntsFromString (char *string, int *vec, int numN, int numC, int shft);

extern void GetFormatsFromString (char *string, int *vec, int numN, int numC);

extern int PrintInts (int *vint, int num);

extern int RemoveInts (int **vint);

extern int CreateDoubles (double **vdouble, int num);

extern int InitDoubles (double *vdouble, int n, double frst, double incr);

extern void GetDoubleFromString (char *string, double *pdbl, int numC);

extern void GetDoublesFromString (char *string, double *vec, int numN, int numC);

extern int PrintDoubles (double *vdouble, int num);

extern int RemoveDoubles (double **vdouble);

#endif

