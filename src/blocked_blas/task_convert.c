#include "task_convert.h"
#include <math.h>
#include <float.h>

void task_d2s(int bm, int bn, int m, int n, double *D, float *S) 
{
	int j;
	int xcase=0;
	for ( j=0; j<bn; ++j) { 
		int i;
		for ( i=0; i<bm; ++i ) { 
			S[j*m+i] = (float) D[j*m+i];
			/* Handling the case of infinity */
//			 switch (isinff(S[j*m+i])) {
//				 case 1:
//					 S[j*m+i] = FLT_MAX;
//					 xcase=1;
//					 break;
//				 case -1:
//					 S[j*m+i] = -FLT_MAX;
//					 xcase=1;
//					 break;
//				 default:
//					 break;
//			 }
		}
	}
//	if (xcase)
//		printf("d2s inf\n");
}

void task_s2d(int bm, int bn, int m, int n, float *S, double *D) 
{
	int j;
	for ( j=0; j<bn; ++j ) { 
		int i;
		for ( i=0; i<bm; ++i ) { 
			D[j*m+i] = (double) S[j*m+i];
		}
	}
}
