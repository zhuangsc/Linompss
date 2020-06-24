#include "chol_check.h"

#include <math.h>
#include <stdlib.h>

#include "fptype.h"
#include "fplapack.h"
#include "densutil.h"



// Perform equivalent to Octave's command:
// norm(triu(abs(M1)-abs(M2)),1)/norm(M1,1);
double __check ( double * M1, double * M2, int m) {
	double *diffmat  = (double *) calloc( (m*m), sizeof(double) );
	if (diffmat == NULL) { 
		perror("Error Allocating diffmat"); 
		exit(1); 
	}
	double *sumvec   = (double *) calloc( (m),   sizeof(double) );
  	if (sumvec == NULL) { 
		perror("Error Allocating sumvec"); 
		exit(1); 
	}
	double *sumvecM1 = (double *) calloc( (m),   sizeof(double) );
  	if (sumvecM1 == NULL) { 
		perror("Error Allocating sumvecM1"); 
		exit(1); 
	}


  	int j;
	for (j=0; j<m; j++) {
		int i;
		for (i=j; i<m; i++) {
      			double tmp = fabs(M1[j*m+i]);
      			sumvecM1[j] += tmp;
      			tmp -= fabs(M2[j*m+i]);
			//printf("dif (%i,%i) %x\n", i, j, &diffmat[j*m+i]);
      			diffmat[j*m+i] = tmp;
	      		sumvec[j] += fabs(tmp);
    		}
  	}

	double norm1=0.0, norm1M1=0.0;
	for (j=0; j<m; j++) {
		norm1= max(norm1,sumvec[j]);
		norm1M1 = max(norm1M1,sumvecM1[j]);
	}
	norm1 /= norm1M1;

	free(sumvec);
	free(sumvecM1);
	free(diffmat);

	return(norm1);
}


int chol_check(int check, int m, int mr, int ts, fp_t **Ah, fp_t *A, fp_t *Rcheck) {
	if( check > 0 ) {
#if USE_PRL
    		fp_t *R = CMH2CM(m, m, ts, ts, Ah[ 0 ], mr);
#else
    		fp_t *R = CMH2CM(m, m, ts, ts, Ah[ 0 ], mr);
#endif
		printf( "calling dpotrf for residual checking...");
		fflush(0);

		int IINFO;
		dpotrf_("Lower", &m, Rcheck, &m, &IINFO);
		printf( "...done (IINFO: %d)\n", IINFO);
		if(IINFO < 0) {
			perror("Error in dpotrf\n");
		}

		fp_t norm_diff = __check( Rcheck, R, m);
		printf( "norm of difference: %22.16g\n", norm_diff);

		free(Rcheck);

		int ret=(norm_diff > 0.0001)?1:0;

		return ret;
	}

	return 0;
}
