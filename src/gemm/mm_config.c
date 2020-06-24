#include "matmul_config.h"


#include <stdio.h>
#include <stdlib.h>
#include "fptype.h"
#include "fpmatr.h"


extern int m;
extern int k;
extern int n;
extern int b;
extern int c;
extern int d;
extern int lda;
extern int ldb;
extern int ldc;
extern int reps;
extern int check;
extern fp_t alpha;
extern fp_t beta;
extern ompssblas_t transa;
extern ompssblas_t transb;

int matmul_config(int argc, char *argv[]) 
{ 
	if (argc < 12 ) {
		fprintf(stderr, "use: %s m n k b d c lda ldb ldc transa transb [alpha] [beta] [reps] [check]\n", argv[0]);
		return 1;
	}

  	m = atoi(argv[1]);
  	n = atoi(argv[2]);
  	k = atoi(argv[3]);
  	b = atoi(argv[4]);
  	d = atoi(argv[5]);
  	c = atoi(argv[6]);
	lda = atoi(argv[7]);
	ldb = atoi(argv[8]);
	ldc = atoi(argv[9]);

	if ( b > m || c > k || d > n ) {
		fprintf(stderr, "error: block dimension exceeds matrix dimension\n");
		return 1;
	}

	int ta = atoi(argv[10]);
	int tb = atoi(argv[11]);
	if ( !ta )
		transa = OMPSSBLAS_NTRANSP;
	else
		transa = OMPSSBLAS_TRANSP;

	if ( !tb )
		transb = OMPSSBLAS_NTRANSP;
	else
		transb = OMPSSBLAS_TRANSP;

	alpha = 1.0;
	beta = 1.0;
	reps = 1;
	check = 0;

	if (argc > 12) {
		alpha = atof(argv[12]);
		if (argc > 13) {
			beta = atof(argv[13]);
			if (argc > 14) {
				reps = atoi(argv[14]);
				if (argc > 15) {
					check = atoi(argv[15]);
				}
			}
		}
	}

	return 0;
}
