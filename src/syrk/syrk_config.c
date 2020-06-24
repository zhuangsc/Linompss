#include "syrk_config.h"


#include <stdio.h>
#include <stdlib.h>
#include "fptype.h"
#include "fpmatr.h"

extern ompssblas_t uplo;
extern ompssblas_t trans;
extern fp_t alpha;
extern fp_t beta;
extern int n;
extern int k;
extern int b;
extern int d;
extern int lda;
extern int ldc;
extern int reps;
extern int check;


int syrk_config(int argc, char *argv[]) { 
	if (argc < 8 ) {
		fprintf(stderr, "use: %s uplo trans n k b lda ldc [alpha] [beta] [reps] [check]\n", argv[0]);
		return 1;
	}

	int u = atoi(argv[1]);
	int t = atoi(argv[2]);
	if ( !u ) {
		uplo = OMPSSBLAS_LOWERTRIANG;
		if ( !t ) {
			trans = OMPSSBLAS_NTRANSP;
		} else {
			trans = OMPSSBLAS_TRANSP;
		}
	} else {
		uplo = OMPSSBLAS_UPPERTRIANG;
		if ( !t ) {
			trans = OMPSSBLAS_NTRANSP;
		} else {
			trans = OMPSSBLAS_TRANSP;
		}
	}
		
  	n = atoi(argv[3]);
  	k = atoi(argv[4]);
  	b = atoi(argv[5]);
	lda = atoi(argv[6]);
	ldc = atoi(argv[7]);

	if ( b > n || b > k ) {
		fprintf(stderr, "use: block size exceeds matrix size\n");
		return 1;
	}

	alpha = 0.0;
	beta = 0.0;
	reps = 1;
	check = 0;
	if (argc > 7) {
		alpha = atof(argv[8]);
		if (argc > 8) {
			beta = atof(argv[9]);
			if (argc > 9) {
				reps = atoi(argv[10]);
				if (argc > 10) {
					check = atoi(argv[11]);
				}
			}
		}
	}

	return 0;
}
