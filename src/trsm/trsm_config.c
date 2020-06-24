#include "trsm_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "fptype.h"
#include "fpblas.h"

extern int m;
extern int n;
extern int b;
extern int lda;
extern int ldb;
extern fp_t alpha;
extern int reps;
extern int check;
extern ompssblas_t uplo;
extern ompssblas_t side;
extern ompssblas_t trans;
extern ompssblas_t diag;

int trsm_config(int argc, char *argv[]) 
{
	if (argc < 8 ) {
		printf("usage: %s side uplo m n b lda ldb [alpha] [trans] [diag] [reps] [check]\n",argv[0]);
		return 1;
	}

	if ( !atoi(argv[1]) )
		side = OMPSSBLAS_LEFT;
	else
		side = OMPSSBLAS_RIGHT;

	if ( !atoi(argv[2]) )
		uplo = OMPSSBLAS_LOWERTRIANG;
	else
		uplo = OMPSSBLAS_UPPERTRIANG;

	m = atoi(argv[3]);
	n = atoi(argv[4]);
	b = atoi(argv[5]);
	lda = atoi(argv[6]);
	ldb = atoi(argv[7]);

	alpha = 1.0;
	trans = OMPSSBLAS_NTRANSP;
	diag = OMPSSBLAS_NDIAGUNIT;
	reps = 1;
	check = 1;

	if (argc > 8 ) {
		alpha = atof(argv[8]);
		if ( argc > 9 ) {
			if (atoi(argv[9])) {
				trans = OMPSSBLAS_TRANSP;
			}
			if ( argc > 10) {
				if ( atoi(argv[10]) ) {
					diag = OMPSSBLAS_DIAGUNIT;
				}
				if ( argc > 11 ) {
					reps = atoi(argv[11]);
					if ( argc > 12 ) {
						check = atoi(argv[12]);
					}
				}
			}
		}
	}

	return 0;
}
