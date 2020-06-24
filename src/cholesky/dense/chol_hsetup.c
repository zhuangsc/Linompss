#include "chol_setup.h"

#include <stdlib.h>
#include <time.h>

#include "fptype.h"
#include "genmat.h"
#include "densutil.h"


extern fp_t **Ah;
extern fp_t *Aorig;
extern fp_t *A;


int chol_setup(int check, int m, int mr, int ts, int bs, int tm, int mleft) { 
  	// Creation of hierarchical matrices. 
	Ah = (fp_t **) malloc( tm * tm * sizeof(fp_t *) );
	if ( Ah == NULL ) {
		return 1;
	}

	// Creation of data matrices
	A = calloc( mr * mr, sizeof(fp_t) );
	if ( A == NULL ) return 1;

	// Initialization of Ah
	int j;
	for ( j=0; j<tm; ++j ) {
    		int i;
		for ( i=0; i<tm; ++i ) {
			Ah[j*tm+i] = (fp_t *) &A[j*ts*mr+i*ts*ts];
		}
  	}

	GENMAT_HYPER_SPD(m, tm, mleft, ts, A);

	Aorig=NULL;
	if ( check ) {
		Aorig = CMH2CM( m, m, ts, ts, Ah[0], mr);
	}

	return 0;
}


void chol_shutdown() {
	free(A);
}
