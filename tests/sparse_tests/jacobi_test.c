#include <time.h>
#include <sys/time.h>
#include "mkl.h"
#include "bsblas_gemv_csr.h"
#include "fptype.h"
#include "hb.h"
#include "hbext.h"
#include "hbconvrt.h"
#include "iohb.h"
#include "genmat.h"
#include "symfac.h"
#include "vector.h"
#include "array.h"
#include "ompss_sparse_jacobi.h"


#ifdef SINGLE_PRECISION

#define fp				float
#define OMPSS_JACOBI	ompss_csr_sjacobi
#define MKL_CSRMV		mkl_scsrmv

#else

#define fp				double
#define OMPSS_JACOBI	ompss_csr_djacobi
#define MKL_CSRMV		mkl_dcsrmv

#endif


hbmat_t* Ahb;
hbmat_t* Ahbh;
fp *v_b; //RHS
fp *v_x; //Result with true value
fp *v_x0;//Result from Jacobi
int format; //CSC/CSR?
int bs; //Block size
int dim; //Dimension
int max_iter;
hbmat_t **diagL; //Factorized diagonal blocks
hbmat_t **A1;
int lookahead;
fp threshold;
char *fname;
int check;
int rep;
fp *work;

fp VECTOR_2NORM(fp *v, int length){
	fp x0 = 0;
	int i;
	for( i = 0; i < length; ++i ) {
		x0 += v[i] * v[i];
	}
	return FP_SQRT(x0);
}

int jacobi_config(int argc, char *argv[]) 
{
	if ( argc < 3 ) {
		printf("Usage %s HBfile b [csc(0)/csr(1)] [iter] [lookhead] [threshold] [check] [rep]\n", argv[0]);
		return 1;
	}

	fname = argv[1];
	bs = atoi(argv[2]);

	format = 1;
	max_iter = 1;
	lookahead = 0;
	threshold = 0;
	check = 0;
	rep = 1;

	if ( argc > 3 ) {
		format = atoi(argv[3]);

		if ( argc > 4 ) {
			max_iter = atoi(argv[4]);

			if ( argc > 5 ) {
				lookahead = atoi(argv[5]);

				if ( argc > 6 ) {
					threshold = atof(argv[6]);

					if ( argc > 7 ) {
						check = atoi(argv[7]);

						if ( argc > 8 ) {
							rep = atoi(argv[8]);
						}
					}
				}
			}			
		}
	}

	return 0;
}

int jacobi_setup (const char *fname) 
{
	Ahb = (hbmat_t*) malloc(sizeof(hbmat_t));
	hbmat_t *A = Ahb;

	readHB_newmat(fname, &(A->m), &(A->n), &(A->elemc), &(A->vptr), &(A->vpos), (fp **)&(A->vval));

	one2zero(A);

	dim = A->n ;
	v_b = malloc(dim * sizeof(fp));
	v_x0 = malloc(2 * dim* sizeof(fp));

	v_x = GENMAT(dim, 1);

	fp alpha = FP_ONE; fp beta = FP_NOUGHT;
	MKL_CSRMV("N", &A->m, &A->n, &alpha, "GLNC", A->vval, A->vpos, A->vptr, (A->vptr)+1, v_x, &beta, v_b);

	work = malloc((1+2*dim) * sizeof(fp));
	return 0;
}

fp jacobi_check(int check, int res_p)
{
	int res = 0;
	if ( check ) {
		fp x0;
		fp_t *v_res = &v_x0[res_p];
		fp *cnorm_x0 = malloc( dim * sizeof(fp));
		fp norm_b = VECTOR_2NORM(v_b, dim);
		int i;
		for ( i = 0; i < dim; ++i )
			cnorm_x0[i] = v_b[i];
		fp alpha = FP_MONE; fp beta = FP_ONE;
		struct timeval start, stop;
		unsigned long elapsed = 0;
		gettimeofday(&start, NULL);
		MKL_CSRMV("N", &(Ahb->m), &(Ahb->n), &alpha, "GLNC", Ahb->vval, Ahb->vpos, Ahb->vptr, (Ahb->vptr)+1, v_res, &beta, cnorm_x0);
		x0 = VECTOR_2NORM(cnorm_x0, dim);
		gettimeofday(&stop, NULL);
		elapsed = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
		printf("2-norm B-A*x0: %e in %lu us\n", x0/norm_b, elapsed);
		res = x0/norm_b;
	}
	return res;
}

int main(int argc, char* argv[])
{

	if ( jacobi_config(argc, argv) != 0 ) {
		return 1;
	}

	if ( jacobi_setup(fname) != 0 ) {
		return 2;
	}

	printf("setup\n");

	int res_p;
	OMPSS_JACOBI(bs, Ahb, v_x0, v_b, max_iter, lookahead, threshold, work, &res_p);

	printf("ompss_jacobi\n");
	jacobi_check(check, res_p);
	printf("check\n");

	return 0;

}
