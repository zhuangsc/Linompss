#include "densutil.h"

#include "fptype.h"
#include "fpblas.h"
#include "fplapack.h"

#ifdef SINGLE_PRECISION

#define __dmat_cp 			dmat_scp
#define __dmat_relerr 		dmat_srelerr
#define __vector_2norm		svector_2norm
#define __mat_normdiff		dmat_snormdiff
#define __mat_1norm			smat_1norm
#define __mat_energy_sym	smat_energy_sym
#define __vec_sanity_check	svec_sanity_check

#endif

#ifdef DOUBLE_PRECISION

#define __dmat_cp 			dmat_dcp
#define __dmat_relerr 		dmat_drelerr
#define __vector_2norm		dvector_2norm
#define __mat_normdiff		dmat_dnormdiff
#define __mat_1norm			dmat_1norm
#define __mat_energy_sym	dmat_energy_sym
#define __vec_sanity_check	dvec_sanity_check

#endif


fp_t * __dmat_cp(int m, int n, fp_t *A, int lda) 
{
	int mn = m * n;
	fp_t *C = calloc(mn, sizeof(fp_t));

	int j;
	for ( j=0; j<n; ++j ) {
		int i;
		for ( i=0; i<m; ++i ) {
			C[j*lda+i] = A[j*lda+i];
		}
	}

	return C;
}


// norm(tril(abs(A)-abs(B)),1) / norm(A,1);
fp_t __dmat_relerr(ompssblas_t trian, int m, int n, fp_t *A, fp_t* B) 
{
	int m2 = m * n;
	fp_t *D  = calloc(m2, sizeof(fp_t));
	if ( D == NULL) { 
		fprintf(stderr, "relerr: failed to allocate D\n");
		return 1.0;
	}

	fp_t *sumA = calloc(n, sizeof(fp_t));
  	if (sumA == NULL) { 
		fprintf(stderr, "relerr: failed to allocate sumA\n");
		return 2.0;
	}

	fp_t *sumD = calloc(n, sizeof(fp_t));
  	if ( sumD == NULL ) { 
		fprintf(stderr, "relerr: failed to allocate sumD\n");
		return 3.0;
	}

  	int j;
	for ( j=0; j<n; ++j ) {
		int start = TEST_LOWER(trian)? j : 0;
		start = TEST_UPPER(trian)? 0 : start;
		int stop = TEST_UPPER(trian)? j : m;

		int i;
		for ( i=start; i<stop; ++i ) {
			fp_t tmp = FP_ABS(A[j*m+i]);
			sumA[j] += tmp;
			tmp -= FP_ABS(B[j*m+i]);
			D[j*m+i] = tmp;
			sumD[j] += FP_ABS(tmp);
		}
  	}

	fp_t norm1 = 0.0, norm1M1 = 0.0;
	for ( j=0; j<n; ++j ) {
		norm1= max(norm1, sumD[j]);
		norm1M1 = max(norm1M1, sumA[j]);
	}
	norm1 /= norm1M1;

	free(sumA);
	free(sumD);
	free(D);

	return(norm1);
}


fp_t __vector_2norm(fp_t *v, int length)
{
	fp_t x0 = 0;
	int i;
	for( i = 0; i < length; ++i ) {
		x0 += v[i] * v[i];
	}
	return FP_SQRT(x0);
}

fp_t __mat_1norm(fp_t *mat, int m, int n)
{

	int nan = 0;
	fp_t maximum = 0.0;
	int j;
	for ( j = 0; j < n; ++j ) {
		fp_t acc = 0.0;
		int i;
		for ( i = 0; i < m; ++i) {
			if (fpclassify(mat[j*m+i]) == FP_NAN)
				nan++;
			acc += FP_ABS(mat[j*m+i]);
		}
		maximum = acc>maximum ? acc : maximum;
	}
	if ( nan )
		fprintf(stderr, "%d NaNs found\n", nan);
	return maximum;

}

fp_t __mat_normdiff(char norm, int m, int n, fp_t *A, int lda, fp_t *A0, int lda0) 
{

	fp_t errA, err0;
	errA = LAPACK_lange(norm, m, n, A, lda);
	err0 = LAPACK_lange(norm, m, n, A0, lda0);
	return (FP_ABS(errA-err0)/err0);

}


fp_t * CMH2CM(int m, int n, int tr, int tc, fp_t *Ah, int ldAh) 
{
    fp_t * A = (fp_t *) malloc( m * n * sizeof(fp_t) );
    if ( A == NULL ) {
        return NULL;
    }

    int i;
    for( i = 0; i < m; ++i ) {
        int j;
        for( j = 0; j < n; ++j ) {
            A[j*m+i] = Ah[ (j/tc)*ldAh*tc + (j%tc)*tr + (i/tr)*tr*tc + (i%tr) ] ;
        }
    }

    return A;
}

fp_t __mat_energy_sym(int n, int lda, int bs, fp_t *mat, fp_t *dist)
{
	int c = 0;
	fp_t max = 0.0;
	fp_t tol = 0.0;
	int i;
	for ( i = 0; i < n; i += bs ) {
		fp_t energy = 0.0;
		fp_t *mptr = &mat[i];
		int rows = n-i < bs ? n-i : bs;
		int cols = n;
		int j;
		for ( j = 0; j < cols; j++ ) {
			int k;
			for ( k = 0; k < rows; k++ ){
				energy += mptr[k*lda+j];
			}
		}
		dist[c++] = energy;
		tol += energy;
		max = isless(max, energy) ? energy : max;
	}

	fp_t avg = tol/c;
	fp_t var = 0.0;
	for ( i = 0; i < c; i++ ) {
		var += FP_POW((dist[i]-avg), 2.0);
	}
	var /= c;
	fprintf(stdout,"Tol: %e, avg: %e, var: %e\n", tol, avg, var);
	return max;
}

void __vec_sanity_check(char *n, fp_t *vec, int len)
{
	int nans; int infs; int zeros; int subs; int norms;
	infs = nans = zeros = subs = norms = 0;
	int extreme=0;

	int i;
	for ( i = 0; i < len; i++ ) {
		if ( isgreater(vec[i],1E37) )
			extreme++;
		int state = fpclassify(vec[i]);
		switch (state) {
			case FP_NAN:
				nans++;
				break;
			case FP_INFINITE:
				infs++;
				break;
			case FP_ZERO:
				zeros++;
				break;
			case FP_SUBNORMAL:
				subs++;
				break;
			case FP_NORMAL:
				norms++;
				break;
			default:
				break;
		}
	}
	fprintf(stdout,"Vec: %s, Total: %d, NaNs: %d, INFs: %d, Zeros: %d, Subs: %d, Normals: %d, X: %d\n", n, len, nans, infs, zeros, subs, norms, extreme);
}
