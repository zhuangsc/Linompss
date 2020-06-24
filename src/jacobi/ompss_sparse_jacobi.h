#ifdef __cplusplus
extern "C" {
#endif


#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define LIBOMPSS_DLL_EXPORTED
#endif


#include "hb.h"

#ifdef SINGLE_PRECISION

#define ompss_csr_jacobi	ompss_csr_sjacobi

#endif

#ifdef DOUBLE_PRECISION

#define ompss_csr_jacobi	ompss_csr_djacobi

#endif

/* 	A: (in)		SPD matrix in HB, order n
	b: (in)		block size
	work: (in)	work array, size n
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED void ompss_csr_sjacobi(int b, hbmat_t *A, float *v_x, float *v_b, int max_iter,\
		int lookahead, float threshold, float *work, int *res_p);

extern LIBOMPSS_DLL_EXPORTED void ompss_csr_djacobi(int b, hbmat_t *A, double *v_x, double *v_b, int max_iter,\
		int lookahead, double threshold, double *work, int *res_p);

#ifdef __cplusplus
}
#endif
