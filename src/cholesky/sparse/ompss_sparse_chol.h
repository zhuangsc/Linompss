#ifdef __cplusplus
extern "C" {
#endif


#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define LIBOMPSS_DLL_EXPORTED
#endif

#ifdef SINGLE_PRECISION

#define ompss_csr_chol_ll	ompss_csr_schol_ll
#define ompss_csc_chol_ll	ompss_csc_schol_ll

#endif

#ifdef DOUBLE_PRECISION

#define ompss_csr_chol_ll	ompss_csr_dchol_ll
#define ompss_csc_chol_ll	ompss_csc_dchol_ll

#endif

#include "hb.h"
#include "hbext.h"
#include "hbconvrt.h"

/* 	A: (in)		SPD matrix in HB, order n
	b: (in)		block size
	work: (in)	work array, size n
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csr_dchol_ll(int b, hbmat_t *A, int *work);  

extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csc_dchol_ll(int b, hbmat_t *A, int *work);

extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csr_schol_ll(int b, hbmat_t *A, int *work);  

extern LIBOMPSS_DLL_EXPORTED hbmat_t* ompss_csc_schol_ll(int b, hbmat_t *A, int *work);

#ifdef __cplusplus
}
#endif
