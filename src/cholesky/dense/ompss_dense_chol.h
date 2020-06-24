#ifdef __cplusplus
extern "C" {
#endif


#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define LIBOMPSS_DLL_EXPORTED
#endif


/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_schol_hll(int n, int b, int t, float **Ah);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_hll(int n, int b, int t, double **Ah);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_schol_hrl(int n, int b, int t, float **Ah);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_hrl(int n, int b, int t, double **Ah);


/* 	A: (in)		SPD matrix 
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_schol_ll(int n, int b, int t, float *A, int lda);

extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_ll(int n, int b, int t, double *A, int lda);

extern LIBOMPSS_DLL_EXPORTED int ompss_schol_rl(int n, int b, int t, float *A, int lda);

extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_rl(int n, int b, int t, double *A, int lda);


extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_nhll(int mt, int b, int t, double **Ah);

extern LIBOMPSS_DLL_EXPORTED int ompss_schol_nhll(int mt, int b, int t, float **Ah);

extern LIBOMPSS_DLL_EXPORTED int ompss_dchol_nhrl(int mt, int b, int t, double **Ah);

extern LIBOMPSS_DLL_EXPORTED int ompss_schol_nhrl(int mt, int b, int t, float **Ah);


#ifdef __cplusplus
}
#endif
