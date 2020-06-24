#ifndef __HB_H__
#define __HB_H__


#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


#define MAT_CSC		0
#define MAT_CSR 	1
#define OFFS_F 		1
#define OFFS_C 		0

#define MAT_SYM		0x80000000
#define MAT_SPD		0x40000000
#define FRM_CSC		0x20000000
#define FRM_CSR		0x0

#define IS_SYM(f)	( f & MAT_SYM )
#define IS_SPD(f)	( f & MAT_SPD )
#define IS_CSC(f)	( (f & FRM_CSC) == FRM_CSC )
#define IS_CSR(f)	( (f & FRM_CSC) == FRM_CSR )



typedef struct strhbmat {
	int m, n;
	int elemc;
	int *vptr;
	int *vpos;
	void *vval;
	int *vdiag;
	int *udiagc;
	int b;
	int type;
	struct strhbmat *orig;
	struct strhbmat *trans;
	struct strhbmat *hyper;
	int orig_row;
	int orig_col;
	int *e_tree;
	int FACT;

	/*
	 * The following for hyper-matrix only
	 */
	int *vptr_pool;
	int *vpos_pool;
	void *vval_pool;
	int vptr_unit;
	int vpos_unit;
	int vval_unit;
	int vptr_pp;
	int vpos_pp;
	int vval_pp;
	pthread_mutex_t* mtx;

} hbmat_t;


static inline void hb_flag(hbmat_t *A, int flag) 
{ 
	A->type |= flag; 
}


static inline int hb_CSR(hbmat_t *A) 
{ 
	return ( (A->type & FRM_CSC) == FRM_CSR ); 
}

static inline int hb_CSC(hbmat_t *A) 
{ 
	return ( (A->type & FRM_CSC) == FRM_CSC ); 
}


#ifdef SINGLE_PRECISION

#define hb_print 				hb_print_single
#define hb_print_CSC			hb_print_CSC_single
#define hb_print_CSC2			hb_print_CSC2_single
#define hb_print_dense			hb_print_dense_single
#define hb_print_struc			hb_print_struc_single
#define hbb_print_dense			hb_hbb_print_dense_single
#define hb_diff					hb_diff_single
#define hb_sanity_check			hb_sanity_check_single
#define hb_init					hb_init_single
#define hb_init_basic			hb_init_basic_single
#define hbh_init				hbh_init_single
#define hb_reset				hb_reset_single
#define get_sdpos				get_sdpos_single
#define hb_cp					hb_cp_single
#define hb_energy_row			hb_energy_row_single
#define hb_sym_diag_block		hb_sym_diag_block_single

#else

#define hb_print 				hb_print_double
#define hb_print_CSC			hb_print_CSC_double
#define hb_print_CSC2			hb_print_CSC2_double
#define hb_print_dense			hb_print_dense_double
#define hb_print_struc			hb_print_struc_double
#define hbb_print_dense			hb_hbb_print_dense_double
#define hb_diff					hb_diff_double
#define hb_sanity_check			hb_sanity_check_double
#define hb_init					hb_init_double
#define hb_init_basic			hb_init_basic_double
#define hbh_init				hbh_init_double
#define hb_reset				hb_reset_double
#define get_sdpos				get_sdpos_double
#define hb_cp					hb_cp_double
#define hb_energy_row			hb_energy_row_double
#define hb_sym_diag_block		hb_sym_diag_block_double

#endif


#define LIBBBLAS_EXPORT __attribute__((__visibility__("default")))


#if 0
These print functions should be in matfprint, if you still need them
extern LIBBBLAS_EXPORT void hb_print_single(FILE *f, const char *name, hbmat_t *A, int full);
extern LIBBBLAS_EXPORT void hb_print_CSC_single(char *fname, hbmat_t *A);
extern LIBBBLAS_EXPORT void hb_print_CSC2_single(char *fname, hbmat_t *A);
extern LIBBBLAS_EXPORT void hb_print_dense_single( FILE* str, char * name, hbmat_t *A, int force );
extern LIBBBLAS_EXPORT void hb_print_struc_single(FILE* f, const char *name, hbmat_t *Ahb);
extern LIBBBLAS_EXPORT void hbb_print_dense_single( FILE* str, char * name, hbmat_t *A );
void hb_print_double(FILE *f, const char *name, hbmat_t *A, int full);
void hb_print_CSC_double(char *fname, hbmat_t *A);
void hb_print_CSC2_double(char *fname, hbmat_t *A);
void hb_print_dense_double( FILE* str, char * name, hbmat_t *A, int force );
void hb_print_struc_double(FILE* f, const char *name, hbmat_t *Ahb);
void hbb_print_dense_double( FILE* str, char * name, hbmat_t *A );
#endif
extern LIBBBLAS_EXPORT void hb_free(hbmat_t *A);
extern LIBBBLAS_EXPORT void hbh_free(hbmat_t *A);
extern LIBBBLAS_EXPORT void hbh_free2(hbmat_t *A);

extern LIBBBLAS_EXPORT int hb_diff_single(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hb_init_single(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hb_init_basic_single(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hbh_init_single(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hb_reset_single(hbmat_t *A);
extern LIBBBLAS_EXPORT void hb_sanity_check_single(const char *N, hbmat_t* A, int is_hbh);
extern LIBBBLAS_EXPORT int* get_sdpos_single(hbmat_t * A);
extern LIBBBLAS_EXPORT hbmat_t* hb_cp_single(hbmat_t *A);
extern LIBBBLAS_EXPORT float* hb_energy_row_single(hbmat_t *A, int bs);
extern LIBBBLAS_EXPORT void hb_sym_diag_block_single(hbmat_t *src_mat, int bsze, hbmat_t *diagb);

/*------------------------double------------------------------------*/
extern LIBBBLAS_EXPORT int hb_diff_double(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hb_init_double(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hb_init_basic_double(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hbh_init_double(hbmat_t *A, hbmat_t *B);
extern LIBBBLAS_EXPORT void hb_reset_double(hbmat_t *A);
extern LIBBBLAS_EXPORT void hb_sanity_check_double(const char *N, hbmat_t* A, int is_hbh);
extern LIBBBLAS_EXPORT int* get_sdpos_double(hbmat_t * A);
extern LIBBBLAS_EXPORT hbmat_t* hb_cp_double(hbmat_t *A);
extern LIBBBLAS_EXPORT double* hb_energy_row_double(hbmat_t *A, int bs);
extern LIBBBLAS_EXPORT void hb_sym_diag_block_double(hbmat_t *src_mat, int bsze, hbmat_t *diagb);


#endif // __HB_H__
