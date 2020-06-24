#ifndef __HBCONVRT_H__
#define __HBCONVRT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "hb.h"
#include "vector.h"

#define CSC2CSR		0
#define CSR2CSC		1

#ifdef SINGLE_PRECISION

#define HBH2HB					shbh2hb
#define	HB2HBH_block			shb2hbh_block
#define hb2hbh_symcsr			hb2hbh_ssymcsr
#define HB_SYM_EXPAND			shb_sym_expand

#define csc2csr					csc2csr_single
#define hb_csrcsc				hb_csrcsc_single
#define hbh_csrcsc				hbh_csrcsc_single
#define hb_transpose			hb_transpose_single
#define hbh_hyper_transpose		hbh_hyper_transpose_single
#define hb2hbb					hb2hbb_single
#define hbb2csrb				hbb2csrb_single
#define hb2hbh					hb2hbh_single
#define set_diag				set_diag_single

#else

#define HBH2HB					dhbh2hb
#define	HB2HBH_block			dhb2hbh_block
#define hb2hbh_symcsr			hb2hbh_dsymcsr
#define HB_SYM_EXPAND			dhb_sym_expand

#define csc2csr					csc2csr_double
#define hb_csrcsc				hb_csrcsc_double
#define hbh_csrcsc				hbh_csrcsc_double
#define hb_transpose			hb_transpose_double
#define hbh_hyper_transpose		hbh_hyper_transpose_double
#define hb2hbb					hb2hbb_double
#define hbb2csrb				hbb2csrb_double
#define hb2hbh					hb2hbh_double
#define set_diag				set_diag_double

#endif


#define LIBBBLAS_EXPORT __attribute__((__visibility__("default")))


extern LIBBBLAS_EXPORT hbmat_t* shbh2hb(hbmat_t *A);
extern LIBBBLAS_EXPORT void* shb2hbh_block(int I, int J, hbmat_t *A, int b, hbmat_t *block);
extern LIBBBLAS_EXPORT hbmat_t* hb2hbh_ssymcsr(hbmat_t *A, int b, int *etree, int alloc);
extern LIBBBLAS_EXPORT void shb_sym_expand(hbmat_t *A, hbmat_t *B);

hbmat_t* csc2csr_single(hbmat_t *A);
void hb_csrcsc_single(hbmat_t *Acsr, hbmat_t *Acsc, int csr2csc);
void hbh_csrcsc_single(hbmat_t *Acsr, hbmat_t *Acsc, int csr2csc);
hbmat_t* hb_transpose_single(hbmat_t *A);
hbmat_t* hbh_hyper_transpose_single(hbmat_t *A);
hbmat_t *hb2hbb_single(hbmat_t *A, int b);
hbmat_t* hbb2csrb_single(hbmat_t *A);
hbmat_t* hb2hbh_single(hbmat_t *A, int b, int is_csr);


extern LIBBBLAS_EXPORT hbmat_t* dhbh2hb(hbmat_t *A);
extern LIBBBLAS_EXPORT void* dhb2hbh_block(int I, int J, hbmat_t *A, int b, hbmat_t *block);
extern LIBBBLAS_EXPORT hbmat_t* hb2hbh_dsymcsr(hbmat_t *A, int b, int *etree, int alloc);
extern LIBBBLAS_EXPORT void dhb_sym_expand(hbmat_t *A, hbmat_t *B);

extern LIBBBLAS_EXPORT hbmat_t* csc2csr_double(hbmat_t *A);
extern LIBBBLAS_EXPORT void hbh_csrcsc_double(hbmat_t *Acsr, hbmat_t *Acsc, int csr2csc);
void hb_csrcsc_double(hbmat_t *Acsr, hbmat_t *Acsc, int csr2csc);
hbmat_t* hb_transpose_double(hbmat_t *A);
hbmat_t* hbh_hyper_transpose_double(hbmat_t *A);
hbmat_t *hb2hbb_double(hbmat_t *A, int b);
hbmat_t* hbb2csrb_double(hbmat_t *A);
hbmat_t* hb2hbh_double(hbmat_t *A, int b, int is_csr);

typedef struct _sparse_nodes {
	int row;
	int col;
	fp_t val;
	struct sparse_node *next;
	struct sparse_node *current;
} _sn_t;

#endif // __HBCONVRT_H__
