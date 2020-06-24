#ifndef _CS_H
#define _CS_H

#include "stdio.h"
#include "stdlib.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 1		    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"	    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    int nzmax ;	    /* maximum number of entries */
    int m ;	    /* number of rows */
    int n ;	    /* number of columns */
    int *p ;	    /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;	    /* row indices, size nzmax */
    float *x ;	    /* numerical values, size nzmax */
    int nz ;	    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_add_single (const cs *A, const cs *B, float alpha, float beta) ;
int cs_cholsol_single (const cs *A, float *b, int order) ;
int cs_dupl_single (cs *A) ;
int cs_entry_single (cs *T, int i, int j, float x) ;
int cs_lusol_single (const cs *A, float *b, int order, float tol) ;
int cs_gaxpy_single (const cs *A, const float *x, float *y) ;
cs *cs_multiply_single (const cs *A, const cs *B) ;
int cs_qrsol_single (const cs *A, float *b, int order) ;
cs *cs_transpose_single (const cs *A, int values) ;
cs *cs_triplet_single (const cs *T) ;
float cs_norm_single (const cs *A) ;
int cs_print_single (const cs *A, int brief) ;
cs *cs_load_single (FILE *f) ;
/* utilities */
void *cs_calloc_single (int n, size_t size) ;
void *cs_free_single (void *p) ;
void *cs_realloc_single (void *p, int n, size_t size, int *ok) ;
cs *cs_spalloc_single (int m, int n, int nzmax, int values, int triplet) ;
cs *cs_spfree_single (cs *A) ;
int cs_sprealloc_single (cs *A, int nzmax) ;
void *cs_malloc_single (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    int *Pinv ;	    /* inverse row perm. for QR, fill red. perm for Chol */
    int *Q ;	    /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;	    /* column pointers for Cholesky, row counts for QR */
    int m2 ;	    /* # of rows for QR, after adding fictitious rows */
    int lnz ;	    /* # entries in L for LU or Cholesky; in V for QR */
    int unz ;	    /* # entries in U for LU; in R for QR */
} css ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs *L ;	    /* L for LU and Cholesky, V for QR */
    cs *U ;	    /* U for LU, R for QR, not used for Cholesky */
    int *Pinv ;	    /* partial pivoting for LU */
    float *B ;	    /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    int *P ;	    /* size m, row permutation */
    int *Q ;	    /* size n, column permutation */
    int *R ;	    /* size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q) */
    int *S ;	    /* size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q) */
    int nb ;	    /* # of blocks in fine dmperm decomposition */
    int rr [5] ;    /* coarse row decomposition */
    int cc [5] ;    /* coarse column decomposition */
} csd ;

int cs_cholsol2 (int n, css *S, csn *N, float *b, float *x);
int *cs_amd_single (const cs *A, int order) ;
csn *cs_chol_single (const cs *A, const css *S) ;
csd *cs_dmperm_single (const cs *A) ;
int cs_droptol_single (cs *A, float tol) ;
int cs_dropzeros_single (cs *A) ;
int cs_happly_single (const cs *V, int i, float beta, float *x) ;
int cs_ipvec_single (int n, const int *P, const float *b, float *x) ;
int cs_lsolve_single (const cs *L, float *x) ;
int cs_ltsolve_single (const cs *L, float *x) ;
csn *cs_lu_single (const cs *A, const css *S, float tol) ;
cs *cs_permute_single (const cs *A, const int *P, const int *Q, int values) ;
int *cs_pinv_single (const int *P, int n) ;
int cs_pvec_single (int n, const int *P, const float *b, float *x) ;
csn *cs_qr_single (const cs *A, const css *S) ;
css *cs_schol_single (const cs *A, int order) ;
css *cs_sqr_single (const cs *A, int order, int qr) ;
cs *cs_symperm_single (const cs *A, const int *Pinv, int values) ;
int cs_usolve_single (const cs *U, float *x) ;
int cs_utsolve_single (const cs *U, float *x) ;
int cs_updown_single (cs *L, int sigma, const cs *C, const int *parent) ;
/* utilities */
css *cs_sfree_single (css *S) ;
csn *cs_nfree_single (csn *N) ;
csd *cs_dfree_single (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
int *cs_counts_single (const cs *A, const int *parent, const int *post, int ata) ;
int cs_cumsum_single (int *p, int *c, int n) ;
int cs_dfs_single (int j, cs *L, int top, int *xi, int *pstack, const int *Pinv) ;
int *cs_etree_single (const cs *A, int ata) ;
int cs_fkeep_single (cs *A, int (*fkeep) (int, int, float, void *), void *other) ;
float cs_house_single (float *x, float *beta, int n) ;
int *cs_maxtrans_single (const cs *A) ;
int *cs_post_single (int n, const int *parent) ;
int cs_reach_single (cs *L, const cs *B, int k, int *xi, const int *Pinv) ;
csd *cs_scc_single (cs *A) ;
int cs_scatter_single (const cs *A, int j, float beta, int *w, float *x, int mark,
    cs *C, int nz) ;
int cs_splsolve_single (cs *L, const cs *B, int k, int *xi, float *x,
    const int *Pinv) ;
int cs_tdfs_single (int j, int k, int *head, const int *next, int *post,
    int *stack) ;
/* utilities */
csd *cs_dalloc_single (int m, int n) ;
cs *cs_done_single (cs *C, void *w, void *x, int ok) ;
int *cs_idone_single (int *p, cs *C, void *w, int ok) ;
csn *cs_ndone_single (csn *N, cs *C, void *w, void *x, int ok) ;
csd *cs_ddone_single (csd *D, cs *C, void *w, int ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)
#endif
