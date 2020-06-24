#ifndef _CS_H
#define _CS_H

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
    double *x ;	    /* numerical values, size nzmax */
    int nz ;	    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_add_double (const cs *A, const cs *B, double alpha, double beta) ;
int cs_cholsol_double (const cs *A, double *b, int order) ;
int cs_dupl_double (cs *A) ;
int cs_entry_double (cs *T, int i, int j, double x) ;
int cs_lusol_double (const cs *A, double *b, int order, double tol) ;
int cs_gaxpy_double (const cs *A, const double *x, double *y) ;
cs *cs_multiply_double (const cs *A, const cs *B) ;
int cs_qrsol_double (const cs *A, double *b, int order) ;
cs *cs_transpose_double (const cs *A, int values) ;
cs *cs_triplet_double (const cs *T) ;
double cs_norm_double (const cs *A) ;
int cs_print_double (const cs *A, int brief) ;
cs *cs_load_double (FILE *f) ;
/* utilities */
void *cs_calloc_double (int n, size_t size) ;
void *cs_free_double (void *p) ;
void *cs_realloc_double (void *p, int n, size_t size, int *ok) ;
cs *cs_spalloc_double (int m, int n, int nzmax, int values, int triplet) ;
cs *cs_spfree_double (cs *A) ;
int cs_sprealloc_double (cs *A, int nzmax) ;
void *cs_malloc_double (int n, size_t size) ;

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
    double *B ;	    /* beta [0..n-1] for QR */
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

int cs_cholsol2 (int n, css *S, csn *N, double *b, double *x);
int *cs_amd_double (const cs *A, int order) ;
csn *cs_chol_double (const cs *A, const css *S) ;
csd *cs_dmperm_double (const cs *A) ;
int cs_droptol_double (cs *A, double tol) ;
int cs_dropzeros_double (cs *A) ;
int cs_happly_double (const cs *V, int i, double beta, double *x) ;
int cs_ipvec_double (int n, const int *P, const double *b, double *x) ;
int cs_lsolve_double (const cs *L, double *x) ;
int cs_ltsolve_double (const cs *L, double *x) ;
csn *cs_lu_double (const cs *A, const css *S, double tol) ;
cs *cs_permute_double (const cs *A, const int *P, const int *Q, int values) ;
int *cs_pinv_double (const int *P, int n) ;
int cs_pvec_double (int n, const int *P, const double *b, double *x) ;
csn *cs_qr_double (const cs *A, const css *S) ;
css *cs_schol_double (const cs *A, int order) ;
css *cs_sqr_double (const cs *A, int order, int qr) ;
cs *cs_symperm_double (const cs *A, const int *Pinv, int values) ;
int cs_usolve_double (const cs *U, double *x) ;
int cs_utsolve_double (const cs *U, double *x) ;
int cs_updown_double (cs *L, int sigma, const cs *C, const int *parent) ;
/* utilities */
css *cs_sfree_double (css *S) ;
csn *cs_nfree_double (csn *N) ;
csd *cs_dfree_double (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
int *cs_counts_double (const cs *A, const int *parent, const int *post, int ata) ;
int cs_cumsum_double (int *p, int *c, int n) ;
int cs_dfs_double (int j, cs *L, int top, int *xi, int *pstack, const int *Pinv) ;
int *cs_etree_double (const cs *A, int ata) ;
int cs_fkeep_double (cs *A, int (*fkeep) (int, int, double, void *), void *other) ;
double cs_house_double (double *x, double *beta, int n) ;
int *cs_maxtrans_double (const cs *A) ;
int *cs_post_double (int n, const int *parent) ;
int cs_reach_double (cs *L, const cs *B, int k, int *xi, const int *Pinv) ;
csd *cs_scc_double (cs *A) ;
int cs_scatter_double (const cs *A, int j, double beta, int *w, double *x, int mark,
    cs *C, int nz) ;
int cs_splsolve_double (cs *L, const cs *B, int k, int *xi, double *x,
    const int *Pinv) ;
int cs_tdfs_double (int j, int k, int *head, const int *next, int *post,
    int *stack) ;
/* utilities */
csd *cs_dalloc_double (int m, int n) ;
cs *cs_done_double (cs *C, void *w, void *x, int ok) ;
int *cs_idone_double (int *p, cs *C, void *w, int ok) ;
csn *cs_ndone_double (csn *N, cs *C, void *w, void *x, int ok) ;
csd *cs_ddone_double (csd *D, cs *C, void *w, int ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)
#endif
