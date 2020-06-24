#ifndef __BBLAS_GEMM_H__
#define __BBLAS_GEMM_H__


#include "fptype.h"
#include "as_man.h"
#include "fpmatr.h"
#include "async_struct.h"

#include "task_gemm.h"
//#include "selfsched.h"



#ifdef SINGLE_PRECISION

#define BBLAS_GEMM 				bblas_sgemm
#define BBLAS_HYPER_GEMM 		bblas_hyper_sgemm
#define BBLAS_GEMM_PRIOR 		bblas_sgemm_prior
#define BBLAS_GEMM_PRIOR_RELEASE 		bblas_sgemm_prior_release
#define BBLAS_GEMM_PRED 		bblas_sgemm_pred
#define BBLAS_GEMM_SCHED		bblas_sgemm_sched
#define BBLAS_GEMM_PRED_SCHED	bblas_sgemm_pred_sched
#define BBLAS_GEMM_PROFPRED		bblas_sgemm_profpred

#else

#define BBLAS_GEMM 				bblas_dgemm
#define BBLAS_HYPER_GEMM 		bblas_hyper_dgemm
#define BBLAS_GEMM_PRIOR 		bblas_dgemm_prior
#define BBLAS_GEMM_PRIOR_RELEASE 		bblas_dgemm_prior_release
#define BBLAS_GEMM_PRED 		bblas_dgemm_pred
#define BBLAS_GEMM_SCHED		bblas_dgemm_sched
#define BBLAS_GEMM_PRED_SCHED	bblas_dgemm_pred_sched
#define BBLAS_GEMM_PROFPRED		bblas_dgemm_profpred

#endif

static inline void __attribute__((always_inline)) bblas_sgemm(int p, ompssblas_t transa, ompssblas_t transb, int bm, int bn, \
		int bk, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC) {
	if ( alpha == 0.0 )
		return;
	if ( TEST_NTRANSP(transb) ) {
		if ( TEST_NTRANSP(transa) ) {
			/* C = alpha * A * B + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;
					int boffs = j * ldimB;

					task_sgemm(transa, transb, c, d, bk, alpha, &A[0*ldimA+i], ldimA, &B[boffs], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_sgemm(transa, transb, c, d, f, alpha, &A[l*ldimA+i], ldimA, &B[boffs+l], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}
		} else {
			/* C = alpha * A^T * B + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;
					int aoffs = i * ldimA;

					task_sgemm(transa, transb, c, d, bk, alpha, &A[aoffs], ldimA, &B[j*ldimB+0], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_sgemm(transa, transb, c, d, f, alpha, &A[aoffs+l], ldimA, &B[j*ldimB+l], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}
		}
	} else {
		if ( TEST_NTRANSP(transa) ) {
			/* C = alpha * A * B^T + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;

					task_sgemm(transa, transb, c, d, bk, alpha, &A[0*ldimA+i], ldimA, &B[0*ldimB+j], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_sgemm(transa, transb, c, d, f, alpha, &A[l*ldimA+i], ldimA, &B[l*ldimB+j], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}

		} else {
			/* C = alpha * A^T * B^T + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;
					int aoffs = i * ldimA;

					task_sgemm(transa, transb, c, d, bk, alpha, &A[aoffs], ldimA, &B[0*ldimB+j], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_sgemm(transa, transb, c, d, f, alpha, &A[aoffs+l], ldimA, &B[l*ldimB+j], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}
		}
	}
}


static inline void __attribute__((always_inline)) bblas_hyper_sgemm(ompssblas_t transa, ompssblas_t transb, int bm, int bn, int bk, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, 
	float beta, float *C, int ldimC) {
	int i;
	for ( i=0; i<m; i+= bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_sgemm(transa, transb, bm, bn, bk, alpha, &A[0*ldimA+i], bm, &B[j*ldimB+0], bk, beta, &C[j*ldimC+i], bm, 1);
			int l;
			for ( l=bk; l<k; l+=bk) {
				task_sgemm(transa, transb, bm, bn, bk, alpha, &A[l*ldimA+i], bm, &B[j*ldimB+l], bk, s_one, &C[j*ldimC+i], bm, 1);
			}
		}
	}
}


static inline void __attribute__((always_inline)) bblas_dgemm(int p, ompssblas_t transa, ompssblas_t transb, int bm, int bn, \
		int bk, int m, int n, int k, double alpha, double *A, int ldimA, double *B, int ldimB, double beta, double *C, int ldimC) {
	if ( alpha == 0.0 )
		return;
	if ( TEST_NTRANSP(transb) ) {
		if ( TEST_NTRANSP(transa) ) {
			/* C = alpha * A * B + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;
					int boffs = j * ldimB;

					task_dgemm(transa, transb, c, d, bk, alpha, &A[0*ldimA+i], ldimA, &B[boffs], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_dgemm(transa, transb, c, d, f, alpha, &A[l*ldimA+i], ldimA, &B[boffs+l], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}
		} else {
			/* C = alpha * A^T * B + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;
					int aoffs = i * ldimA;

					task_dgemm(transa, transb, c, d, bk, alpha, &A[aoffs], ldimA, &B[j*ldimB+0], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_dgemm(transa, transb, c, d, f, alpha, &A[aoffs+l], ldimA, &B[j*ldimB+l], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}
		}
	} else {
		if ( TEST_NTRANSP(transa) ) {
			/* C = alpha * A * B^T + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;

					task_dgemm(transa, transb, c, d, bk, alpha, &A[0*ldimA+i], ldimA, &B[0*ldimB+j], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_dgemm(transa, transb, c, d, f, alpha, &A[l*ldimA+i], ldimA, &B[l*ldimB+j], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}

		} else {
			/* C = alpha * A^T * B^T + beta * C */
			int i;
			for ( i=0; i<m; i+= bm ) {
				int cs = m - i;
				int c = cs < bm ? cs : bm;

				int j;
				for ( j=0; j<n; j+=bn ) {
					int ds = n - j;
					int d = ds < bn ? ds : bn;
					int coffs = j * ldimC + i;
					int aoffs = i * ldimA;

					task_dgemm(transa, transb, c, d, bk, alpha, &A[aoffs], ldimA, &B[0*ldimB+j], ldimB, beta, &C[coffs], ldimC, p);

					int l;
					for ( l=bk; l<k; l+=bk) {
						int fs = k - l;
						int f = fs < bk ? fs : bk;

						task_dgemm(transa, transb, c, d, f, alpha, &A[aoffs+l], ldimA, &B[l*ldimB+j], ldimB, FP_ONE, &C[coffs], ldimC, p);
					}
				}
			}
		}
	}
}


static inline void __attribute__((always_inline)) bblas_hyper_dgemm(ompssblas_t transa, ompssblas_t transb, int bm, int bn, int bk, int m, int n, int k, 
	double alpha, double *A, int ldimA, double *B, int ldimB, double beta, double *C, int ldimC) {
	int i;
	for ( i=0; i<m; i+= bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_dgemm(transa, transb, bm, bn, bk, alpha, &A[0*ldimA+i*bk], bm, &B[j*ldimB+0], bk, beta, &C[j*ldimC+i*bn], bm, 1);
			int l;
			for ( l=bk; l<k; l+=bk) {
				task_dgemm(transa, transb, bm, bn, bk, alpha, &A[l*ldimA+i*bk], bm, &B[j*ldimB+l*bn], bk, d_one, &C[j*ldimC+i*bn], bm, 1);
			}
		}
	}
}


static inline async_stat_t __attribute__((always_inline)) bblas_sgemm_pred(int bm, int bn, int bk, int m, int n, int k, float alpha, float *A, float *B, float beta, float *C, 
	int it, int est, asman_t *asman, async_t *cond, float **alpha1, double tol) {
	async_stat_t status = STAT_AHEAD;

	int i;
	for ( i=0; i<m; i+= bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			if ( status == STAT_AHEAD ) {
//				status = asman_sbreak(it, est, asman, alpha1, cond, 0, tol);
				printf("bblas_sgemm_pred error!\n");
				if ( status == STAT_BROKEN || status == STAT_CONVERGED ) {
					printf("broken in sGEMM_PRED\n");
					return status;
				}
			}

			task_sgemm(OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, 9999);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_sgemm(OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, FP_ONE, &C[j*m+i], m, 9999);
			}
		}
	}

	return status;
}


static inline async_stat_t __attribute__((always_inline)) bblas_dgemm_pred(int bm, int bn, int bk, int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, 
	int it, int est, asman_t *asman, async_t *cond, double **alpha1, double tol) {
	async_stat_t status = STAT_AHEAD;

	int i;
	for ( i=0; i<m; i+= bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			if ( status == STAT_AHEAD ) {
//				status = asman_dbreak(it, est, asman, alpha1, cond, 0, tol, 0x0);
				printf("bblas_dgemm_pred error!\n");
				if ( status == STAT_BROKEN || status == STAT_CONVERGED ) {
					return status;
				}
			}

			task_dgemm(OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, 9999);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_dgemm(OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, FP_ONE, &C[j*m+i], m, 9999);
			}
		}
	}

	return status;
}


static inline async_stat_t __attribute__((always_inline)) bblas_sgemm_profpred(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, float alpha, float *A, float *B, float beta, float *C, int previt, int est, asman_t *asman, async_t *cond, float **alpha1, double tol) 
{
	async_stat_t status = STAT_AHEAD;

	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			if ( status == STAT_AHEAD ) {
//				status = asman_sbreak(previt, est, asman, alpha1, cond, 0, tol);
				printf("bblas_sgemm_pred error!\n");
				if ( status == STAT_BROKEN || status == STAT_CONVERGED ) {
					printf("broken in sGEMM_PRED_SCHED\n");
					return status;
				}
			}

			int childp = p + prof_priority(&sync->prof, idx);
//			printf("childp: %d, idx: %d\n", childp, idx);

			task_sgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_sgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1);
			}
		}
	}

	return status;
}

static inline async_stat_t __attribute__((always_inline)) bblas_dgemm_profpred(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, int previt, int est, asman_t *asman, async_t *cond, double **alpha1, double tol) 
{
	async_stat_t status = STAT_AHEAD;

	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			if ( status == STAT_AHEAD ) {
//				status = asman_dbreak(previt, est, asman, alpha1, cond, 0, tol);
				printf("bblas_dgemm_pred error!\n");
				if ( status == STAT_BROKEN || status == STAT_CONVERGED ) {
					printf("broken in sGEMM_PRED_SCHED\n");
					return status;
				}
			}

			int childp = p + prof_priority(&sync->prof, idx);
//			printf("childp: %d, idx: %d\n", childp, idx);

			task_dgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_dgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1);
			}
		}
	}

	return status;
}


static inline void __attribute__((always_inline)) bblas_sgemm_prior(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, float alpha, float *A, float *B, float beta, float *C, float *mat_energy) 
{
	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int resolution = sync->prof.resolution;
			int bs = (m+bm-1)/bm;
			float M = mat_energy[bs];
			float scal = (1 / M * resolution);
			float len = fabs(mat_energy[i] * scal);
			int childp = p + (int) len;
//			int childp = p + prof_priority(&sync->prof, idx);
//			printf("p: %d, childp: %d, idx: %d\n", p, childp, idx);

			task_sgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_sgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1);
			}
		}
	}
}

static inline void __attribute__((always_inline)) bblas_dgemm_prior(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, float *mat_energy) 
{
	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int resolution = sync->prof.resolution;
			int bs = (m+bm-1)/bm;
			float M = mat_energy[bs];
			float scal = (1 / M * resolution);
			float len = fabs(mat_energy[i] * scal);
			int childp = p + (int) len;

			task_dgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_dgemmcg(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1);
			}
		}
	}
}


static inline void __attribute__((always_inline)) bblas_sgemm_prior_release(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, float alpha, float *A, float *B, float beta, float *C, int release, float *mat_energy, int *bitmap) 
{
	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int resolution = sync->prof.resolution;
			int bs = (m+bm-1)/bm;
			float M = mat_energy[bs];
			float scal = (1 / M * resolution);
			float len = fabs(mat_energy[idx] * scal);
			int childp = p + (int) len;

			task_sgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx, release, &bitmap[idx]);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_sgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, FP_ONE, &C[j*m+i], m, childp, -1, release, &bitmap[idx]);
			}
		}
	}
}

static inline void __attribute__((always_inline)) bblas_dgemm_prior_release(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, int release, float *mat_energy, int *bitmap) 
{
	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int resolution = sync->prof.resolution;
			int bs = (m+bm-1)/bm;
			float M = mat_energy[bs];
			float scal = (1 / M * resolution);
			float len = fabs(mat_energy[i] * scal);
			int childp = p + (int) len;
//			int childp = p + prof_priority(&sync->prof, idx);

			task_dgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx, release, &bitmap[idx]);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_dgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1, release, &bitmap[idx]);
			}
		}
	}
}


static inline void __attribute__((always_inline)) bblas_sgemm_prior_switch(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, float alpha, float *A, float *B, float beta, float *C, int release_id, float *mat_energy, int *bitmap) 
{
	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int resolution = sync->prof.resolution;
			int bs = (m+bm-1)/bm;
			float M = mat_energy[bs];
			float scal = (1 / M * resolution);
			float len = fabs(mat_energy[idx] * scal);
			int childp = p + (int) len;

			int on = ( idx == release_id ) ? 0 : 1;
			task_sgemmcg_switch(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx, release_id, &bitmap[idx], on);
//			task_sgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx, release, &bitmap[idx]);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

			task_sgemmcg_switch(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, -1, release_id, &bitmap[idx], on);
//				task_sgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, FP_ONE, &C[j*m+i], m, childp, -1, release, &bitmap[idx]);
			}
		}
	}
}

static inline void __attribute__((always_inline)) bblas_dgemm_prior_switch(int it, async_t *sync, async_t *sync2, int p, int bm, int bn, int bk, int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, int release_id, float *mat_energy, int *bitmap) 
{
	int idx;
	int i;
	for ( i=0, idx=0; i<m; i+= bm, ++idx ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int resolution = sync->prof.resolution;
			int bs = (m+bm-1)/bm;
			float M = mat_energy[bs];
			float scal = (1 / M * resolution);
			float len = fabs(mat_energy[i] * scal);
			int childp = p + (int) len;
//			int childp = p + prof_priority(&sync->prof, idx);

			int on = ( idx == release_id ) ? 0 : 1;
			task_dgemmcg_switch(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx, release_id, &bitmap[idx], on);
//			task_dgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, bk, alpha, &A[0*m+i], m, &B[j*k+0], k, beta, &C[j*m+i], m, childp, idx, release, &bitmap[idx]);

			int l;
			for ( l=bk; l<k; l+=bk) {
				int fs = k - l;
				int f = fs < bk ? fs : bk;

				task_dgemmcg_switch(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1, release_id, &bitmap[idx], on);
//				task_dgemmcg_release(it, sync2, OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, c, d, f, alpha, &A[l*m+i], m, &B[j*k+l], k, 1.0, &C[j*m+i], m, childp, -1, release, &bitmap[idx]);
			}
		}
	}
}


#endif // __BBLAS_GEMM_H__
