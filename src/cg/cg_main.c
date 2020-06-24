#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "cg_main.h"

void *Ahbh;
void *Ahb;
void *Acsr;
int n;
int bm;
int cgit;
double prec;
int correction;
//int is_precond;
int rep;
char *aname;
char *rhsfname;
double orth_fac;
int cglog_level;
int cg_ver;

double *rhs;
double *x;

void *preconditioner;
css **S;
csn **N;
cs *Acs;

int main(int argc, char *argv[])
{
	if ( cg_config(argc, argv) ) {
		return 1;
	}

	if ( cg_setup(n, bm, &x, &rhs) ) {
		return 2;
	}

	for ( int i = 0; i < rep; i++ ) {
		switch ( cg_ver ) {
			case 0: 
				printf("ALG1 PCG\n");
				CG_ALG1(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 1 original PCG
				break;
			case 1:
				printf("ALG3 Chronopoulos\n");
				CG_ALG3(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 3 Chronopoulos
				break;
			case 2:
				printf("ALG4 Pipelined\n");
				CG_ALG4(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 4 Pipelined
				break;
			case 3:
				printf("ALG7 Gropp\n");
				CG_ALG7(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 7 Gropp
				break;
			default:
				printf("No algorithm selected\n");
				break;
		}
	}

	return 0;
}


/* Standard PCG */
int CG_ALG1(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(4 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * p_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	double  *s = malloc(n * sizeof(double)); // s = Ap

	double *alpha1 = calloc(2, sizeof(double));
	double *alpha2 = calloc(2, sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
	/* u[0] = M^-1 * r[0] */
	BSBLAS_CHOLSOLV2(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = z[0] */
	BBLAS_COPY(1, bm, 1, n, 1, u[i], p[i]);
	/* alpha1[0] = <r[0], z[0]> */
	BBLAS_DOT_PURE(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* s = A * p[i] */
		HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, p[iprev], FP_NOUGHT, s); 
		/* alpha2[i] = <s, p[i]> */
		BBLAS_DOT_PURE(1, bm, 1, n, 1, s, p[iprev], &alpha2[i]);
		BBLAS_CPAXPY_COMB(bm, 1, n, 1, FP_MONE, &alpha1[iprev], &alpha2[i], s, p[iprev], r[iprev], x[iprev], r[i], x[i]);

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
			HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
		BSBLAS_CHOLSOLV2(1, bm, n, S, N, r[i], u[i]);
		/* alpha1[i+1] = <r, u> */
		BBLAS_DOT_PURE(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);
		/* p[i+1] = u[i+1] + transpose(beta[i]) * p[i] */
		BBLAS_EXTM_AXPY(1, bm, 1, n, 1, &alpha1[i], &alpha1[iprev], p[iprev], u[i], p[i]); 	

		#pragma omp taskwait

		stop_timer(&elapses[k]);
		alpha1[iprev] = alpha2[iprev] = (double) 0;
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}
	memcpy(solution, x[i], n * sizeof(double));
	dump_info("cg_alg1.log", k, residuals, elapses);
	free(pool);
	free(alpha1);
	free(alpha2);
	free(residuals);
	free(elapses);
	return 0;
}

/* Chronopoulos PCG */
int CG_ALG3(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(6 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * s_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w = A * u_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i+1 = w_i+1 + beta * s_i

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;

	BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
	/* u[0] = M^-1 * r[0] */
	BSBLAS_CHOLSOLV2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, u[i], FP_NOUGHT, w[i]);
	/* gamma[0] = <r[0], u[0]> */
	BBLAS_DOT_PURE(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
	/* alpha[0] = gamma[0]/<w[0], u[0]> */
	BBLAS_DOT_PURE(1, bm, 1, n, 1, w[i], u[i], &alpha[i]);

	#pragma omp taskwait on (gamma[i], alpha[i])

	alpha[i] = gamma[i] / alpha[i];

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *wp0 = &(w[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *up1 = &(u[i])[j];
			double *sp1 = &(s[i])[j];
			double *wp1 = &(w[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];

			#pragma omp task out([c]pp1, [c]sp1, [c]xp1, [c]rp1)  priority(1) label(alg3_fuse)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				cblas_dcopy(c, up0, 1, pp1, 1);
				cblas_daxpy(c, beta[iprev], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				cblas_dcopy(c, wp0, 1, sp1, 1);
				cblas_daxpy(c, beta[iprev], sp0, 1, sp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				cblas_dcopy(c, xp0, 1, xp1, 1);
				cblas_daxpy(c, alpha[iprev], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				cblas_dcopy(c, rp0, 1, rp1, 1);
				cblas_daxpy(c, -1*alpha[iprev], sp1, 1, rp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
			HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
		BSBLAS_CHOLSOLV2(1, bm, n, S, N, r[i], u[i]);
		/* w[i+1] = A * u[i+1] */
		HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, u[i], FP_NOUGHT, w[i]); 

		CG_DOT2(1, bm, 1, n, 1, r[i], u[i], &gamma[i], w[i], u[i], &delta);

		#pragma omp taskwait

		beta[i] = gamma[i] / gamma[iprev];
		alpha[i] = gamma[i]/(delta - beta[i] * gamma[i] / alpha[iprev]);

		stop_timer(&elapses[k]);
		gamma[iprev] = delta = (double) 0;

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	dump_info("cg_alg3.log", k, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Pipelined PCG */
int CG_ALG4(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
	/* u[0] = M^-1 * r[0] */
	BSBLAS_CHOLSOLV2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, u[i], FP_NOUGHT, w[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		CG_DOT2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);
		/* m[i] = M^-1 * w[i] */
		BSBLAS_CHOLSOLV2(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, m[i], FP_NOUGHT, n0[i]); 

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				cblas_dcopy(c, np1, 1, zp1, 1);
				cblas_daxpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				cblas_dcopy(c, mp1, 1, qp1, 1);
				cblas_daxpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				cblas_dcopy(c, wp0, 1, sp1, 1);
				cblas_daxpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				cblas_dcopy(c, up0, 1, pp1, 1);
				cblas_daxpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				cblas_dcopy(c, xp0, 1, xp1, 1);
				cblas_daxpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				cblas_dcopy(c, rp0, 1, rp1, 1);
				cblas_daxpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				cblas_dcopy(c, up0, 1, up1, 1);
				cblas_daxpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				cblas_dcopy(c, wp0, 1, wp1, 1);
				cblas_daxpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
			HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
		}

		#pragma omp taskwait

		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	dump_info("cg_alg4.log", k, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Gropp PCG */
int CG_ALG7(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(7 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
	/* u[0] = M^-1 * r[0] */
	BSBLAS_CHOLSOLV2(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = u[0] */
	BBLAS_COPY(1, bm, 1, n, 1, u[i], p[i]);
	/* s[0] = A * p[0] */
	HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, p[i], FP_NOUGHT, s[i]);
	/* gamma[0] = <r[0], u[0]> */
	BBLAS_DOT_PURE(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* delta = <p[i], s[i]> */
		BBLAS_DOT_PURE(1, bm, 1, n, 1, p[iprev], s[iprev], &delta);
		/* q[i] = M^-1 * s[i] */
		BSBLAS_CHOLSOLV2(1, bm, n, S, N, s[iprev], q[i]);

		#pragma omp taskwait on(delta)
		alpha[i] = gamma[iprev]/delta;

		/* Axpy fuse x,r,u */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *xp0 = &(x[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *qp0 = &(q[iprev])[j];

			double *xp1 = &(x[i])[j];
			double *pp1 = &(p[i])[j];
			double *rp1 = &(r[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *qp1 = &(q[i])[j];

			#pragma omp task in([c]qp1) out([c]xp1, [c]rp1, [c]up1) priority(1) label(alg7_fuse0)
//			#pragma omp task in([c]qp1) out([c]xp1, [c]rp1, [c]up1) concurrent([n](u[i])) priority(1) label(alg7_fuse0)
			{
				/* x_i+1 = x_i + alpha_i * p_i */
				cblas_dcopy(c, xp0, 1, xp1, 1);
				cblas_daxpy(c, alpha[i], pp0, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				cblas_dcopy(c, rp0, 1, rp1, 1);
				cblas_daxpy(c, -1*alpha[i], sp0, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				cblas_dcopy(c, up0, 1, up1, 1);
				cblas_daxpy(c, -1*alpha[i], qp1, 1, up1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait			
			BBLAS_COPY(1, bm, 1, n, 1, b, r[i]);
			HBSBLAS_CSRMV(1, bm, FP_MONE, Ahbh, x[i], FP_ONE, r[i]);
		}

		/* gamma[i+1] = <r[i+1], u[i+1]> */
		BBLAS_DOT_PURE(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
		/* w[i+1] = A * u[i+1] */
		HBSBLAS_CSRMV(1, bm, FP_ONE, Ahbh, u[i], FP_NOUGHT, w[i]);

		#pragma omp taskwait on(gamma[i])
		beta[i] = gamma[i]/gamma[iprev];

		/* Axpy fuse p,s */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			#pragma omp task in([c]up1, [c]pp0, [c]wp1, [c]sp0) priority(1) label(alg7_fuse1)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				cblas_dcopy(c, up1, 1, pp1, 1);
				cblas_daxpy(c, beta[i], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				cblas_dcopy(c, wp1, 1, sp1, 1);
				cblas_daxpy(c, beta[i], sp0, 1, sp1, 1);
			}
		}

		#pragma omp taskwait

		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	dump_info("cg_alg7.log", k, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}
