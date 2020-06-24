#include "ompss_cg.h"

//#include "cg_main.h"
//#include "cgmod1_main.h"
#include "cgprof_main.h"


int ompss_scgprof(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, int steps, 
				void *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog, int release, 
				double orth_fac, float *mat_energy, int is_precond, void *zbuff, cs *Acs, css **S, csn **N, int interval, int corrections) {
	return scgprof(bm, bn, n, A, s, b, x, offs, tol, steps, work, works, lookahead, async, profile, warmup, cglog, release, orth_fac, mat_energy, is_precond, zbuff, Acs, S, N, interval, corrections);
}

int ompss_dcgprof(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, int steps, 
		void *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog, int release, 
		double orth_fac, float *mat_energy, int is_precond, void *zbuff, cs *Acs, css **S, csn **N, int interval, int corrections) {
	return dcgprof(bm, bn, n, A, s, b, x, offs, tol, steps, work, works, lookahead, async, profile, warmup, cglog, release, orth_fac, mat_energy, is_precond, zbuff, Acs, S, N, interval, corrections);
}
