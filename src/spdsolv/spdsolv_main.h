#ifndef __SPDSOLV_MAIN_H__
#define __SPDSOLV_MAIN_H__


#ifdef SINGLE_PRECISION

#define SPDSOLV_MAIN	sspdsolv_main
#define OMPSS_CHOL_LL   ompss_schol_ll
#define OMPSS_CHOL_RL   ompss_schol_rl

#else 

#define SPDSOLV_MAIN	dspdsolv_main
#define OMPSS_CHOL_LL   ompss_dchol_ll
#define OMPSS_CHOL_RL   ompss_dchol_rl

#endif

int sspdsolv_main(int m, int n, int b, float *A, float *B); 
int dspdsolv_main(int m, int n, int b, double *A, double *B); 


#endif // __SPDSOLV_MAIN_H__
