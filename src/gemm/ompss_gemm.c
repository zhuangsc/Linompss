#include "ompss_gemm.h"

#include "bblas_gemm.h"


int ompss_hyper_sgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC) {
	bblas_hyper_sgemm(b, d, c, m, n, k, alpha, A, ldimA, B, ldimB, beta, C, ldimC);
}

int ompss_hyper_dgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC) {
	bblas_hyper_sgemm(b, d, c, m, n, k, alpha, A, ldimA, B, ldimB, beta, C, ldimC);
}

int ompss_sgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC) {
	bblas_sgemm(b, d, c, m, n, k, alpha, A, ldimA, B, ldimB, beta, C, ldimC);
}

int ompss_dgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC) {
	bblas_sgemm(b, d, c, m, n, k, alpha, A, ldimA, B, ldimB, beta, C, ldimC);
}
