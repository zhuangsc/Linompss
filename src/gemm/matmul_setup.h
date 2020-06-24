#ifndef __MATMUL_SETUP_H__
#define __MATMUL_SETUP_H__


int matmul_setup(const char *fname, int m, int n, int k, int b, int d, int c, void **A, void **B, void **C);
int matmul_shutdown(int m, int n, int k, int b, int d, int c, void *A, void *B, void *C);


#endif // __MATMUL_SETUP_H__
