#ifndef __TASKS_POTRF_CSR_H__
#define __TASKS_POTRF_CSR_H__


#include "fptype.h"
#include "hb.h"


#define LIBSBBLAS_EXPORT __attribute__((__visibility__("default")))


#ifdef SINGLE_PRECISION

#define TASK_POTRF_CSR 		task_spotrf_csr

#else

#define TASK_POTRF_CSR		task_dpotrf_csr

#endif


extern LIBSBBLAS_EXPORT void task_dpotrf_csr(hbmat_t *A);
extern LIBSBBLAS_EXPORT void task_spotrf_csr(hbmat_t *A);


#endif // __TASKS_POTRF_CSR_H__
