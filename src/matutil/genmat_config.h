#ifndef __GENMAT_CONFIG_H__
#define __GENMAT_CONFIG_H__


#include "ompss_mat.h"


#if SINGLE_PRECISION
#define GENMAT_CONFIG		sgenmat_config
#else
#define GENMAT_CONFIG		dgenmat_config
#endif


#define LIBBBLAS_EXPORT __attribute__((__visibility__("default")))


extern LIBBBLAS_EXPORT void* sgenmat_config(char *app, char *id, int m, int n, int ldim, ompssmat_t deflt); 
extern LIBBBLAS_EXPORT void* dgenmat_config(char *app, char *id, int m, int n, int ldim, ompssmat_t deflt); 


#endif // __GENMAT_CONFIG_H__
