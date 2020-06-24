#ifndef __HBDEBUG_H__
#define __HBDEBUG_H__


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>


#include "hb.h"
#include "hbconvrt.h"
#include "vector.h"


#ifdef SINGLE_PRECISION

#define hb_dense_print 				hb_dense_print_single
#define hb_struc_diff				hb_struc_diff_single
#define csrb_dense_printf			csrb_dense_printf_single
#define hbb_dense_printf			hbb_dense_printf_single
#define hb_elemat					hb_elemat_single
#define hb_vecat					hb_vecat_single
#define	print_etree					print_etree_single
#define print_matrix				print_matrix_single
#define print_address				print_address_single

#else

#define hb_dense_print 				hb_dense_print_double
#define hb_struc_diff				hb_struc_diff_double
#define csrb_dense_printf			csrb_dense_printf_double
#define hbb_dense_printf			hbb_dense_printf_double
#define hb_elemat					hb_elemat_double
#define hb_vecat					hb_vecat_double
#define	print_etree					print_etree_double
#define print_matrix				print_matrix_double
#define print_address				print_address_double

#endif


void hb_dense_print(const char *name, hbmat_t *Ahb);
int hb_struc_diff(hbmat_t *A, hbmat_t *B);
void csrb_dense_printf(FILE *f, const char *name, hbmat_t *Acsrb);
void hbb_dense_printf(FILE *f, const char *name, hbmat_t *Ahbb, int struc, int detail);
void hb_elemat(hbmat_t *A, int i, int j);
void hb_vecat(hbmat_t *A, int i);
void print_etree(int* etree_ptr, int col);
void print_matrix(const hbmat_t* matrix_info, int h, char* name);
void print_address(const hbmat_t* A, int h, char* name);


#endif // __HBDEBUG_H__

