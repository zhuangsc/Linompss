#ifndef __IOHB_H__
#define __IOHB_H__

#include	<stdio.h>
#include	<stdlib.h>
#include	<malloc.h>
#include 	"fptype.h"

#ifdef SINGLE_PRECISION

#define readHB_info			readHB_info_single
#define readHB_header		readHB_header_single
#define readHB_mat			readHB_mat_single
#define readHB_newmat		readHB_newmat_single
#define readHB_aux			readHB_aux_single
#define readHB_newaux		readHB_newaux_single
#define writeHB_mat			writeHB_mat_single
#define readHB_mat_char		readHB_mat_char_single
#define readHB_newmat_char	readHB_newmat_char_single
#define readHB_aux_char		readHB_aux_char_single
#define readHB_newaux_char	readHB_newaux_char_single
#define writeHB_mat_char	writeHB_mat_char_single
#define ParseIfmt			ParseIfmt_single
#define ParseRfmt			ParseRfmt_single
#define IOHBTerminate		IOHBTerminate_single
#define upcase				upcase_single
#define substr				substr_single

#else

#define readHB_info			readHB_info_double
#define readHB_header		readHB_header_double
#define readHB_mat			readHB_mat_double
#define readHB_newmat		readHB_newmat_double
#define readHB_aux			readHB_aux_double
#define readHB_newaux		readHB_newaux_double
#define writeHB_mat			writeHB_mat_double
#define readHB_mat_char		readHB_mat_char_double
#define readHB_newmat_char	readHB_newmat_char_double
#define readHB_aux_char		readHB_aux_char_double
#define readHB_newaux_char	readHB_newaux_char_double
#define writeHB_mat_char	writeHB_mat_char_double
#define ParseIfmt			ParseIfmt_double
#define ParseRfmt			ParseRfmt_double
#define IOHBTerminate		IOHBTerminate_double
#define upcase				upcase_double
#define substr				substr_double

#endif



#ifdef __cplusplus
extern "C" {
#endif
int readHB_info_single(const char* filename, int* M, int* N, int* nz, char** Type, 
                                                      int* Nrhs);

int readHB_header_single(FILE* in_file, char* Title, char* Key, char* Type, 
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt, 
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd, 
                    char *Rhstype);

int readHB_mat_single(const char* filename, int colptr[], int rowind[], 
                                                                 fp_t val[]);

int readHB_newmat_single(const char* filename, int* M, int* N, int* nonzeros, 
                         int** colptr, int** rowind, fp_t** val);

int readHB_aux_single(const char* filename, const char AuxType, fp_t b[]);

int readHB_newaux_single(const char* filename, const char AuxType, fp_t** b);

int writeHB_mat_single(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const fp_t val[], int Nrhs, const fp_t rhs[], 
                        const fp_t guess[], const fp_t exact[],
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int readHB_mat_char_single(const char* filename, int colptr[], int rowind[], 
                                           char val[], char* Valfmt);

int readHB_newmat_char_single(const char* filename, int* M, int* N, int* nonzeros, int** colptr, 
                          int** rowind, char** val, char** Valfmt);

int readHB_aux_char_single(const char* filename, const char AuxType, char b[]);

int readHB_newaux_char_single(const char* filename, const char AuxType, char** b, char** Rhsfmt);

int writeHB_mat_char_single(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const char val[], int Nrhs, const char rhs[], 
                        const char guess[], const char exact[], 
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int ParseIfmt_single(char* fmt, int* perline, int* width);

int ParseRfmt_single(char* fmt, int* perline, int* width, int* prec, int* flag);

void IOHBTerminate_single(char* message);
void upcase_single(char* S);
char* substr_single(const char* S, const int pos, const int len);


////////////////////////////////////
int readHB_info_double(const char* filename, int* M, int* N, int* nz, char** Type, 
                                                      int* Nrhs);

int readHB_header_double(FILE* in_file, char* Title, char* Key, char* Type, 
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt, 
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd, 
                    char *Rhstype);

int readHB_mat_double(const char* filename, int colptr[], int rowind[], 
                                                                 fp_t val[]);

int readHB_newmat_double(const char* filename, int* M, int* N, int* nonzeros, 
                         int** colptr, int** rowind, fp_t** val);

int readHB_aux_double(const char* filename, const char AuxType, fp_t b[]);

int readHB_newaux_double(const char* filename, const char AuxType, fp_t** b);

int writeHB_mat_double(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const fp_t val[], int Nrhs, const fp_t rhs[], 
                        const fp_t guess[], const fp_t exact[],
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int readHB_mat_char_double(const char* filename, int colptr[], int rowind[], 
                                           char val[], char* Valfmt);

int readHB_newmat_char_double(const char* filename, int* M, int* N, int* nonzeros, int** colptr, 
                          int** rowind, char** val, char** Valfmt);

int readHB_aux_char_double(const char* filename, const char AuxType, char b[]);

int readHB_newaux_char_double(const char* filename, const char AuxType, char** b, char** Rhsfmt);

int writeHB_mat_char_double(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const char val[], int Nrhs, const char rhs[], 
                        const char guess[], const char exact[], 
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int ParseIfmt_double(char* fmt, int* perline, int* width);

int ParseRfmt_double(char* fmt, int* perline, int* width, int* prec, int* flag);

void IOHBTerminate_double(char* message);
void upcase_double(char* S);
char* substr_double(const char* S, const int pos, const int len);

#ifdef __cplusplus
}
#endif

#endif
