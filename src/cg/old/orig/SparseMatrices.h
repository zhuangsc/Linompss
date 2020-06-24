#ifndef SparseMatrixTip

#define SparseMatrixTip 1

typedef struct
	{
		int dim1, dim2;
		int *vptr;
		int *vpos;
		double *vval;
	} SparseMatrix, *ptr_SparseMatrix;

extern FILE *OpenFile (char *name, char *attr);

extern void ReadStringFile (FILE *file, char *string, int length);

extern int CreateSparseMatrix (ptr_SparseMatrix spr, int numR, int numC, 
															 int numE, int msr);

extern int RemoveSparseMatrix (ptr_SparseMatrix spr);

extern int PrintSparseMatrix (SparseMatrix spr, int CorF);

extern int ProdSparseMatrixVector (SparseMatrix spr, double *vec, double *res);

extern int ProdSymSparseMatrixVector (SparseMatrix spr, double *vec, double *res);

extern void CreateSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC);

#endif
