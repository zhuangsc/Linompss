#include <stdio.h>
#include <stdlib.h>
#include "ScalarVectors.h"
#include "SparseMatrices.h"

#define length1 82
#define length2 82

FILE *OpenFile (char *name, char *attr)
	{
    FILE *fich;
    if ((fich = fopen (name, attr)) == NULL)
      { printf ("File %s not exists \n", name); exit(1); }
		return fich;
	}

void ReadStringFile (FILE *file, char *string, int length)
	{
		char *s = NULL;
		if ((s = fgets (string, length, file)) == NULL)
			{ printf ("Error en lectura \n"); exit (1); }
	}

int CreateSparseMatrix (ptr_SparseMatrix spr, int numR, int numC, int numE, 
												int msr)
	{
		spr->dim1 = numR; spr->dim2 = numC; 
		CreateInts (&(spr->vptr), numE+numR+1);
		*(spr->vptr) = ((msr)? (numR+1): 0);
		spr->vpos = spr->vptr + ((msr)? 0: (numR+1));
		CreateDoubles (&(spr->vval), numE+(numR+1)*msr);
	}

int PrintSparseMatrix (SparseMatrix spr, int CorF)
	{
		int i, j;

		if (spr.vptr == spr.vpos)
			{
				printf ("Diagonals : \n ");
				for (i=0; i<spr.dim1; i++) printf ("%f ", spr.vval[i]); printf ("\n");
			}

		printf ("Pointers: \n ");
		if (spr.dim1 > 0)
			for (i=0; i<=spr.dim1; i++) printf ("%d ", spr.vptr[i]); printf ("\n");

		printf ("Values: \n");
		for (i=0; i<spr.dim1; i++)
			{ printf (" Row %d --> ", i+CorF);
			for (j=(spr.vptr[i]-CorF); j<(spr.vptr[i+1]-CorF); j++)
				printf ("(%d,%f) ", spr.vpos[j], spr.vval[j]); 
			printf ("\n");	}
		printf ("\n");
	}

int RemoveSparseMatrix (ptr_SparseMatrix spr)
	{
		spr->dim1 = -1; spr->dim2 = -1; 
		RemoveInts (&(spr->vptr));
		RemoveDoubles (&(spr->vval)); 
	}

int ProdSparseMatrixVector (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		double aux;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++)
					aux += spr.vval[j] * vec[spr.vpos[j]];
				res[i] += aux;
			}
	}

int ProdSymSparseMatrixVector (SparseMatrix spr, double *vec, double *res)
	{
		int i, j, k;
		double aux, val;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++) {
					k = spr.vpos[j]; val = spr.vval[j];
					aux += val * vec[k];
					if (k != i) res[k] += (val * vec[i]);
				}
				res[i] += aux;
			}
	}

void CreateSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC)
	{
		FILE *file;
		char string[length1], *s = NULL;
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int *vptr = NULL, *vpos = NULL;
		double *vval = NULL; 
		int lines[5], dim[4], formats[10];

		file = OpenFile (nameFile, "r");
		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0); 
		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		CreateSparseMatrix (spr, dim[0], dim[1], dim[2], 0);
		vptr = spr->vptr; vpos = spr->vpos; vval = spr->vval; 

		ReadStringFile (file, string, length1);
		GetFormatsFromString (string, formats, 2, 16);
		GetFormatsFromString ((string+32), (formats+4), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) ReadStringFile (file, string, length1);

		j = 0;
		for (i = 0; i < lines[1]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = ((dim[0] + 1) - j);
				if (k > formats[0]) k = formats[0];
				GetIntsFromString (string, (vptr+j), k, formats[1], shft);
				j+=formats[0];
			}

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[2]) k = formats[2];
				GetIntsFromString (string, (vpos+j), k, formats[3], shft);
				j+=formats[2];
			}

		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetDoublesFromString (string, (vval+j), k, formats[5]);
				j+=formats[4];
			}

		fclose (file);
	}

