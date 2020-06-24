#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ScalarVectors.h"

int CreateInts (int **vint, int num)
	{
		if ((*vint = (int *) malloc (sizeof(int)*num)) == NULL)
			{ printf ("Memory Error (CreateInts(%d))\n", num); exit (1); }
	}

int InitInts (int *vint, int n, int frst, int incr) 
	{
		int i, *p1 = vint, num = frst;

		for (i=0; i<n; i++) 
			{ *(p1++) = num; num += incr; }
	}
		
int CopyShiftInts (int *src, int *dst, int n, int shft)
	{
		int i, *p1 = src, *p2 = dst;

		for (i=0; i<n; i++)
			*(p2++) = *(p1++) + shft;
	}
	
void GetIntFromString (char *string, int *pnum, int numC, int shft)
	{
		int j = 0, num = 0, neg = 0;
		char *pchar = string;

		while ((j < numC) && ((*pchar < '0') || (*pchar > '9')) &&
           (*pchar != '+') && (*pchar != '-')) { j++; pchar++; }
		if (j < numC)
			{
				if ((*pchar == '+') || (*pchar == '-'))
					{ neg = (*pchar == '-'); j++; pchar++; }
				while ((j < numC) && (*pchar >= '0') && (*pchar <= '9'))
					{ num = num * 10 + (*pchar - 48); j++; pchar++; }
			}
		if (neg) num = -num;
		*pnum = num + shft; 
	}

void GetIntsFromString (char *string, int *vec, int numN, int numC, int shft)
	{
		int i, *pint = vec;
		char *pchar = string;

		for (i=0; i<numN; i++)
			{ GetIntFromString (pchar, (pint++), numC, shft); pchar += numC; }
	}

void GetFormatsFromString (char *string, int *vec, int numN, int numC)
	{
		int i, k = 0;
		int *pint = vec;
		char *pchar = string, *pch = NULL, c = ' ', c2;

		for (i=0; i<numN; i++)
			{ 
				pch = pchar; 
				while (*pch == ' ') pch++;
				sscanf (pch, "(%i%c", pint, &c);
				if ((c == 'P') || (c == 'p')) 
					{ 
						sscanf (pch, "(%i%c%i%c%i.%i)", &k, &c2, pint, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else if ((c == 'E') || (c == 'e') || (c == 'D') || (c == 'd')	||
								 (c == 'F') || (c == 'f') || (c == 'G') || (c == 'g'))
					{ 
						sscanf (pch, "(%i%c%i.%i)", pint, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else 
					{ sscanf (pch, "(%i%c%i)", pint, &c, pint+1); pint += 2; }
				pchar += numC;
			}
	}

int PrintInts (int *vint, int num)
	{
		int i, *pi = vint;

		for (i=0; i<num; i++) printf ("%d ", *(pi++));
		printf ("\n");
	}

int RemoveInts (int **vint)
	{ free (*vint); *vint = NULL; }

int CreateDoubles (double **vdouble, int num)
	{
		if ((*vdouble = (double *) malloc (sizeof(double)*num)) == NULL)
			{ printf ("Memory Error (CreateDoubles(%d))\n", num); exit (1); }
	}

int InitDoubles (double *vdouble, int n, double frst, double incr) 
	{
		int i; 
		double *p1 = vdouble, num = frst;

		for (i=0; i<n; i++) 
			{ *(p1++) = num; num += incr; }
	}
		
void GetDoubleFromString (char *string, double *pdbl, int numC)
	{
		int j, k, exp, neg;
		double num, frac;
		char *pchar = string;

		j = 0; exp = 0; neg = 0; num = 0.0; frac = 1.0;
		while ((j < numC) && ((*pchar < '0') || (*pchar > '9')) && 
					 (*pchar != '+') && (*pchar != '-') && (*pchar != '.')) { j++; pchar++; }
		if (j < numC)
			{
				if ((*pchar == '+') || (*pchar == '-'))
					{ neg = (*pchar == '-'); j++; pchar++; }
				if (j < numC)
					{
						if (*pchar != '.')
							while ((j < numC) && (*pchar >= '0') && (*pchar <= '9'))
								{ num = num * 10 + (*pchar - 48); j++; pchar++; }
						if (j < numC)
							{
								if (*pchar == '.')
									{
										j++; pchar++; 
										while ((j < numC) && (*pchar >= '0') && (*pchar <= '9'))
											{ frac /= 10; num += (*pchar-48) * frac; j++; pchar++; }
									}
								if (neg) num = -num;
								if (j < numC)
									{
										if ((*pchar == 'e') || (*pchar == 'E') || (*pchar == 'd') || (*pchar == 'D'))
											{
												neg = 0; j++; pchar++; 
												if (j < numC)
													{
														if ((*pchar == '+') || (*pchar == '-'))
															{ neg = (*pchar == '-'); j++; pchar++; }
														if (j < numC)
															{
																while ((j < numC) && (*pchar >= '0') && 
																		 (*pchar <= '9'))
																	{ exp = exp*10 + (*pchar-48); j++; pchar++; }
																if (neg) exp = -exp;
																for (k=0; k<exp; k++) num *= 10;
																for (k=0; k>exp; k--) num /= 10;
															}
													}
											}
									}
							}
						else
 							if (neg) num = -num;
					}
			}
		*pdbl = num; 
	}

void GetDoublesFromString (char *string, double *vec, int numN, int numC)
	{
		int i;
		double *paux = vec;
		char *pchar = string;

		for (i=0; i<numN; i++)
			{ GetDoubleFromString (pchar, (paux++), numC); pchar += numC; }
	}

int PrintDoubles (double *vdouble, int num)
	{
		int i;
		double *pd = vdouble;

		for (i=0; i<num; i++) printf ("%e ", *(pd++));
		printf ("\n");
	}

int RemoveDoubles (double **vdouble)
	{ free (*vdouble); *vdouble = NULL; }

