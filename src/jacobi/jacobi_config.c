#include "jacobi_config.h"


#include <stdio.h>
#include <stdlib.h>
#include "fptype.h"


extern int bs;
extern int format;
extern int max_iter;
extern int lookahead;
extern fp_t threshold;
extern char *fname;
extern int check;
extern int rep;


int jacobi_config(int argc, char *argv[]) 
{
	if ( argc < 3 ) {
		printf("Usage %s HBfile b [csc(0)/csr(1)] [iter] [lookhead] [threshold] [check] [rep]\n", argv[0]);
		return 1;
	}

	fname = argv[1];
	bs = atoi(argv[2]);

	format = 1;
	max_iter = 1;
	lookahead = 0;
	threshold = 0;
	check = 0;
	rep = 1;

	if ( argc > 3 ) {
		format = atoi(argv[3]);

		if ( argc > 4 ) {
			max_iter = atoi(argv[4]);

			if ( argc > 5 ) {
				lookahead = atoi(argv[5]);

				if ( argc > 6 ) {
					threshold = atof(argv[6]);

					if ( argc > 7 ) {
						check = atoi(argv[7]);

						if ( argc > 8 ) {
							rep = atoi(argv[8]);
						}
					}
				}
			}			
		}
	}

	return 0;
}
