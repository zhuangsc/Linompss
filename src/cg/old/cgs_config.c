#include "cg_config.h"

#include "fptype.h"
#include "hb.h"
#include "iohb.h"
#include "matfprint.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


extern int n;
extern int bm;
extern int bn;
extern int s;
extern int it;
extern fp_t prec;
extern int lookahead; 
extern int async;
extern int rep;
extern char *rhsfname;
extern void *A;
extern fp_t profile;


int cg_config(int argc, char *argv[]) 
{
	if ( argc < 5 ) {
		fprintf(stderr, "use: %s Afile bm bn s [it] [prec] [lookahead] [async] [profile] [rep] [B]\n", argv[0]);
		return 1;
	}

	struct stat buf;
	if ( stat(argv[1], &buf) != 0 ) {
		fprintf(stderr, "config: could not open %s for reading\n", argv[1]);
		return 2;
	}

	hbmat_t *Ahb = (hbmat_t*) malloc(sizeof(hbmat_t));
	hb_reset(Ahb);
	if ( readHB_newmat(argv[1], &(Ahb->m), &(Ahb->n), &(Ahb->elemc), &(Ahb->vptr), &(Ahb->vpos), (fp_t **)&(Ahb->vval)) ) {
		fprintf(stderr, "config: could not read %s\n", argv[1]);
		return 3;
	}

	printf("warn: triang form of A not yet supported\n");

	n = Ahb->m;
	A = Ahb;

	bm = atoi(argv[2]);
	bn = atoi(argv[3]);
	s = atoi(argv[4]);
	if ( bn > s ) {
		bn = s;
		fprintf(stderr, "warn: bn too large, resetting to %i\n", s);
	}

	it = n;
	prec = -1.0;
	rep = 1;
	lookahead = 0;
	async = 0;
	profile = 1.0;
	rhsfname = NULL;

	if ( argc > 5 ) {
		it = atoi(argv[5]);
		if ( it > n ) {
			fprintf(stderr, "warn: setting it to %i\n", n);
			it = n;
		}
		if ( argc > 6 ) {
			prec = atof(argv[6]);
			if ( argc > 7 ) {
				lookahead = atoi(argv[7]);
				if ( argc > 8 ) {
					async = atoi(argv[8]);

					if ( argc > 9 ) {
						profile  = atof(argv[9]);
						if ( profile <= 0.0 ) {
							fprintf(stderr, "warn: setting profile to 1.0\n");
							profile = 1.0;
						}
			
						if ( argc > 10 ) {
							rep  = atoi(argv[10]);

							if ( argc > 11 ) {
								rhsfname = argv[11];
							}
						}
					}
				}
			}
		}
	}

	return 0;
}
