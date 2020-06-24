#include "cg_config.h"

#include "fptype.h"

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
extern double prec;
extern int lookahead; 
extern int async;
extern int rep;
extern char *rhsfname;
extern double profile;
extern int warmup;


int cg_config(int argc, char *argv[]) { 
	if ( argc < 5) {
		fprintf(stderr, "use: %s n bm bn s [it] [prec] [lookahead] [async] [profile] [rep] [warmup] [B]\n", argv[0]);
		return 1;
	}

	n = atoi(argv[1]);
	bm = atoi(argv[2]);
	bn = atoi(argv[3]);
	s = atoi(argv[4]);
	if ( bn > s ) {
		bn = s;
		fprintf(stderr, "warn: bn too large, resetting to %i\n", s);
	}

#if 0
	if ( n % bm ) {
		n = (( n + bm - 1 ) / bm ) * bm;
		fprintf(stderr, "warn: n mod bm != 0, padded n to %i\n", n);
	}
	if ( s % bn ) {
		s = (( s + bn - 1 ) / bn ) * bn ;
		fprintf(stderr, "warn: s mod bn != 0, padded s to %i\n", s);
	}
#endif 

	it = n;
	prec = -1.0;
	rep = 1;
	lookahead = 0;
	async = 0;
	profile = 1.0;
	warmup = 0;
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

						if ( profile != 1.0 && async == 0 ) {
							fprintf(stderr, "warn: profiling requires async\n");
							profile = 1.0;
						}
			
						if ( argc > 10 ) {
							rep  = atoi(argv[10]);

							if ( argc > 11 ) {
								warmup = atoi(argv[11]);

								if ( argc > 12 ) {
									rhsfname = argv[12];
								}
							}
						}
					}
				}
			}
		}
	}

	printf("n %i bm %i bn %i s %i it %i prec %.16e lookahead %i async %i profile %e rep %i warmup %i B %s\n", n, bm, bn, s, it, prec, lookahead, async, profile, rep, warmup, rhsfname);

	if ( profile != 1.0 || async == 1 ) {
		printf("warn: don't forget schedule=bf and schedule-priority!\n");
	}

	return 0;
}
