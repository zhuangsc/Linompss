#include "fptype.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#if USE_SPARSE
#include "hb.h"
#include "iohb.h"
#include "hbconvrt.h"
#include "scsparse.h"
extern void *Ahb;
extern void *Acsr;
#endif

extern int n;
extern int bm;
extern int bn;
extern int s;
extern int it;
extern double prec;
extern int is_precond;
extern int async;
extern int rep;
extern char *rhsfname;
extern char *aname;
extern double profile;
extern int warmup;

extern int refit;
extern int release;
extern double icriteria;
extern double ii_distance;
extern double orth_fac;
extern int cglog_level;

extern int interval;
extern int corrections;

int itref_config(int argc, char *argv[]) 
{

#if USE_SPARSE
	if ( argc < 5 ) {
		fprintf(stderr, "use: %s [bm] [interval] [corrections] [refit] [it] [precision] [is_precond] [async] [profile] [rep] [warmup] [release_id] [icriteria] [ii_distance] [orth_fac] [HB_MAT] [FULL?] [B] [cglog]\n", argv[0]);
		return 1;
	}

	int is_full = 1;
	bm = atoi(argv[1]);
	interval = atoi(argv[2]);
	corrections = atoi(argv[3]);
	refit = atoi(argv[4]);

	if ( interval < corrections ) {
		fprintf(stderr, "Warning: interval < corrections, changing interval = corrections+1\n");
		corrections = interval - 1;
	}
	it = refit;
	prec = -1.0;
	rep = 1;
	is_precond = 0;
	async = 0;
	profile = 1.0;
	warmup = 0;
	release = 1;
	icriteria = 1E-5;
	ii_distance = 1E-2;
	rhsfname = NULL;
	aname = NULL;
	orth_fac = 1.0;
	cglog_level = 1;

	if ( argc > 5 ) {
		it = atoi(argv[5]);
		if ( argc > 6 ) {
			prec = atof(argv[6]);
			if ( argc > 7 ) {
				is_precond = atoi(argv[7]);
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
									release = atoi(argv[12]);
									if ( argc > 13 ) {
										icriteria = atof(argv[13]);
										if ( argc > 14 ) {
											ii_distance = atof(argv[14]);
											if ( argc > 15 ) {
												orth_fac = atof(argv[15]);
												if ( argc > 16) {
													aname = argv[16];
													if ( argc > 17 ) {
														is_full = atoi(argv[17]);
														if ( argc > 18 ) {
																rhsfname = argv[18];
															if ( argc > 19 ) {
																cglog_level = atoi(argv[19]);
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

#else
	if ( argc < 6) {
		fprintf(stderr, "use: %s n bm [interval] [corrections] refit [it] [prec] [is_precond] [async] [profile] [rep] [warmup] [release] [icriteria] [ii_distance] [orth_fac] [A] [B] [cglog_level]\n", argv[0]);
		return 1;
	}

	n = atoi(argv[1]);
	bm = atoi(argv[2]);
	interval = atoi(argv[3]);
	corrections = atoi(argv[4]);
	if ( interval < corrections ) {
		fprintf(stderr, "Warning: interval < corrections, changing interval = corrections+1\n");
		corrections = interval - 1;
	}

//	if ( bn > s ) {
//		bn = s;
//		fprintf(stderr, "warn: bn too large, resetting to %i\n", s);
//	}
	refit = atoi(argv[5]);

	it = n;
	prec = -1.0;
	rep = 1;
	is_precond = 0;
	async = 0;
	profile = 1.0;
	warmup = 0;
	release = 1;
	icriteria = 1E-5;
	ii_distance = 1E-2;
	rhsfname = NULL;
	aname = NULL;
	orth_fac = 1.0;
	cglog_level = 1;

	if ( argc > 6 ) {
		it = atoi(argv[6]);
		if ( it > n ) {
			fprintf(stderr, "warn: setting it to %i\n", n);
			it = n;
		}
		if ( argc > 7 ) {
			prec = atof(argv[7]);
			if ( argc > 8 ) {
				is_precond = atoi(argv[8]);
				if ( argc > 9 ) {
					async = atoi(argv[9]);

					if ( argc > 10 ) {
						profile  = atof(argv[10]);
						if ( profile <= 0.0 ) {
							fprintf(stderr, "warn: setting profile to 1.0\n");
							profile = 1.0;
						}

						if ( profile != 1.0 && async == 0 ) {
							fprintf(stderr, "warn: profiling requires async\n");
							profile = 1.0;
						}
			
						if ( argc > 11 ) {
							rep  = atoi(argv[11]);
							if ( argc > 12 ) {
								warmup = atoi(argv[12]);
								if ( argc > 13 ) {
									release = atoi(argv[13]);
									if ( argc > 14 ) {
										icriteria = atof(argv[14]);
										if ( argc > 15 ) {
											ii_distance = atof(argv[15]);
											if ( argc > 16 ) {
												orth_fac = atof(argv[16]);
												if ( argc>17) {
													aname = argv[17];
													if ( argc > 18 ) {
														rhsfname = argv[18];
														if ( argc > 19 ) {
															cglog_level = atoi(argv[19]);
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
#endif

//	if ( profile != 1.0 || async == 1 ) {
//		printf("warn: don't forget schedule=bf and schedule-priority!\n");
//	}
#if USE_SPARSE
	hbmat_t *tmp = malloc(sizeof(hbmat_t));
	hb_reset(tmp);
	readHB_newmat_double(aname, &(tmp->m), &(tmp->n), &(tmp->elemc), &(tmp->vptr), &(tmp->vpos), (double **)&(tmp->vval));
	tmp->b = bm;
	one2zero(tmp);
	Acsr = tmp;
	if ( !is_full ) {
		Ahb = malloc(sizeof(hbmat_t));
		dhb_sym_expand(Ahb, tmp);
		hb_free(tmp);
		Acsr = Ahb;
	}
	n = ((hbmat_t*)Acsr)->m;
	if ( it > n )
		it = n;
#endif

	printf("n %i bm %i interval %i corrections %i refit %i it %i prec %.16e is_precond %i async %i profile %f rep %i \
			warmup %i release_id %d icriteria %E ii_distance %E orth_fac: %E A %s B: %s, loglevel: %d\n", \
			n, bm, interval, corrections, refit, it, prec, is_precond, async, profile, rep, warmup, \
			release, icriteria, ii_distance, orth_fac, aname, rhsfname, cglog_level);

	return 0;
}
