#include "qr_config.h"

#include <stdlib.h>
#include <stdio.h>


extern int m;
extern int n;
extern int tr;
extern int tc;
extern int bs;
extern int reps;
extern int check;
extern int hierfactor;


int qr_config(int argc, char *argv[]) 
{
	if (argc < 5 ) {
		fprintf(stderr, "usage: %s m n tr tc [bs] [reps] [check] [hierfactor]\n",argv[0]);
		return 1;
  	}

	m = atoi(argv[1]);
	n = atoi(argv[2]);
	tr = atoi(argv[3]);
	tc = atoi(argv[4]);


	bs = tc;
	reps = 1;
	check = 0;

	if (tc > tr) {
		fprintf(stderr, "tr must be greater than or equal to tc\n");
		return 2;
	}
	
	if (tr % tc) {
		fprintf(stderr, "tr must be an integral multiple of tc\n");
		return 3;
	}

	if (argc > 5) {
		bs = atoi(argv[5]);
		if (tc % bs) {
			fprintf(stderr, "bs must be an integer multiple of tc\n");
			return 4;
		}
		if (argc > 6) {
			reps=atoi(argv[6]);

			if (argc > 7) {
  				check=atoi(argv[7]);
      		}
		}
	}

	return 0;
}
