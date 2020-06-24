#include "gensolv_config.h"

#include <stdlib.h>
#include <stdio.h>


extern int m;
extern int n;
extern int b;
extern int reps;
extern int check;


int gensolv_config(int argc, char *argv[]) {
	if (argc < 4 ) {
		printf("usage: %s m n b [reps] [check]\n",argv[0]);
		return 1;
	}

	m = atoi(argv[1]);
	n = atoi(argv[2]);
	b = atoi(argv[3]);

#if 0
	if ( m % b ) {
		m = ( ( m + b - 1 ) / b ) * b;
		printf("warn: padding m to %i\n", m);
	}
#endif

	reps = 1;
	check = 0;

	if ( argc > 4 ) {
		reps = atoi(argv[4]);
		if ( argc > 5 ) {
			check = atoi(argv[5]);
		}
  	}

	return 0;
}
