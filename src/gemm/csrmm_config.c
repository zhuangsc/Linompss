#include "matmul_config.h"


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


extern int m;
extern int n;
extern int k;
extern int b;
extern int reps;
extern int check;
extern char *fname;


int matmul_config(int argc, char *argv[]) 
{ 
	if (argc < 4 ) {
		fprintf(stderr, "use: %s HBfile n b [reps] [check]\n", argv[0]);
		return 1;
	}

  	fname = argv[1];
	struct stat buf;
	if ( stat(fname, &buf) == -1 ) {
		fprintf(stderr, "config: file %s does not exist\n", fname);
		return 2;
	}
	
  	n = atoi(argv[2]);
  	b = atoi(argv[3]);

	reps = 1;
	check = 0;
	if (argc > 4) {
		reps = atoi(argv[4]);
		if (argc > 5) {
			check = atoi(argv[5]);
		}
	}

	return 0;
}
