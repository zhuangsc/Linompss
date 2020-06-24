#include "lu_config.h"

#include <stdio.h>
#include <stdlib.h>

extern int m;
extern int n;
extern int t;
extern int reps;
extern int check;

int lu_config(int argc, char* argv[])
{
	if (argc < 3 ) {
		printf("%s rows columns tilesize [reps] [check]\n",argv[0]);
		return 1;
  	}

  	m = atoi(argv[1]);
  	n = atoi(argv[2]);
  	t = atoi(argv[3]);
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
