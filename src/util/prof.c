#include "prof.h"


void prof_print(prof_t *prof)
{
#if 0
	int s = prof->s;
	int r = prof->read;
	float *p = prof->p + r * s;
	
	int k;
	for ( k=0; k<s; ++k ) {
		printf("%e %i\n", p[k], prior[k]);
	}
#endif
}


