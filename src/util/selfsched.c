#include "selfsched.h"

#include <stdlib.h>



int sched_ipmalloc(int initid, int c, selfsched_t *p) 
{
	if ( c > 64 ) {
		return 1; 
	}

	p->ready = 0;
	p->id = initid;
	p->c = c;
	p->fshed = 0;

	pthread_mutex_init(&p->mutex, NULL);
	pthread_cond_init(&p->cond, NULL);

	return 0;
}


void * sched_malloc(int initid, int c) 
{
	selfsched_t *s = malloc(sizeof(selfsched_t));
	
	if ( sched_ipmalloc(initid, c, s) ) {
		free(s);
		return NULL;
	}

	return s;
}


void sched_free(selfsched_t *sched) 
{
	pthread_mutex_destroy(&sched->mutex);
	pthread_cond_destroy(&sched->cond);

	free(sched);
}


void sched_alldone(int n, selfsched_t *sched, int k) 
{
	int c;
	for ( c=0; c<n; ++c ) {
		selfsched_t *csched = &sched[c];

		int partc = csched->c;
		int i;
		for ( i=0; i<partc; ++i ) {
			sched_done(csched, k, i);
		}
	}
}
