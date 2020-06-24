#ifndef __SELFSCHED_H__
#define __SELFSCHED_H__


#include <pthread.h>
#include <stdint.h>
#include <stdio.h>

typedef struct selfsched { 
	pthread_cond_t cond;
	pthread_mutex_t mutex;
	uint64_t ready;
	uint64_t blocked;
	int id;
	int fshed;
	int c;
} selfsched_t;


void * sched_malloc(int initid, int c);
int sched_ipmalloc(int initid, int c, selfsched_t *p);
void sched_free(selfsched_t *sched);
void sched_alldone(int n, selfsched_t *sched, int k);


/* k : identifier for current iteration 
   i : identifier for block */
static inline __attribute__((always_inline)) void sched_done(selfsched_t *sched, int k, int i) 
{
	uint64_t flag = 1 << i;
	pthread_mutex_t *mutex = &sched->mutex;
	int reset;
	int ok;
	int fshed;

	pthread_mutex_lock(mutex);
	
	int bc = sched->c;
	int id = sched->id;
	uint64_t ready = sched->ready;	
	reset = k!=id;
	fshed = sched->fshed;
	ok = fshed == bc * bc;
	ready = reset ? 0 : ready;
	sched->id += reset ? 1 : 0;
	sched->fshed = reset ? 0 : fshed;
	sched->ready = ready | flag;	

	//printf("sched done k %i id %i i %i reset %i\n", k, id, i, reset);

	uint64_t waiting = sched->blocked & flag;

	if ( waiting ) { /* wake up waiting tasks */
		sched->blocked &= ~flag;
		pthread_cond_broadcast(&sched->cond);
	}

	pthread_mutex_unlock(mutex);

	if ( !ok && reset ) {
		printf("selfsched broke %i %i\n", k, fshed);
	}
}


static inline __attribute__((always_inline)) void sched_wait(selfsched_t *sched, int i) 
{
	uint64_t flag = 1 << i;
	pthread_mutex_t *mutex = &sched->mutex;

	pthread_mutex_lock(mutex);
	
	uint64_t ready = sched->ready & flag;	
	flag &= ready==0 ? ~0 : 0; 
	//printf("wait: %i ready %i\n", i, ready);
	sched->blocked |= flag;
	sched->fshed += 1;
	while ( ready == 0 ) {
		//printf("blocked!\n");
		pthread_cond_wait(&sched->cond, mutex);
		ready = sched->ready & flag;
	} 

	pthread_mutex_unlock(mutex);
}


#endif // __SELFSCHED_H__
