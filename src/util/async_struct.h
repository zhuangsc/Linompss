#ifndef __ASYNC_STRUCT_H__
#define __ASYNC_STRUCT_H__


#include <pthread.h>
#include <stdio.h>

#include "prof.h"


typedef enum 
{
	STAT_AHEAD, 
	STAT_BROKEN, 
	STAT_SYNC, 
	STAT_CONVERGED, 
	STAT_ERROR
} async_stat_t;


typedef struct asynct 
{
	int id; /* -1 if not initialized yet */
	int create;
	int dot_control; 
	int ccnt;
	int pcnt;
	int pcompl;
	int consume;
	int wait;
	int flags;
	volatile int ready;
	void *log;
	FILE *logf;
	int logc;
	int logs;
	prof_t prof;
	pthread_cond_t cond;
	pthread_mutex_t mutex;
} async_t;


#define LIBBBLAS_EXPORT __attribute__((__visibility__("default")))


static inline __attribute__((always_inline)) void async_setup(async_t *sync) {
	sync->id = -1;
}


extern LIBBBLAS_EXPORT async_t * async_init(async_t *sync, int c, int n, int b, int initc, int log);
extern LIBBBLAS_EXPORT void async_profile(async_t *sync, float prof);
extern LIBBBLAS_EXPORT void async_fini(async_t *sync, int c);
extern LIBBBLAS_EXPORT void async_log(async_t *sync, unsigned int size, const char *id);


#endif // __ASYNC_STRUCT_H__
