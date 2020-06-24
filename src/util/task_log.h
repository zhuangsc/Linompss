#ifndef __TASK_LOG_H__
#define __TASK_LOG_H__


#include <stdio.h>


struct asynct;


typedef enum {
	/*0*/ EVENT_ASYNC_FRACTION, 
	/*1*/ EVENT_PRIORITY, 
	/*2*/ EVENT_1NORM, 
	/*3*/ EVENT_PARTICIPATE, 
	/*4*/ EVENT_REFUSED, 
	/*5*/ EVENT_RESIDUAL, 
	/*6*/ EVENT_ELAPSED
} event_t;


typedef struct logt {
	int id;
	int type;
	int i;
	float f;
} log_t;


extern void log_record(struct asynct* sync, int id, event_t e, int i, float f);
extern void log_locknrecord(struct asynct *sync, int id, event_t e, int i, float f);
extern void log_flush(FILE *logf, log_t *log, int logc);  


#endif // __TASK_LOG_H__
