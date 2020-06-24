#include "task_log.h"


#include "async_struct.h"


void log_record(struct asynct *sync, int id, event_t e, int i, float f) 
{
	int logs = sync->logs;
	if ( logs > 0 ) {
		log_t *log = sync->log;
		int logc = sync->logc;

		log[logc].id = id;
		log[logc].type = e;
		log[logc].i = i;
		log[logc++].f = f;

		if ( logc == logs ) {
			log_flush(sync->logf, log, logc);
			logc = 0;
		}

		sync->logc = logc;
	}
}

void log_locknrecord(struct asynct *sync, int id, event_t e, int i, float f) 
{
	int logs = sync->logs;
	if ( logs > 0 ) {
		log_t *log = sync->log;

		pthread_mutex_t *mutex = &sync->mutex;	
		pthread_mutex_lock(mutex);

		int logc = sync->logc;

		log[logc].id = id;
		log[logc].type = e;
		log[logc].i = i;
		log[logc++].f = f;

		if ( logc == logs ) {
			log_flush(sync->logf, log, logc);
			logc = 0;
		}

		sync->logc = logc;

		pthread_mutex_unlock(mutex);
	}
}

void log_flush(FILE *logf, log_t *log, int logc) 
{
	int cid = -2;
	int i;	
	for ( i=0; i<logc; ++i ) {
		int id = log[i].id; 
		int nok = cid==-2 ? 0 : abs(id-cid)>1;
		fprintf(logf, "%i: %i %.16e %i\n", id, log[i].type, log[i].f, log[i].i);

		if ( nok ) {
			fprintf(stderr, "err in logf\n");
			fprintf(stderr, "%i %i\n", id, cid);
		}
		
		cid = id;
	}
}
