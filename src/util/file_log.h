#ifndef __FILE_LOG_H__
#define __FILE_LOG_H__

#include <sys/time.h>
#include <sys/utsname.h>

typedef struct logf {
	char log_name[32];
	unsigned long elapsed;
	struct timeval start;
	struct timeval end;
	int rep;
	int num_threads;
	struct utsname info;

} log_f;

void logf_init(log_f *log, int rep, char *name);
void logf_record_stime(log_f *log);
void logf_record_edimt(log_f *log);
void logf_record_elapse(log_f *log);
void logf_append(log_f *log, char *info);
void logf_dump(log_f *log);


#endif //__FILE_LOG_H__
