#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "file_log.h"


void logf_init(log_f *log, int rep, char *name)
{

	uname(&log->info);
	snprintf(log->log_name, 32, "ompss-%s.stats", name);
	log->elapsed = 0;
	log->rep = rep;
	/* Zhuang, this should be checked / determined by the environment that invokes the app, not by the app itself */
//	log->num_threads = omp_get_max_threads();

}

void logf_record_stime(log_f *log)
{
	gettimeofday(&log->start, NULL);
}

void logf_record_etime(log_f *log)
{
	gettimeofday(&log->end, NULL);
	logf_record_elapse(log);
}

void logf_record_elapse(log_f *log)
{
	log->elapsed += (log->end.tv_sec - log->start.tv_sec) * 1e6 + (log->end.tv_usec - log->start.tv_usec);
}

void logf_append(log_f *log, char *info)
{
	FILE *logf = fopen(log->log_name, "a");
	fprintf(logf, "%s\n", info);
	fclose(logf);
}

void logf_dump(log_f *log)
{
	FILE *logf = fopen(log->log_name, "w");
	fprintf(logf, "%s %s %s %s\n", log->info.nodename, log->info.sysname, log->info.release, log->info.machine);
	fprintf(logf, "time: %.5f ms\n", log->elapsed/((double)log->rep * 1e3));
	fprintf(logf, "reps: %d\n", log->rep);
	fprintf(logf, "threads: %d\n", log->num_threads);
	fclose(logf);
}

