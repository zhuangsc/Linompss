lib_LTLIBRARIES = libasync.la

include_HEADERS = as_man.h async_struct.h selfsched.h prof.h task_log.h

libasync_la_SOURCES 	= async_struct.c async_struct.h as_man.c as_man.h selfsched.c selfsched.h task_log.c task_log.h 
libasync_la_CFLAGS 		= $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -fvisibility-inlines-hidden -I$(top_srcdir)/src/fp 
