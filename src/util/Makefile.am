lib_LTLIBRARIES = libompssutil.la 
noinst_LTLIBRARIES = libsutil.la libdutil.la 

UTIL_SOURCES = array.h vector.h \
				selfsched.h selfsched.c \
				async_struct.h async_struct.c \
				prof.c prof.h \
				file_log.c file_log.h \
				task_log.c task_log.h \
				strutil.c strutil.h \
				ompss_apps.c ompss_app.h \
				ompss_options.c ompss_options.h

libsutil_la_SOURCES = as_man.c as_man.h
libsutil_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/fp
libdutil_la_SOURCES = as_man.c as_man.h
libdutil_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp


libompssutil_la_SOURCES 	= $(UTIL_SOURCES)
libompssutil_la_LIBADD 		= -lpthread libsutil.la libdutil.la
libompssutil_la_CFLAGS 		= $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/fp 
