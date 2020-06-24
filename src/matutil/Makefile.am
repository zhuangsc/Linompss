noinst_LTLIBRARIES = libsmatutil.la libdmatutil.la
lib_LTLIBRARIES = libmatutil.la

include_HEADERS = hb.h

MATUTILFP = matfprint.c matfprint.h \
			matfread.c matfread.h \
			genmat.c genmat.h \
			densutil.h densutil.c \
			hbdebug.h hbdebug.c \
			hbconvrt.h hbconvrt.c \
			symfac.h symfac.c \
			hbpad.h hbpad.c \
			iohb.h iohb.c \
			hb.h hb.c \
			genmat_config.h genmat_config.c 


libsmatutil_la_SOURCES = $(MATUTILFP)
libsmatutil_la_CPPFLAGS  = -I$(top_srcdir)/src -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas
libsmatutil_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DSINGLE_PRECISION

libdmatutil_la_SOURCES = $(MATUTILFP)
libdmatutil_la_CPPFLAGS  = -I$(top_srcdir)/src -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas
libdmatutil_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DDOUBLE_PRECISION


libmatutil_la_SOURCES   =   convertmat.c convertmat.h hbmalloc.c hbext.c hbext.h 
libmatutil_la_CPPFLAGS  =   -I$(top_srcdir)/src -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas
libmatutil_la_CFLAGS    =   $(OMPSS_CFLAGS) $(LAPACK_CFLAGS)
libmatutil_la_LIBADD    =   libsmatutil.la libdmatutil.la $(top_builddir)/src/fp/libfp.la
