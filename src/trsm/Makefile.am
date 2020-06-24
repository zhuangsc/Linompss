trsm_sources = trsm.c trsm_main.c trsm_config.c trsm_setup.c trsm_check.c
trsm_headers = trsm_main.h trsm_config.h trsm_setup.h trsm_check.h

include_HEADERS = ompss_trsm.h

noinst_LTLIBRARIES = libstrsm.la libdtrsm.la
#
libstrsm_la_SOURCES = trsm_main.c
libstrsm_la_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/util

libdtrsm_la_SOURCES = trsm_main.c
libdtrsm_la_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/util

lib_LTLIBRARIES = libompsstrsm.la
libompsstrsm_la_SOURCES = ompss_trsm.c ompss_trsm.h
libompsstrsm_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompsstrsm_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
libompsstrsm_la_LIBADD = $(builddir)/libstrsm.la $(builddir)/libdtrsm.la $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/blocked_blas/libbblas.la


bin_PROGRAMS = strsm dtrsm

strsm_SOURCES = $(trsm_sources) $(trsm_headers)
strsm_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
strsm_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil \
				   -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/util
strsm_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la \
				   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la

dtrsm_SOURCES = $(trsm_sources) $(trsm_headers)
dtrsm_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
dtrsm_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil \
				   -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/util
dtrsm_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la \
				   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la
