lib_LTLIBRARIES = libsbblas.la 
include_HEADERS = tasks_trsm_csr.h tasks_gemv_csr.h tasks_potrf_csr.h bsblas_gemv_csr.h bsblas_csrmmb.h bsblas_csrmm.h async.h

noinst_LTLIBRARIES = libssbblas.la libdsbblas.la 

BBLAS_SOURCES = tasks_trsm_csr.c tasks_gemv_csr.c tasks_potrf_csr.c task_csrmmb.c bsblas_gemv_csr.h async.c \
				bsblas_csrmm.h task_csrmm.h task_csrmm.c

libssbblas_la_SOURCES = $(BBLAS_SOURCES)
libssbblas_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -DSINGLE_PRECISION -I$(top_srcdir)/src -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/csparse

libdsbblas_la_SOURCES = $(BBLAS_SOURCES)
libdsbblas_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -DDOUBLE_PRECISION -I$(top_srcdir)/src -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/csparse

libsbblas_la_SOURCES 	= #$(BBLAS_SOURCES)
libsbblas_la_LIBADD 	= libssbblas.la libdsbblas.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la 
libsbblas_la_CFLAGS 	= $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -I$(top_srcdir)/src -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/csparse
