jacobi_sources = jacobi.c jacobi_config.c jacobi_setup.c jacobi_main.c jacobi_check.c
jacobi_headers = jacobi_config.h jacobi_setup.h jacobi_main.h jacobi_check.h

bin_PROGRAMS = 
if USE_MKL
bin_PROGRAMS += djacobi sjacobi
endif

sjacobi_SOURCES = $(jacobi_sources) $(jacobi_headers)
sjacobi_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
sjacobi_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
sjacobi_LDADD = $(top_builddir)/src/matutil/libmatutil.la \
				$(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/fp/libfp.la \
				$(top_builddir)/src/blocked_blas/libbblas.la \
				$(top_builddir)/src/util/libompssutil.la

djacobi_SOURCES = $(jacobi_sources) $(jacobi_headers)
djacobi_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
djacobi_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
djacobi_LDADD = $(top_builddir)/src/matutil/libmatutil.la \
				$(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/fp/libfp.la \
				$(top_builddir)/src/blocked_blas/libbblas.la \
				$(top_builddir)/src/util/libompssutil.la


noinst_LTLIBRARIES = libsjacobi.la libdjacobi.la


libsjacobi_la_SOURCES = ompss_sparse_jacobi.c ompss_sparse_jacobi.h jacobi_main.c jacobi_main.h
libsjacobi_la_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libdjacobi_la_SOURCES = jacobi_main.c jacobi_main.h ompss_sparse_jacobi.c ompss_sparse_jacobi.h 
libdjacobi_la_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


lib_LTLIBRARIES = libompssjacobi.la

libompssjacobi_la_SOURCES = 
libompssjacobi_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompssjacobi_la_CFLAGS  = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING  -I$(top_srcdir)/src/matutil \
							-I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas \
							-I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
libompssjacobi_la_LIBADD  = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/blocked_sblas/libsbblas.la \
							$(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la \
							$(top_builddir)/src/util/libompssutil.la $(builddir)/libsjacobi.la $(builddir)/libdjacobi.la
