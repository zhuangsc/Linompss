gauss_seidel_sources = gauss-seidel.c gauss-seidel_config.c gauss-seidel_setup.c gauss-seidel_main.c gauss-seidel_check.c
gauss_seidel_headers = gauss-seidel_config.h gauss-seidel_setup.h gauss-seidel_main.h gauss-seidel_check.h

bin_PROGRAMS = 
if USE_MKL
bin_PROGRAMS += sgauss-seidel dgauss-seidel
endif

sgauss_seidel_SOURCES = $(gauss_seidel_sources) $(gauss_seidel_headers)
sgauss_seidel_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
sgauss_seidel_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
sgauss_seidel_LDADD = $(top_builddir)/src/matutil/libmatutil.la \
				$(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/fp/libfp.la \
				$(top_builddir)/src/blocked_blas/libbblas.la

dgauss_seidel_SOURCES = $(gauss_seidel_sources) $(gauss_seidel_headers)
dgauss_seidel_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
dgauss_seidel_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dgauss_seidel_LDADD = $(top_builddir)/src/matutil/libmatutil.la \
				$(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/fp/libfp.la \
				$(top_builddir)/src/blocked_blas/libbblas.la



noinst_LTLIBRARIES = libsgs.la libdgs.la

libsgs_la_SOURCES = gauss-seidel_main.c gauss-seidel_main.h ompss_sparse_gs.c ompss_sparse_gs.h
libsgs_la_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libdgs_la_SOURCES = gauss-seidel_main.c gauss-seidel_main.h ompss_sparse_gs.c ompss_sparse_gs.h
libdgs_la_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
				 -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


lib_LTLIBRARIES = libompssgs.la

libompssgs_la_SOURCES =
libompssgs_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompssgs_la_CFLAGS  = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING  -I$(top_srcdir)/src/matutil \
							-I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/blocked_blas \
							-I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
libompssgs_la_LIBADD  = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/blocked_sblas/libsbblas.la \
							$(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la \
							$(top_builddir)/src/util/libompssutil.la $(builddir)/libsgs.la $(builddir)/libdgs.la
