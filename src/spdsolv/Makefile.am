spdsolv_sources = spdsolv.c spdsolv_main.c spdsolv_config.c spdsolv_setup.c spdsolv_check.c
spdsolv_headers = spdsolv_main.h spdsolv_config.h spdsolv_setup.h spdsolv_check.h

bin_PROGRAMS = sspdsolv dspdsolv

sspdsolv_SOURCES = $(spdsolv_sources) $(spdsolv_headers)
sspdsolv_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
sspdsolv_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil \
				   -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas \
				   -I$(top_srcdir)/src/cholesky/dense
sspdsolv_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la \
				   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/cholesky/dense/libompssdchol.la \
				   $(top_builddir)/src/util/libompssutil.la

dspdsolv_SOURCES = $(spdsolv_sources) $(spdsolv_headers)
dspdsolv_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
dspdsolv_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil \
				   -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas \
				   -I$(top_srcdir)/src/cholesky/dense
dspdsolv_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la \
				   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/cholesky/dense/libompssdchol.la \
				   $(top_builddir)/src/util/libompssutil.la
