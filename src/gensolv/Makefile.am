gensolv_sources = gensolv.c gensolv_main.c gensolv_config.c gensolv_setup.c gensolv_check.c
gensolv_headers = gensolv_main.h gensolv_config.h gensolv_setup.h gensolv_check.h

bin_PROGRAMS = sgensolv dgensolv

sgensolv_SOURCES = $(gensolv_sources) $(gensolv_headers)
sgensolv_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
sgensolv_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil \
				   -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/lu
sgensolv_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la \
				   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/lu/libompsslu.la

dgensolv_SOURCES = $(gensolv_sources) $(gensolv_headers)
dgensolv_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
dgensolv_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil \
				   -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/lu
dgensolv_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la \
				   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/lu/libompsslu.la
