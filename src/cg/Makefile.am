bin_PROGRAMS = scg


scg_SOURCES = cg_main.c cg_config.c cg_setup.c 
scg_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm -lpthread
scg_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas\
				-I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse
scg_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/blocked_sblas/libsbblas.la\
				$(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/csparse/libdcsparse.la

