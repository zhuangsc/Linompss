gen_sources = geqr.c qr_check.h geqr_setup.h qr_config.h
geqr_sources = qr_check.c qr_config.c geqr_setup.c 


noinst_LTLIBRARIES = libsgeqrf.la libdgeqrf.la 


libdgeqrf_la_SOURCES = geqrf_main.c 
libdgeqrf_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libsgeqrf_la_SOURCES = geqrf_main.c 
libsgeqrf_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


bin_PROGRAMS = dgeqrf sgeqrf


dgeqrf_SOURCES = $(gen_sources) $(geqr_sources)
dgeqrf_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_GEQRT -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dgeqrf_LDADD =  $(builddir)/libdgeqrf.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS)

sgeqrf_SOURCES = $(gen_sources) $(geqr_sources)
sgeqrf_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_GEQRT -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
sgeqrf_LDADD =  $(builddir)/libsgeqrf.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS)


#lib_LTLIBRARIES = libompssmm.la
#libompssmm_la_SOURCES = ompss_mm.c ompss_mm.h 
#libompssmm_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#libompssmm_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas
#libompssmm_la_LIBADD = $(builddir)/libdgemm.la $(builddir)/libsgemm.la $(builddir)/libdcsrmmb.la $(builddir)/libscsrmmb.la \
#						$(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS)


