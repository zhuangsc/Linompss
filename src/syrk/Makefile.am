gen_sources = syrk.c syrk_check.h syrk_setup.h 
syrk_sources = syrk_check.c syrk_config.c syrk_setup.c syrk_main.h syrk_main.c

include_HEADERS = ompss_syrk.h

noinst_LTLIBRARIES = libssyrk.la libdsyrk.la

libdsyrk_la_SOURCES = syrk_main.c 
libdsyrk_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_SYRK -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libssyrk_la_SOURCES = syrk_main.c 
libssyrk_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_SYRK -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

bin_PROGRAMS = dsyrk ssyrk

dsyrk_SOURCES = $(gen_sources) $(syrk_sources)
dsyrk_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_SYRK -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dsyrk_LDADD =  $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS)

ssyrk_SOURCES = $(gen_sources) $(syrk_sources)
ssyrk_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_SYRK -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
ssyrk_LDADD =  $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS)


lib_LTLIBRARIES = libompsssyrk.la
libompsssyrk_la_SOURCES = ompss_syrk.c ompss_syrk.h 
libompsssyrk_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompsssyrk_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
libompsssyrk_la_LIBADD = $(builddir)/libssyrk.la $(builddir)/libdsyrk.la $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/blocked_blas/libbblas.la


