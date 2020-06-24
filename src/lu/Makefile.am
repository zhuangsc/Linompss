base_sources = lu.c lu_kernels.c lu_check.c lu_setup.c lu_config.c
base_headers = lu_kernels.h lu_check.h lu_setup.h lu_config.h

include_HEADERS = ompss_lu.h

noinst_LTLIBRARIES = libompssslu.la libompssdlu.la

libompssslu_la_SOURCES = lu_main.c lu_kernels.c lull_main.c
libompssslu_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/fp

libompssdlu_la_SOURCES = lu_main.c lu_kernels.c lull_main.c
libompssdlu_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp

lib_LTLIBRARIES = libompsslu.la
libompsslu_la_SOURCES = ompss_lu.c ompss_lu.h
libompsslu_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompsslu_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAGS_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/util
libompsslu_la_LIBADD = $(builddir)/libompssslu.la $(builddir)/libompssdlu.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la

bin_PROGRAMS = slu dlu dlull slull


slu_SOURCES = lu_main.c lu_main.h $(base_sources) $(base_headers)
slu_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
slu_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_LAPACK_STYLE -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matutil
slu_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS) $(top_builddir)/src/matutil/libmatutil.la

dlu_SOURCES = lu_main.c lu_main.h $(base_sources) $(base_headers)
dlu_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dlu_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_LAPACK_STYLE -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matutil
dlu_LDADD   = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS) $(top_builddir)/src/matutil/libmatutil.la

slull_SOURCES = lull_main.c lull_main.h $(base_sources) $(base_headers)
slull_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
slull_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_LL -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matutil
slull_LDADD  = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS) $(top_builddir)/src/matutil/libmatutil.la

dlull_SOURCES = lull_main.c lull_main.h $(base_sources) $(base_headers)
dlull_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dlull_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_LL -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matutil
dlull_LDADD  = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS) $(top_builddir)/src/matutil/libmatutil.la

#lut_SOURCES = lu_main.c lu_main.h $(base_sources) $(base_headers)
#lut_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#lut_CFLAGS = $(OMPSS_CFLAGS) --instrument -DUSE_LAPACK_STYLE

#slurecurs_SOURCES = lurecurs_main.c lurecurs_main.h $(base_sources) $(base_headers)
#slurecurs_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#slurecurs_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_RECURS -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util
#slurecurs_LDADD  = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS)
#
#dlurecurs_SOURCES = lurecurs_main.c lurecurs_main.h $(base_sources) $(base_headers)
#dlurecurs_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#dlurecurs_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_RECURS -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util
#dlurecurs_LDADD  = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS)

#lurecurst_SOURCES = lurecurs_main.c lurecurs_main.h $(base_sources) $(base_headers)
#lurecurst_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#lurecurst_CFLAGS = $(OMPSS_CFLAGS) --instrument -DUSE_RECURS

#lullt_SOURCES = lull_main.c lull_main.h $(base_sources) $(base_headers)
#lullt_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#lullt_CFLAGS = $(OMPSS_CFLAGS) --instrument -DUSE_LL
