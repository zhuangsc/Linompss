gen_sources = matmul.c matmul_check.h matmul_setup.h 
gemm_sources = mm_check.c mm_config.c mm_setup.c mm_main.h
csrmmb_sources = csrmmb_config.c csrmmb_setup.c csrmmb_check.c 
csrmm_sources = csrmm_config.c csrmm_setup.c csrmm_check.c 

include_HEADERS = ompss_mm.h

noinst_LTLIBRARIES = libsgemm.la libdgemm.la \
					libscsrmmb.la libdcsrmmb.la \
					libscsrmm.la libdcsrmm.la 


libdgemm_la_SOURCES = mm_main.c 
libdgemm_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libsgemm_la_SOURCES = mm_main.c 
libsgemm_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_GEMM -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libdcsrmmb_la_SOURCES = csrmmb_main.c
libdcsrmmb_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_CSRMMB -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libscsrmmb_la_SOURCES = csrmmb_main.c
libscsrmmb_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_CSRMMB -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libdcsrmm_la_SOURCES = csrmm_main.c
libdcsrmm_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_CSRMM -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libscsrmm_la_SOURCES = csrmm_main.c
libscsrmm_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_CSRMM -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


bin_PROGRAMS = dmatmul smatmul 
if USE_MKL
bin_PROGRAMS += dcsrmmb scsrmmb scsrmm
endif


dmatmul_SOURCES = $(gen_sources) $(gemm_sources)
dmatmul_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_GEMM -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dmatmul_LDADD =  $(builddir)/libdgemm.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la

smatmul_SOURCES = $(gen_sources) $(gemm_sources)
smatmul_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_GEMM -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
smatmul_LDADD = $(builddir)/libsgemm.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la

dcsrmmb_SOURCES =  $(gen_sources) $(csrmmb_sources)
dcsrmmb_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_CSRMMB -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dcsrmmb_LDADD =  $(builddir)/libdcsrmmb.la $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la 

scsrmmb_SOURCES = $(gen_sources) $(csrmmb_sources)
scsrmmb_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_CSRMMB -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
scsrmmb_LDADD =  $(builddir)/libscsrmmb.la $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la

dcsrmm_SOURCES =  $(gen_sources) $(csrmm_sources)
dcsrmm_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_CSRMM -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dcsrmm_LDADD =  $(builddir)/libdcsrmm.la $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la

scsrmm_SOURCES = $(gen_sources) $(csrmm_sources)
scsrmm_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_CSRMM -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
scsrmm_LDADD =  $(builddir)/libscsrmm.la $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la




lib_LTLIBRARIES = libompssmm.la
libompssmm_la_SOURCES = ompss_mm.c ompss_mm.h 
libompssmm_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompssmm_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
libompssmm_la_LIBADD = $(builddir)/libdgemm.la $(builddir)/libsgemm.la $(builddir)/libdcsrmmb.la $(builddir)/libscsrmmb.la  $(top_builddir)/src/blocked_sblas/libsbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/blocked_blas/libbblas.la


