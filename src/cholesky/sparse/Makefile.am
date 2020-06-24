SUBDIRS = kernels


chol_sources = chol.c 
chol_headers = chol_setup.h chol_config.h

chols_ll_sources = chols_llmain.c chols_config.c chols_setup.c chols_warm.c
chols_ll_headers = chols_llmain.h chol_config.h chols_setup.h chols_warm.h


bin_PROGRAMS = dchols_ll schols_ll
include_HEADERS = ompss_sparse_chol.h


# convenience libraries for the main functions
noinst_LTLIBRARIES = libdsll.la libssll.la


# sparse double left-looking 
libssll_la_SOURCES = ompss_sparse_chol.c ompss_sparse_chol.h chols_llmain.c chols_llmain.h 
libssll_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(srcdir)/kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

libdsll_la_SOURCES = ompss_sparse_chol.c ompss_sparse_chol.h chols_llmain.c chols_llmain.h 
libdsll_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(srcdir)/kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


lib_LTLIBRARIES = libompssschol.la
libompssschol_la_SOURCES = 
libompssschol_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompssschol_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(srcdir)/kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp
libompssschol_la_LIBADD = $(builddir)/libdsll.la $(builddir)/libssll.la $(builddir)/kernels/libskernels.la $(top_builddir)/src/matutil/libmatutil.la

# Sparse, left-looking
dchols_ll_SOURCES = $(chol_sources) $(chol_headers) $(chols_ll_sources) $(chols_ll_headers)
dchols_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
dchols_ll_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_SLL -I$(srcdir)/kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/fp
dchols_ll_LDADD = $(builddir)/kernels/libskernels.la $(top_builddir)/src/matutil/libmatutil.la \
				  $(top_builddir)/src/util/libompssutil.la

schols_ll_SOURCES = $(chol_sources) $(chol_headers) $(chols_ll_sources) $(chols_ll_headers)
schols_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
schols_ll_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_SLL -I$(srcdir)/kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/fp
schols_ll_LDADD = $(builddir)/kernels/libskernels.la $(top_builddir)/src/matutil/libmatutil.la \
				  $(top_builddir)/src/util/libompssutil.la

