
chol_sources = chol.c chol_hcheck.c 
chol_headers = chol_setup.h chol_check.h 

chol_hll_sources = chol_hllmain.c chol_hsetup.c chol_config.c 
chol_hll_headers = chol_llmain.h chol_config.h 
chol_hrl_sources = chol_hrlmain.c chol_hsetup.c chol_config.c 
chol_hrl_headers = chol_rlmain.h chol_config.h 

#chol_ll_sources = chol_llmain.c chol_setup.c chol_config.c 
#chol_ll_headers = chol_llmain.h chol_config.h 
chol_rl_sources = chol_rlmain.c chol_setup.c chol_config.c 
chol_rl_headers = chol_rlmain.h chol_config.h 

prlchol_sources = chol_prlmain.c chol_psetup.c chol_config.c 
prlchol_headers = chol_prlmain.h chol_config.h 
chol_nll_sources = chol_nllmain.c chol_hsetup.c chol_config.c 
chol_nll_headers = chol_nllmain.h chol_config.h 
chol_nrl_sources = chol_nrlmain.c chol_hsetup.c chol_config.c 
chol_nrl_headers = chol_nrlmain.h chol_config.h 


bin_PROGRAMS = dchol_hll schol_hll dchol_ll schol_ll dchol_rl schol_rl dchol_nll schol_nll dchol_nrl schol_nrl dchol_hrl schol_hrl
include_HEADERS = ompss_dense_chol.h


# convenience libraries for the main functions
noinst_LTLIBRARIES = libdhrl.la libshrl.la libdhll.la libshll.la \
			libdrl.la libsrl.la libdll.la libsll.la \
			libdnll.la libsnll.la libdnrl.la libsnrl.la

# dense hypermatrix double left-looking 
libdhll_la_SOURCES = chol_hllmain.c chol_hllmain.h
libdhll_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
# dense hypermatrix float left-looking 
libshll_la_SOURCES = chol_hllmain.c chol_hllmain.h
libshll_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
# dense hypermatrix double right-looking 
libdhrl_la_SOURCES = chol_hrlmain.c chol_rlmain.h
libdhrl_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
# dense hypermatrix float right-looking 
libshrl_la_SOURCES = chol_hrlmain.c chol_rlmain.h
libshrl_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

# dense column-major matrix double left-looking 
libdll_la_SOURCES = chol_llmain.c chol_llmain.h
libdll_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
# dense column-major matrix float left-looking 
libsll_la_SOURCES = chol_llmain.c chol_llmain.h
libsll_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
# dense column-major matrix double right-looking 
libdrl_la_SOURCES = chol_rlmain.c chol_rlmain.h
libdrl_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
# dense column-major matrix float right-looking 
libsrl_la_SOURCES = chol_rlmain.c chol_rlmain.h
libsrl_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas

# dense double nested right-looking 
libdnrl_la_SOURCES = chol_nrlmain.c chol_nrlmain.h
libdnrl_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
# dense float nested right-looking 
libsnrl_la_SOURCES = chol_nrlmain.c chol_nrlmain.h
libsnrl_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
# dense double nested left-looking 
libdnll_la_SOURCES = chol_nllmain.c chol_nllmain.h
libdnll_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
# dense float nested left-looking 
libsnll_la_SOURCES = chol_nllmain.c chol_nllmain.h
libsnll_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
libdsll_la_CFLAGS = $(OMPSS_CFLAGS) -I$(srcdir)/sparse_kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


lib_LTLIBRARIES = libompssdchol.la
libompssdchol_la_SOURCES = ompss_dense_chol.c ompss_dense_chol.h 
libompssdchol_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompssdchol_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util
libompssdchol_la_LIBADD = $(builddir)/libdhll.la $(builddir)/libshll.la $(builddir)/libdhrl.la $(builddir)/libshrl.la \
						$(builddir)/libdll.la $(builddir)/libsll.la $(builddir)/libdrl.la $(builddir)/libsrl.la \
						$(builddir)/libdnll.la $(builddir)/libsnll.la $(builddir)/libdnrl.la $(builddir)/libsnrl.la \
						$(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

# Dense hypermatrix left-looking
dchol_hll_SOURCES = $(chol_sources) $(chol_headers) $(chol_hll_sources) $(chol_hll_headers)
dchol_hll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_hll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_HLL -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dchol_hll_LDADD = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la \
				  $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

schol_hll_SOURCES = $(chol_sources) $(chol_headers) $(chol_hll_sources) $(chol_hll_headers)
schol_hll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_hll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_HLL -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
schol_hll_LDADD = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/fp/libfp.la \
				  $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

# Dense column-major matrix left-looking
dchol_ll_SOURCES = chol.c chol_config.c chol_setup.c chol_check.c $(chol_headers) 
dchol_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) 
dchol_ll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_LL -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dchol_ll_LDADD = $(builddir)/libdll.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

schol_ll_SOURCES = chol.c chol_config.c chol_setup.c chol_check.c $(chol_headers)  
schol_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_ll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_LL -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas
schol_ll_LDADD = $(builddir)/libsll.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

# Dense column-major matrix right-looking
dchol_rl_SOURCES = chol.c chol_check.c chol_config.c chol_setup.c chol_rlmain.c 
dchol_rl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_rl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_RL -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dchol_rl_LDADD = $(builddir)/libdll.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

schol_rl_SOURCES = chol.c chol_check.c chol_config.c chol_setup.c chol_rlmain.c 
schol_rl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_rl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_RL -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
schol_rl_LDADD = $(builddir)/libsll.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la


# Dense hypermatrix right-looking
dchol_hrl_SOURCES = $(chol_sources) $(chol_headers) $(chol_hrl_sources) $(chol_hrl_headers)
dchol_hrl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_hrl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_HRL -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
dchol_hrl_LDADD = $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la

schol_hrl_SOURCES = $(chol_sources) $(chol_headers) $(chol_hrl_sources) $(chol_hrl_headers)
schol_hrl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_hrl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_HRL -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
schol_hrl_LDADD = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/util/libompssutil.la



# Dense, nested, left-looking
dchol_nll_SOURCES = $(chol_sources) $(chol_headers) $(chol_nll_sources) $(chol_nll_headers)
dchol_nll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_nll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NLL -DDOUBLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
dchol_nll_LDADD =  $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la

schol_nll_SOURCES = $(chol_sources) $(chol_headers) $(chol_nll_sources) $(chol_nll_headers)
schol_nll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_nll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NLL -DSINGLE_PRECISION -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
schol_nll_LDADD = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la


# Dense, nested, right-looking
dchol_nrl_SOURCES = $(chol_sources) $(chol_headers) $(chol_nrl_sources) $(chol_nrl_headers)
dchol_nrl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_nrl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NRL -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas
dchol_nrl_LDADD = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la

schol_nrl_SOURCES = $(chol_sources) $(chol_headers) $(chol_nrl_sources) $(chol_nrl_headers)
schol_nrl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_nrl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NRL -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/blocked_blas
schol_nrl_LDADD = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la


# Panels, instead of blocks, right-looking 
#prlchol_SOURCES = $(chol_sources) $(chol_headers) $(prlchol_sources) $(prlchol_headers)
#prlchol_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
#prlchol_CFLAGS = $(OMPSS_CFLAGS) -DUSE_PRL -I$(srcdir)/kernels -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/blas
#prlchol_LDADD = $(builddir)/kernels/libdkernels.la
