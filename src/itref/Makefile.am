bin_PROGRAMS = itref sitref


itref_SOURCES = itref.c itref_config.c itref_setup.c 
itref_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm -lpthread
itref_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil \
                -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse
itref_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la \
				$(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/csparse/libdcsparse.la

sitref_SOURCES = itref.c itref_config.c itref_setup.c 
sitref_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm -lpthread
sitref_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_SPARSE -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas\
				-I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse
sitref_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/blocked_sblas/libsbblas.la\
				$(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/csparse/libdcsparse.la


#cgbaust_SOURCES = cg.c cg_bau.c $(cg_sources_common) $(cg_headers_common) 
#cgbaust_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#cgbaust_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) $(OMPSS_TRACE_FLAGS) -DUSE_CGBAU -DOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil 
#cgbaust_LDADD   = $(srcdir)/dense_kernels/libdbblas.a $(srcdir)/matutil/libmatutil.a 

#noinst_LTLIBRARIES = libscg.la 
#libscg_la_SOURCES = cg_main.c cg_main.h blas.h cgas_workspace.h as_man.c as_man.h
#libscg_la_CFLAGS = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/fp
#libscg_la_LIBADD = $(builddir)/fp/libfp.la $(builddir)/dense_kernels/libsbblas.la

