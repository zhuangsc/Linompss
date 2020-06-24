cg_sources_common = cg.c
cg_headers_common = cg_setup.h cg_config.h 


bin_PROGRAMS = #scgprof dcgprof scgpipe dcgpipe scgmod1 dcgmod1 
#bin_PROGRAMS = scg dcg scgs dcgs scgsched scgsched2 scgmod1 dcgmod1 scgmod1s dcgmod1s


noinst_LTLIBRARIES = libscgprof.la libsscgprof.la libsdcgprof.la libdcgprof.la
#					libscgpipe.la libdcgpipe.la \
#					libscgmod1.la libdcgmod1.la \
#					libscg.la libdcg.la 


### DENSE
libscgprof_la_SOURCES = cgprof_main.c cgprof_main.h cgas_workspace
libscgprof_la_CFLAGS =  $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse
libdcgprof_la_SOURCES = cgprof_main.c cgprof_main.h cgas_workspace
libdcgprof_la_CFLAGS =  $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse

### SPARSE
libsscgprof_la_SOURCES = cgprof_main.c cgprof_main.h cgas_workspace
libsscgprof_la_CFLAGS =  $(OMPSS_CFLAGS) -DUSE_SPARSE -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse
libsdcgprof_la_SOURCES = cgprof_main.c cgprof_main.h cgas_workspace
libsdcgprof_la_CFLAGS =  $(OMPSS_CFLAGS) -DUSE_SPARSE -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/csparse


#libscgmod1_la_SOURCES = cgmod1_main.c cgmod1_main.h cgmod1_workspace
#libscgmod1_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
#libdcgmod1_la_SOURCES = cgmod1_main.c cgmod1_main.h cgmod1_workspace
#libdcgmod1_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
#
#libscgpipe_la_SOURCES = cgpipe_main.c cgpipe_main.h cgpipe_workspace
#libscgpipe_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
#libdcgpipe_la_SOURCES = cgpipe_main.c cgpipe_main.h cgpipe_workspace
#libdcgpipe_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
#
#libscg_la_SOURCES = cg_main.c cg_main.h cgas_workspace.h
#libscg_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -DUSE_DENSE -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
#libdcg_la_SOURCES = cg_main.c cg_main.h cgas_workspace
#libdcg_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -DUSE_DENSE -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas


# dense
#dcg_SOURCES = cg_config.c cg_setup.c cg_main.c cgas_workspace.h $(cg_sources_common) $(cg_headers_common) 
#dcg_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#dcg_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CG -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp \
#				-I$(top_srcdir)/src/blas
#dcg_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la
#
#scg_SOURCES = cg_config.c cg_setup.c cg_main.c cgas_workspace.h $(cg_sources_common) $(cg_headers_common) 
#scg_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scg_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CG -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
#			-I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas
#scg_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la


# dense, self-scheduling
#scgsched_SOURCES = cg_config.c cg_setup.c cgsched_main.c cgsched_main.h cgas_workspace.h $(cg_sources_common) $(cg_headers_common) 
#scgsched_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scgsched_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGSCHED -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas 
#scgsched_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la
#
## dense, profiling
#scgprof_SOURCES = cg_config.c cg_setup.c cgprof_main.c cgprof_main.h cgas_workspace.h $(cg_sources_common) $(cg_headers_common) 
#scgprof_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scgprof_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGPROF -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas 
#scgprof_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la
#
#dcgprof_SOURCES = cg_config.c cg_setup.c cgprof_main.c cgprof_main.h cgas_workspace.h $(cg_sources_common) $(cg_headers_common) 
#dcgprof_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#dcgprof_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGPROF -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas 
#dcgprof_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la




#scgpipe_SOURCES = cg_config.c cg_setup.c cgpipe_main.c cgpipe_main.h cgpipe_workspace.h $(cg_sources_common) $(cg_headers_common) 
#scgpipe_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scgpipe_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGPIPE -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas 
#scgpipe_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la
#
#dcgpipe_SOURCES = cg_config.c cg_setup.c cgpipe_main.c cgpipe_main.h cgpipe_workspace.h $(cg_sources_common) $(cg_headers_common) 
#dcgpipe_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#dcgpipe_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGPIPE -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas 
#dcgpipe_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la




#scgmod1_SOURCES = cg_config.c cg_setup.c $(cg_sources_common) $(cg_headers_common) 
#scgmod1_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scgmod1_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGMOD1 -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas \
#				-I$(top_srcdir)/src/blocked_blas 
#scgmod1_LDADD   = $(builddir)/libscgmod1.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la

#dcgmod1_SOURCES = cg_config.c cg_setup.c $(cg_sources_common) $(cg_headers_common) 
#dcgmod1_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#dcgmod1_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGMOD1 -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas \
#                -I$(top_srcdir)/src/blocked_blas
#dcgmod1_LDADD   = $(builddir)/libdcgmod1.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la


# sparse
#scgs_SOURCES = cg.c cgs_config.c cg_setup.c cg_main.h cg_main.c cgas_workspace.h
#scgs_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CG -DUSE_SPARSE -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil \
#				-I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas 
#scgs_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scgs_LDADD   = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/blocked_sblas/libsbblas.la
#
#dcgs_SOURCES = cg.c cgs_config.c cg_setup.c cg_main.h cg_main.c cgas_workspace.h
#dcgs_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CG -DUSE_SPARSE -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/matutil \
#				-I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas 
#dcgs_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#dcgs_LDADD   = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/blocked_sblas/libsbblas.la


# sparse, different version
#scgmod1s_SOURCES = cgs_config.c cg_setup.c cgmod1_main.c cgmod1_main.h cgmod1_workspace.h $(cg_sources_common) $(cg_headers_common) 
#scgmod1s_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#scgmod1s_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGMOD1 -DUSE_SPARSE -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas \
#				-I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas
#scgmod1s_LDADD   = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/blocked_sblas/libsbblas.la 
#
#dcgmod1s_SOURCES = cg_config.c cg_setup.c cgmod1_main.c cgmod1_main.h cgmod1_workspace.h $(cg_sources_common) $(cg_headers_common) 
#dcgmod1s_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
#dcgmod1s_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DUSE_CGMOD1 -DUSE_SPARSE -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas \
#                -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas
#dcgmod1s_LDADD   = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/blocked_sblas/libsbblas.la 


#itref_SOURCES = itref.c itref_config.c itref_setup.c cg_setup.h cgmod1_main.c cgmod1_main.h cg_main.c cg_main.h cgprof_main.c cgprof_main.h cgas_workspace.h
#itref_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm -lpthread
#itref_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil \
#                -I$(top_srcdir)/src/util -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blas 
#itref_LDADD   = $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la 


#libscg_la_SOURCES = cg_main.c cg_main.h cg_workspace.h \
#					cgmod1_main.c cgmod1_main.h cgmod1_workspace.h 
#libscg_la_CFLAGS = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util
#libscg_la_LIBADD = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/blocked_sblas/libsbblas.la 
#libdcg_la_SOURCES = cg_main.c cg_main.h cg_workspace.h \
#					cgmod1_main.c cgmod1_main.h cgmod1_workspace.h 
#libdcg_la_CFLAGS = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util
#libdcg_la_LIBADD = $(top_builddir)/src/matutil/libmatutil.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/blocked_sblas/libsbblas.la 


lib_LTLIBRARIES = libompsscg.la libompssscg.la
libompsscg_la_SOURCES = ompss_cg.c ompss_cg.h
libompsscg_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompsscg_la_CFLAGS =  $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
					   -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/csparse
libompsscg_la_LIBADD = libscgprof.la libdcgprof.la \
					   $(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) \
					   $(top_builddir)/src/csparse/libscsparse.la

libompssscg_la_SOURCES = ompss_cg.c ompss_cg.h
libompssscg_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompssscg_la_CFLAGS =  $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DUSE_SPARSE -DLIBOMPSS_BUILDING -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util \
						-I$(top_srcdir)/src/fp  -I$(top_srcdir)/src/blocked_blas -I$(top_srcdir)/src/blocked_sblas -I$(top_srcdir)/src/csparse
libompssscg_la_LIBADD = libsscgprof.la libsdcgprof.la \
					   $(top_builddir)/src/blocked_blas/libbblas.la  $(top_builddir)/src/blocked_sblas/libsbblas.la \
					   $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS) \
					   $(top_builddir)/src/csparse/libscsparse.la


#libompsscg_la_LIBADD = libscgprof.la libdcgprof.la libscgpipe.la libdcgpipe.la libscgmod1.la libdcgmod1.la libscg.la libdcg.la \
#						$(top_builddir)/src/blocked_blas/libbblas.la $(top_builddir)/src/util/libompssutil.la $(top_builddir)/src/fp/libfp.la $(LAPACK_LIBS)

