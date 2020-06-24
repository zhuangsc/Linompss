lib_LTLIBRARIES = libbblas.la 
include_HEADERS = bblas_dot.h bblas_axpy.h bblas_gemm.h task_gemm.h tasks_potrf.h tasks_trsm.h tasks_syrk.h tasks_nested_gemm.h tasks_nested_potrf.h tasks_nested_trsm.h tasks_nested_syrk.h

noinst_LTLIBRARIES = libsbblas.la libdbblas.la 

BBLAS_SOURCES = bblas_gemv.h task_gemv.c task_gemv.h \
				bblas_axpy.h task_axpy.h task_axpy.c bblas_axpy_async.h task_axpy_async.c \
				bblas_dot.h task_dot.c bblas_dot_async.h task_dot_async.h task_dot_async.c \
				bblas_copy.h task_copy.h task_copy.c \
				bblas_trsm.h tasks_trsm.c tasks_trsm.h tasks_potrf.c tasks_potrf.h \
				bblas_gemm.h task_gemm.h task_gemm.c tasks_syrk.h tasks_syrk.c bblas_syrk.h tasks_nested_gemm.c tasks_nested_potrf.c tasks_nested_trsm.c tasks_nested_syrk.c tasks_nested_gemm.h tasks_nested_potrf.h tasks_nested_trsm.h tasks_nested_syrk.h \
				task_geqrf.h task_geqrf.c task_larfb.h task_larfb.c task_tsqrt.h task_tsqrt.c qr_assist.h qr_assist.c task_tsmqr.h task_tsmqr.c \
				task_util.c task_util.h bblas_util.h


BUTIL_SOURCES = task_convert.c bblas_convert.h 

libsbblas_la_SOURCES = $(BBLAS_SOURCES)
libsbblas_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -DSINGLE_PRECISION -I$(top_srcdir)/src -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util

libdbblas_la_SOURCES = $(BBLAS_SOURCES)
libdbblas_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -DDOUBLE_PRECISION -I$(top_srcdir)/src -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util

#libbutil_la_SOURCES = $(BUTIL_SOURCES)
#libbutil_la_CFLAGS  = $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/src/fp

libbblas_la_SOURCES 	= $(BUTIL_SOURCES) 
libbblas_la_LIBADD 		= libsbblas.la libdbblas.la $(top_builddir)/src/fp/libfp.la $(top_builddir)/src/util/libompssutil.la $(LAPACK_LIBS) 
libbblas_la_CFLAGS 		= $(OMPSS_CFLAGS) $(LAPACK_CFLAGS) -fvisibility=hidden -I$(top_srcdir)/src -I$(top_srcdir)/src/blas -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/fp -I$(top_srcdir)/src/util
