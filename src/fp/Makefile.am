noinst_LTLIBRARIES = libfp.la 

include_HEADERS = fptype.h fpmatr.h fpblas.h fpsblas.h

FP_SOURCES = fpblas.h fplapack.h fptype.c fptype.h fpsblas.h fpmatr.c

libfp_la_SOURCES = $(FP_SOURCES)
libfp_la_CPPFLAGS = -I$(top_srcdir)/src/blas
libfp_la_CFLAGS  = $(OMPSS_CFLAGS)
