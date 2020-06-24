noinst_LTLIBRARIES = libskernels.la libsskernels.la libdskernels.la

libsskernels_la_SOURCES = chols_kernels_csc.c chols_kernels_csr.c chols_kernels.h
libsskernels_la_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util  -I$(top_srcdir)/src/fp

libdskernels_la_SOURCES = chols_kernels_csc.c chols_kernels_csr.c chols_kernels.h
libdskernels_la_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util  -I$(top_srcdir)/src/fp

libskernels_la_SOURCES = #chols_kernels.c chols_kernels_upper.c chols_kernels.h
libskernels_la_CFLAGS  = $(OMPSS_CFLAGS) -I$(top_srcdir)/src/matutil -I$(top_srcdir)/src/util  -I$(top_srcdir)/src/fp
libskernels_la_LIBADD  = libsskernels.la libdskernels.la
