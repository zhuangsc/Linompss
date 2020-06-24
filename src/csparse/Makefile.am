#noinst_LTLIBRARIES = libscsparse.la libdcsparse.la
lib_LTLIBRARIES = libscsparse.la libdcsparse.la

include_HEADERS = scsparse.h dcsparse.h

libscsparse_la_SOURCES = scsparse.c scsparse.h
libscsparse_la_CPPFLAGS = 
libscsparse_la_CFLAGS  = $(OMPSS_CFLAGS)

libdcsparse_la_SOURCES = dcsparse.c dcsparse.h
libdcsparse_la_CPPFLAGS = 
libdcsparse_la_CFLAGS  = $(OMPSS_CFLAGS)
