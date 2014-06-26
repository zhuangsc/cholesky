noinst_LTLIBRARIES = libmatutil.la libsgenmat.la libdgenmat.la

include_HEADERS = hb.h

libsgenmat_la_SOURCES = genmat.c genmat.h
libsgenmat_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src

libdgenmat_la_SOURCES = genmat.c genmat.h
libdgenmat_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src

libmatutil_la_SOURCES = hb.c  hbconvrt.c  hbconvrt.h  hbdebug.c  hbdebug.h  hbext.c  hbext.h  hb.h  iohb.c  symfac.c symfac.h iohb.h matfprint.c matfprint.h hbpad.c hbpad.h
libmatutil_la_CFLAGS  = $(OMPSS_CFLAGS) -I$(top_srcdir)/src/util -I$(top_srcdir)/src
