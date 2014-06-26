noinst_LTLIBRARIES = libddkernels.la libsdkernels.la libdndkernels.la libsndkernels.la 


# Single precision, Dense
libsdkernels_la_SOURCES = chol_kernels.c chol_kernels.h
libsdkernels_la_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src

# Double precision, Dense
libddkernels_la_SOURCES = chol_kernels.c chol_kernels.h
libddkernels_la_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src


# Single precision, Nested, Dense
libsndkernels_la_SOURCES = chol_kernels_nested.c chol_kernels_nested.h
libsndkernels_la_CFLAGS  = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(top_srcdir)/src

# Double precision, Nested, Dense
libdndkernels_la_SOURCES = chol_kernels_nested.c chol_kernels_nested.h
libdndkernels_la_CFLAGS  = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(top_srcdir)/src
