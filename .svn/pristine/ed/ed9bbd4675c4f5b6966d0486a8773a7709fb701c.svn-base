SUBDIRS = dense_kernels sparse_kernels matutil


chol_sources = chol.c chol_utils.c chol_check.c fptype.c
chol_headers = chol_setup.h chol_utils.h chol_check.h fptype.h
chol_ll_sources = chol_llmain.c chol_setup.c chol_config.c 
chol_ll_headers = chol_llmain.h chol_config.h 
chol_rl_sources = chol_rlmain.c chol_setup.c chol_config.c 
chol_rl_headers = chol_rlmain.h chol_config.h 
prlchol_sources = chol_prlmain.c chol_psetup.c chol_config.c 
prlchol_headers = chol_prlmain.h chol_config.h 
chol_nll_sources = chol_nllmain.c chol_setup.c chol_config.c 
chol_nll_headers = chol_nllmain.h chol_config.h 
chol_nrl_sources = chol_nrlmain.c chol_setup.c chol_config.c 
chol_nrl_headers = chol_nrlmain.h chol_config.h 
chols_ll_sources = chols_llmain.c chols_config.c chols_setup.c chols_warm.c
chols_ll_headers = chols_llmain.h chol_config.h chols_setup.h chols_warm.h


bin_PROGRAMS = dchol_ll schol_ll dchol_rl schol_rl dchol_nll schol_nll dchol_nrl schol_nrl chols_ll

include_HEADERS = ompss_cholesky.h


# convenience libraries for the main functions
noinst_LTLIBRARIES = libdll.la libsll.la libdrl.la libsrl.la \
			libdnll.la libsnll.la libdnrl.la libsnrl.la \
			libdsll.la
# dense double left-looking 
libdll_la_SOURCES = chol_llmain.c chol_llmain.h
libdll_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense float left-looking 
libsll_la_SOURCES = chol_llmain.c chol_llmain.h
libsll_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense double right-looking 
libdrl_la_SOURCES = chol_rlmain.c chol_rlmain.h
libdrl_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense float right-looking 
libsrl_la_SOURCES = chol_rlmain.c chol_rlmain.h
libsrl_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense double nested right-looking 
libdnrl_la_SOURCES = chol_nrlmain.c chol_nrlmain.h
libdnrl_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense float nested right-looking 
libsnrl_la_SOURCES = chol_nrlmain.c chol_nrlmain.h
libsnrl_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense double nested left-looking 
libdnll_la_SOURCES = chol_nllmain.c chol_nllmain.h
libdnll_la_CFLAGS = $(OMPSS_CFLAGS) -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# dense float nested left-looking 
libsnll_la_SOURCES = chol_nllmain.c chol_nllmain.h
libsnll_la_CFLAGS = $(OMPSS_CFLAGS) -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
# sparse double left-looking 
libdsll_la_SOURCES = chols_llmain.c chols_llmain.h
libdsll_la_CFLAGS = $(OMPSS_CFLAGS) -I$(srcdir)/sparse_kernels -I$(srcdir)/matutil -I$(srcdir)/util


lib_LTLIBRARIES = libompsschol.la
libompsschol_la_SOURCES = ompss_cholesky.c ompss_cholesky.h fptype.c
libompsschol_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
libompsschol_la_CFLAGS = $(OMPSS_CFLAGS) $(CFLAG_VISIBILITY) -DLIBOMPSS_BUILDING -I$(srcdir)/dense_kernels -I$(srcdir)/sparse_kernels -I$(srcdir)/matutil -I$(srcdir)/util
libompsschol_la_LIBADD = $(builddir)/libdll.la $(builddir)/libsll.la $(builddir)/libdrl.la $(builddir)/libsrl.la \
			$(builddir)/libdnll.la $(builddir)/libsnll.la $(builddir)/libdnrl.la $(builddir)/libsnrl.la \
			$(builddir)/libdsll.la \
			$(builddir)/dense_kernels/libdndkernels.la $(builddir)/dense_kernels/libsndkernels.la \
			$(builddir)/dense_kernels/libddkernels.la $(builddir)/dense_kernels/libsdkernels.la \
			$(builddir)/sparse_kernels/libskernels.la \
			$(builddir)/matutil/libmatutil.la

# Dense, left-looking
dchol_ll_SOURCES = $(chol_sources) $(chol_headers) $(chol_ll_sources) $(chol_ll_headers)
dchol_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_ll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_LL -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
dchol_ll_LDADD = $(builddir)/dense_kernels/libddkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libdgenmat.la

schol_ll_SOURCES = $(chol_sources) $(chol_headers) $(chol_ll_sources) $(chol_ll_headers)
schol_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_ll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_LL -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
schol_ll_LDADD = $(builddir)/dense_kernels/libsdkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libsgenmat.la


# Dense, right-looking
dchol_rl_SOURCES = $(chol_sources) $(chol_headers) $(chol_rl_sources) $(chol_rl_headers)
dchol_rl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_rl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_RL -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
dchol_rl_LDADD = $(builddir)/dense_kernels/libddkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libdgenmat.la

schol_rl_SOURCES = $(chol_sources) $(chol_headers) $(chol_rl_sources) $(chol_rl_headers)
schol_rl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_rl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_RL -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
schol_rl_LDADD = $(builddir)/dense_kernels/libsdkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libsgenmat.la


# Dense, nested, left-looking
dchol_nll_SOURCES = $(chol_sources) $(chol_headers) $(chol_nll_sources) $(chol_nll_headers)
dchol_nll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_nll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NLL -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
dchol_nll_LDADD = $(builddir)/dense_kernels/libdndkernels.la $(builddir)/dense_kernels/libddkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libdgenmat.la

schol_nll_SOURCES = $(chol_sources) $(chol_headers) $(chol_nll_sources) $(chol_nll_headers)
schol_nll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_nll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NLL -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
schol_nll_LDADD = $(builddir)/dense_kernels/libsndkernels.la $(builddir)/dense_kernels/libsdkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libsgenmat.la


# Dense, nested, right-looking
dchol_nrl_SOURCES = $(chol_sources) $(chol_headers) $(chol_nrl_sources) $(chol_nrl_headers)
dchol_nrl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
dchol_nrl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NRL -DDOUBLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
dchol_nrl_LDADD = $(builddir)/dense_kernels/libdndkernels.la $(builddir)/dense_kernels/libddkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libdgenmat.la

schol_nrl_SOURCES = $(chol_sources) $(chol_headers) $(chol_nrl_sources) $(chol_nrl_headers)
schol_nrl_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
schol_nrl_CFLAGS = $(OMPSS_CFLAGS) -DUSE_NRL -DSINGLE_PRECISION -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
schol_nrl_LDADD = $(builddir)/dense_kernels/libsndkernels.la $(builddir)/dense_kernels/libsdkernels.la $(builddir)/matutil/libmatutil.la $(builddir)/matutil/libsgenmat.la


# Panels, instead of blocks, right-looking 
prlchol_SOURCES = $(chol_sources) $(chol_headers) $(prlchol_sources) $(prlchol_headers)
prlchol_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)
prlchol_CFLAGS = $(OMPSS_CFLAGS) -DUSE_PRL -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util
prlchol_LDADD = $(builddir)/dense_kernels/libdkernels.la

# Sparse, left-looking
chols_ll_SOURCES = $(chol_sources) $(chol_headers) $(chols_ll_sources) $(chols_ll_headers)
chols_ll_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) #-L/home/szhuang/mnt/Experiments/algebra/mm/build/lib -lompssgemm
chols_ll_CFLAGS = $(OMPSS_CFLAGS) -DUSE_SLL -I$(srcdir)/sparse_kernels -I$(srcdir)/dense_kernels -I$(srcdir)/matutil -I$(srcdir)/util #-I/home/szhuang/mnt/Experiments/algebra/mm/build/include
chols_ll_LDADD = $(builddir)/sparse_kernels/libskernels.la $(builddir)/matutil/libmatutil.la
