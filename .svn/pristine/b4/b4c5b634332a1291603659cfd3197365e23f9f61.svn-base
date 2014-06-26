#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([cholesky],[1.0],[pieter.bellens@bsc.es])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/chol.c])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_FILES([Makefile 
		src/Makefile
		src/dense_kernels/Makefile
		src/sparse_kernels/Makefile
		src/matutil/Makefile])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC([mcc cc])

# Checks for libraries.
m4_include([m4/m4_ax_blas.m4])
AX_BLAS([])
m4_include([m4/m4_ax_lapack.m4])
AX_LAPACK([])
LT_INIT([])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_TIME
gl_VISIBILITY

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday])

AM_PROG_CC_C_O 

case "x$CC" in
xmcc)
  AC_DEFINE([USE_OMPSS], 1, [Should we build with OMPSs.])
  AC_SUBST([OMPSS_CFLAGS],[--ompss])
;;
xcc)
;;
*)
  AC_MSG_ERROR([Cannot find a suitable compiler])
esac

AC_OUTPUT
