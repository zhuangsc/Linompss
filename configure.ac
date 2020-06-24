#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([linompss],[1.0],[pieter.bellens@bsc.es, sicong.zhuang@bsc.es])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/cholesky/dense/chol.c])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile 
		src/Makefile
		src/csparse/Makefile
		src/util/Makefile
		src/jacobi/Makefile
		src/gauss_seidel/Makefile
		src/cholesky/Makefile
		src/cholesky/dense/Makefile
		src/cholesky/dense/kernels/Makefile
		src/cholesky/sparse/Makefile
		src/cholesky/sparse/kernels/Makefile
		src/blocked_blas/Makefile
		src/blocked_sblas/Makefile
		src/trsm/Makefile
		src/spdsolv/Makefile
		src/cg/Makefile
		src/gemm/Makefile
		src/syrk/Makefile
		src/fp/Makefile
		src/lu/Makefile
		src/gensolv/Makefile
		src/geqr/Makefile
		src/itref/Makefile
		src/matutil/Makefile])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC([smpcc cc])

# Checks for libraries.
m4_include([m4/m4_ax_blas.m4])
AX_BLAS([])
m4_include([m4/m4_ax_lapack.m4])
AX_LAPACK([])
AC_CHECK_LIB(config, config_lookup, [], [], [])
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
AC_CHECK_FUNCS([gettimeofday])

AM_PROG_CC_C_O 

case "x$CC" in
xsmpcc)
  AC_DEFINE([USE_OMPSS], 1, [Should we build with OMPSs.])
  AC_SUBST([OMPSS_CFLAGS],[--ompss])
;;
xcc)
;;
*)
  AC_MSG_ERROR([Cannot find a suitable compiler])
esac

AS_IF([test "x$ax_mkl" = xno], [AC_DEFINE([USE_MKL], [0], [Do we use Intel MKL.])], [AC_DEFINE([USE_MKL], [1], [Do we use Intel MKL.])])
AM_CONDITIONAL([USE_MKL], [test "x$ax_mkl" = xyes])

AC_OUTPUT
