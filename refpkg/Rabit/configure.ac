dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Rabit, 1.0)


AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()

AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

#GSL check
AC_CHECK_HEADERS([gsl/gsl_blas.h gsl/gsl_cdf.h], ,AC_MSG_ERROR([Cannot find gsl headers]))
AC_CHECK_LIB(gslcblas, cblas_dgemm, ,AC_MSG_ERROR([Cannot find cblas_dgemm]))
AC_CHECK_LIB(gsl, gsl_cdf_gaussian_Pinv, ,AC_MSG_ERROR([Cannot find gsl_cdf_gaussian_Pinv]))
AC_CHECK_LIB(gsl, gsl_linalg_cholesky_decomp, ,AC_MSG_ERROR([Cannot find gsl_linalg_cholesky_decomp]))
AC_CHECK_LIB(gsl, gsl_linalg_cholesky_invert, ,AC_MSG_ERROR([Cannot find gsl_linalg_cholesky_invert]))

# LAPACK check
AC_CHECK_LIB(lapack, dgeqrf_, , AC_MSG_ERROR([Cannot find lapack dgeqrf_]), [-lblas -lgfortran])

# disable default LIBS
AC_SUBST(LIBS, "")

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT
