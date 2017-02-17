#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(nvxs, 1.0.2, nagadomi@nurs.or.jp)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([.])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
#AC_PROG_CXX
#AC_PROG_CC
#AC_PROG_CPP
#AC_PROG_INSTALL
#AC_PROG_LN_S
#AC_PROG_MAKE_SET

# Checks for libraries.
#AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stdlib.h string.h strings.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
#AC_FUNC_SETVBUF_REVERSED
AC_CHECK_FUNCS([memset])

# OpenMP
AC_ARG_ENABLE(openmp,
        AC_HELP_STRING([--enable-openmp],
                [fast parallel scan with OpenMP (default: disable)]),
        [CFLAGS+=" -fopenmp";LDFLAGS+=" -lgomp"],
        [])

#AC_CONFIG_FILES([nv_core/Makefile nv_num/Makefile nv_ip/Makefile nv_ml/Makefile nv_face/Makefile Makefile])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
