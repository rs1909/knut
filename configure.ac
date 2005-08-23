AC_INIT([PDDE-CONT], [0.9.16])
AC_DEFINE(DEFS,"")
AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

if test "$GXX" == yes; then
  NOWARN_CXXFLAGS="$CXXFLAGS -fno-inline"
  CXXFLAGS="$NOWARN_CXXFLAGS -W -Wall"
  AC_SUBST(CXXFLAGS NOWARN_CXXFLAGS)
fi

AC_ARG_WITH(no-optimization-flags,
[  --with-no-optimization-flags=FLAGS       Choose an optimization flag to give togive to the compilers (default: -O0)],
[NOOPT_CFLAGS=$withval],[NOOPT_CFLAGS="-O0"])

AC_SUBST(NOOPT_CFLAGS)

# AC_DEFINE([DEBUG])

AC_PROG_INSTALL

AC_CONFIG_FILES([Makefile  
	src/Makefile
	bin/Makefile
	bin/pcompile
	libf2c/Makefile
	laarpack/Makefile
	UMFPACK/Makefile
	ATLAS/Makefile
	doc/Makefile
	doc/fig/Makefile
	doc/html/Makefile
	doc/html/index.html
	demo/Makefile 
	demo/glass/Makefile 
	demo/logistic/Makefile 
	demo/milloneup/Makefile
	demo/milltwodown/Makefile 
	demo/milltwofull/Makefile 
	demo/milltwoup/Makefile 
	demo/tq/Makefile 
	demo/traffic/Makefile 
	demo/turn/Makefile 
	demo/laser/Makefile ])
	
AC_OUTPUT
