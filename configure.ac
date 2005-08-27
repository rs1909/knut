AC_INIT([PDDE-CONT], [0.9.17-rc2])
AC_DEFINE(DEFS,"")
AM_INIT_AUTOMAKE

AC_PROG_CC
AC_C_INLINE
AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG(C++)
AC_MSG_CHECKING([whether -fvisibility-inlines-hidden is supported by $CXX])
TMP_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -fvisibility-inlines-hidden"
AC_COMPILE_IFELSE([int main() { return 0; }],[echo "yes"],
[ CXXFLAGS=$TMP_CXXFLAGS
echo "no"])

AC_MSG_CHECKING([whether -fPIC --shared is supported by $CXX])
TMP_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -fPIC --shared"
AC_COMPILE_IFELSE([int main() { return 0; }],
[ CXXFLAGS=$TMP_CXXFLAGS; MOD_CXXFLAGS="-fPIC --shared"; echo "yes"],
[ CXXFLAGS=$TMP_CXXFLAGS; MOD_CXXFLAGS=""; echo "no"])

if test "$GXX" == yes; then
  NOWARN_CXXFLAGS="$CXXFLAGS -fno-rtti -fno-exceptions"
  CXXFLAGS="$NOWARN_CXXFLAGS -W -Wall -Woverloaded-virtual"
  MOD_CXXFLAGS="$NOWARN_CXXFLAGS $MOD_CXXFLAGS"  
  AC_SUBST(CXXFLAGS NOWARN_CXXFLAGS MOD_CXXFLAGS)
fi

AC_ARG_WITH(no-optimization-flags,
[  --with-no-optimization-flags=FLAGS       Choose an optimization flag for not optimizing some parts of the code (default: -O0)],
[NOOPT_CFLAGS=$withval],[NOOPT_CFLAGS="-O0"])

AC_SUBST(NOOPT_CFLAGS)

AC_DEFINE([DEBUG])

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
