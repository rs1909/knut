AC_INIT([PDDE-CONT], [1.0.0])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE

AC_PROG_CC
AC_C_INLINE
AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG(C++)

REVISION="`cat REVISION`"
AC_SUBST(REVISION)
AC_DEFINE([PKG_REV],["@REVISION@"])

UNSUPPORTED_OS="PDDE-CONT is currently unsupported on your platform.
If you are interested in making it work on your platform, you are
more than *welcome*. Contact szalai@mm.bme.hu for details."

MISSING_DLFCN_H="Missing dlfcn.h. The software could not be built."
MISSING_WINDOWS_H="Missing windows.h. The software could not be built."
MISSING_DLOPEN="Missing dlopen. The software could not be built."
NO_SHARED="Cannot find a method to create dynamically loadable objects."

AC_DEFUN([PDDE_DISP_ARCH], [
  echo ============================================================
  echo  Setting up build environment for ${target_cpu}-${target_os}
  echo ============================================================
])

AC_DEFUN([PDDE_CHECK_UNIX], [
  AC_CHECK_HEADERS( dlfcn.h, [ ], [ AC_MSG_ERROR($MISSING_DLFCN_H) ] )
  AC_CHECK_LIB( dl, dlopen, 
    [ DL_LIB="-ldl" ], [ 
      AC_CHECK_LIB( ,dlopen, [ DL_LIB=""], [ AC_MSG_ERROR($MISSING_DLOPEN) ]) 
    ])
  AC_SUBST(DL_LIB)
])

AC_DEFUN([PDDE_CHECK_SHARED],
[ AC_MSG_CHECKING([whether --shared -fPIC is supported by $CXX])
  TMP_CXXFLAGS=$CXXFLAGS
  CXXFLAGS="$TMP_CXXFLAGS --shared -fPIC"
  AC_LINK_IFELSE([int main() { return 0; }],
  [ MOD_CXXFLAGS=$CXXFLAGS; CXXFLAGS=$TMP_CXXFLAGS; echo "yes"], 
  [ echo "no"
    AC_MSG_CHECKING([whether -G -KPIC is supported by $CXX])
    CXXFLAGS="$TMP_CXXFLAGS -G -KPIC"
    AC_LINK_IFELSE([int main() { return 0; }], 
      [ MOD_CXXFLAGS=$CXXFLAGS; CXXFLAGS=$TMP_CXXFLAGS; echo "yes"],
      [ AC_MSG_ERROR($NO_SHARED) ])
  ])
])

case "${target_os}" in
  linux*)
    PDDE_DISP_ARCH()
    PDDE_CHECK_SHARED()
    PDDE_CHECK_UNIX()
    ;;
  mingw32*)
    PDDE_DISP_ARCH()    
    PDDE_CHECK_SHARED()
    AC_CHECK_HEADERS( windows.h, [ echo "yes" ], [ AC_MSG_ERROR($MISSING_WINDOWS_H) ])
    ;;
  solaris*)
    PDDE_DISP_ARCH()
    PDDE_CHECK_SHARED()
    PDDE_CHECK_UNIX()
    ;;
  irix*)
    PDDE_DISP_ARCH()
    PDDE_CHECK_SHARED()
    PDDE_CHECK_UNIX()
    ;;
  *)
    AC_MSG_ERROR($UNSUPPORTED_OS)
    ;;
esac

dnl setting up warning message handling
dnl if GNU C++ is used

case "$CXX" in
  *g++*)
    echo "Setting up the warning flags for $CXX"
    NOWARN_CXXFLAGS="-fno-rtti -fno-exceptions"
    CXXFLAGS="$CXXFLAGS -W -Wall $NOWARN_CXXFLAGS"
    MOD_CXXFLAGS="$MOD_CXXFLAGS $NOWARN_CXXFLAGS"  
    AC_SUBST(CXXFLAGS NOWARN_CXXFLAGS MOD_CXXFLAGS)
    
    AC_MSG_CHECKING([whether -fvisibility-inlines-hidden is supported by $CXX])
    TMP_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS -fvisibility-inlines-hidden"
    AC_RUN_IFELSE([int main() { return 0; }],[echo "yes"],
    [ CXXFLAGS=$TMP_CXXFLAGS
    echo "no"])
    ;;
  *)
    echo "Not using g++ and not setting warning and visibility flags for $CXX"
    ;;
esac

dnl Setting up the cflags for not optimizing
dnl the dlamch.c in laarpack.
dnl This is needed if the architecture has wider
dnl floating point registers than double's stored size

case "$CC" in
  gcc*) NOOPT_CFLAGS="-O0" ;;
  *) NOOPT_CFLAGS=$CFLAGS ;;
esac

echo "Flag for not optimizing with $CC: $NOOPT_CFLAGS"

dnl This also can be specified by a command line argument

AC_ARG_WITH(no-optimization-flags,
[  --with-no-optimization-flags=FLAGS       Choose an optimization flag for not optimizing some parts of the code (default: -O0)],
[NOOPT_CFLAGS=$withval], )

AC_SUBST(NOOPT_CFLAGS)

dnl specfy the locatio of the ATLAS libraries

AC_ARG_WITH(atlas,
[  --with-atlas=<directory>       specify the location of the atlas libraries (default: -O0)],
[ATLAS_DIR=$withval], [ATLAS_DIR="`pwd`/ATLAS/lib"])
AC_SUBST(ATLAS_DIR)

dnl remove "#" if the software crashes

# AC_DEFINE([DEBUG])
# CFLAGS="$CFLAGS -fno-inline"
# CXXFLAGS="$CXXFLAGS -fno-inline"

AC_CHECK_PROG(SED, sed, sed, no)

AC_PROG_INSTALL

AC_CONFIG_FILES([Makefile  
	src/Makefile
	bin/Makefile
	bin/pcompile
	libf2c/Makefile
	laarpack/Makefile
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
	demo/turntwo/Makefile 
	demo/laser/Makefile 
	UMFPACK/Makefile 
	UMFPACK/AMD/Makefile 
	UMFPACK/AMD/Source/Makefile 
	UMFPACK/UMFPACK/Makefile
	UMFPACK/UMFPACK/Source/Makefile
	pdde.spec 
	pdde-sse.spec 
	pdde-sse2.spec 
	pdde-3dnow.spec 
	])
	
AC_OUTPUT
