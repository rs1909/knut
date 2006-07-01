PROJECT( PDDE )

#------------------------------------------------------------------------------
# SETTING UP THE VARIABLES
#------------------------------------------------------------------------------

SET( PACKAGE_NAME "PDDE-CONT" )
SET( PACKAGE_VERSION "1.0.4" )
FILE( READ REVISION REV )
STRING( REPLACE "\\n" "" PACKAGE_REVISION ${REV} )
CONFIGURE_FILE( src/config.h.cmake src/config.h )

FIND_FILE( HAVE_CONFIG config.h ${PDDE_SOURCE_DIR}/src )
IF( EXISTS ${HAVE_CONFIG} )
        ADD_DEFINITIONS( -DHAVE_CONFIG_H )
ENDIF( EXISTS ${HAVE_CONFIG} )

FIND_PROGRAM( SED NAMES sed )
IF( NOT EXISTS ${SED} )
	MESSAGE( SEND_ERROR "The sed utility is mot found" )
ENDIF( NOT EXISTS ${SED} )

FIND_PROGRAM( SHELL NAMES sh bash )
IF( NOT EXISTS ${SHELL} )
        MESSAGE( SEND_ERROR "Cannot find a shell program" )
ENDIF( NOT EXISTS ${SHELL} )

# finding GCC or the Intel compiler

SET( C_WARNINGS "-W -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Wstrict-prototypes -Wpointer-arith" )
SET( CXX_WARNINGS "-W -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Wpointer-arith" )

FIND_PROGRAM( C_COMPILER NAMES gcc icc )
IF( EXISTS ${C_COMPILER} )
	MESSAGE( "A C compiler has been found: ${C_COMPILER}" )
	SET( CMAKE_C_COMPILER "${C_COMPILER}" CACHE STRING "C compiler" FORCE )
	SET( CMAKE_C_FLAGS "-O2 ${C_WARNINGS}" CACHE STRING "CFLAGS" FORCE )
ENDIF( EXISTS ${C_COMPILER} )

FIND_PROGRAM( CXX_COMPILER NAMES g++ icpc )
IF( EXISTS ${CXX_COMPILER} )
	MESSAGE( "A C++ compiler has been found: ${CXX_COMPILER}" )
	SET( CMAKE_CXX_COMPILER ${CXX_COMPILER} CACHE STRING "C++ compiler" FORCE )
        SET( CMAKE_CXX_FLAGS "-O2 ${CXX_WARNINGS}" CACHE STRING "CXX_FLAGS" FORCE )
ENDIF( EXISTS ${CXX_COMPILER} )

CONFIGURE_FILE( bin/pcompile.cmake bin/pcompile )
INSTALL_PROGRAMS( /bin FILES bin/pcompile )

#------------------------------------------------------------------------------
# DISTRIBUTING WITH SVN
#------------------------------------------------------------------------------

SET( DIST_DIR pdde-cont-${PACKAGE_VERSION} )
ADD_CUSTOM_TARGET( dist 
	COMMAND ln -s . ${DIST_DIR}\; 
		svn log -r HEAD -q | sed -e s/' |.*.|'// -e s/' \(.*\)'// -e /'---'/d | awk -- '/r/ { printf \"%s\", $$0 }' >REVISION\;
		svn list -r HEAD -R --non-interactive | awk -- '/.*./ { print \"${DIST_DIR}/\"$$0}' >tmp_filelist\;
		echo >>tmp_filelist \"${DIST_DIR}/REVISION\"\;
		tar cjf ${DIST_DIR}.tar.bz2 --no-recursion -h -T tmp_filelist\;
		rm tmp_filelist ${DIST_DIR}\; )

#------------------------------------------------------------------------------
# ENTERING SUBDIRECTORIES
#------------------------------------------------------------------------------

# fortran to C support routines
# this is required by LAARPACK
ADD_SUBDIRECTORY(libf2c)

# LAPACK and ARPACK
# this produces a library "liblaarpack.a"
ADD_SUBDIRECTORY(laarpack)

# these are the UMFPACK routines
ADD_SUBDIRECTORY(UMFPACK)

# This is where the whole computation is
# this produces a library "libpdde.a"
ADD_SUBDIRECTORY(src)

# The command line interface
ADD_SUBDIRECTORY(cli)

# Documentation
ADD_SUBDIRECTORY(doc)

# demo programs
ADD_SUBDIRECTORY(demo)

