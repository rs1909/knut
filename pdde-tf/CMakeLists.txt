PROJECT( PDDE )

#------------------------------------------------------------------------------
# SETTING UP THE VARIABLES
#------------------------------------------------------------------------------

SET( PACKAGE_NAME "PDDE-CONT" )
SET( PACKAGE_VERSION "1.0.4" )
FILE( READ REVISION REV )
STRING( REPLACE "\\n" "" PACKAGE_REVISION ${REV} )
CONFIGURE_FILE( src/config.h.cmake src/config.h )
CONFIGURE_FILE( pdde.spec.cmake pdde.spec )

FIND_FILE( HAVE_CONFIG config.h ${PDDE_SOURCE_DIR}/src )
IF( EXISTS ${HAVE_CONFIG} )
        ADD_DEFINITIONS( -DHAVE_CONFIG_H )
ENDIF( EXISTS ${HAVE_CONFIG} )

FIND_PROGRAM( SED NAMES sed )
IF( NOT EXISTS ${SED} )
	MESSAGE( SEND_ERROR "The sed utility is mot found" )
ENDIF( NOT EXISTS ${SED} )

FIND_PROGRAM( SHELL NAMES sh bash )
IF( EXISTS ${SHELL} )
	MESSAGE( "The shell is found: ${SHELL}" )
ELSE( EXISTS ${SHELL} )
	MESSAGE( SEND_ERROR "Cannot find a shell program" )
ENDIF( EXISTS ${SHELL} )


# finding GCC or the Intel compiler

SET( GCC_C_WARN "-W -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Wstrict-prototypes -Wpointer-arith" )
SET( GCC_CXX_WARN "-W -Wall -Wno-unused-parameter -Wno-unknown-pragmas -Wpointer-arith" )
SET( INTEL_C_WARN  )
SET( INTEL_CXX_WARN  )

SET( GCC_CFLAGS "-O2" )
SET( GCC_NOOPT_CFLAGS "-O0" )
SET( GCC_CXXFLAGS "-O2" )

SET( INTEL_GCC_COMPATIBILITY "340" )
SET( INTEL_CFLAGS "-gcc-version=${INTEL_GCC_COMPATIBILITY} -O2" )
SET( INTEL_NOOPT_CFLAGS "-gcc-version=${INTEL_GCC_COMPATIBILITY} -O0" )
SET( INTEL_CXXFLAGS "-gcc-version=${INTEL_GCC_COMPATIBILITY} -cxxlib-icc -xB -O2" )

FIND_PROGRAM( C_COMPILER NAMES gcc icc )
IF( EXISTS ${C_COMPILER} )
	MESSAGE( "A C compiler has been found: ${C_COMPILER}" )
	SET( CMAKE_C_COMPILER ${C_COMPILER} CACHE STRING "C compiler" FORCE )
	GET_FILENAME_COMPONENT( CC_NAME ${C_COMPILER} NAME )
	IF( ${CC_NAME} STREQUAL icc )
		SET( CMAKE_C_FLAGS "${INTEL_CFLAGS} ${INTEL_C_WARN}" CACHE STRING "CFLAGS" FORCE )
		SET( CMAKE_C_NOOPT_FLAGS "${INTEL_NOOPT_CFLAGS} ${INTEL_C_WARN}" )
	ELSE( ${CC_NAME} STREQUAL icc )
		SET( CMAKE_C_FLAGS "${GCC_CFLAGS} ${GCC_C_WARN}" CACHE STRING "CFLAGS" FORCE )
		SET( CMAKE_C_NOOPT_FLAGS "${GCC_NOOPT_CFLAGS} ${GCC_C_WARN}" )
	ENDIF( ${CC_NAME} STREQUAL icc )
ENDIF( EXISTS ${C_COMPILER} )

FIND_PROGRAM( CXX_COMPILER NAMES g++ icpc )
IF( EXISTS ${CXX_COMPILER} )
	MESSAGE( "A C++ compiler has been found: ${CXX_COMPILER}" )
	SET( CMAKE_CXX_COMPILER ${CXX_COMPILER} CACHE STRING "C++ compiler" FORCE )
	GET_FILENAME_COMPONENT( CXX_NAME ${CXX_COMPILER} NAME )
	IF( ${CXX_NAME} STREQUAL icpc )
		SET( CMAKE_CXX_FLAGS "${INTEL_CXXFLAGS} ${INTEL_CXX_WARN}" CACHE STRING "CFLAGS" FORCE )
	ELSE( ${CXX_NAME} STREQUAL icpc )
		SET( CMAKE_CXX_FLAGS "${GCC_CXXFLAGS} ${GCC_CXX_WARN}" CACHE STRING "CFLAGS" FORCE )
	ENDIF( ${CXX_NAME} STREQUAL icpc )
ENDIF( EXISTS ${CXX_COMPILER} )

CONFIGURE_FILE( bin/pcompile.cmake bin/pcompile )
INSTALL_PROGRAMS( /bin FILES bin/pcompile )

#------------------------------------------------------------------------------
# Finding ATLAS
#------------------------------------------------------------------------------

SET( ATLAS_SEARCH_PATHS ${PDDE_SOURCE_DIR}/ATLAS/lib /usr/lib /usr/lib/atlas
	/usr/lib/sse /usr/lib/sse2 /usr/lib/3dnow )

FIND_LIBRARY( ATLAS_CBLAS NAMES cblas PATHS ${ATLAS_SEARCH_PATHS} )

FIND_LIBRARY( ATLAS_ATLAS NAMES atlas PATHS ${ATLAS_SEARCH_PATHS} )

IF( EXISTS ${ATLAS_CBLAS} AND EXISTS ${ATLAS_ATLAS} )
	GET_FILENAME_COMPONENT( CBLAS_PATH ${ATLAS_CBLAS} PATH )
	GET_FILENAME_COMPONENT( ATLAS_PATH ${ATLAS_ATLAS} PATH )
	MESSAGE( "ATLAS was found in ${CBLAS_PATH} and ${ATLAS_PATH}" )
ELSE( EXISTS ${ATLAS_CBLAS} AND EXISTS ${ATLAS_ATLAS} )
        MESSAGE( SEND_ERROR "ATLAS was NOT found!" )
ENDIF( EXISTS ${ATLAS_CBLAS} AND EXISTS ${ATLAS_ATLAS} )

#------------------------------------------------------------------------------
# DISTRIBUTING WITH SVN
#------------------------------------------------------------------------------

SET( DIST_DIR pdde-cont-${PACKAGE_VERSION} )
ADD_CUSTOM_TARGET( dist 
	COMMAND ln -s . ${DIST_DIR}\; 
		svn log -r HEAD -q | sed -e s/' |.*.|'// -e s/' \(.*\)'// -e /'---'/d | awk -- '/r/ { printf \"%s\", $$0 }' >REVISION\;
		svn list -r HEAD -R --non-interactive | awk -- '/.*./ { print \"${DIST_DIR}/\"$$0}' >tmp_filelist\;
		echo >>tmp_filelist \"${DIST_DIR}/REVISION\"\;
		echo >>tmp_filelist \"${DIST_DIR}/pdde.spec\"\;
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

