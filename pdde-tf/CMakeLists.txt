PROJECT( PDDE )

#------------------------------------------------------------------------------
# SETTING UP THE VARIABLES
#------------------------------------------------------------------------------

SET( PACKAGE_NAME "PDDE-CONT" )
SET( PACKAGE_VERSION "1.0.4" )
FILE( READ REVISION PACKAGE_REVISION )

SET( PDDE_GUI TRUE )

IF( NOT WIN32 )
  FIND_PROGRAM( SED NAMES sed )
  IF( NOT EXISTS ${SED} )
    MESSAGE( SEND_ERROR "The sed utility is not found" )
  ENDIF( NOT EXISTS ${SED} )

  FIND_PROGRAM( SHELL NAMES sh bash )
  IF( EXISTS ${SHELL} )
    MESSAGE( "The shell is found: ${SHELL}" )
  ELSE( EXISTS ${SHELL} )
    MESSAGE( SEND_ERROR "Cannot find a shell program" )
  ENDIF( EXISTS ${SHELL} )
  SET( PDDE_DEMO_DIR share/${PACKAGE_NAME}/demo/ CACHE INTERNAL "demo dir" )
  SET( PDDE_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
ELSE( NOT WIN32 )
  SET( PDDE_DEMO_DIR demo/ CACHE INTERNAL "demo dir" )
  SET( PDDE_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
ENDIF( NOT WIN32 )

#------------------------------------------------------------------------------
# Writing the generated files
#------------------------------------------------------------------------------

CONFIGURE_FILE( src/config.h.cmake src/config.h )
CONFIGURE_FILE( pdde.spec.cmake pdde.spec )
CONFIGURE_FILE( doc/html/index.html.in doc/html/index.html )
CONFIGURE_FILE( bin/pcompile.cmake bin/pcompile )

IF( NOT WIN32 )
  INSTALL( PROGRAMS bin/pcompile DESTINATION ${PDDE_BIN_DIR} )
ENDIF( NOT WIN32 )

ADD_DEFINITIONS( -DHAVE_CONFIG_H )
# ADD_DEFINITIONS( -DHAVE_CONFIG_H -DDEBUG )

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
	MESSAGE( STATUS "ATLAS was found in ${CBLAS_PATH} and ${ATLAS_PATH}" )
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

# graphical user interface
ADD_SUBDIRECTORY(gui)

# demo programs
ADD_SUBDIRECTORY(demo)
