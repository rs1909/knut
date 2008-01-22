PROJECT( PDDE )

#------------------------------------------------------------------------------
# SETTING UP THE VARIABLES
#------------------------------------------------------------------------------

SET( PACKAGE_NAME "pdde-cont" )
SET( PACKAGE_VERSION "1.2.0" )
FILE( READ ${PDDE_SOURCE_DIR}/REVISION PACKAGE_REVISION )

SET( PDDE_GUI TRUE )

IF( NOT WIN32 )
  FIND_PROGRAM( SED NAMES sed )
  IF( NOT EXISTS ${SED} )
    MESSAGE( SEND_ERROR "The sed utility is not found" )
  ENDIF( NOT EXISTS ${SED} )

  FIND_PROGRAM( SHELL NAMES sh bash )
  IF( EXISTS ${SHELL} )
    MESSAGE( STATUS "The shell is found: ${SHELL}" )
  ELSE( EXISTS ${SHELL} )
    MESSAGE( SEND_ERROR "Cannot find a shell program" )
  ENDIF( EXISTS ${SHELL} )
  SET( PDDE_DEMO_DIR share/${PACKAGE_NAME}/demo/ CACHE INTERNAL "demo dir" )
  SET( PDDE_MATLAB_DIR share/${PACKAGE_NAME}/matlab/ CACHE INTERNAL "matlab dir" )
  SET( PDDE_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
ELSE( NOT WIN32 )
  SET( PDDE_DEMO_DIR demo/ CACHE INTERNAL "demo dir" )
  SET( PDDE_MATLAB_DIR matlab/ CACHE INTERNAL "matlab dir" )
  SET( PDDE_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
ENDIF( NOT WIN32 )

IF( MSVC )
  ADD_DEFINITIONS(-Dinline=_inline -D_USE_MATH_DEFINES)
ENDIF( MSVC )

#------------------------------------------------------------------------------
# finding Qt4
#------------------------------------------------------------------------------

SET( QT_MIN_VERSION "4.3.0" )
FIND_PACKAGE( Qt4 )

IF(WIN32)
  FIND_LIBRARY( PDDE_QTMAIN qtmain ${QT_LIBRARY_DIR} )
  MESSAGE( STATUS "QtMain as ${PDDE_QTMAIN} was found" )
ENDIF(WIN32)

#------------------------------------------------------------------------------
# Writing the generated files
#------------------------------------------------------------------------------

CONFIGURE_FILE( src/config.h.cmake src/config.h )
CONFIGURE_FILE( pdde.spec.cmake pdde.spec )
CONFIGURE_FILE( doc/html/index.html.in doc/html/index.html )
CONFIGURE_FILE( bin/pcompile.cmake bin/pcompile )
IF(WIN32)
  CONFIGURE_FILE( wininstaller.nsi.cmake wininstaller.nsi )
ENDIF(WIN32)

IF( NOT WIN32 )
  INSTALL( PROGRAMS ${PDDE_BINARY_DIR}/bin/pcompile DESTINATION ${PDDE_BIN_DIR} )
ENDIF( NOT WIN32 )

# no debugging with release version
IF( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )
  ADD_DEFINITIONS( -DHAVE_CONFIG_H -DDEBUG )
  SET(PDDE_RELEASE FALSE)
ELSE( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )
  ADD_DEFINITIONS( -DHAVE_CONFIG_H )
  SET(PDDE_RELEASE TRUE)
ENDIF( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )

#------------------------------------------------------------------------------
# Finding ATLAS
#------------------------------------------------------------------------------

SET( ATLAS_SEARCH_PATHS ${PDDE_SOURCE_DIR}/ATLAS/lib /usr/lib64/atlas /usr/lib /usr/lib/atlas
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

SET( DIST_DIR ${PACKAGE_NAME}-${PACKAGE_VERSION} )
ADD_CUSTOM_TARGET( dist 
	COMMAND ln -s ${PDDE_SOURCE_DIR} ${DIST_DIR}\;
		ln -s ${PDDE_BINARY_DIR}/pdde.spec ${DIST_DIR}/pdde.spec\;
		pushd ${PDDE_SOURCE_DIR}\;
                git log -n1 HEAD | grep commit | sed -e s/'commit '// | awk -- '// { printf \"%s\", $$0 }' >${PDDE_SOURCE_DIR}/REVISION\;
                git log -n1 HEAD | grep 'Date' | sed -e s/'Date: '/','/ | awk -- '// { printf \"%s\", $$0 }' >>${PDDE_SOURCE_DIR}/REVISION\;
                git-ls-files | awk -- '/.*./ { print \"${DIST_DIR}/\"$$0}' >${PDDE_BINARY_DIR}/tmp_filelist\;
                popd\;
		svn log -r HEAD -q ${PDDE_SOURCE_DIR} | sed -e s/' |.*.|'// -e s/' \(.*\)'// -e /'---'/d | awk -- '/r/ { printf \"%s\", $$0 }' >${PDDE_SOURCE_DIR}/REVISION\;
		svn list -r HEAD -R --non-interactive ${PDDE_SOURCE_DIR} | awk -- '/.*./ { print \"${DIST_DIR}/\"$$0}' >tmp_filelist\;
		echo >>tmp_filelist \"${DIST_DIR}/REVISION\"\;
		echo >>tmp_filelist \"${DIST_DIR}/pdde.spec\"\;
		tar cjf ${DIST_DIR}.tar.bz2 --exclude *.svg --no-recursion -h -T tmp_filelist\;
		rm tmp_filelist ${DIST_DIR}/pdde.spec ${DIST_DIR}\; )

#------------------------------------------------------------------------------
# ENTERING SUBDIRECTORIES
#------------------------------------------------------------------------------

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

# graphical user interface
IF(QT4_QMAKE_FOUND)
  IF(PDDE_QTMAIN)
    SET(BUILDGUI TRUE)
  ENDIF(PDDE_QTMAIN)
  IF(NOT WIN32)
    SET(BUILDGUI TRUE)
  ENDIF(NOT WIN32)
ENDIF(QT4_QMAKE_FOUND)

IF(BUILDGUI)
  ADD_SUBDIRECTORY(gui)
  MESSAGE( STATUS "GUI will be built" )
ELSE(BUILDGUI)  
  MESSAGE( STATUS "GUI will not be built" )
ENDIF(BUILDGUI)

# demo programs
ADD_SUBDIRECTORY(demo)

# MATLAB scripts
ADD_SUBDIRECTORY(matlab)
