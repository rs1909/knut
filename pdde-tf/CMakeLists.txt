cmake_minimum_required(VERSION 2.6)
PROJECT( KNUT C CXX Fortran )

#------------------------------------------------------------------------------
# SETTING UP THE VARIABLES
#------------------------------------------------------------------------------

SET( PACKAGE_NAME "Knut" )
SET( PACKAGE_VERSION "2.0.0" )
FILE( READ ${KNUT_SOURCE_DIR}/REVISION PACKAGE_REVISION )

SET( KNUT_GUI TRUE )

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
  IF(APPLE)
    SET( KNUT_DEMO_DIR demo/ CACHE INTERNAL "demo dir" )
    SET( KNUT_MATLAB_DIR matlab/ CACHE INTERNAL "matlab dir" )
    SET( KNUT_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
    SET( KNUT_INCLUDE_DIR Headers CACHE INTERNAL "include dir" )
  ELSE(APPLE)
    SET( KNUT_DEMO_DIR share/${PACKAGE_NAME}/demo/ CACHE INTERNAL "demo dir" )
    SET( KNUT_MATLAB_DIR share/${PACKAGE_NAME}/matlab/ CACHE INTERNAL "matlab dir" )
    SET( KNUT_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
    SET( KNUT_INCLUDE_DIR include CACHE INTERNAL "include dir" )
  ENDIF(APPLE)
ELSE( NOT WIN32 )
  SET( KNUT_DEMO_DIR demo/ CACHE INTERNAL "demo dir" )
  SET( KNUT_MATLAB_DIR matlab/ CACHE INTERNAL "matlab dir" )
  SET( KNUT_BIN_DIR bin/ CACHE INTERNAL "bin dir" )
ENDIF( NOT WIN32 )

IF( MSVC )
  ADD_DEFINITIONS(-Dinline=_inline -D_USE_MATH_DEFINES)
ENDIF( MSVC )

#------------------------------------------------------------------------------
# setting up Mini Xml
#------------------------------------------------------------------------------

SET( MXML_INCLUDE_DIRS ${KNUT_SOURCE_DIR}/mxml )

#------------------------------------------------------------------------------
# finding Qt4
#------------------------------------------------------------------------------

SET( QT_MIN_VERSION "4.4.0" )
FIND_PACKAGE( Qt4 )

IF(WIN32)
  FIND_LIBRARY( KNUT_QTMAIN qtmain ${QT_LIBRARY_DIR} )
  MESSAGE( STATUS "QtMain as ${KNUT_QTMAIN} was found" )
ENDIF(WIN32)

#------------------------------------------------------------------------------
# Writing the generated files
#------------------------------------------------------------------------------

CONFIGURE_FILE( src/config.h.cmake src/config.h )
CONFIGURE_FILE( knut.spec.cmake knut.spec )
CONFIGURE_FILE( doc/html/index.html.cmake doc/html/index.html )
CONFIGURE_FILE( bin/pcompile.cmake bin/pcompile )
IF(WIN32)
  CONFIGURE_FILE( wininstaller.nsi.cmake wininstaller.nsi )
ENDIF(WIN32)

IF( NOT WIN32 )
  INSTALL( PROGRAMS ${KNUT_BINARY_DIR}/bin/pcompile DESTINATION ${KNUT_BIN_DIR} )
ENDIF( NOT WIN32 )

# no debugging with release version
IF( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )
  ADD_DEFINITIONS( -DHAVE_CONFIG_H -DDEBUG )
  SET(KNUT_RELEASE FALSE)
ELSE( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )
  ADD_DEFINITIONS( -DHAVE_CONFIG_H )
  SET(KNUT_RELEASE TRUE)
ENDIF( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )

#------------------------------------------------------------------------------
# Finding ATLAS
#------------------------------------------------------------------------------

IF(APPLE)
  # We need to search for ATLAS in MacPorts' location, too: /opt/local/lib
  SET( ATLAS_SEARCH_PATHS ${KNUT_SOURCE_DIR}/ATLAS/lib /opt/local/lib /usr/lib )
  FIND_LIBRARY( ATLAS_CBLAS NAMES libf77blas.a cblas PATHS ${ATLAS_SEARCH_PATHS} )
  FIND_LIBRARY( ATLAS_ATLAS NAMES libatlas.a atlas PATHS ${ATLAS_SEARCH_PATHS} )
ELSE(APPLE)
  # For Linux distributions ATLAS can be located in many places
  SET( ATLAS_SEARCH_PATHS ${KNUT_SOURCE_DIR}/ATLAS/lib /usr/lib64/atlas /usr/lib /usr/lib/atlas
    /usr/lib/sse /usr/lib/sse2 /usr/lib/3dnow )
  FIND_LIBRARY( ATLAS_CBLAS NAMES f77blas PATHS ${ATLAS_SEARCH_PATHS} )
  FIND_LIBRARY( ATLAS_ATLAS NAMES atlas PATHS ${ATLAS_SEARCH_PATHS} )
ENDIF(APPLE)

IF( EXISTS ${ATLAS_CBLAS} AND EXISTS ${ATLAS_ATLAS} )
	MESSAGE( STATUS "ATLAS was found as ${ATLAS_CBLAS} and ${ATLAS_ATLAS}" )
ELSE( EXISTS ${ATLAS_CBLAS} AND EXISTS ${ATLAS_ATLAS} )
        MESSAGE( SEND_ERROR "ATLAS was NOT found!" )
ENDIF( EXISTS ${ATLAS_CBLAS} AND EXISTS ${ATLAS_ATLAS} )

#------------------------------------------------------------------------------
# DISTRIBUTING WITH SVN
#------------------------------------------------------------------------------

SET( DIST_DIR ${PACKAGE_NAME}-${PACKAGE_VERSION} )
ADD_CUSTOM_TARGET( dist 
	COMMAND ln -s ${KNUT_SOURCE_DIR} ${DIST_DIR}\;
		ln -s ${KNUT_BINARY_DIR}/knut.spec ${DIST_DIR}/knut.spec\;
		pushd ${KNUT_SOURCE_DIR}\;
                git log -n1 HEAD | grep commit | sed -e s/'commit '// | awk -- '// { printf \"%s\", $$0 }' >${KNUT_SOURCE_DIR}/REVISION\;
                git log -n1 HEAD | grep 'Date' | sed -e s/'Date: '/','/ | awk -- '// { printf \"%s\", $$0 }' >>${KNUT_SOURCE_DIR}/REVISION\;
                git ls-files | awk -- '/.*./ { print \"${DIST_DIR}/\"$$0}' >${KNUT_BINARY_DIR}/tmp_filelist\;
                popd\;
		echo >>tmp_filelist \"${DIST_DIR}/REVISION\"\;
		echo >>tmp_filelist \"${DIST_DIR}/knut.spec\"\;
		tar cjf ${DIST_DIR}.tar.bz2 --exclude *.svg --no-recursion -h -T tmp_filelist\;
		rm tmp_filelist ${DIST_DIR}/knut.spec ${DIST_DIR}\; )

#------------------------------------------------------------------------------
# CREATING A MAC OSX DISK IMAGE
#------------------------------------------------------------------------------

IF(APPLE)
  # FOR DETERMINING THE APPLICATION FILE PATH
  FIND_LIBRARY(COREFOUNDATION_LIBRARY CoreFoundation )
  # CREATING A MAC OSX DISK IMAGE  
  ADD_CUSTOM_TARGET( dmg
    COMMAND pushd ${CMAKE_INSTALL_PREFIX}\;
    rm -f "${KNUT_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}.dmg"\;
    hdiutil create -volname "Knut" -srcfolder Knut.app/ -srcfolder matlab/ -srcfolder demo/ "${KNUT_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}.dmg"\;
    popd\;)
  # CREATING A MAC OSX BUNDLE
  ADD_CUSTOM_TARGET( frameworks
    COMMAND rm -rf "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks"\;
    mkdir -p "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks"\;
    rm -rf "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Plugins"\;
    mkdir -p "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Plugins"\;
    cp -R ${QT_QTCORE_LIBRARY_RELEASE} "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks"\;
    cp -R ${QT_QTGUI_LIBRARY_RELEASE} "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks"\;
    cp -R ${QT_QTSVG_LIBRARY_RELEASE} "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks"\;
    cp -R ${QT_PLUGINS_DIR}/imageformats "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Plugins"\;
    "LIBNAMES=(QtCore.framework/Versions/4/QtCore" "QtGui.framework/Versions/4/QtGui" "QtSvg.framework/Versions/4/QtSvg)"\;
    rm `find "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks" -type f | grep Headers`\;
    for ix in \$\${!LIBNAMES[*]}\; do
      LOCNAMES[\$\$ix]=`otool -L gui/Knut.app/Contents/MacOS/Knut | grep \$\${LIBNAMES[\$\$ix]} | sed -e s/\\\(.*\\\)//g -e s/\\ //g`\;
      # CHANGING THE LIBRARY IDS
      install_name_tool -id "\@executable_path/../Frameworks/\$\${LIBNAMES[\$\$ix]}"
"${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks/\$\${LIBNAMES[\$\$ix]}"\;
      # CHANGING THE DEPENDENCIES in KNUT
      install_name_tool -change "\$\${LOCNAMES[\$\$ix]}" "\@executable_path/../Frameworks/\$\${LIBNAMES[\$\$ix]}"
"${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/MacOS/Knut"\;
      # REMOVING UNNECESSARY ARCHS
      lipo "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks/\$\${LIBNAMES[\$\$ix]}" -thin ${CMAKE_OSX_ARCHITECTURES} -output
"${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks/\$\${LIBNAMES[\$\$ix]}"\;
    done\;
    for ix in \$\${!LIBNAMES[*]}\; do
      for iy in \$\${!LIBNAMES[*]}\; do
        # CHANGING THE DEPENDENCIES IN THE LIBRARIES
        install_name_tool -change "\$\${LOCNAMES[\$\$ix]}" "\@executable_path/../Frameworks/\$\${LIBNAMES[\$\$ix]}"
"${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Frameworks/\$\${LIBNAMES[\$\$iy]}"\;
      done\;
      # COPYING THE PLUGINS AND CHANGING THE DEPENDENCIES
      for plg in "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Plugins/imageformats/*.dylib"\; do
        install_name_tool -change "\$\${LOCNAMES[\$\$ix]}" "\@executable_path/../Frameworks/\$\${LIBNAMES[\$\$ix]}" "\$\$plg"\;
      done\;
    done\;
    #
    for plg in "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Plugins/imageformats/*.dylib"\; do
      lipo "\$\$plg" -thin ${CMAKE_OSX_ARCHITECTURES} -output "\$\$plg"\;
    done\;
    # SETTING UP THE PLUGIN LOCATIONS
    rm -f "${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Resources/qt.conf"\;
    echo "[Paths]" >"${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Resources/qt.conf"\;
    echo "Plugins = Plugins" >"${CMAKE_INSTALL_PREFIX}/Knut.app/Contents/Resources/qt.conf"\;)
ENDIF(APPLE)

#------------------------------------------------------------------------------
# ENTERING SUBDIRECTORIES
#------------------------------------------------------------------------------

# LAPACK and ARPACK
# this produces a library "liblaarpack.a"
ADD_SUBDIRECTORY(laarpack)

# these are the UMFPACK routines
ADD_SUBDIRECTORY(UMFPACK)

# The Mini-Xml library
ADD_SUBDIRECTORY(mxml)

# This is where the whole computation is
# this produces a library "libknut.a"
ADD_SUBDIRECTORY(src)

# The command line interface
# IF(NOT APPLE)
  # NO NEED FOR THE CLI
  ADD_SUBDIRECTORY(cli)
# ENDIF(NOT APPLE)

# graphical user interface
IF(QT4_QMAKE_FOUND)
  IF(KNUT_QTMAIN)
    SET(BUILDGUI TRUE)
  ENDIF(KNUT_QTMAIN)
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
