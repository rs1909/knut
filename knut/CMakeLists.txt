cmake_minimum_required(VERSION 3.9.6)
PROJECT( KNUT )

#------------------------------------------------------------------------------
# SETTING UP THE VARIABLES
#------------------------------------------------------------------------------

#IF(WIN32)
#    SET(BUILD_SHARED_LIBRARIES OFF)
#    SET(CMAKE_EXE_LINKER_FLAGS "-static-libgcc")
#    SET(CMAKE_LINK_SEARCH_START_STATIC TRUE)
#    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#endif()

SET( PACKAGE_NAME "Knut" )
SET( PACKAGE_VERSION "11" )
FILE( READ ${KNUT_SOURCE_DIR}/REVISION PACKAGE_REVISION )
SET( PACKAGE_COPYRIGHT "Copyright (c) 2002-2017, Robert Szalai")
SET( PACKAGE_URL "http://seis.bris.ac.uk/~rs1909/knut")

get_filename_component(CMAKE_FIND_ROOT_TMP ${CMAKE_CXX_COMPILER} PATH)
string(REPLACE "/" "\\" CMAKE_FIND_ROOT_PATH ${CMAKE_FIND_ROOT_TMP} )

enable_language (Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

STRING(TOLOWER "${PACKAGE_NAME}-cli" CLI_NAME)
IF(APPLE)
  SET(GUI_NAME ${PACKAGE_NAME})
ELSE(APPLE)
  STRING(TOLOWER "${PACKAGE_NAME}" GUI_NAME)
ENDIF(APPLE)
SET(VFKNUT_NAME "vfknut")
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
  SET( KNUT_INCLUDE_DIR include CACHE INTERNAL "include dir" )
ENDIF( NOT WIN32 )

IF( MSVC )
  ADD_DEFINITIONS(-Dinline=_inline -D_USE_MATH_DEFINES)
ENDIF( MSVC )

IF(WIN32)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
# ELSE(WIN32)
#     SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
ENDIF()

#------------------------------------------------------------------------------
# setting up Mini Xml
#------------------------------------------------------------------------------

# SET( MXML_INCLUDE_DIRS ${KNUT_SOURCE_DIR}/mxml )

#------------------------------------------------------------------------------
# finding libstdc++
#------------------------------------------------------------------------------

IF(${CMAKE_CXX_FLAGS} MATCHES ".*-m64.*")
  SET(OSX_ARCH "x86_64" CACHE INTERNAL "The build architecture" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_CXX_COMPILER} -m64 -print-file-name=libstdc++.a
    OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE )
ELSE(${CMAKE_CXX_FLAGS} MATCHES ".*-m64.*")
  SET(OSX_ARCH "i386" CACHE INTERNAL "The build architecture" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.a
    OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE )
ENDIF(${CMAKE_CXX_FLAGS} MATCHES ".*-m64.*")

MESSAGE( STATUS "The C++ runtime library is ${RES}" )
SET( KNUT_CXX_RUNTIME ${RES} CACHE INTERNAL "The fortran runtime library" )
SET( KNUT_NO_EXCEPTIONS "-fno-rtti -fno-exceptions" CACHE INTERNAL "No exceptions")

#------------------------------------------------------------------------------
# finding libgfortran
#------------------------------------------------------------------------------

IF(${CMAKE_Fortran_FLAGS} MATCHES ".*-m64.*")
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -m64 -print-file-name=libgfortran.a
    OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE )
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -m64 -print-file-name=libgcc.a
    OUTPUT_VARIABLE RESGCC OUTPUT_STRIP_TRAILING_WHITESPACE )
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -m64 -print-file-name=libquadmath.a
    OUTPUT_VARIABLE RESQUAD OUTPUT_STRIP_TRAILING_WHITESPACE )
ELSE(${CMAKE_Fortran_FLAGS} MATCHES ".*-m64.*")
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.a
    OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE )
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgcc.a
    OUTPUT_VARIABLE RESGCC OUTPUT_STRIP_TRAILING_WHITESPACE )
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libquadmath.a
    OUTPUT_VARIABLE RESQUAD OUTPUT_STRIP_TRAILING_WHITESPACE )
ENDIF(${CMAKE_Fortran_FLAGS} MATCHES ".*-m64.*")

MESSAGE( STATUS "The FORTRAN runtime library is ${RES}, ${RESGCC} and ${RESQUAD}" )
SET( KNUT_FORTRAN_RUNTIME ${RES} ${RESGCC} ${RESQUAD} CACHE INTERNAL "The fortran runtime library" )

# IF(${CMAKE_Fortran_FLAGS} MATCHES ".*-m64.*")
#   EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -m64 -print-file-name=libgfortran.so
#    OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE )
# ELSE(${CMAKE_Fortran_FLAGS} MATCHES ".*-m64.*")
#  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.so
#    OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE )
# ENDIF(${CMAKE_Fortran_FLAGS} MATCHES ".*-m64.*")

MESSAGE( STATUS "The DYNAMIC FORTRAN runtime library is ${RES}" )
SET( KNUT_FORTRAN_DYN_RUNTIME ${RES} CACHE INTERNAL "The fortran runtime library" )
get_filename_component(KNUT_FORTRAN_DYN_RUNTIME_PATH ${KNUT_FORTRAN_DYN_RUNTIME} PATH)

#------------------------------------------------------------------------------
# Find MATLAB mex compiler
#------------------------------------------------------------------------------

FIND_PROGRAM(MEX_EXECUTABLE mex)
IF( EXISTS ${MEX_EXECUTABLE} )
  IF(NOT WIN32)
    ADD_DEFINITIONS( -fPIC )
  ENDIF(NOT WIN32)
  get_filename_component(MATLAB_BIN_PATH ${MEX_EXECUTABLE} PATH CACHE)
#  SET( CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} -fPIC )
ENDIF( EXISTS ${MEX_EXECUTABLE} )

#------------------------------------------------------------------------------
# Finding openBLAS
#------------------------------------------------------------------------------

ENABLE_LANGUAGE( Fortran )

IF(NOT WIN32)
  SET(THREADLIB -pthread)
ENDIF(NOT WIN32)

FIND_LIBRARY( OPENBLAS64 NAMES openblas64 PATHS ${OPENBLAS_SEARCH_PATHS} ${KNUT_SOURCE_DIR} )
IF( EXISTS ${OPENBLAS64} )
  MESSAGE( STATUS "openBLAS64 was found as ${OPENBLAS64}." )
  SET( BLAS_LIBS ${OPENBLAS64} ${THREADLIB} )
  SET( LONGBLAS TRUE )
ELSE( EXISTS ${OPENBLAS64} )
  FIND_LIBRARY( OPENBLAS NAMES openblas PATHS ${OPENBLAS_SEARCH_PATHS} ${KNUT_SOURCE_DIR} )
  IF( EXISTS ${OPENBLAS} )
    MESSAGE( STATUS "openBLAS was found as ${OPENBLAS}." )
    SET( BLAS_LIBS ${OPENBLAS} ${THREADLIB} )
  ELSE( EXISTS ${OPENBLAS} )
    SET( USE_REFBLAS TRUE )
    SET( BLAS_LIBS -L${KNUT_BINARY_DIR}/laarpack -lREF_BLAS -lAUX_BLAS )
    MESSAGE( WARNING "Using reference BLAS. Please provide an optimised version." )
  ENDIF( EXISTS ${OPENBLAS} )
ENDIF( EXISTS ${OPENBLAS64} )

#------------------------------------------------------------------------------
# Finding system libraries
#------------------------------------------------------------------------------

# only works if normal BLAS is used
IF( EXISTS ${OPENBLAS} )
  FIND_LIBRARY( UMFPACK NAMES umfpack )
  IF( EXISTS ${UMFPACK} )
    FIND_LIBRARY( AMD NAMES amd )
    IF( EXISTS ${AMD} )
      FIND_PATH(UMF_INCLUDE NAMES umfpack.h
        PATHS CMAKE_INCLUDE_PATH
          /usr/include
          /usr/include/suitesparse/
          ${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include)
      FIND_PATH(AMD_INCLUDE NAMES amd.h
        PATHS CMAKE_INCLUDE_PATH
          /usr/include
          /usr/include/suitesparse/
          ${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include)
    ENDIF( EXISTS ${AMD} )
    MESSAGE( STATUS "Using system UMFPACK and AMD: ${UMFPACK} and ${AMD}." )
    MESSAGE( STATUS "The includes are at ${UMF_INCLUDE} and ${AMD_INCLUDE}." )
  ENDIF( EXISTS ${UMFPACK} )
  FIND_LIBRARY( LAPACK NAMES lapack )
  FIND_LIBRARY( ARPACK NAMES arpack )
  IF((EXISTS ${LAPACK}) AND (EXISTS ${ARPACK}))
    SET(LAARPACK ${LAPACK} ${ARPACK})
    MESSAGE( STATUS "Using system LAPACK and ARPACK ${LAARPACK}." )
  ENDIF((EXISTS ${LAPACK}) AND (EXISTS ${ARPACK}))
ENDIF( EXISTS ${OPENBLAS} )

IF( NOT EXISTS ${AMD} )
  SET( UMF_INCLUDE ${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include UMF_INCLUDE ${KNUT_SOURCE_DIR}/UMFPACK/UFconfig )
  SET( AMD_INCLUDE ${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include )
  SET( UMFPACK ${KNUT_BINARY_DIR}/UMFPACK/UMFPACK/Source/libumfpack.a )
  SET( AMD ${KNUT_BINARY_DIR}/UMFPACK/AMD/Source/libamd.a )
ENDIF( NOT EXISTS ${AMD} )

IF((NOT EXISTS ${LAPACK}) OR (NOT EXISTS ${ARPACK}))
  SET( LAARPACK ${KNUT_BINARY_DIR}/laarpack/liblaarpack.a )
ENDIF((NOT EXISTS ${LAPACK}) OR (NOT EXISTS ${ARPACK}))

# uncomment this once upstream fixes the MXML bug
# FIND_LIBRARY( MXML_LIB NAMES mxml )
IF( EXISTS ${MXML_LIB} )
  FIND_PATH( MXML_INC NAMES mxml.h CMAKE_INCLUDE_PATH /usr/include )
  IF( EXISTS ${MXML_INC} )
    SET( KNUT_MXML_LIB ${MXML_LIB} )
    SET( KNUT_MXML_INC ${MXML_INC} )
  ENDIF( EXISTS ${MXML_INC} )
ENDIF( EXISTS ${MXML_LIB} )
IF( NOT EXISTS ${MXML_INC} )
  SET( KNUT_MXML_LIB mxml )
  SET( KNUT_MXML_INC ${KNUT_SOURCE_DIR}/mxml )
ENDIF( NOT EXISTS ${MXML_INC} )

MESSAGE( STATUS "The mxml library is ${KNUT_MXML_LIB}." )
MESSAGE( STATUS "The mxml include directory is ${KNUT_MXML_INC}." )

#------------------------------------------------------------------------------
# Finding Qt
#------------------------------------------------------------------------------

find_package(Qt5Core)
if(${Qt5Core_FOUND})
  message(STATUS "--> Qt5Core found.")
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  # Instruct CMake to run moc automatically when needed.
  find_package(Qt5Widgets REQUIRED)
  find_package(Qt5PrintSupport REQUIRED)
  find_package(Qt5Svg REQUIRED)
#  find_package(Qt5Charts REQUIRED)

  if(WIN32)
    get_target_property(QtCore_location Qt5::Core LOCATION)
    get_target_property(QtGui_location Qt5::Gui LOCATION)
    get_target_property(QtWidgets_location Qt5::Widgets LOCATION)
    get_target_property(QtPrintSupport_location Qt5::PrintSupport LOCATION)
    get_target_property(QtSvg_location Qt5::Svg LOCATION)
#    get_target_property(QtCharts_location Qt5::Charts LOCATION)
    string(REPLACE "/" "\\" QtCoreLoc ${QtCore_location} )
    string(REPLACE "/" "\\" QtGuiLoc ${QtGui_location} )
    string(REPLACE "/" "\\" QtWidgetsLoc ${QtWidgets_location} )
    string(REPLACE "/" "\\" QtPrintSupportLoc ${QtPrintSupport_location} )
    string(REPLACE "/" "\\" QtSvgLoc ${QtSvg_location} )
#    string(REPLACE "/" "\\" QtChartsLoc ${QtCharts_location} )
    get_filename_component(QtCorePath_TMP ${QtCore_location} PATH)
    string(REPLACE "/" "\\" QtCorePath ${QtCorePath_TMP} )
  endif(WIN32)
endif(${Qt5Core_FOUND})


#------------------------------------------------------------------------------
# Writing the generated files
#------------------------------------------------------------------------------

CONFIGURE_FILE( src/config.h.cmake src/config.h )
CONFIGURE_FILE( knut.spec.cmake knut.spec )
CONFIGURE_FILE( gui/MacOS/Info.plist.cmake gui/MacOS/Info.plist )
CONFIGURE_FILE( doc/html/index.html.cmake doc/html/index.html )
CONFIGURE_FILE( bin/pcompile.cmake bin/pcompile )
IF(WIN32)
  CONFIGURE_FILE( wininstaller.nsi.cmake wininstaller.nsi )
#  CONFIGURE_FILE( w64installer.nsi.cmake w64installer.nsi )
ENDIF(WIN32)

IF( NOT WIN32 )
  INSTALL( PROGRAMS ${KNUT_BINARY_DIR}/bin/pcompile DESTINATION ${KNUT_BIN_DIR} )
ENDIF( NOT WIN32 )

# no debugging with release version
IF( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )
  ADD_DEFINITIONS( -DHAVE_CONFIG_H -DDEBUG )
  SET(KNUT_RELEASE FALSE)
  SET(CLI_NAME "${CLI_NAME}-debug")
  SET(GUI_NAME "${GUI_NAME}-debug")
  SET(VFKNUT_NAME "${VFKNUT_NAME}-debug")
ELSE( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )
  ADD_DEFINITIONS( -DHAVE_CONFIG_H )
  SET(KNUT_RELEASE TRUE)
ENDIF( CMAKE_BUILD_TYPE MATCHES "[D,d][E,e][B,b][U,u][G,g]" )

#------------------------------------------------------------------------------
# DISTRIBUTING WITH GIT
#------------------------------------------------------------------------------

SET( DIST_DIR ${PACKAGE_NAME}-${PACKAGE_VERSION} )
ADD_CUSTOM_TARGET( dist
	COMMAND ln -s ${KNUT_SOURCE_DIR} ${DIST_DIR}\;
		ln -s ${KNUT_BINARY_DIR}/knut.spec ${DIST_DIR}/knut.spec\;
		pushd ${KNUT_SOURCE_DIR}\;
                ./mkrevision.sh\;
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
    rm -f "${KNUT_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}.dmg" \;
    hdiutil create -volname "Knut" -srcfolder Knut.app
      -srcfolder demo/
      "${KNUT_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}.dmg" \;
    popd\;)
ENDIF(APPLE)

#------------------------------------------------------------------------------
# ENTERING SUBDIRECTORIES
#------------------------------------------------------------------------------

# LAPACK and ARPACK
# this produces a library "liblaarpack.a"
IF((NOT EXISTS ${LAPACK}) OR (NOT EXISTS ${ARPACK}))
  ADD_SUBDIRECTORY(laarpack)
ENDIF((NOT EXISTS ${LAPACK}) OR (NOT EXISTS ${ARPACK}))

# these are the UMFPACK routines
IF(NOT EXISTS ${AMD})
  ADD_SUBDIRECTORY(UMFPACK)
ENDIF(NOT EXISTS ${AMD})

# The Mini-Xml library
IF( NOT EXISTS ${MXML_INC} )
  ADD_SUBDIRECTORY(mxml)
ENDIF( NOT EXISTS ${MXML_INC} )

# This is where the whole computation is
# this produces a library "libknut.a"
ADD_SUBDIRECTORY(src)

# The command line interface
ADD_SUBDIRECTORY(cli)

# graphical user interface
if(${Qt5Core_FOUND})
ADD_SUBDIRECTORY(gui)
endif(${Qt5Core_FOUND})

# demo programs
ADD_SUBDIRECTORY(demo)

# MATLAB scripts (disable, it does not work)
# IF( EXISTS ${MEX_EXECUTABLE} )
#   ADD_SUBDIRECTORY(matlab)
# ENDIF( EXISTS ${MEX_EXECUTABLE} )
