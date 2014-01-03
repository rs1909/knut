LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src
  ${KNUT_BINARY_DIR}/src
  ${KNUT_SOURCE_DIR}/laarpack
  ${KNUT_MXML_INC}
  ${UMF_INCLUDE}
  ${AMD_INCLUDE})

SET( MEX_INCLUDES
  -I${KNUT_SOURCE_DIR}/src
  -I${KNUT_BINARY_DIR}/src
  -I${KNUT_SOURCE_DIR}/laarpack
  -I${KNUT_MXML_INC}
  -I${UMF_INCLUDE}
  -I${AMD_INCLUDE})

SET(INTERNAL_LIBS -L${KNUT_BINARY_DIR}/src -L${KNUT_BINARY_DIR}/laarpack
-L${KNUT_BINARY_DIR}/UMFPACK/AMD/Source -L${KNUT_BINARY_DIR}/UMFPACK/UMFPACK/Source -lknutlib -llaarpack -lumfpack -lamd -lREF_BLAS -lAUX_BLAS -l${KNUT_MXML_LIB})
SET(EXTERNAL_LIBS ${GINAC_LDFLAGS} ${BLAS_LIBS}
	-l${CMAKE_DL_LIBS})

IF(${CMAKE_DL_LIBS} MATCHES ".*lib.*")
	SET(DYNLOAD_LIB -l${CMAKE_DL_LIBS})
ENDIF(${CMAKE_DL_LIBS} MATCHES ".*lib.*")

# string(REPLACE " " ";" KNUT_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})

IF(APPLE)
  SET(MEX_SUFF maci64)
  SET(DIR_SUFF maci64)
  SET(CFFRAME LDFLAGS=\"-framework CoreFoundation -undefined dynamic_lookup -bundle\")
  SET(STDLIB "-stdlib=libc++")
ELSE(APPLE)
  SET(MEX_SUFF a64)
  SET(DIR_SUFF glnxa64)
ENDIF(APPLE)

IF( LONGBLAS )
  SET(MWBLAS -largeArrayDims)
  MESSAGE(STATUS "Using LONG BLAS")
ELSE( LONGBLAS )
  SET(MWBLAS -compatibleArrayDims)
  MESSAGE(STATUS "Using regular BLAS")
ENDIF( LONGBLAS )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src
        ${KNUT_BINARY_DIR}/src
        ${KNUT_SOURCE_DIR}/laarpack
        ${MATLAB_BIN_PATH}/../extern/include/ )

ADD_LIBRARY( mxknut STATIC knut.cpp mxdata.cpp mxdata.h )
SET_TARGET_PROPERTIES( mxknut PROPERTIES COMPILE_FLAGS
"-DMATLAB_MEX_FILE")

# add_custom_target(knut.o COMMAND
# ${CMAKE_CXX_COMPILER} -DMATLAB_MEX_FILE -std=c++11 ${STDLIB} ${KNUT_CXX_FLAGS} -fPIC -c ${KNUT_SOURCE_DIR}/matlab/knut.cpp -I${KNUT_SOURCE_DIR}/src/ -I${KNUT_SOURCE_DIR}/laarpack/ -I${KNUT_BINARY_DIR}/src -I${MATLAB_BIN_PATH}/../extern/include/ -o knut.o
# DEPENDS ${CLI_NAME})

add_custom_target(knut.mex COMMAND
${CMAKE_AR} x libmxknut.a \;
mex -cxx ${MWBLAS} -output knut CXX='${CMAKE_CXX_COMPILER}' CFLAGS='-std=c++11 ${STDLIB} ${KNUT_CXX_FLAGS}'
knut.cpp.o mxdata.cpp.o
${KNUT_BINARY_DIR}/src/libknutlib.a ${UMFPACK} ${AMD} ${KNUT_MXML_LIB} ${LAARPACK}
-lmwblas -lut -lgfortran -ldl ${CFFRAME}
DEPENDS mxknut)

# add_custom_target(knut.mex COMMAND
# ${CMAKE_CXX_COMPILER} -std=c++11 ${KNUT_CXX_FLAGS} -shared knut.o -Wl,-rpath-link,${MATLAB_BIN_PATH}/${MEX_SUFF}/ -L${MATLAB_BIN_PATH}/${DIR_SUFF}/ -lmx -lmex -lut ${KNUT_BINARY_DIR}/src/libknutlib.a ${KNUT_BINARY_DIR}/vfgen/libvfgenlib.a ${GINAC_LDFLAGS} ${KNUT_BINARY_DIR}/UMFPACK/UMFPACK/Source/libumfpack.a ${KNUT_BINARY_DIR}/UMFPACK/AMD/Source/libamd.a ${KNUT_BINARY_DIR}/mxml/libmxml.a ${KNUT_BINARY_DIR}/laarpack/liblaarpack.a
# ${KNUT_BINARY_DIR}/laarpack/libREF_BLAS.a ${KNUT_BINARY_DIR}/laarpack/libAUX_BLAS.a
# -lgfortran -lc++ ${CFFRAME} -ldl -o knut.mex${MEX_SUFF}
# DEPENDS knut.o)

# ${KNUT_BINARY_DIR}/laarpack/libREF_BLAS.a ${KNUT_BINARY_DIR}/laarpack/libAUX_BLAS.a
# ${BLAS_LIBS}
INSTALL( FILES branch2knmat.m knmat2branch.m dde2knut.m knut2dde.m torplot.m
         DESTINATION ${KNUT_MATLAB_DIR} )
