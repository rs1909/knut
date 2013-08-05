LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )

SET( MEX_INCLUDES 
	-I${KNUT_SOURCE_DIR}/src
	-I${KNUT_BINARY_DIR}/src
	-I${KNUT_SOURCE_DIR}/vfgen
	-I${GINAC_INCLUDE_DIRS}
	-I${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include
	-I${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include
	-I${KNUT_SOURCE_DIR}/UMFPACK/UFconfig
	-I${KNUT_SOURCE_DIR}/ATLAS/include
	-I${KNUT_SOURCE_DIR}/laarpack )

SET(INTERNAL_LIBS -L${KNUT_BINARY_DIR}/src -L${KNUT_BINARY_DIR}/laarpack
-L${KNUT_BINARY_DIR}/vfgen -L${KNUT_BINARY_DIR}/UMFPACK/AMD/Source -L${KNUT_BINARY_DIR}/UMFPACK/UMFPACK/Source -L${KNUT_BINARY_DIR}/mxml -lknutlib -llaarpack -lumfpack -lamd -lREF_BLAS -lAUX_BLAS -l${VFGEN_LIBS} -lmxml)
SET(EXTERNAL_LIBS ${GINAC_LDFLAGS} ${BLAS_LIBS}
	-l${CMAKE_DL_LIBS})

IF(${CMAKE_DL_LIBS} MATCHES ".*lib.*")
	SET(DYNLOAD_LIB -l${CMAKE_DL_LIBS})
ENDIF(${CMAKE_DL_LIBS} MATCHES ".*lib.*")

string(REPLACE " " ";" KNUT_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})

add_custom_target(knut.o COMMAND
${CMAKE_CXX_COMPILER} -std=c++11 -stdlib=libc++ ${KNUT_CXX_FLAGS} -fPIC -c ${KNUT_SOURCE_DIR}/matlab/knut.cpp -I${KNUT_SOURCE_DIR}/src/ -I${KNUT_SOURCE_DIR}/laarpack/ -I${KNUT_BINARY_DIR}/src -I${MATLAB_BIN_PATH}/../extern/include/ -o knut.o
DEPENDS ${CLI_NAME})

IF(APPLE)
  SET(MEX_SUFF maci64)
  SET(DIR_SUFF maci64)
  SET(CFFRAME LDFLAGS=\"-framework CoreFoundation -undefined dynamic_lookup -bundle\")
ELSE(APPLE)
  SET(MEX_SUFF a64)
  SET(DIR_SUFF glnxa64)
ENDIF(APPLE)

IF( LONGBLAS )
  SET(MWBLAS -lmwblas)
ELSE( LONGBLAS )
  SET(MWBLAS -lmwblascompat32)
ENDIF( LONGBLAS )

add_custom_target(knut.mex COMMAND mex -cxx CXX='clang++' CFLAGS='-std=c++11 -stdlib=libc++ ${KNUT_CXX_FLAGS}' knut.o 
${KNUT_BINARY_DIR}/src/libknutlib.a ${KNUT_BINARY_DIR}/UMFPACK/UMFPACK/Source/libumfpack.a ${KNUT_BINARY_DIR}/UMFPACK/AMD/Source/libamd.a ${KNUT_BINARY_DIR}/mxml/libmxml.a ${KNUT_BINARY_DIR}/laarpack/liblaarpack.a 
-lmwblas -lut -lgfortran -ldl ${CFFRAME}
DEPENDS knut.o)

# add_custom_target(knut.mex COMMAND 
# ${CMAKE_CXX_COMPILER} -std=c++11 ${KNUT_CXX_FLAGS} -shared knut.o -Wl,-rpath-link,${MATLAB_BIN_PATH}/${MEX_SUFF}/ -L${MATLAB_BIN_PATH}/${DIR_SUFF}/ -lmx -lmex -lut ${KNUT_BINARY_DIR}/src/libknutlib.a ${KNUT_BINARY_DIR}/vfgen/libvfgenlib.a ${GINAC_LDFLAGS} ${KNUT_BINARY_DIR}/UMFPACK/UMFPACK/Source/libumfpack.a ${KNUT_BINARY_DIR}/UMFPACK/AMD/Source/libamd.a ${KNUT_BINARY_DIR}/mxml/libmxml.a ${KNUT_BINARY_DIR}/laarpack/liblaarpack.a 
# ${KNUT_BINARY_DIR}/laarpack/libREF_BLAS.a ${KNUT_BINARY_DIR}/laarpack/libAUX_BLAS.a
# -lgfortran -lc++ ${CFFRAME} -ldl -o knut.mex${MEX_SUFF}
# DEPENDS knut.o)

# ${KNUT_BINARY_DIR}/laarpack/libREF_BLAS.a ${KNUT_BINARY_DIR}/laarpack/libAUX_BLAS.a 
# ${BLAS_LIBS}
INSTALL( FILES branch2knmat.m knmat2branch.m dde2knut.m knut2dde.m torplot.m
         DESTINATION ${KNUT_MATLAB_DIR} )
