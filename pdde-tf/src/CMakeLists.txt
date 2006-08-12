ADD_LIBRARY( pddelib plot.cpp ncolloc.cpp torcolloc.cpp point.cpp torpoint.cpp 
	matrix.cpp cspblas.c cspblas.h spmatrix.cpp testfunct.cpp
	hypermatrix.h matrix.h ncolloc.h torcolloc.h
	pddesys.h plot.h point.h torpoint.h spmatrix.h testfunct.h
	system.cpp system.h pointtype.h pderror.h pderror.cpp
	mat4data.cpp mat4data.h constants.cpp constants.h )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/UMFPACK/AMD/Include
	${PDDE_SOURCE_DIR}/UMFPACK/UMFPACK/Include
	${PDDE_SOURCE_DIR}/ATLAS/include
	${PDDE_SOURCE_DIR}/src )

INSTALL_FILES( /include FILES matrix.h pddesys.h )
