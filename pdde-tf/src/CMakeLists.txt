ADD_LIBRARY( pddelib 
	basecolloc.h basepoint.h basepoint.cpp
	constants.cpp constants.h cspblas.c cspblas.h hypermatrix.h
	mat4data.cpp mat4data.h matrix.cpp matrix.h
	ncolloc.cpp ncolloc.h pddesys.h pderror.h
	plot.cpp plot.h point.cpp point.h pointtype.h
	polynomial.cpp polynomial.h spmatrix.cpp spmatrix.h
	system.cpp system.h testfunct.cpp testfunct.h
	torcolloc.cpp torcolloc.h torpoint.cpp torpoint.h 
	multipliers.cpp multipliers.h)

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/UMFPACK/AMD/Include
	${PDDE_SOURCE_DIR}/UMFPACK/UMFPACK/Include
	${PDDE_SOURCE_DIR}/UMFPACK/UFconfig
	${PDDE_SOURCE_DIR}/ATLAS/include
	${PDDE_SOURCE_DIR}/laarpack
	${PDDE_SOURCE_DIR}/src 
	${PDDE_BINARY_DIR}/src )

INSTALL( FILES matrix.h pddesys.h DESTINATION include/ )
