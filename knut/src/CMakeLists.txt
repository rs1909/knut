ADD_LIBRARY( knutlib 
	basecolloc.h basecolloc.cpp
	basepoint.h basepoint.cpp
	basecomp.h basecomp.cpp
	constants.cpp constants.h cspblas.c cspblas.h hypermatrix.h
	mat4data.cpp mat4data.h matrix.cpp matrix.h
	ncolloc.cpp ncolloc.h knutsys.h knerror.h knerror.cpp
	plot.cpp plot.h point.cpp point.h pointtype.h
	polynomial.cpp polynomial.h spmatrix.cpp spmatrix.h
	system.cpp system.h testfunct.cpp testfunct.h
	torcolloc.cpp torcolloc.h torpoint.cpp torpoint.h 
	multipliers.cpp multipliers.h
	odecolloc.cpp odecolloc.h
	odepoint.cpp odepoint.h )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include
	${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include
	${KNUT_SOURCE_DIR}/UMFPACK/UFconfig
	${KNUT_SOURCE_DIR}/ATLAS/include
	${KNUT_SOURCE_DIR}/laarpack
	${KNUT_SOURCE_DIR}/vfgen
	${KNUT_SOURCE_DIR}/src 
	${KNUT_BINARY_DIR}/src 
	${MXML_INCLUDE_DIRS})

IF(NOT APPLE)
  INSTALL( FILES matrix.h knutsys.h pddesys.h DESTINATION include/ )
ENDIF(NOT APPLE)
