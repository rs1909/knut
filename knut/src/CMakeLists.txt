ADD_LIBRARY( knutlib STATIC
	base58.h base58.cpp
	City.h City.cpp
	basecolloc.h basecolloc.cpp
	basepoint.h basepoint.cpp
	basecomp.h basecomp.cpp
	constants.cpp constants_compat.cpp constants.h cspblas.c cspblas.h
	hypermatrix.h mat4data.cpp mat4data.h matrix.cpp matrix.h
	ncolloc.cpp ncolloc.h knutsys.h knerror.h knerror.cpp
	plot.cpp plot.h point.cpp point.h pointtype.h pointtype.cpp
	polynomial.cpp polynomial.h spmatrix.cpp spmatrix.h
#	system.cpp system.h
	expr.cpp expr.h
	exprsystem.cpp exprsystem.h
        testfunct.cpp testfunct.h
	torcolloc.cpp torcolloc.h torpoint.cpp torpoint.h
	multipliers.cpp multipliers.h
	odecolloc.cpp odecolloc.h
	odepoint.cpp odepoint.h
	stcolloc.cpp stcolloc.h
        stpoint.cpp stpoint.h )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src
  ${KNUT_BINARY_DIR}/src
  ${KNUT_SOURCE_DIR}/laarpack
  ${KNUT_MXML_INC}
  ${UMF_INCLUDE}
  ${AMD_INCLUDE})

IF(NOT APPLE)
  INSTALL( FILES matrix.h knutsys.h pddesys.h DESTINATION include/ )
ENDIF(NOT APPLE)
