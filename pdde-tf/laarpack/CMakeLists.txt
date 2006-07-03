SET( BLASINC  ${PDDE_SOURCE_DIR}/ATLAS/include )
SET( F2CINC ${PDDE_SOURCE_DIR}/libf2c )
INCLUDE_DIRECTORIES( ${BLASINC} ${F2CINC} )

ADD_LIBRARY( laarpack laarpack.c dlamch.o cblaswrap.h dlamch.c )
SET_SOURCE_FILES_PROPERTIES( dlamch.c PROPERTIES COMPILE_FLAGS -ffloat-store )
