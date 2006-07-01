SET( BLASINC  ${PDDE_SOURCE_DIR}/ATLAS/include )
INCLUDE_DIRECTORIES( ${BLASINC} )
ADD_LIBRARY( laarpack laarpack.c dlamch.o cblaswrap.h f2c.h )

SET_SOURCE_FILES_PROPERTIES( dlamch.o PROPERTIES GENERATED TRUE )
ADD_CUSTOM_TARGET( dlamch.o
    COMMAND ${CMAKE_C_COMPILER} -O0 -I${BLASINC} -c dlamch.c -o dlamch.o )
ADD_DEPENDENCIES( laarpack dlamch.o )
