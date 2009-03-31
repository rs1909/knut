SET( NAME glass )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${KNUT_SOURCE_DIR}/demo/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp 
    cfile-start.knut cfile-nssw.knut cfile-pdsw.knut cfile-pd.knut
    DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND knut -c ${CFILEDIR}/cfile-start.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND knut -c ${CFILEDIR}/cfile-nssw.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND knut -c ${CFILEDIR}/cfile-pdsw.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-4 
    COMMAND knut -c ${CFILEDIR}/cfile-pd.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
