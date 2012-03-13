SET( NAME turn )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp ${NAME}.vf
	cfile-start.knut cfile-ns.knut cfile-nssw.knut DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-2
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-ns.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-3
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-nssw.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
