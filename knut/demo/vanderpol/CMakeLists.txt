SET( NAME vanderpol )

SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
INSTALL( FILES ${NAME}.vf cfile-start.knut DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${CLI_NAME} )
