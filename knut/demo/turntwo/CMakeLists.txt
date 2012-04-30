SET( NAME turntwo )

SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
INSTALL( FILES ${NAME}.vf
	cfile-start.knut cfile-nssw.knut cfile-ns-1.knut
	cfile-ns-2.knut cfile-lp-1.knut cfile-lp-2.knut DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-2
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-nssw.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-3
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-ns-1.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-4
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-ns-2.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-5
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-1.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-6
	COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-2.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
