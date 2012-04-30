SET( NAME traffic3cars )

SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
INSTALL( FILES ${NAME}.vf
	cfile-start-a.knut cfile-start-b.knut
	cfile-ns-a.knut cfile-ns-b.knut
	cfile-nssw-a.knut cfile-nssw-b.knut
	cfile-lp-a.knut cfile-lp-b.knut
	README sys-trafficncars.cpp DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start-b.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-nssw-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-nssw-b.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-b.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
