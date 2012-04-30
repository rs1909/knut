SET( NAME tq )

SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
INSTALL( FILES ${NAME}.vf
        cfile-lp-1-a.knut cfile-lp-2-a.knut cfile-ns.knut cfile-pdsw-2.knut cfile-tr.knut
        cfile-lp-1-b.knut cfile-lp-2-b.knut cfile-pdsw-1.knut cfile-start.knut DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-pdsw-1.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-pdsw-2.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-1-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-1-b.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-2-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-7
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-lp-2-b.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-8
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-ns.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-9
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-tr.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
