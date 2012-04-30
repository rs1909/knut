SET( NAME glass )

SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
INSTALL( FILES mackey-glass.vf
    cfile-start.knut cfile-nssw.knut cfile-pdsw.knut cfile-pd.knut
    DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-start.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-nssw.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-pdsw.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-4 
    COMMAND exec ${CLI_NAME} -c ${CFILEDIR}/cfile-pd.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
