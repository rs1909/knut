SET( NAME laser )

SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
INSTALL( FILES ${NAME}.vf
    cfile-start.knut cfile-hopfsw.knut cfile-psol.knut cfile-ns.knut cfile-pd.knut cfile-lp.knut
    DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND exec knut -c ${CFILEDIR}/cfile-start.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND exec knut -c ${CFILEDIR}/cfile-hopfsw.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND exec knut -c ${CFILEDIR}/cfile-psol.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND exec knut -c ${CFILEDIR}/cfile-ns.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND exec knut -c ${CFILEDIR}/cfile-pd.knut )
ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND exec knut -c ${CFILEDIR}/cfile-lp.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
