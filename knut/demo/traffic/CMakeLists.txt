SET( NAME traffic3cars )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp ${NAME}.vf
	cfile-start-a.knut cfile-start-b.knut
	cfile-ns-a.knut cfile-ns-b.knut
	cfile-nssw-a.knut cfile-nssw-b.knut
	cfile-lp-a.knut cfile-lp-b.knut
	README sys-trafficncars.cpp DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND exec knut -c ${CFILEDIR}/cfile-start-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND exec knut -c ${CFILEDIR}/cfile-start-b.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND exec knut -c ${CFILEDIR}/cfile-nssw-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND exec knut -c ${CFILEDIR}/cfile-nssw-b.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND exec knut -c ${CFILEDIR}/cfile-lp-a.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND exec knut -c ${CFILEDIR}/cfile-lp-b.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
