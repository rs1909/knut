SET( NAME turntwo )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${CMAKE_INSTALL_PREFIX}/${KNUT_DEMO_DIR}/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
	cfile-start.knut cfile-nssw.knut cfile-ns-1.knut
	cfile-ns-2.knut cfile-lp-1.knut cfile-lp-2.knut DESTINATION ${KNUT_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
	COMMAND exec knut -c ${CFILEDIR}/cfile-start.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-2
	COMMAND exec knut -c ${CFILEDIR}/cfile-nssw.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-3
	COMMAND exec knut -c ${CFILEDIR}/cfile-ns-1.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-4
	COMMAND exec knut -c ${CFILEDIR}/cfile-ns-2.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-5
	COMMAND exec knut -c ${CFILEDIR}/cfile-lp-1.knut )

ADD_CUSTOM_TARGET( ${NAME}-test-6
	COMMAND exec knut -c ${CFILEDIR}/cfile-lp-2.knut )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
