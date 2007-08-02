SET( NAME tq )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${PDDE_SOURCE_DIR}/demo/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
	cfile-lp-1-a cfile-lp-2-a cfile-ns cfile-pdsw-2 cfile-tr
	cfile-lp-1-b cfile-lp-2-b cfile-pdsw-1 cfile-start
        cfile-lp-1-a.xml cfile-lp-2-a.xml cfile-ns.xml cfile-pdsw-2.xml cfile-tr.xml
        cfile-lp-1-b.xml cfile-lp-2-b.xml cfile-pdsw-1.xml cfile-start.xml DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND pdde -c ${CFILEDIR}/cfile-start -o out-start.mat -b out-start.br )

ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND pdde -c ${CFILEDIR}/cfile-pdsw-1 -i out-start.mat -o out-pdsw-1.mat -b out-pdsw-1.br )

ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND pdde -c ${CFILEDIR}/cfile-pdsw-2 -i out-start.mat -o out-pdsw-2.mat -b out-pdsw-2.br )

ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND pdde -c ${CFILEDIR}/cfile-lp-1-a -i out-pdsw-1.mat -o out-lp-1-a.mat -b out-lp-1-a.br )

ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND pdde -c ${CFILEDIR}/cfile-lp-1-b -i out-pdsw-1.mat -o out-lp-1-b.mat -b out-lp-1-b.br )

ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND pdde -c ${CFILEDIR}/cfile-lp-2-a -i out-pdsw-1.mat -o out-lp-2-a.mat -b out-lp-2-a.br )

ADD_CUSTOM_TARGET( ${NAME}-test-7
    COMMAND pdde -c ${CFILEDIR}/cfile-lp-2-b -i out-pdsw-1.mat -o out-lp-2-b.mat -b out-lp-2-b.br )

ADD_CUSTOM_TARGET( ${NAME}-test-8
    COMMAND pdde -c ${CFILEDIR}/cfile-ns -i out-start.mat -o out-ns.mat -b out-ns.br )

ADD_CUSTOM_TARGET( ${NAME}-test-9
    COMMAND pdde -c ${CFILEDIR}/cfile-tr -i out-start.mat -o out-tr.mat -b out-tr.br )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )