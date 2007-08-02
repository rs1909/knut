SET( NAME traffic3cars )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${PDDE_SOURCE_DIR}/demo/traffic )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
	cfile-start-a cfile-start-b
	cfile-ns-a cfile-ns-b
	cfile-nssw-a cfile-nssw-b
	cfile-lp-a cfile-lp-b
	cfile-start-a.xml cfile-start-b.xml
	cfile-ns-a.xml cfile-ns-b.xml
	cfile-nssw-a.xml cfile-nssw-b.xml
	cfile-lp-a.xml cfile-lp-b.xml
	README sys-trafficncars.cpp DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND pdde -c ${CFILEDIR}/cfile-start-a -o out-start-a.mat -b out-start-a.br )

ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND pdde -c ${CFILEDIR}/cfile-start-b -o out-start-b.mat -b out-start-b.br )

ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND pdde -c ${CFILEDIR}/cfile-nssw-a -i out-start-a.mat -o out-nssw-a.mat -b out-nssw-a.br )

ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND pdde -c ${CFILEDIR}/cfile-nssw-b -i out-start-b.mat -o out-nssw-b.mat -b out-nssw-b.br )

ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND pdde -c ${CFILEDIR}/cfile-lp-a -i out-nssw-a.mat -o out-lp-a.mat -b out-lp-a.br )

ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND pdde -c ${CFILEDIR}/cfile-lp-b -i out-nssw-b.mat -o out-lp-b.mat -b out-lp-b.br )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
