SET( NAME turntwo )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${PDDE_SOURCE_DIR}/demo/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
	cfile-start cfile-nssw cfile-ns-1 cfile-ns-2 cfile-lp-1 cfile-lp-2 
	cfile-start.xml cfile-nssw.xml cfile-ns-1.xml
	cfile-ns-2.xml cfile-lp-1.xml cfile-lp-2.xml DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
	COMMAND pdde -c ${CFILEDIR}/cfile-start -o out-start.mat -b out-start.br )

ADD_CUSTOM_TARGET( ${NAME}-test-2
	COMMAND pdde -c ${CFILEDIR}/cfile-nssw -i out-start.mat -o out-nssw.mat -b out-nssw.br )

ADD_CUSTOM_TARGET( ${NAME}-test-3
	COMMAND pdde -c ${CFILEDIR}/cfile-ns-1 -i out-start.mat -o out-ns-1.mat -b out-ns-1.br )

ADD_CUSTOM_TARGET( ${NAME}-test-4
	COMMAND pdde -c ${CFILEDIR}/cfile-ns-2 -i out-start.mat -o out-ns-2.mat -b out-ns-2.br )

ADD_CUSTOM_TARGET( ${NAME}-test-5
	COMMAND pdde -c ${CFILEDIR}/cfile-lp-1 -i out-nssw.mat -o out-lp-1.mat -b out-lp-1.br )

ADD_CUSTOM_TARGET( ${NAME}-test-6
	COMMAND pdde -c ${CFILEDIR}/cfile-lp-2 -i out-nssw.mat -o out-lp-2.mat -b out-lp-2.br )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
