SET( NAME glass )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${PDDE_SOURCE_DIR}/demo/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp 
    cfile-start     cfile-nssw     cfile-pdsw     cfile-pd 
    cfile-start.xml cfile-nssw.xml cfile-pdsw.xml cfile-pd.xml
    DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND pdde -c ${CFILEDIR}/cfile-start -o out-start.mat -b out-start.br )
ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND pdde -c ${CFILEDIR}/cfile-nssw -i out-start.mat -o out-nssw.mat -b out-nssw.br )
ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND pdde -c ${CFILEDIR}/cfile-pdsw -i out-nssw.mat -o out-pdsw.mat -b out-pdsw.br )
ADD_CUSTOM_TARGET( ${NAME}-test-4 
    COMMAND pdde -c ${CFILEDIR}/cfile-pd -i out-nssw.mat -o out-pd.mat -b out-pd.br )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
