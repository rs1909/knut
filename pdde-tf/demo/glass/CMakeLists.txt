SET( NAME glass )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp 
    cfile-start cfile-nssw cfile-pdsw cfile-pd 
    cfile-start.xml cfile-nssw.xml cfile-pdsw.xml cfile-pd.xml
    DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test
    COMMAND pdde -c cfile-start -o out-start.mat -b br-start\;
	pdde -c cfile-nssw -i out-start.mat -o out-nssw.mat -b br-nssw\;
	pdde -c cfile-pdsw -i out-nssw.mat -o out-pdsw.mat -b br-pdsw )
	
ADD_CUSTOM_TARGET( ${NAME}-test-pd 
    COMMAND pdde -c cfile-pd -i out-nssw.mat -o out-pd.mat -b br-pd )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-pd ${NAME} )
