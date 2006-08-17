SET( NAME glass )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
    FILES sys-${NAME}.cpp 
    cfile-start cfile-pdsw cfile-pd 
    cfile-start.xml cfile-pdsw.xml cfile-pd.xml )

ADD_CUSTOM_TARGET( ${NAME}-test
    COMMAND pdde -c cfile-start -o in-start -b br-start\;
	pdde -c cfile-pdsw -i in-start -o in-pdsw -b br-pdsw )
	
ADD_CUSTOM_TARGET( ${NAME}-test-pd 
    COMMAND pdde -c cfile-pd -i in-start -o in-pd -b br-pd )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-pd ${NAME} )
