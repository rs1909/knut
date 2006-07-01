SET( NAME glass )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
    FILES sys-${NAME}.cpp 
    cfile-start cfile-pdsw cfile-pd )

ADD_CUSTOM_TARGET( test
    COMMAND pdde -c cfile-start -o in-start -b br-start\;
	pdde -c cfile-pdsw -i in-start -o in-pdsw -b br-pdsw )
	
ADD_CUSTOM_TARGET( test-pd 
    COMMAND pdde -c cfile-pd -i in-start -o in-pd -b br-pd )

ADD_DEPENDENCIES( test ${NAME} )
ADD_DEPENDENCIES( test-pd ${NAME} )
