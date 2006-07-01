SET( NAME milloneup )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
    FILES sys-${NAME}.cpp 
    cfile-start cfile-ns cfile-tr )

ADD_CUSTOM_TARGET( test
    COMMAND pdde -c cfile-start -o out-start -b br-start \;
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns )
	
ADD_CUSTOM_TARGET( test-tr
    COMMAND pdde -c cfile-tr -i out-start -o out-tr -b br-tr )

ADD_DEPENDENCIES( test ${NAME} )
ADD_DEPENDENCIES( test-tr ${NAME} )
