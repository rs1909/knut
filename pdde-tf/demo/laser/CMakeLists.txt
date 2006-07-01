SET( NAME laser )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
    FILES sys-${NAME}.cpp
    cfile-start cfile-ns cfile-pd cfile-lp in-start )

ADD_CUSTOM_TARGET( test
    COMMAND pdde -c cfile-start -i in-start -o out-start -b br-start )
	
ADD_CUSTOM_TARGET( test-ns
    COMMAND pdde -c cfile-ns -i out-start -o out-ns -b branch-ns )

ADD_CUSTOM_TARGET( test-pd
    COMMAND pdde -c cfile-pd -i out-start -o out-pd -b branch-pd )

ADD_CUSTOM_TARGET( test-lp
    COMMAND pdde -c cfile-lp -i out-start -o out-lp -b branch-lp )

ADD_DEPENDENCIES( test ${NAME} )
ADD_DEPENDENCIES( test-ns ${NAME} )
ADD_DEPENDENCIES( test-pd ${NAME} )
ADD_DEPENDENCIES( test-lp ${NAME} )
