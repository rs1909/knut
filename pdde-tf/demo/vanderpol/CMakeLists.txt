SET( NAME vanderpol )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
	FILES sys-${NAME}.cpp 
	cfile-start )

ADD_CUSTOM_TARGET( test
	COMMAND pdde -c cfile-start -o out-start -b br-start )

ADD_DEPENDENCIES( test ${NAME} )
