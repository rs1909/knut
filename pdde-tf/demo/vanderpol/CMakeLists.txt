SET( NAME vanderpol )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
	FILES sys-${NAME}.cpp 
	cfile-start cfile-start.xml )

ADD_CUSTOM_TARGET( ${NAME}-test
	COMMAND pdde -c cfile-start -o out-start -b br-start )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
