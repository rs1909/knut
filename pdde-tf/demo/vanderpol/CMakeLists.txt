SET( NAME vanderpol )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${PDDE_SOURCE_DIR}/demo/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp 
	cfile-start cfile-start.xml DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
	COMMAND pdde -c ${CFILEDIR}/cfile-start -o out-start.mat -b br-start )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
