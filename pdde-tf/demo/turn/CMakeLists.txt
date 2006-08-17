SET( NAME turn )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
	FILES sys-${NAME}.cpp 
	cfile-start cfile-ns cfile-nssw
	cfile-start.xml cfile-ns.xml cfile-nssw.xml )

ADD_CUSTOM_TARGET( ${NAME}-test
	COMMAND pdde -c cfile-start -o out-start -b br-start \;
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns \;
	pdde -c cfile-nssw -i out-start -o out-nssw -b br-nssw )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
