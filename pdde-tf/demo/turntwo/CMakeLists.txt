SET( NAME turntwo )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
	FILES sys-${NAME}.cpp
	cfile-start cfile-nssw cfile-ns-1 cfile-ns-2 cfile-lp-1 cfile-lp-2 )

ADD_CUSTOM_TARGET( test
	COMMAND pdde -c cfile-start -o out-start -b br-start )
ADD_CUSTOM_TARGET( test-nssw
	COMMAND pdde -c cfile-nssw -i out-start -o out-nssw -b br-nssw )
ADD_CUSTOM_TARGET( test-ns
	COMMAND pdde -c cfile-ns-1 -i out-start -o out-ns-1 -b br-ns-1 \;
	pdde -c cfile-ns-2 -i out-start -o out-ns-2 -b br-ns-2 )
ADD_CUSTOM_TARGET( test-lp
	COMMAND pdde -c cfile-lp-1 -i out-nssw -o out-lp-1 -b br-lp-1 \;
	pdde -c cfile-lp-2 -i out-nssw -o out-lp-2 -b br-lp-2 )

ADD_DEPENDENCIES( test ${NAME} )
ADD_DEPENDENCIES( test-nssw ${NAME} )
ADD_DEPENDENCIES( test-ns ${NAME} )
ADD_DEPENDENCIES( test-lp ${NAME} )
