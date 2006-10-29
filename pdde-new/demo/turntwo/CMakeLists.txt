SET( NAME turntwo )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
	cfile-start cfile-nssw cfile-ns-1 cfile-ns-2 cfile-lp-1 cfile-lp-2 
	cfile-start.xml cfile-nssw.xml cfile-ns-1.xml
	cfile-ns-2.xml cfile-lp-1.xml cfile-lp-2.xml DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test
	COMMAND pdde -c cfile-start -o out-start.mat -b br-start )
ADD_CUSTOM_TARGET( ${NAME}-test-nssw
	COMMAND pdde -c cfile-nssw -i out-start.mat -o out-nssw.mat -b br-nssw )
ADD_CUSTOM_TARGET( ${NAME}-test-ns
	COMMAND pdde -c cfile-ns-1 -i out-start.mat -o out-ns-1.mat -b br-ns-1 \;
	pdde -c cfile-ns-2 -i out-start.mat -o out-ns-2.mat -b br-ns-2 )
ADD_CUSTOM_TARGET( ${NAME}-test-lp
	COMMAND pdde -c cfile-lp-1 -i out-nssw.mat -o out-lp-1.mat -b br-lp-1 \;
	pdde -c cfile-lp-2 -i out-nssw.mat -o out-lp-2.mat -b br-lp-2 )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-nssw ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-ns ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-lp ${NAME} )
