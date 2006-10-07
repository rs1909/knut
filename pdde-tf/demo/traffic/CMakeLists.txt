SET( NAME traffic3cars )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
	cfile-start-a cfile-start-b
	cfile-ns-a cfile-ns-b
	cfile-nssw-a cfile-nssw-b
	cfile-lp-a cfile-lp-b
	cfile-start-a.xml cfile-start-b.xml
	cfile-ns-a.xml cfile-ns-b.xml
	cfile-nssw-a.xml cfile-nssw-b.xml
	cfile-lp-a.xml cfile-lp-b.xml
	README sys-trafficncars.cpp DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test1a
    COMMAND pdde -c cfile-start-a -o out-start-a.mat -b br-start-a )

ADD_CUSTOM_TARGET( ${NAME}-test1b
    COMMAND pdde -c cfile-start-b -o out-start-b.mat -b br-start-b )

ADD_CUSTOM_TARGET( ${NAME}-test2a
    COMMAND pdde -c cfile-nssw-a -i out-start-a.mat -o out-nssw-a.mat -b br-nssw-a )

ADD_CUSTOM_TARGET( ${NAME}-test2b
    COMMAND pdde -c cfile-nssw-b -i out-start-b.mat -o out-nssw-b.mat -b br-nssw-b )

ADD_CUSTOM_TARGET( ${NAME}-test3a
    COMMAND pdde -c cfile-lp-a -i out-nssw-a.mat -o out-lp-a.mat -b br-lp-a )

ADD_CUSTOM_TARGET( ${NAME}-test3b
    COMMAND pdde -c cfile-lp-b -i out-nssw-b.mat -o out-lp-b.mat -b br-lp-b )

ADD_DEPENDENCIES( ${NAME}-test1a ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test1b ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test2a ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test2b ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test3a ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test3b ${NAME} )
