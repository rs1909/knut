SET( NAME traffic3cars )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
	FILES sys-${NAME}.cpp
	cfile-start-a cfile-start-b
	cfile-ns-a cfile-ns-b
	cfile-nssw-a cfile-nssw-b
	cfile-lp-a cfile-lp-b
	README sys-trafficncars.cpp )

ADD_CUSTOM_TARGET( ${NAME}-test1a
    COMMAND pdde -c cfile-start-a -o out-start-a -b br-start-a )

ADD_CUSTOM_TARGET( ${NAME}-test1b
    COMMAND pdde -c cfile-start-b -o out-start-b -b br-start-b )

ADD_CUSTOM_TARGET( ${NAME}-test2a
    COMMAND pdde -c cfile-nssw-a -i out-start-a -o out-nssw-a -b br-nssw-a )

ADD_CUSTOM_TARGET( ${NAME}-test2b
    COMMAND pdde -c cfile-nssw-b -i out-start-b -o out-nssw-b -b br-nssw-b )

ADD_CUSTOM_TARGET( ${NAME}-test3a
    COMMAND pdde -c cfile-lp-a -i out-nssw-a -o out-lp-a -b br-lp-a )

ADD_CUSTOM_TARGET( ${NAME}-test3b
    COMMAND pdde -c cfile-lp-b -i out-nssw-b -o out-lp-b -b br-lp-b )

ADD_DEPENDENCIES( ${NAME}-test1a ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test1b ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test2a ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test2b ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test3a ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test3b ${NAME} )
