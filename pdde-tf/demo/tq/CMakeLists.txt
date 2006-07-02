SET( NAME tq )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL_FILES( /share/${PACKAGE_NAME}/demo/${NAME}
	FILES sys-${NAME}.cpp
	cfile-lp-1-a cfile-lp-2-a cfile-ns cfile-pdsw-2 cfile-tr
	cfile-lp-1-b cfile-lp-2-b cfile-pdsw-1 cfile-start )

ADD_CUSTOM_TARGET( test
    COMMAND pdde -c cfile-start -o out-start -b br-start )
	
ADD_CUSTOM_TARGET( test-tr
    COMMAND pdde -c cfile-tr -i out-start -o out-tr -b br-tr )

ADD_CUSTOM_TARGET( test-pdsw
    COMMAND pdde -c cfile-pdsw-1 -i out-start -o out-pdsw-1 -b br-pdsw-1 \;
	pdde -c cfile-pdsw-2 -i out-start -o out-pdsw-2 -b br-pdsw-2 )

ADD_CUSTOM_TARGET( test-lp-1-a
    COMMAND pdde -c cfile-lp-1-a -i out-pdsw-1 -o out-lp-1-a -b br-lp-1-a )

ADD_CUSTOM_TARGET( test-lp-1-b
    COMMAND pdde -c cfile-lp-1-b -i out-pdsw-1 -o out-lp-1-b -b br-lp-1-b )

ADD_CUSTOM_TARGET( test-lp-2-a
    COMMAND pdde -c cfile-lp-2-a -i out-pdsw-1 -o out-lp-2-a -b br-lp-2-a )

ADD_CUSTOM_TARGET( test-lp-2-b
    COMMAND pdde -c cfile-lp-2-b -i out-pdsw-1 -o out-lp-2-b -b br-lp-2-b )

ADD_CUSTOM_TARGET( test-ns
    COMMAND pdde -c cfile-ns -i out-start -o out-ns -b br-ns )

ADD_DEPENDENCIES( test ${NAME} )
ADD_DEPENDENCIES( test-tr ${NAME} )
ADD_DEPENDENCIES( test-pdsw ${NAME} )
ADD_DEPENDENCIES( test-lp-1-a ${NAME} )
ADD_DEPENDENCIES( test-lp-1-b ${NAME} )
ADD_DEPENDENCIES( test-lp-2-a ${NAME} )
ADD_DEPENDENCIES( test-lp-2-b ${NAME} )
ADD_DEPENDENCIES( test-ns ${NAME} )
