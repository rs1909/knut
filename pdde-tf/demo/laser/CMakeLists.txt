SET( NAME laser )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
    cfile-start cfile-ns cfile-pd cfile-lp in-start from_biftool.mat
    cfile-start.xml cfile-ns.xml cfile-pd.xml cfile-lp.xml DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test
    COMMAND pdde -c cfile-start -i in-start -o out-start -b br-start )
	
ADD_CUSTOM_TARGET( ${NAME}-test-ns
    COMMAND pdde -c cfile-ns -i out-start -o out-ns -b branch-ns )

ADD_CUSTOM_TARGET( ${NAME}-test-pd
    COMMAND pdde -c cfile-pd -i out-start -o out-pd -b branch-pd )

ADD_CUSTOM_TARGET( ${NAME}-test-lp
    COMMAND pdde -c cfile-lp -i out-start -o out-lp -b branch-lp )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-ns ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-pd ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-lp ${NAME} )
