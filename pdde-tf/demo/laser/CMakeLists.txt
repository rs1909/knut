SET( NAME laser )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
    cfile-start cfile-hopfsw cfile-psol cfile-ns cfile-pd cfile-lp
    cfile-start.xml cfile-hopfsw.xml cfile-psol.xml cfile-ns.xml cfile-pd.xml cfile-lp.xml
    DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test
    COMMAND pdde -c cfile-start -o out-start.mat -b br-start\;
    pdde -c cfile-hopfsw -i out-start.mat -o out-hopfsw.mat -b br-hopfsw\;
    pdde -c cfile-psol -i out-hopfsw.mat -o out-psol.mat -b br-psol\; )
	
ADD_CUSTOM_TARGET( ${NAME}-test-ns
    COMMAND pdde -c cfile-ns -i out-psol.mat -o out-ns.mat -b branch-ns )

ADD_CUSTOM_TARGET( ${NAME}-test-pd
    COMMAND pdde -c cfile-pd -i out-psol.mat -o out-pd.mat -b branch-pd )

ADD_CUSTOM_TARGET( ${NAME}-test-lp
    COMMAND pdde -c cfile-lp -i out-psol.mat -o out-lp.mat -b branch-lp )

ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-ns ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-pd ${NAME} )
ADD_DEPENDENCIES( ${NAME}-test-lp ${NAME} )
