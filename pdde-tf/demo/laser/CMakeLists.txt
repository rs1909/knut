SET( NAME laser )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET( CFILEDIR ${PDDE_SOURCE_DIR}/demo/${NAME} )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
    cfile-start     cfile-hopfsw     cfile-psol     cfile-ns     cfile-pd     cfile-lp
    cfile-start.xml cfile-hopfsw.xml cfile-psol.xml cfile-ns.xml cfile-pd.xml cfile-lp.xml
    DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test-1
    COMMAND pdde -c ${CFILEDIR}/cfile-start -o out-start.mat -b out-start.br )
ADD_CUSTOM_TARGET( ${NAME}-test-2
    COMMAND pdde -c ${CFILEDIR}/cfile-hopfsw -i out-start.mat -o out-hopfsw.mat -b out-hopfsw.br)
ADD_CUSTOM_TARGET( ${NAME}-test-3
    COMMAND pdde -c ${CFILEDIR}/cfile-psol -i out-hopfsw.mat -o out-psol.mat -b out-psol.br )
ADD_CUSTOM_TARGET( ${NAME}-test-4
    COMMAND pdde -c ${CFILEDIR}/cfile-ns -i out-psol.mat -o out-ns.mat -b out-ns.br )
ADD_CUSTOM_TARGET( ${NAME}-test-5
    COMMAND pdde -c ${CFILEDIR}/cfile-pd -i out-psol.mat -o out-pd.mat -b out-pd.br )
ADD_CUSTOM_TARGET( ${NAME}-test-6
    COMMAND pdde -c ${CFILEDIR}/cfile-lp -i out-psol.mat -o out-lp.mat -b out-lp.br )

ADD_DEPENDENCIES( ${NAME}-test-1 ${NAME} )
