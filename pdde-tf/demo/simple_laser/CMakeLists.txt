SET( NAME simple_laser )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/src )
ADD_LIBRARY( ${NAME} MODULE sys-${NAME}.cpp )
SET_TARGET_PROPERTIES( ${NAME} PROPERTIES PREFIX sys- SUFFIX .so )
INSTALL( FILES sys-${NAME}.cpp
    cfile-start cfile-hbsw cfile-hbct
    DESTINATION ${PDDE_DEMO_DIR}/${NAME} )

ADD_CUSTOM_TARGET( ${NAME}-test
    COMMAND pushd ${CMAKE_INSTALL_PREFIX}/${PDDE_DEMO_DIR}/${NAME}\;
    pcompile sys-${NAME}.cpp\;
    pdde -c cfile-start -o out-start.mat -b out-start.br\;
    pdde -c cfile-hbsw -i out-start.mat -o out-hbsw.mat -b out-hbsw.br\;
    pdde -c cfile-hbct -i out-start.mat -o out-hbct.mat -b out-hbct.br\; 
    popd\;)
	
ADD_DEPENDENCIES( ${NAME}-test ${NAME} )
