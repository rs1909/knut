SET( ALL_SOURCES amd_aat.c amd_1.c amd_2.c amd_dump.c
	amd_postorder.c amd_post_tree.c amd_defaults.c
	amd_order.c amd_control.c amd_info.c amd_valid.c amd_preprocess.c
	../Include/amd_internal.h amd_global.c )

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include
                     ${KNUT_SOURCE_DIR}/UMFPACK/AMD/Source
                     ${KNUT_SOURCE_DIR}/UMFPACK/UFconfig
                     ${KNUT_BINARY_DIR}/src/ )

SET( FLT DINT ) # from DINT DLONG

ADD_LIBRARY( amd STATIC ${ALL_SOURCES} )
SET_SOURCE_FILES_PROPERTIES( ${ALL_SOURCES} PROPERTIES DEFINE_SYMBOL ${FLT} )
