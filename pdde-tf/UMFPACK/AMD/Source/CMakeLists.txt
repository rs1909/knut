SET( ALL_SOURCES amd_aat.c amd_1.c amd_2.c amd_dump.c
	amd_postorder.c amd_post_tree.c amd_defaults.c
	amd_order.c amd_control.c amd_info.c amd_valid.c amd_preprocess.c
	amd_internal.h amd_global.c )

INCLUDE_DIRECTORIES( ${PDDE_SOURCE_DIR}/UMFPACK/AMD/Include ${PDDE_SOURCE_DIR}/UMFPACK/AMD/Source )

SET( FLT DINT ) # from DINT DLONG

ADD_LIBRARY( amd ${ALL_SOURCES} )
SET_SOURCE_FILES_PROPERTIES( ${ALL_SOURCES} PROPERTIES DEFINE_SYMBOL ${FLT} )
