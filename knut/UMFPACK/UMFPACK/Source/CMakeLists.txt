SET( UMFINT_C umf_analyze.c umf_apply_order.c umf_colamd.c umf_free.c umf_fsize.c
	umf_is_permutation.c umf_malloc.c umf_realloc.c umf_report_perm.c
	umf_singletons.c )

SET( UMFINT_H umf_analyze.h umf_apply_order.h umf_colamd.h umf_free.h umf_fsize.h
	umf_is_permutation.h umf_malloc.h umf_realloc.h umf_report_perm.h
	umf_singletons.h )

SET( UMF_C umf_assemble.c umf_blas3_update.c umf_build_tuples.c umf_create_element.c
	umf_dump.c umf_extend_front.c umf_garbage_collection.c umf_get_memory.c
	umf_init_front.c umf_kernel.c umf_kernel_init.c umf_kernel_wrapup.c
	umf_local_search.c umf_lsolve.c umf_ltsolve.c umf_mem_alloc_element.c
	umf_mem_alloc_head_block.c umf_mem_alloc_tail_block.c
	umf_mem_free_tail_block.c umf_mem_init_memoryspace.c
	umf_report_vector.c umf_row_search.c umf_scale_column.c
	umf_set_stats.c umf_solve.c umf_symbolic_usage.c umf_transpose.c
	umf_tuple_lengths.c umf_usolve.c umf_utsolve.c umf_valid_numeric.c
	umf_valid_symbolic.c umf_grow_front.c umf_start_front.c umf_2by2.c
	umf_store_lu.c umf_scale.c )

SET( UMF_H umf_assemble.h umf_blas3_update.h umf_build_tuples.h umf_create_element.h
	umf_dump.h umf_extend_front.h umf_garbage_collection.h umf_get_memory.h
	umf_init_front.h umf_kernel.h umf_kernel_init.h umf_kernel_wrapup.h
	umf_local_search.h umf_lsolve.h umf_ltsolve.h umf_mem_alloc_element.h
	umf_mem_alloc_head_block.h umf_mem_alloc_tail_block.h
	umf_mem_free_tail_block.h umf_mem_init_memoryspace.h
	umf_report_vector.h umf_row_search.h umf_scale_column.h
	umf_set_stats.h umf_solve.h umf_symbolic_usage.h umf_transpose.h
	umf_tuple_lengths.h umf_usolve.h umf_utsolve.h umf_valid_numeric.h
	umf_valid_symbolic.h umf_grow_front.h umf_start_front.h umf_2by2.h
	umf_store_lu.h umf_scale.h
	umf_internal.h umf_config.h umf_version.h )

SET( UMFPACK_C umfpack_col_to_triplet.c umfpack_defaults.c umfpack_free_numeric.c
	umfpack_free_symbolic.c umfpack_get_numeric.c umfpack_get_lunz.c
	umfpack_get_symbolic.c umfpack_get_determinant.c umfpack_numeric.c
	umfpack_qsymbolic.c umfpack_report_control.c umfpack_report_info.c
	umfpack_report_matrix.c umfpack_report_numeric.c umfpack_report_perm.c
	umfpack_report_status.c umfpack_report_symbolic.c umfpack_report_triplet.c
	umfpack_report_vector.c umfpack_solve.c umfpack_symbolic.c
	umfpack_transpose.c umfpack_triplet_to_col.c umfpack_scale.c
#	umfpack_load_numeric.c umfpack_save_numeric.c
#	umfpack_load_symbolic.c umfpack_save_symbolic.c
)

SET( GENERIC_C umfpack_timer.c umfpack_tictoc.c umfpack_global.c cholmod_blas.h )

SET( UMF_GENERATED_C umf_ltsolve_conj.c umf_utsolve_conj.c
	umf_triplet_map.c umf_triplet_values.c umf_triplet_map_values.c
	umf_assemble_fixq.c
	umf_store_lu_drop.c
	umfpack_wsolve.c )

#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

INCLUDE_DIRECTORIES( 
	${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include
	${KNUT_SOURCE_DIR}/UMFPACK/AMD/Source
	${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include
	${KNUT_SOURCE_DIR}/UMFPACK/UFconfig 
	${KNUT_SOURCE_DIR}/ATLAS/include )

SET( FLT DINT ) # form DINT DLONG ZINT ZLONG 

ADD_DEFINITIONS( -DNO_TIMER )

ADD_LIBRARY( umfpack
	${UMFINT_C} ${UMFINT_H}
	${UMF_C} ${UMF_H} ${UMFPACK_C}
	${GENERIC_C} ${UMF_GENERATED_C} )

SET_SOURCE_FILES_PROPERTIES( ${UMFINT_C} ${UMFINT_H}
	${UMF_C} ${UMF_H} ${UMFPACK_C}
	${GENERIC_C} ${UMF_GENERATED_C}
	PROPERTIES COMPILE_FLAGS "-D${FLT}" )
