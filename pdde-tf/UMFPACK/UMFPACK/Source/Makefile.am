
UMF_C = umf_assemble.c umf_blas3_update.c umf_build_tuples.c umf_create_element.c \
	umf_dump.c umf_extend_front.c umf_garbage_collection.c umf_get_memory.c \
	umf_init_front.c umf_kernel.c umf_kernel_init.c umf_kernel_wrapup.c \
	umf_local_search.c umf_lsolve.c umf_ltsolve.c umf_mem_alloc_element.c \
	umf_mem_alloc_head_block.c umf_mem_alloc_tail_block.c \
	umf_mem_free_tail_block.c umf_mem_init_memoryspace.c \
	umf_report_vector.c umf_row_search.c umf_scale_column.c \
	umf_set_stats.c umf_solve.c umf_symbolic_usage.c umf_transpose.c \
	umf_tuple_lengths.c umf_usolve.c umf_utsolve.c umf_valid_numeric.c \
	umf_valid_symbolic.c umf_grow_front.c umf_start_front.c umf_2by2.c \
	umf_store_lu.c umf_scale.c

UMFINT_C = umf_analyze.c umf_apply_order.c umf_colamd.c umf_free.c umf_fsize.c \
	umf_is_permutation.c umf_malloc.c umf_realloc.c umf_report_perm.c \
	umf_singletons.c

UMFPACK_C = umfpack_col_to_triplet.c umfpack_defaults.c umfpack_free_numeric.c \
	umfpack_free_symbolic.c umfpack_get_numeric.c umfpack_get_lunz.c \
	umfpack_get_symbolic.c umfpack_get_determinant.c umfpack_numeric.c \
	umfpack_qsymbolic.c umfpack_report_control.c umfpack_report_info.c \
	umfpack_report_matrix.c umfpack_report_numeric.c umfpack_report_perm.c \
	umfpack_report_status.c umfpack_report_symbolic.c umfpack_report_triplet.c \
	umfpack_report_vector.c umfpack_solve.c umfpack_symbolic.c \
	umfpack_transpose.c umfpack_triplet_to_col.c umfpack_scale.c \
	umfpack_load_numeric.c umfpack_save_numeric.c \
	umfpack_load_symbolic.c umfpack_save_symbolic.c

GENERIC_C = umfpack_timer.c umfpack_tictoc.c umfpack_global.c

UMF_H = umf_assemble.c umf_blas3_update.h umf_build_tuples.h umf_create_element.h \
	umf_dump.h umf_extend_front.h umf_garbage_collection.h umf_get_memory.h \
	umf_init_front.h umf_kernel.h umf_kernel_init.h umf_kernel_wrapup.h \
	umf_local_search.h umf_lsolve.h umf_ltsolve.h umf_mem_alloc_element.h \
	umf_mem_alloc_head_block.h umf_mem_alloc_tail_block.h \
	umf_mem_free_tail_block.h umf_mem_init_memoryspace.h \
	umf_report_vector.h umf_row_search.h umf_scale_column.h \
	umf_set_stats.h umf_solve.h umf_symbolic_usage.h umf_transpose.h \
	umf_tuple_lengths.h umf_usolve.h umf_utsolve.h umf_valid_numeric.h \
	umf_valid_symbolic.h umf_grow_front.h umf_start_front.h umf_2by2.h \
	umf_store_lu.h umf_scale.h \
	umf_internal.h umf_config.h umf_version.h

UMFINT_H = umf_analyze.h umf_apply_order.h umf_colamd.h umf_free.h umf_fsize.h \
	umf_is_permutation.h umf_malloc.h umf_realloc.h umf_report_perm.h \
	umf_singletons.h

#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

INCLUDES = -I../../AMD/Include -I../../AMD/Source -I../Include -I../../../ATLAS/include

AM_CFLAGS = $(CFLAGS) -DCBLAS

#-------------------------------------------------------------------------------
# archive files for each version
#-------------------------------------------------------------------------------

LIB_I = lib_I.a
LIB_L = # lib_L.a
LIB_DI = lib_DI.a lib_DI_SH.a lib_DI_T_NX.a lib_DI_T_MN.a lib_DI_T_MX.a lib_DI_AF.a lib_DI_LUD.a lib_DI_WS.a
LIB_DL = # lib_DL.a lib_DL_SH.a lib_DL_T_NX.a lib_DL_T_MN.a lib_DL_T_MX.a lib_DL_AF.a lib_DL_LUD.a lib_DL_WS.a
LIB_ZI = # lib_ZI.a lib_ZI_SH.a lib_ZI_T_NX.a lib_ZI_T_MN.a lib_ZI_T_MX.a lib_ZI_AF.a lib_ZI_LUD.a lib_ZI_WS.a
LIB_ZL = # lib_ZL.a lib_ZL_SH.a lib_ZL_T_NX.a lib_ZL_T_MN.a lib_ZL_T_MX.a lib_ZL_AF.a lib_ZL_LUD.a lib_ZL_WS.a
LIB_GN = lib_GN.a

noinst_LIBRARIES = ../Lib/libumfpack.a $(LIB_I) $(LIB_L) $(LIB_DI) $(LIB_DL) $(LIB_ZI) $(LIB_ZL) $(LIB_GN)

___Lib_libumfpack_a_SOURCES = 

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

lib_I_a_SOURCES = $(UMFINT_C) $(UMFINT_H)
lib_I_a_CFLAGS = $(AM_CFLAGS) -DDINT

# lib_L_a_SOURCES = $(UMFINT_C) $(UMFINT_H)
# lib_L_a_CFLAGS = $(AM_CFLAGS) -DLONG

#-------------------------------------------------------------------------------
# compile each routine in the DI version
#-------------------------------------------------------------------------------

lib_DI_a_SOURCES = $(UMF_C) $(UMFPACK_C) $(UMF_H)
lib_DI_a_CFLAGS = $(AM_CFLAGS) -DDINT

lib_DI_SH_a_SOURCES = umf_ltsolve.c umf_utsolve.c umf_ltsolve.c umf_utsolve.c 
lib_DI_SH_a_CFLAGS = $(AM_CFLAGS) -DDINT -DCONJUGATE_SOLVE

lib_DI_T_NX_a_SOURCES = umf_triplet.c umf_triplet.h
lib_DI_T_NX_a_CFLAGS = $(AM_CFLAGS) -DDINT -DDO_VALUES

lib_DI_T_MN_a_SOURCES = umf_triplet.c umf_triplet.h
lib_DI_T_MN_a_CFLAGS = $(AM_CFLAGS) -DDINT -DDO_MAP

lib_DI_T_MX_a_SOURCES = umf_triplet.c umf_triplet.c
lib_DI_T_MX_a_CFLAGS = $(AM_CFLAGS) -DDINT -DDO_MAP -DDO_VALUES

lib_DI_AF_a_SOURCES = umf_assemble.c umf_assemble.h
lib_DI_AF_a_CFLAGS = $(AM_CFLAGS) -DDINT -DFIXQ

lib_DI_LUD_a_SOURCES = umf_store_lu.c umf_store_lu.h
lib_DI_LUD_a_CFLAGS = $(AM_CFLAGS) -DDINT -DDROP

lib_DI_WS_a_SOURCES = umfpack_solve.c
lib_DI_WS_a_CFLAGS = $(AM_CFLAGS) -DDINT -DWSOLVE

#-------------------------------------------------------------------------------
# compile each routine in the DL version
#-------------------------------------------------------------------------------

# lib_DL_a_SOURCES = $(UMF_C) $(UMFPACK_C) $(UMF_H)
# lib_DL_a_CFLAGS = $(AM_CFLAGS) -DDLONG
# 
# lib_DL_SH_a_SOURCES = umf_ltsolve.c umf_utsolve.c umf_ltsolve.c umf_utsolve.c 
# lib_DL_SH_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DCONJUGATE_SOLVE
# 
# lib_DL_T_NX_a_SOURCES = umf_triplet.c umf_triplet.h
# lib_DL_T_NX_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DDO_VALUES
# 
# lib_DL_T_MN_a_SOURCES = umf_triplet.c umf_triplet.h
# lib_DL_T_MN_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DDO_MAP
# 
# lib_DL_T_MX_a_SOURCES = umf_triplet.c umf_triplet.c
# lib_DL_T_MX_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DDO_MAP -DDO_VALUES
# 
# lib_DL_AF_a_SOURCES = umf_assemble.c umf_assemble.h
# lib_DL_AF_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DFIXQ
# 
# lib_DL_LUD_a_SOURCES = umf_store_lu.c umf_store_lu.h
# lib_DL_LUD_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DDROP
# 
# lib_DL_WS_a_SOURCES = umfpack_solve.c
# lib_DL_WS_a_CFLAGS = $(AM_CFLAGS) -DDLONG -DWSOLVE

#-------------------------------------------------------------------------------
# compile each routine in the ZI version
#-------------------------------------------------------------------------------

# lib_ZI_a_SOURCES = $(UMF_C) $(UMFPACK_C) $(UMF_H)
# lib_ZI_a_CFLAGS = $(AM_CFLAGS) -DZINT
# 
# lib_ZI_SH_a_SOURCES = umf_ltsolve.c umf_utsolve.c umf_ltsolve.c umf_utsolve.c 
# lib_ZI_SH_a_CFLAGS = $(AM_CFLAGS) -DZINT -DCONJUGATE_SOLVE
# 
# lib_ZI_T_NX_a_SOURCES = umf_triplet.c umf_triplet.h
# lib_ZI_T_NX_a_CFLAGS = $(AM_CFLAGS) -DZINT -DDO_VALUES
# 
# lib_ZI_T_MN_a_SOURCES = umf_triplet.c umf_triplet.h
# lib_ZI_T_MN_a_CFLAGS = $(AM_CFLAGS) -DZINT -DDO_MAP
# 
# lib_ZI_T_MX_a_SOURCES = umf_triplet.c umf_triplet.c
# lib_ZI_T_MX_a_CFLAGS = $(AM_CFLAGS) -DZINT -DDO_MAP -DDO_VALUES
# 
# lib_ZI_AF_a_SOURCES = umf_assemble.c umf_assemble.h
# lib_ZI_AF_a_CFLAGS = $(AM_CFLAGS) -DZINT -DFIXQ
# 
# lib_ZI_LUD_a_SOURCES = umf_store_lu.c umf_store_lu.h
# lib_ZI_LUD_a_CFLAGS = $(AM_CFLAGS) -DZINT -DDROP
# 
# lib_ZI_WS_a_SOURCES = umfpack_solve.c
# lib_ZI_WS_a_CFLAGS = $(AM_CFLAGS) -DZINT -DWSOLVE

#-------------------------------------------------------------------------------
# compile each routine in the DL version
#-------------------------------------------------------------------------------

# lib_ZL_a_SOURCES = $(UMF_C) $(UMFPACK_C) $(UMF_H)
# lib_ZL_a_CFLAGS = $(AM_CFLAGS) -DZLONG
# 
# lib_ZL_SH_a_SOURCES = umf_ltsolve.c umf_utsolve.c umf_ltsolve.c umf_utsolve.c 
# lib_ZL_SH_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DCONJUGATE_SOLVE
# 
# lib_ZL_T_NX_a_SOURCES = umf_triplet.c umf_triplet.h
# lib_ZL_T_NX_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DDO_VALUES
# 
# lib_ZL_T_MN_a_SOURCES = umf_triplet.c umf_triplet.h
# lib_ZL_T_MN_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DDO_MAP
# 
# lib_ZL_T_MX_a_SOURCES = umf_triplet.c umf_triplet.c
# lib_ZL_T_MX_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DDO_MAP -DDO_VALUES
# 
# lib_ZL_AF_a_SOURCES = umf_assemble.c umf_assemble.h
# lib_ZL_AF_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DFIXQ
# 
# lib_ZL_LUD_a_SOURCES = umf_store_lu.c umf_store_lu.h
# lib_ZL_LUD_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DDROP
# 
# lib_ZL_WS_a_SOURCES = umfpack_solve.c
# lib_ZL_WS_a_CFLAGS = $(AM_CFLAGS) -DZLONG -DWSOLVE

#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

lib_GN_a_SOURCES = $(GENERIC_C)

#-------------------------------------------------------------------------------
# Create the libumfpack.a library, but don't create the intermediate libraries
#-------------------------------------------------------------------------------

ALL_DI_OBJECTS = $(lib_DI_a_OBJECTS) $(lib_DI_SH_a_OBJECTS) \
	$(lib_DI_T_NX_a_OBJECTS) $(lib_DI_T_MN_a_OBJECTS) $(lib_DI_T_MX_a_OBJECTS) \
	$(lib_DI_AF_a_OBJECTS) $(lib_DI_LUD_a_OBJECTS) $(lib_DI_WS_a_OBJECTS)

ALL_DL_OBJECTS = $(lib_DL_a_OBJECTS) $(lib_DL_SH_a_OBJECTS) \
	$(lib_DL_T_NX_a_OBJECTS) $(lib_DL_T_MN_a_OBJECTS) $(lib_DL_T_MX_a_OBJECTS) \
	$(lib_DL_AF_a_OBJECTS) $(lib_DL_LUD_a_OBJECTS) $(lib_DL_WS_a_OBJECTS)

ALL_ZI_OBJECTS = $(lib_ZI_a_OBJECTS) $(lib_ZI_SH_a_OBJECTS) \
	$(lib_ZI_T_NX_a_OBJECTS) $(lib_ZI_T_MN_a_OBJECTS) $(lib_ZI_T_MX_a_OBJECTS) \
	$(lib_ZI_AF_a_OBJECTS) $(lib_ZI_LUD_a_OBJECTS) $(lib_ZI_WS_a_OBJECTS)

ALL_ZL_OBJECTS = $(lib_ZL_a_OBJECTS) $(lib_ZL_SH_a_OBJECTS) \
	$(lib_ZL_T_NX_a_OBJECTS) $(lib_ZL_T_MN_a_OBJECTS) $(lib_ZL_T_MX_a_OBJECTS) \
	$(lib_ZL_AF_a_OBJECTS) $(lib_ZL_LUD_a_OBJECTS) $(lib_ZL_WS_a_OBJECTS)

ALL_OBJECTS = $(lib_I_a_OBJECTS) $(lib_L_a_OBJECTS) $(lib_GN_a_OBJECTS) \
	$(ALL_DI_OBJECTS) $(ALL_DL_OBJECTS) $(ALL_ZI_OBJECTS) $(ALL_ZL_OBJECTS)

../Lib/libumfpack.a: $(ALL_OBJECTS)
	$(AR) $(ARFLAGS) ../Lib/libumfpack.a $(ALL_OBJECTS)
