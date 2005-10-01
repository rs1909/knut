LIB_I = lib_I.a
LIB_L = # lib_L.a

noinst_LIBRARIES = $(LIB_I) $(LIB_L) ../Lib/libamd.a

___Lib_libamd_a_SOURCES =

lib_I_a_SOURCES = amd_aat.c amd_1.c amd_2.c amd_dump.c \
	amd_postorder.c amd_post_tree.c amd_defaults.c \
	amd_order.c amd_control.c amd_info.c amd_valid.c amd_preprocess.c \
	amd_internal.h

# lib_L_a_SOURCES = amd_aat.c amd_1.c amd_2.c amd_dump.c \
# 	amd_postorder.c amd_post_tree.c amd_defaults.c \
# 	amd_order.c amd_control.c amd_info.c amd_valid.c amd_preprocess.c \
# 	amd_internal.h

lib_I_a_CFLAGS = $(AM_CFLAGS) -DDINT
# lib_L_a_CFLAGS = $(AM_CFLAGS) -DDLONG

INCLUDES = -I../Include

../Lib/libamd.a: $(lib_I_a_OBJECTS) $(lib_L_a_OBJECTS)
	$(AR) $(ARFLAGS) ../Lib/libamd.a $(lib_I_a_OBJECTS) $(lib_L_a_OBJECTS)
