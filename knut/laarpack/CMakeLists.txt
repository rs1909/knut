ADD_LIBRARY( laarpack laarpack.f dlamch.f debug.h stat.h )
SET_SOURCE_FILES_PROPERTIES( dlamch.f PROPERTIES COMPILE_FLAGS -ffloat-store )
ADD_LIBRARY( AUX_BLAS xerbla.f )
