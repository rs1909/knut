noinst_LIBRARIES = libpdde.a

libpdde_a_SOURCES =	plot.cpp ncolloc.cpp torcolloc.cpp point.cpp torpoint.cpp \
			matrix.cpp cspblas.c cspblas.h spmatrix.cpp testfunct.cpp \
			hypermatrix.h matrix.h ncolloc.h torcolloc.h \
			pddesys.h plot.h point.h torpoint.h spmatrix.h testfunct.h \
			system.cpp system.h pointtype.h pderror.h pderror.cpp \
			parameters.cpp parameters.h
			
incdir = $(includedir)
inc_HEADERS = pddesys.h matrix.h
INCLUDES = -I$(top_srcdir)/UMFPACK/AMD/Include -I$(top_srcdir)/UMFPACK/UMFPACK/Include -I$(top_srcdir)/ATLAS/include
