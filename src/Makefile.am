bin_PROGRAMS = ../bin/pdde ../bin/cfile09to10

___bin_pdde_SOURCES =	main.cpp charmat.cpp plot.cpp ncolloc.cpp torcolloc.cpp point.cpp torpoint.cpp \
			matrix.cpp cspblas.c cspblas.h spmatrix.cpp hypermatrix.cpp \
			charmat.h hypermatrix.h matrix.h ncolloc.h torcolloc.h pddesys.h plot.h point.h torpoint.h spmatrix.h \
			system.cpp system.h pointtype.h error.h error.cpp
			
___bin_cfile09to10_SOURCES = cfile09to10.cpp pointtype.h error.h error.cpp

incdir = $(includedir)
inc_HEADERS = pddesys.h matrix.h
INCLUDES = -I$(top_srcdir)/UMFPACK/UMFPACK/Include -I$(top_srcdir)/ATLAS/include
AM_LDFLAGS = -L$(top_srcdir)/UMFPACK/UMFPACK/Lib \
	-L$(top_srcdir)/UMFPACK/AMD/Lib \
	-L$(top_srcdir)/ATLAS/lib \
	-L$(top_srcdir)/libf2c \
	-L$(top_srcdir)/laarpack -Wl,--as-needed

LDADD = -llaarpack -lumfpack -lamd -lcblas -latlas -lpddef2c -ldl 
