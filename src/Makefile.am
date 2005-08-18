bin_PROGRAMS = ../bin/pdde
___bin_pdde_SOURCES =	main.cpp charmat.cpp plot.cpp hypermatrix.cpp ncolloc.cpp torcolloc.cpp point.cpp torpoint.cpp spmatrix.cpp matrix.cpp \
			charmat.h hypermatrix.h matrix.h ncolloc.h torcolloc.h pddesys.h plot.h point.h torpoint.h spmatrix.h \
			system.cpp system.h pointtype.h
incdir = $(includedir)
inc_HEADERS = pddesys.h matrix.h
INCLUDES = -I$(top_srcdir)/UMFPACK/UMFPACK/Include
AM_LDFLAGS = -L$(top_srcdir)/UMFPACK/UMFPACK/Lib \
	-L$(top_srcdir)/UMFPACK/AMD/Lib \
	-L$(top_srcdir)/ATLAS/lib \
	-L$(top_srcdir)/libf2c \
	-L$(top_srcdir)/laarpack

LDADD = -llaarpack -lumfpack -lamd -lcblas -latlas -lpddef2c -ldl
