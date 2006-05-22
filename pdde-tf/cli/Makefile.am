bin_PROGRAMS = ../bin/pdde

___bin_pdde_SOURCES = main.cpp

INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/UMFPACK/AMD/Include \
	-I$(top_srcdir)/UMFPACK/UMFPACK/Include \
	-I$(top_srcdir)/ATLAS/include

___bin_pdde_LDFLAGS = -L$(top_srcdir)/src \
	-L$(top_srcdir)/UMFPACK/UMFPACK/Lib \
	-L$(top_srcdir)/UMFPACK/AMD/Lib \
	-L$(ATLAS_DIR) \
	-L$(top_srcdir)/libf2c \
	-L$(top_srcdir)/laarpack

___bin_pdde_LDADD = -lpdde -llaarpack -lumfpack -lamd -lcblas -latlas -lpddef2c @DL_LIB@
