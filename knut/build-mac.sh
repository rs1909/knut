#!/bin/sh
SCRIPTDIR=${0/build-mac.sh/}
cmake -G "Unix Makefiles"\
	-D CMAKE_CXX_COMPILER="g++-4.2"\
	-D CMAKE_C_COMPILER="gcc-4.2"\
	-D CMAKE_Fortran_COMPILER="gfortran-mp-4.3"\
	-D CMAKE_C_FLAGS="-m64 -g -O2 -W -Wall -Wno-unused -Wno-unknown-pragmas"\
	-D CMAKE_CXX_FLAGS="-m64 -g -O2 -W -Wall -Wno-unused -Wconversion"\
	-D CMAKE_Fortran_FLAGS="-m64 -g -O3"\
	-D CMAKE_INSTALL_PREFIX=$HOME/Apps/Knut\
	-D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
make install
