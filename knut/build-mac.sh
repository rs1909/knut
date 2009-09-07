#!/bin/sh
SCRIPTDIR=${0/build-mac.sh/}
cmake -G "Unix Makefiles"\
	-D CMAKE_CXX_COMPILER="g++"\
	-D CMAKE_C_COMPILER="gcc"\
	-D CMAKE_Fortran_COMPILER="gfortran"\
	-D CMAKE_C_FLAGS="-g -O2 -W -Wall -Wno-unused -Wno-unknown-pragmas"\
	-D CMAKE_CXX_FLAGS="-g -O2 -W -Wall -Wno-unused -Wconversion"\
	-D CMAKE_Fortran_FLAGS="-g -O3"\
	-D CMAKE_INSTALL_PREFIX=$HOME/Apps/Knut\
	-D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
make install
