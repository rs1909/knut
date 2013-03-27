#!/bin/sh
SCRIPTDIR=${0/build-mac-debug.sh/}
pushd $SCRIPTDIR
if test -d ../.git ; then 
	echo 'It is a git repository, updating REVISION' 
	./mkrevision.sh
	fi
popd  
cmake -G "Unix Makefiles"\
	-D CMAKE_CXX_COMPILER="g++"\
	-D CMAKE_C_COMPILER="gcc"\
	-D CMAKE_Fortran_COMPILER="gfortran"\
	-D CMAKE_C_FLAGS_DEBUG="-m64 -g -O0 -W -Wall -Wno-unused -Wno-unknown-pragmas -fno-inline"\
	-D CMAKE_CXX_FLAGS_DEBUG="-m64 -g -O0 -std=c++11 -W -Wall -Wno-unused -Wconversion -fno-inline"\
	-D CMAKE_Fortran_FLAGS_DEBUG="-m64 -g -O0"\
	-D CMAKE_INSTALL_PREFIX=$HOME/Applications/Knut\
	-D CMAKE_BUILD_TYPE=DEBUG $SCRIPTDIR
make
make install
