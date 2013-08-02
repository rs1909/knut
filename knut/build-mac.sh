#!/bin/sh
# CMAKE_OSX_ARCHITECTURES is necessary to remove other architectures
SCRIPTDIR=${0/build-mac.sh/}
pushd $SCRIPTDIR
if test -d ../.git ; then 
	echo 'It is a git repository, updating REVISION' 
	./mkrevision.sh
	fi
popd
cmake -G "Unix Makefiles"\
	-D CMAKE_CXX_COMPILER="clang++"\
	-D CMAKE_C_COMPILER="clang"\
	-D CMAKE_Fortran_COMPILER="gfortran"\
	-D CMAKE_C_FLAGS="-m64 -g -O2 -W -Wall -Wno-unused-parameter -Wno-unknown-pragmas"\
	-D CMAKE_CXX_FLAGS="-m64 -g -O2 -std=c++11 -stdlib=libc++ -W -Wall -Wno-unused-parameter"\
	-D CMAKE_Fortran_FLAGS="-m64 -g -O3"\
	-D CMAKE_INSTALL_PREFIX=$HOME/Applications/Knut\
	-D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
