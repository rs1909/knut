#!/bin/sh
SCRIPTDIR=${0/build-mac-debug.sh/}
pushd $SCRIPTDIR
if test -d ../.git ; then 
	echo 'It is a git repository, updating REVISION' 
	./mkrevision.sh
	fi
popd  
cmake -G "Unix Makefiles"\
	-D CMAKE_CXX_COMPILER="g++-4.2"\
	-D CMAKE_C_COMPILER="gcc-4.2"\
	-D CMAKE_Fortran_COMPILER="/opt/local/bin/gfortran"\
	-D CMAKE_C_FLAGS="-m64 -g -O0 -fno-inline -W -Wall -Wno-unused -Wno-unknown-pragmas"\
	-D CMAKE_CXX_FLAGS="-m64 -g -O0 -fno-inline -W -Wall -Wno-unused"\
	-D CMAKE_Fortran_FLAGS="-m64 -g -O0"\
	-D CMAKE_INSTALL_PREFIX=$HOME/Applications/Knut\
	-D CMAKE_BUILD_TYPE=DEBUG $SCRIPTDIR
make
make install
