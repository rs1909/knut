#!/bin/sh
SCRIPTDIR=${0/build-mac.sh/}
cmake -G "Unix Makefiles"\
	-D CMAKE_CXX_COMPILER="/usr/bin/g++-4.2"\
	-D CMAKE_C_COMPILER="/usr/bin/gcc-4.2"\
	-D CMAKE_OSX_ARCHITECTURES="x86_64"\
	-D CMAKE_INSTALL_PREFIX=$HOME/Apps/Knut\
	-D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
make install
