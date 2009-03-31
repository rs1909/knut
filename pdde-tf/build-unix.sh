#!/bin/sh
SCRIPTDIR=${0/build-unix.sh/}
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=$HOME/knut -D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
make install
