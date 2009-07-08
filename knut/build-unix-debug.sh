#!/bin/sh
SCRIPTDIR=${0/build-unix-debug.sh/}
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=$HOME/knut -D CMAKE_BUILD_TYPE=DEBUG $SCRIPTDIR
make
make install
