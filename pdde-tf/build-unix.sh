#!/bin/sh
# SCRIPTDIR=${0/build-unix.sh/}
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=$HOME/pdde-cont -D CMAKE_BUILD_TYPE=RELEASE `pwd`
make
make install
