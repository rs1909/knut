#!/bin/sh
SCRIPTDIR=${0/configure-unix.sh/}
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=$HOME/pdde-cont -D CMAKE_BUILD_TYPE=RELEASE $SCRIPTDIR
make -C $SCRIPTDIR
make install -C $SCRIPTDIR
