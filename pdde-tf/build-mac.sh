#!/bin/sh
SCRIPTDIR=${0/build-mac.sh/}
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=$HOME/Apps/pdde-cont -D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
make install
