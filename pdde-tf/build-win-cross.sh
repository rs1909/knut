#!/bin/sh
SCRIPTDIR=${0/build-win-cross.sh/}
rm -f CMakeCache.txt
cmake -DCMAKE_TOOLCHAIN_FILE=$SCRIPTDIR/wintoolchain -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=$HOME/WinKnut/ $SCRIPTDIR
make
make install
