#!/bin/sh
SCRIPTDIR=${0/build-win-cross.sh/}
rm -f CMakeCache.txt
cmake -DCMAKE_TOOLCHAIN_FILE=$SCRIPTDIR/wintoolchain \
	-DPKG_CONFIG_EXECUTABLE=`which mingw32-pkg-config` \
	-DCMAKE_BUILD_TYPE=RELEASE \
	-DCMAKE_INSTALL_PREFIX=$HOME/WinKnut/ $SCRIPTDIR
make
make install
