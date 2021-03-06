#!/bin/sh
SCRIPTDIR=${0/build-win-cross.sh/}
pushd $SCRIPTDIR
if test -d ../.git ; then 
	echo 'It is a git repository, updating REVISION' 
	./mkrevision.sh
	fi
popd
rm -f CMakeCache.txt
cmake -DCMAKE_TOOLCHAIN_FILE=$SCRIPTDIR/wintoolchain \
	-DPKG_CONFIG_EXECUTABLE=`which mingw32-pkg-config` \
	-DCMAKE_BUILD_TYPE=RELEASE \
	-DCMAKE_C_FLAGS="-g -O2 -W -Wall -Wno-unused -Wno-unknown-pragmas" \
	-DCMAKE_CXX_FLAGS="-g -O2 -std=c++11 -W -Wall -Wextra -Wno-unused" \
	-DCMAKE_Fortran_FLAGS="-g -O3" \
	-DCMAKE_INSTALL_PREFIX=$HOME/WinKnut/ $SCRIPTDIR
#make
#make install
#makensis wininstaller.nsi
