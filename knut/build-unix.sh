#!/bin/bash
SCRIPTDIR=${0/build-unix.sh/}
pushd $SCRIPTDIR
if test -d ../.git ; then 
	echo 'It is a git repository, updating REVISION' 
	./mkrevision.sh
	fi
popd
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=$HOME/knut -D CMAKE_BUILD_TYPE=RELWITHDEBINFO $SCRIPTDIR
make
