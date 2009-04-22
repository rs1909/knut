#! /bin/sh
git describe --abbrev=100 --tags | sed -e s/.*-g//g | awk -- '// { printf "%s", $0 }' >REVISION
git log -n1 HEAD | grep 'Date' | sed -e s/'Date:\ *'/', '/ | awk -- '// { printf "%s", $0 }' >>REVISION
