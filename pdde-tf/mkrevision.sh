#! /bin/sh
git log -n1 HEAD | grep commit | sed -e s/'commit '// | awk -- '// { printf "%s", $0 }' >REVISION
git log -n1 HEAD | grep 'Date' | sed -e s/'Date: '/','/ | awk -- '// { printf "%s", $0 }' >>REVISION
