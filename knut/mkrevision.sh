#! /bin/sh
git log | grep '^commit' | wc -l | sed -e s/'^\ *'// | awk -- '// { printf "%s", $0 }' >REVISION
git show HEAD | head -n3 | grep 'Date' | sed -e s/'Date:\ *'/', '/ | awk -- '// { printf "%s", $0 }' >>REVISION
git show HEAD | head -n3 | grep 'commit' | sed -e s/'commit\ *'/', '/ | awk -- '// { printf "%s", $0 }' >>REVISION