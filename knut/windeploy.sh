#!/bin/bash
# first argument is the binary name
# second argument is the directory of deployment
WINEDEBUG=+loaddll wine $1 2> dll.log

grep Loaded dll.log | grep -v 'system32\|:load_builtin_dll' \
  | awk -F'"' '{print $2}' \
  | sed -e 's@\\\\@/@g' -e 's/^[A-Z]://' \
  | sort > dll.lst

mkdir -p "$2"/{imageformats,platforms}
for i in imageformats platforms ; do
	grep "/plugins/$i" dll.lst | xargs -r cp -t "$2"/$i
done
grep -v '/plugins/' dll.lst | xargs -r cp -t "$2"
