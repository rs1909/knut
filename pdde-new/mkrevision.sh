#! /bin/sh
svn log -r HEAD -q | sed -e s/' |.*.|'// -e s/' (.*)'// -e /'---'/d | awk -- '/r/ { printf "%s", $0 }' >REVISION
