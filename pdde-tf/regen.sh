#! /bin/sh
if ! `test -e ChangeLog`; then touch ChangeLog ; fi
svn log -r HEAD -q | sed -e s/' |.*.|'// -e s/' (.*)'// -e /'---'/d >REVISION
touch configure.ac
aclocal
automake -a -c
autoconf
autoheader