#! /bin/sh
if ! `test -e ChangeLog`; then touch ChangeLog ; fi
svn log -r HEAD -q | sed -e s/' |.*.|'// -e s/' (.*)'// -e s/' '/'\\ '/g -e /'---'/d >REVISION
aclocal
automake -a -c
autoconf
