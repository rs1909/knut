#! /bin/sh
if ! `test -e ChangeLog`; then touch ChangeLog ; fi
aclocal
automake -a -c
autoconf
