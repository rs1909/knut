#! /bin/sh
if ! `test -e ChangeLog`; then touch ChangeLog ; fi
./mkrevision.sh
touch configure.ac
aclocal
automake -a -c
autoconf
autoheader