#!/bin/sh
astyle -s2 -C -S -b --min-conditional-indent=0 -p -U -V -M80 --mode=c $1
