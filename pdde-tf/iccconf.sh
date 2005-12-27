#!/bin/sh
./configure --prefix=$HOME/pdde-cont/ "CC=icc" "CFLAGS=-xB -O2 -g" "CXX=icpc" "CXXFLAGS=-cxxlib-icc -xB -O2 -g"