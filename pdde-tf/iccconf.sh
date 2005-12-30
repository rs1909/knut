#!/bin/sh
./configure --prefix=$HOME/pdde-cont/ 'CC=icc' 'CFLAGS=-gcc-version=340 -g -O2 -xB' 'CXX=icpc' 'CXXFLAGS=-gcc-version=340 -cxxlib-icc -g -O2 -xB'