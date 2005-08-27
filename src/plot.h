// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef PLOT_H
#define PLOT_H

#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>

using namespace std;

class GnuPlot{
  FILE*          prg;
  ostringstream  pl;
  vector<string> cmd;
  vector<string> data;
public:
  inline GnuPlot() : cmd(), data() { prg = popen("gnuplot","w"); }
  inline ~GnuPlot() { pclose(prg); }
  void SetTerminal( char *term );
  void SetOutput( char *term );
  void SetPointSize( double ps );
  void SetNoKey();
  void SetXLabel( char *lb );
  void SetYLabel( char *lb );
  void Plot( int n, char *style = "" );
  void AddData( int n, double x, double y );
  void AddSData( int n, double x, double y );
  void Show();
  void Save( char * filename );
  void Clear();
};

#endif
