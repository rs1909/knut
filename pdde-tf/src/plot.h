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

#ifdef _MSC_VER
#define popen _popen
#define pclose _pclose
#endif

#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>

class GnuPlot
{
    FILE*               prg;
    std::ostringstream  pl;
    std::vector<std::string> cmd;
    std::vector<std::string> data;
  public:
    inline GnuPlot() : cmd(), data()
    {
      prg = popen("gnuplot", "w");
    }
    inline ~GnuPlot()
    {
      pclose(prg);
    }
    void SetTerminal(const char *term);
    void SetOutput(const char *term);
    void SetPointSize(double ps);
    void SetNoKey();
    void SetXLabel(const char *lb);
    void SetYLabel(const char *lb);
    void Plot(int n, const char *style = "");
    void AddData(int n, double x, double y);
    void AddSData(int n, double x, double y);
    void Show();
    void Save(const char * filename);
    void Clear();
};

#endif
