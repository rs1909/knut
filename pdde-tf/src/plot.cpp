// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "plot.h"
#include <iostream>

void GnuPlot::SetTerminal( const char *term )
{
  pl<<"set terminal "<<term<<'\n';
}

void GnuPlot::SetOutput( const char *term )
{
  pl<<"set output "<<term<<'\n';
}

void GnuPlot::SetPointSize( double ps )
{
  pl<<"set pointsize "<<ps<<'\n';
}
void GnuPlot::SetNoKey()
{
  pl<<"set nokey\n";
}
void GnuPlot::SetXLabel( const char *lb )
{
  pl<<"set xlabel \""<<lb<<"\"\n";
}
void GnuPlot::SetYLabel( const char *lb )
{
  pl<<"set ylabel \""<<lb<<"\"\n";
}

void GnuPlot::Plot( int n, const char *style )
{
  if( n >= (int)cmd.size() ) { 
    for( int i = cmd.size(); i <= n; i++ ) {
      cmd.insert(cmd.end(), std::string() );
      data.insert(data.end(), std::string() );
    }
  }
  cmd[n] = style;
}

void GnuPlot::AddData( int n, double x, double y )
{
  std::ostringstream tmp;
  if( n >= (int)cmd.size() ){ std::cout<<"Gnuplot::AddData: Bad index\n"; exit(-1); }
  
  tmp<<x<<'\t'<<y<<'\n'; 
  data[n] += tmp.str();
}

void GnuPlot::AddSData( int n, double x, double y )
{
  std::ostringstream tmp;
  if( n >= (int)cmd.size() ){ std::cout<<"Gnuplot::AddData: Bad index\n"; exit(-1); }
  
  tmp<<"\n\n"<<x<<'\t'<<y<<'\n'; 
  data[n] += tmp.str();
}
  
void GnuPlot::Show()
{
  if( cmd.size() != 0 ) pl<<"plot ";
  for( int i = 0; i < (int)cmd.size(); i++) {
    if( i == 0 ) pl<<" '-' "<<cmd[i];
    else pl<<", '-' "<<cmd[i];
  }
  pl<<'\n';
  for( int i = 0; i < (int)cmd.size(); i++) {
    pl<<data[i]<<"e\n";
  }
  pl.flush();
  fputs( (pl.str()).c_str(), prg );
  fflush( prg );
}
void GnuPlot::Save( const char * filename )
{
  pl<<"set terminal postscript eps enh 22\n"<<"set output \""<<filename<<"\"\n";
  if( cmd.size() != 0 ) pl<<"plot ";
  for( int i = 0; i < (int)cmd.size(); i++) {
    if( i == 0 ) pl<<" '-' "<<cmd[i];
    else pl<<", '-' "<<cmd[i];
  }
  pl<<'\n';
  for( int i = 0; i < (int)cmd.size(); i++) {
    pl<<data[i]<<"e\n";
  }
  pl.flush();
  fputs( (pl.str()).c_str(), prg );
  fflush( prg );
}
void GnuPlot::Clear()
{
  fputs( "clear\nreset\n", prg );
  fflush( prg );
  
  cmd.erase(cmd.begin(), cmd.end());
  data.erase(data.begin(), data.end());
  pl.str("");
}

