/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <cmath>
#include "knutsys.h"

// system definition
// par(0) = T the period length 
// par(1) = L/N
// par(2) = the delay
// par(3) = alpha
// par(4) = v_0
// ...
// x[0]    = v_0
// x[1]    = v_1
// ...
// x[N-1]  = v_N 
// x[N]    = h_0
// x[N+1]  = h_1
// ...
// x[2N-2] = h_{N-1}
// dim = 2N-1 

extern "C"
{

int sys_ndim() { return 5; }
int sys_npar() { return 5; }
int sys_ntau() { return 2; }
int sys_nderi() { return 0; }

void sys_p_tau( Array2D<double>& out, const Array1D<double>& time, const Array1D<double>& par )
{
  for (int idx = 0; idx < time.size(); ++idx)
  {
    out(0,idx) = 0.0;
    out(1,idx) = par(2);
  }
}

void sys_p_dtau( Array2D<double>& out, const Array1D<double>& time, const Array1D<double>& par, int vp )
{
  for (int idx = 0; idx < time.size(); ++idx)
  {
    out(0,idx) = 0.0;
    if( vp == 2 ) out(1,idx) = 1.0;
    else out(1,idx) = 0.0;
  }
}

static inline double V( double x )
{
  return pow( ( x - 1 ) / 1, 3.0 ) / ( pow( ( x - 1 ) / 1, 3.0 ) + 1);
}

void sys_p_rhs( Array2D<double>& out, const Array1D<double>& time, const Array3D<double>& yy, const Array1D<double>& par, int sel )
{
#define xx(i,j) yy(i-1,j-1,idx)
#define f(i,j) out(i-1,idx)
  for (int idx = 0; idx < time.size(); ++idx)
  {
    if ( xx(4,2) > 1 ){
      f(1,1) = par(3) * ( par(4) * V(xx(4,2)) - xx(1,1) ); 
    }else{
      f(1,1) = - par(3) * xx(1,1);
    }
  }

  for (int idx = 0; idx < time.size(); ++idx)
  {
    if( xx(5,2) > 1 ){
      f(2,1) = par(3) * ( par(4) * V(xx(5,2)) - xx(2,1) ); 
    }else{
      f(2,1) = - par(3) * xx(2,1);
    }
  }

  for (int idx = 0; idx < time.size(); ++idx)
  {
    if( ( 3*par(1) - xx(5,2) - xx(4,2) ) > 1 ){
      f(3,1) = par(3) * ( par(4) * V( 3*par(1) - xx(5,2) - xx(4,2) ) - xx(3,1) );
    }else{
      f(3,1) = - par(3) * xx(3,1);
    }
  }

  for (int idx = 0; idx < time.size(); ++idx)
  {
    f(4,1) = xx(2,1) - xx(1,1);
    f(5,1) = xx(3,1) - xx(2,1);
  }
}

void sys_p_deri( Array3D<double>&, const Array1D<double>& time, const Array3D<double>& yy, const Array1D<double>& par, int sel, 
                 int nx, const int* vx, int np, const int* vp, const Array3D<double>& vv )
{

}

void sys_stpar( Vector& par )
{
  par(0) = 1.1; // T the period length 
  par(1) = 2.0; // L/N
  par(2) = 1.0; // the delay
  par(3) = 1.0; // alpha
  par(4) = 1.0; // v_0
}

void sys_stsol( Vector& out, double t )
{
  out(0) = 0.5;// + 1.3*sin(2*M_PI*t);
  out(1) = 0.5;// + 1.3*(2*M_PI*t+2*M_PI/3);
  out(2) = 0.5;// + 1.3*sin(2*M_PI*t+4*M_PI/3);
  out(3) = 2.0;// - 0.5*cos(2*M_PI*t+2*M_PI/3);
  out(4) = 2.0;// - 0.5*cos(2*M_PI*t+4*M_PI/3);
}

}
