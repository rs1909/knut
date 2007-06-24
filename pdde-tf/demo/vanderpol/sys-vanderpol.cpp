/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/

#include <cmath>
#include "pddesys.h"

extern "C"
{

int sys_ndim(){ return 2; }
int sys_npar(){ return 2; }
int sys_ntau(){ return 1; }
int sys_nderi(){ return 0; }


void sys_tau( Vector& out, double t, const Vector& par )
{
  out(0) = 0.0;
}

void sys_dtau( Vector& out, double t, const Vector& par, int vp )
{
	out(0) = 0.0;
}

void sys_rhs( Vector& out, double t, const Matrix& x, const Vector& par )
{
	out(0) = x(1,0);
	out(1) = -x(0,0) - par(1)*(x(0,0)*x(0,0) - 1.0)*x(1,0);
}

void sys_deri( Matrix &out, double t, const Matrix& x, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{

}

void sys_stpar( Vector& par )
{
	par(0) = 2.0*M_PI;
	par(1) = 0.1;
}

void sys_stsol( Vector& out, double t )
{
	out(0) = sin(2.0*M_PI*t);
	out(1) = cos(2.0*M_PI*t);
}

} // extern "C"

