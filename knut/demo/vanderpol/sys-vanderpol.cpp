/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/

#include <cmath>
#include "knutsys.h"

extern "C"
{

int sys_ndim(){ return 2; }
int sys_npar(){ return 2; }
int sys_ntau(){ return 1; }
int sys_nderi(){ return 0; }


void sys_tau( KNVector& out, double t, const KNVector& par )
{
  out(0) = 0.0;
}

void sys_dtau( KNVector& out, double t, const KNVector& par, int vp )
{
	out(0) = 0.0;
}

void sys_rhs( KNVector& out, double t, const KNMatrix& x, const KNVector& par )
{
	out(0) = x(1,0);
	out(1) = -x(0,0) - par(1)*(x(0,0)*x(0,0) - 1.0)*x(1,0);
}

void sys_deri( KNMatrix &out, double t, const KNMatrix& x, const KNVector& par, 
	       int nx, const int* vx, int np, const int* vp, const KNMatrix& vv )
{

}

void sys_stpar( KNVector& par )
{
	par(0) = 2.0*M_PI;
	par(1) = 0.1;
}

void sys_stsol( KNVector& out, double t )
{
	out(0) = sin(2.0*M_PI*t);
	out(1) = cos(2.0*M_PI*t);
}

} // extern "C"

