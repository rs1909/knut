/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <cmath>
#include "knutsys.h"

// system definition for the turning process
// par(0) = T the period length 
// par(1) = k1 cutting coefficient
// par(2) = the delay
// ...
// x[0] = xdot(t)
// x[1] = x(t)
// x[2] = x(t-tau_1)
// ...

#define ZETA (0.015)
#define H0 (1.0)

static inline double cf( double x )
{
	if( x > 0 ) return pow(x,0.75)/0.75;
	else
	{
		std::cout<<".";
		return 0.0;
	}
}

extern "C"
{

int sys_ndim(){ return 2; }
int sys_npar(){ return 3; }
int sys_ntau(){ return 2; }
int sys_nderi(){ return 0; }


void sys_tau( KNVector& out, double t, const KNVector& par )
{
  out(0) = 0.0;
  out(1) = par(2);
}

void sys_dtau( KNVector& out, double t, const KNVector& par, int vp )
{
	out(0) = 0.0;
	if( vp == 2 ) out(1) = 1.0;
	else out(1) = 0.0;
}

void sys_rhs( KNVector& out, double t, const KNMatrix& x, const KNVector& par )
{
	out(0) = x(1,0);
	out(1) = -x(0,0) - 2*ZETA*x(1,0) + par(1)* (cf(H0+x(0,1)-x(0,0))-cf(H0));
}

void sys_deri( KNMatrix &out, double t, const KNMatrix& x, const KNVector& par, 
	       int nx, const int* vx, int np, const int* vp, const KNMatrix& vv )
{

}

void sys_stpar( KNVector& par )
{
	par(0) = 0.4;
	par(1) = 0.0;
	par(2) = 1.5*M_PI;
}

void sys_stsol( KNVector& out, double t )
{
	out(0) = 0.0;
	out(1) = 0.0;
}

} // extern "C"

