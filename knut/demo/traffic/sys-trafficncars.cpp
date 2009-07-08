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

#define NCARS 17

extern "C"
{

int sys_ndim() { return 2*NCARS-1; }
int sys_npar() { return 5; }
int sys_ntau() { return 2; }
int sys_nderi() { return 0; }

void sys_tau( Vector& out, double t, const Vector& par )
{
	out(0) = 0.0;
	out(1) = par(2);
}

void sys_dtau( Vector& out, double t, const Vector& par, int vp )
{
	out(0) = 0.0;
	if( vp == 2 ) out(1) = 1.0;
	else out(1) = 0.0;
}

static inline double V( double x )
{
	return pow( ( x - 1 ) / 1, 3.0 ) / ( pow( ( x - 1 ) / 1, 3.0 ) + 1);
}

void sys_rhs( Vector& out, double t, const Matrix& yy, const Vector& par )
{
#define xx(i,j) yy(i-1,j-1)
#define f(i,j) out(i-1)

	for( int i = 1; i < NCARS; i++ )
	{
		if ( xx(i+NCARS,2) > 1 )
		{
			f(i,1) = par(3) * ( par(4) * V(xx(i+NCARS,2)) - xx(i,1) ); 
		}else
		{
			f(i,1) = -par(3) * xx(i,1);
		}
	}

	double sumi = NCARS*par(1);
	for( int i = 1; i < NCARS; i++ )
	{
		sumi = sumi - xx(i+NCARS,2);
	}
	
	if( sumi > 1 )
	{
		f(NCARS,1) = par(3) * ( par(4) * V(sumi) - xx(NCARS,1) );
	}else
	{
		f(NCARS,1) = - par(3) * xx(NCARS,1);
	}

	for( int i = 1; i < NCARS; i++ )
	{
		f(i+NCARS,1) = xx(i+1,1) - xx(i,1);
	}
	
#undef f
#undef xx
}

void sys_deri( Matrix &out, double t, const Matrix& x, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{

}

void sys_stpar( Vector& par )
{
	par(0) = 1.1; // T the period length 
	par(1) = 2.0; // L/NCARS
	par(2) = 1.0; // the delay
	par(3) = 1.0; // alpha
	par(4) = 1.0; // v_0
}

void sys_stsol( Vector& out, double t )
{
#define f(i) out(i-1)
	for( int i=1; i<NCARS+1; i++ ){
	    f(i) = 0.5;
	}
	for( int i=1; i<NCARS; i++ ){
	    f(i+NCARS) = 2.0;
    }
#undef f
}

}
