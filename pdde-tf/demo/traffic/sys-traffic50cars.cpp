/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <cmath>
#include "pddesys.h"

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

int sys_ndim() { return 99; }
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

	for( int i = 1; i < 50; i++ )
	{
		if ( xx(i+50,2) > 1 )
		{
			f(i,1) = par(3) * ( par(4) * V(xx(i+50,2)) - xx(i,1) ); 
		}else
		{
			f(i,1) = -par(3) * xx(i,1);
		}
	}

	double sumi = 50*par(1);
	for( int i = 1; i < 50; i++ )
	{
		sumi = sumi - xx(i+50,2);
	}
	
	if( sumi > 1 )
	{
		f(50,1) = par(3) * ( par(4) * V(sumi) - xx(50,1) );
	}else
	{
		f(50,1) = - par(3) * xx(50,1);
	}

	for( int i = 1; i < 50; i++ )
	{
		f(i+50,1) = xx(i+1,1) - xx(i,1);
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
	par(1) = 2.0; // L/N
	par(2) = 1.0; // the delay
	par(3) = 1.0; // alpha
	par(4) = 1.0; // v_0
}

void sys_stsol( Vector& out, double t )
{
#define f(i) out(i-1)
	for( int i=1; i<51; i++ ){
	    f(i) = 0.5;
	}
	for( int i=1; i<50; i++ ){
	    f(i+50) = 2.0;
    }
#undef f
}

}
