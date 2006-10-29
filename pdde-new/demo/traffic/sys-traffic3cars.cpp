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

int sys_ndim() { return 5; }
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

	if ( xx(4,2) > 1 ){
		f(1,1) = par(3) * ( par(4) * V(xx(4,2)) - xx(1,1) ); 
	}else{
		f(1,1) = - par(3) * xx(1,1);
	}
	
	if( xx(5,2) > 1 ){
		f(2,1) = par(3) * ( par(4) * V(xx(5,2)) - xx(2,1) ); 
	}else{
		f(2,1) = - par(3) * xx(2,1);
	}
	
	if( ( 3*par(1) - xx(5,2) - xx(4,2) ) > 1 ){
		f(3,1) = par(3) * ( par(4) * V( 3*par(1) - xx(5,2) - xx(4,2) ) - xx(3,1) );
	}else{
		f(3,1) = - par(3) * xx(3,1);
	}
	
	f(4,1) = xx(2,1) - xx(1,1);
	f(5,1) = xx(3,1) - xx(2,1);
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
	out(0) = 0.5;// + 1.3*sin(2*M_PI*t);
	out(1) = 0.5;// + 1.3*(2*M_PI*t+2*M_PI/3);
	out(2) = 0.5;// + 1.3*sin(2*M_PI*t+4*M_PI/3);
	out(3) = 2.0;// - 0.5*cos(2*M_PI*t+2*M_PI/3);
	out(4) = 2.0;// - 0.5*cos(2*M_PI*t+4*M_PI/3);
}

}
