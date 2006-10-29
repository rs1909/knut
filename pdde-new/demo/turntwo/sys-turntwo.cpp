/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <cmath>
#include "pddesys.h"

using namespace std;

// MODEL OF TURNING WITH FLYOVER
//
// system definition
// par(0) = the period length
// par(1) = cutting coefficient
// par(2) = the delay
// ...

#define ZETA (0.01)
#define H0 (1.0)

// #define cf_a (2.13727)
// #define cf_b (-1.28236/H0)
// #define cf_c (0.467529/H0/H0)

#define cf_a (270.0/161.0)
#define cf_b (-72.0/161.0)
#define cf_c (5.0/69.0)

// #define cf_a (120.0/59.0)
// #define cf_b (-48.0/59.0/H0)
// #define cf_c (35.0/177.0/H0/H0)

static inline double cf_p( double x )
{
	return x*( cf_a + cf_b*x + cf_c*x*x );
}

static inline double d_cf_p( double x )
{
	return cf_a + 2*cf_b*x + 3*cf_c*x*x;
}

static inline double dd_cf_p( double x )
{
	return 2*cf_b + 6*cf_c*x;
}

// WITH THE THREE-QURTER RULE
//
// static inline double cf_p( double x )
// {
// 	if( x > 0 ) return pow(x,0.75)/0.75;
// 	else return 0.0;
// }
// 
// static inline double d_cf_p( double x )
// {
// 	if( x > 0 ) return pow(x,-0.25);
// 	else return 0.0;
// }
// 
// static inline double dd_cf_p( double x )
// {
// 	if( x > 0 ) return -0.25*pow(x,-1.25);
// 	else return 0.0;
// }

// SMOOTH RHS

#define CC (15.0)

static inline double cg2( double x )
{
	return (1.0+tanh(CC*x))/2.0;
}

static inline double d_cg2( double x )
{
	return CC*(1.0/pow(cosh(CC*x),2.0))/2.0;
}

static inline double dd_cg2( double x )
{
	return -CC*CC*(1.0/pow(cosh(CC*x),2.0))*tanh(CC*x);
}

// NONSMOOTH RHS
//
// static inline double cg2( double x )
// {
// 	if( x > 0 ) return 1.0;
// 	else return 0.0;
// }
// 
// static inline double d_cg2( double x )
// {
//   return 0.0;
// }
// 
// static inline double dd_cg2( double x )
// {
//   return 0.0;
// }

static inline double cf( double x )
{
	return cg2(x)*cf_p(x);
}

static inline double d_cf( double x )
{
	return d_cg2(x)*cf_p(x) + cg2(x)*d_cf_p(x);
}

static inline double dd_cf( double x )
{
	return dd_cg2(x)*cf_p(x) + 2*d_cg2(x)*d_cf_p(x) + cg2(x)*dd_cf_p(x);
}

// CUTTING FORCE FUNCTION

#define AD (H0+x2-x1)
#define A1 (H0+x1-x0)
#define A2 (2*H0+x2-x0)

static inline double scf( double x0, double x1, double x2 )
{
	return  cf(A1)*cg2(AD) + cf(A2)*cg2(-AD);
}

// THIS IS THE SYSTEM DEFINITION

extern "C"
{

int sys_ndim(){ return 2; }
int sys_npar(){ return 3; }
int sys_ntau(){ return 3; }
int sys_nderi(){ return 0; }

void sys_tau( Vector& out, double t, const Vector& par )
{
	out(0) = 0.0;
	out(1) = 1.0*par(2);
	out(2) = 2.0*par(2);
}

void sys_dtau( Vector& out, double t, const Vector& par, int vp )
{
	switch( vp )
	{
		case 0:
		case 1:
			out(0) = 0.0;
			out(1) = 0.0;
			out(2) = 0.0;
			break;
		case 2:
			out(0) = 0.0;
			out(1) = 1.0;
			out(2) = 0.0;
			break;
		default:
			cout << "dtau: not implemented\n";
			break;
	}
}

void sys_rhs( Vector& out, double t, const Matrix& x, const Vector& par )
{
	out(0) = x(1,0);
	out(1) = -x(0,0) - 2*ZETA*x(1,0) + par(1)*scf( x(0,0), x(0,1), x(0,2) );
}

void sys_deri( Matrix &out, double t, const Matrix& x, const Vector& par,
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{
}

void sys_stpar( Vector& par )
{
	par(0) = 0.4;
	par(1) = 0.0;
	par(2) = 1.5*M_PI;
}

void sys_stsol( Vector& out, double t )
{
	out(0) = 0.0;
	out(1) = 0.0;
}

} // extern "C"
