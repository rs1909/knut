/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <cmath>
#include "pddesys.h"

using namespace std;

// system definition
// par(0) = T the period length (=tau)
// par(1) = k1 cutting coefficient
// ...
// x[0] = xdot(t)
// x[1] = x(t)
// x[2] = x(t-tau_1)
// ...

#define ZETA (0.0331166)
#define H0 (1.0)
#define TAU2 (0.0631373)

#define cf_a (2.13727)
#define cf_b (-1.28236/H0)
#define cf_c (0.467529/H0/H0)

static inline double cf_p( double x )
{
	if( x > 0 ) return pow(x,0.75)/0.75;
	else return 0.0;
//  return x*( cf_a + cf_b*x + cf_c*x*x );
}

static inline double d_cf_p( double x )
{
	if( x > 0 ) return pow(x,-0.25);
	else return 0.0;
//  return cf_a + 2*cf_b*x + 3*cf_c*x*x;
}

static inline double dd_cf_p( double x )
{
	if( x > 0 ) return -0.25*pow(x,-1.25);
	else return 0.0;
//  return 2*cf_b + 6*cf_c*x;
}

// SMOOTHING THE PER 2 THINGS
#define CC (10.0)

static inline double cg2( double x )
{
	if( x > 0 ) return 1.0;
	else return 0.0;
//  return (1.0+tanh(CC*x))/2.0;
}

static inline double d_cg2( double x )
{
	return 0.0;
//  return CC*(1.0/pow(cosh(CC*x),2.0))/2.0;
}

static inline double dd_cg2( double x )
{
	return 0.0;
//  return -CC*CC*(1.0/pow(cosh(CC*x),2.0))*tanh(CC*x);
}

static inline double cf( double x )
{
  return cg2(x)*cf_p(x); ///AAA*pow(BBB*x,2.7)/(1.0+BBB*BBB*x*x);
}

static inline double d_cf( double x )
{
  return d_cg2(x)*cf_p(x) + cg2(x)*d_cf_p(x); // BBB*AAA*(2.7*pow(BBB*x,1.7)+0.7*pow(BBB*x,3.7))/pow(1.0+BBB*BBB*x*x,2.0);
}

static inline double dd_cf( double x )
{
  return dd_cg2(x)*cf_p(x) + 2*d_cg2(x)*d_cf_p(x) + cg2(x)*dd_cf_p(x);
}

static inline double scf( double x0, double x1, double x2 )
{
  return  cg2(H0+x2-x1)*cf(H0+x1-x0) + cg2(-H0+x1-x2)*cf(2*H0+x2-x0);
}

static inline double d_scf_x0( double x0, double x1, double x2 )
{
  return -(cg2(H0 - x1 + x2)*d_cf(H0 - x0 + x1)) - 
   cg2(-H0 + x1 - x2)*d_cf(2*H0 - x0 + x2);
}

static inline double d_scf_x1( double x0, double x1, double x2 )
{
  return  cg2(H0 - x1 + x2)*d_cf(H0 - x0 + x1) + 
   cf(2*H0 - x0 + x2)*d_cg2(-H0 + x1 - x2) - 
   cf(H0 - x0 + x1)*d_cg2(H0 - x1 + x2);
}

static inline double d_scf_x2( double x0, double x1, double x2 )
{
  return  cg2(-H0 + x1 - x2)*d_cf(2*H0 - x0 + x2) - 
   cf(2*H0 - x0 + x2)*d_cg2(-H0 + x1 - x2) + 
   cf(H0 - x0 + x1)*d_cg2(H0 - x1 + x2);
}

static inline double d_scf_x0x0( double x0, double x1, double x2 )
{
  return  cg2(H0 - x1 + x2)*dd_cf(H0 - x0 + x1) + 
   cg2(-H0 + x1 - x2)*dd_cf(2*H0 - x0 + x2);
}

static inline double d_scf_x1x1( double x0, double x1, double x2 )
{
  return  -2*d_cf(H0 - x0 + x1)*d_cg2(H0 - x1 + x2) + 
   cg2(H0 - x1 + x2)*dd_cf(H0 - x0 + x1) + 
   cf(2*H0 - x0 + x2)*dd_cg2(-H0 + x1 - x2) + 
   cf(H0 - x0 + x1)*dd_cg2(H0 - x1 + x2);
}

static inline double d_scf_x2x2( double x0, double x1, double x2 )
{
  return -2*d_cf(2*H0 - x0 + x2)*d_cg2(-H0 + x1 - x2) + 
   cg2(-H0 + x1 - x2)*dd_cf(2*H0 - x0 + x2) + 
   cf(2*H0 - x0 + x2)*dd_cg2(-H0 + x1 - x2) + 
   cf(H0 - x0 + x1)*dd_cg2(H0 - x1 + x2);
}

static inline double d_scf_x0x1( double x0, double x1, double x2 )
{
  return -d_cf(2*H0 - x0 + x2)*d_cg2(-H0 + x1 - x2) + d_cf(H0 - x0 + x1)*d_cg2(H0 - x1 + x2) - cg2(H0 - x1 + x2)*dd_cf(H0 - x0 + x1);
}

static inline double d_scf_x0x2( double x0, double x1, double x2 )
{
  return d_cf(2*H0 - x0 + x2)*d_cg2(-H0 + x1 - x2) - 
   d_cf(H0 - x0 + x1)*d_cg2(H0 - x1 + x2) - 
   cg2(-H0 + x1 - x2)*dd_cf(2*H0 - x0 + x2);
}

static inline double d_scf_x1x2( double x0, double x1, double x2 )
{
  return d_cf(2*H0 - x0 + x2)*d_cg2(-H0 + x1 - x2) + d_cf(H0 - x0 + x1)*d_cg2(H0 - x1 + x2) - 
   cf(2*H0 - x0 + x2)*dd_cg2(-H0 + x1 - x2) - cf(H0 - x0 + x1)*dd_cg2(H0 - x1 + x2);
}

int Sys::ndim(){ return 2; }
int Sys::npar(){ return 2; }
int Sys::ntau(){ return 3; }

void Sys::tau( Vector& out, double t, const Vector& par )
{
	out(0) = 0.0;
	out(1) = 0.5*par(0);
	out(2) = par(0);
}

void Sys::dtau( Vector& out, double t, const Vector& par, int vp )
{
	switch( vp )
	{
		case 0:
			out(0) = 0.0;
			out(1) = 0.5;
			out(2) = 1.0;
			break;
		case 1:
			out(0) = 0.0;
			out(1) = 0.0;
			out(2) = 0.0;
			break;
		default:
			cout << "dtau: not implemented\n";
			break;
	}
}

void Sys::rhs( Vector& out, double t, const Matrix& x, const Vector& par )
{
	double g;
	if( (t < 0)||(t > 1) ) cout << "rhs: t is not element of the interval\n";
	double tt = 2*t - floor(2*t);
	if( tt <= TAU2 )
	{
		g = 1.0;
	}else
	{
		g = 0.0;
	}
   
  out(0) = x(1,0);
  out(1) = -x(0,0) - 2*ZETA*x(1,0) + g*par(1)*scf( x(0,0), x(0,1), x(0,2) );
}

void Sys::deri( Matrix &out, double t, const Matrix& x, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{
	double g;
	if( (t < 0)||(t > 1) ) cout << "deri: t is not element of the interval\n";
	double tt = 2*t - floor(2*t);
	if( tt <= TAU2 ){
		g = 1.0;
	}else{
	g = 0.0;
  }
  
  // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
  if( (nx == 1) && (np == 0) ){
    switch( vx[0] ){
    case 0:
      out(0,0) = 0.0; 
      out(0,1) = 1.0;
      out(1,0) = -1.0 + g*par(1)*d_scf_x0( x(0,0), x(0,1), x(0,2) );
      out(1,1) = -2.0*ZETA;
      break;
    case 1:
      out(0,0) = 0.0; 
      out(0,1) = 0.0;
      out(1,0) = g*par(1)*d_scf_x1( x(0,0), x(0,1), x(0,2) );
      out(1,1) = 0.0;
      break;
    case 2:
      out(0,0) = 0.0; 
      out(0,1) = 0.0;
      out(1,0) = g*par(1)*d_scf_x2( x(0,0), x(0,1), x(0,2) );
      out(1,1) = 0.0;
      break;
    default:
      cout << "deri: not implemented\n";
      break;
    }		
  }     
  // derivatives w.r.t. the parameters, purely, so this results a vector
  if( (nx == 0) && (np == 1) ){
    switch( vp[0] ){
    case 0: //  T period length
      cout << "deri w.r.t. T: not implemented\n";
      break;
    case 1:
      out(0) = 0.0;
      out(1) = g*scf( x(0,0), x(0,1), x(0,2) );
      break;
    default:
      cout << "deri: not implemented\n";
      break;
    }
  }
  // second derivatives wrt. xx
  if( (nx == 2) && (np == 0) ){
    switch( vx[0] ){ // phi(.,vx[0])
    case 0:
      switch( vx[1] ){
      case 0:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x0x0( x(0,0), x(0,1), x(0,2) )*vv(0,0);
	out(1,1) = 0.0;
	break;
      case 1:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x0x1( x(0,0), x(0,1), x(0,2) )*vv(0,0);
	out(1,1) = 0.0;
	break;
      case 2:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x0x2( x(0,0), x(0,1), x(0,2) )*vv(0,0);
	out(1,1) = 0.0;
	break;
      default:
	cout << "deri: not implemented\n";
	break;
      }
      break;
    case 1:
      switch( vx[1] ){
      case 0:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x0x1( x(0,0), x(0,1), x(0,2) )*vv(0,1);
	out(1,1) = 0.0;
	break;
      case 1:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x1x1( x(0,0), x(0,1), x(0,2) )*vv(0,1);
	out(1,1) = 0.0;
	break;
      case 2:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x1x2( x(0,0), x(0,1), x(0,2) )*vv(0,1);
	out(1,1) = 0.0;
	break;
      default:
	cout << "deri: not implemented\n";
	break;
      }
      break;
    case 2:
      switch( vx[1] ){
      case 0:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x0x2( x(0,0), x(0,1), x(0,2) )*vv(0,2);
	out(1,1) = 0.0;
	break;
      case 1:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x1x2( x(0,0), x(0,1), x(0,2) )*vv(0,2);
	out(1,1) = 0.0;
	break;
      case 2:
	out(0,0) = 0.0;
	out(0,1) = 0.0;
	out(1,0) = g*par(1)*d_scf_x2x2( x(0,0), x(0,1), x(0,2) )*vv(0,2);
	out(1,1) = 0.0;
	break;
      default:
	cout << "deri: not implemented\n";
	break;
      }
      break;
    default:
      cout << "deri: not implemented\n";
      break;
    }
  }
  // mixed derivative wrt to xx and par
  if( (nx == 1) && (np == 1) ){
    switch( vp[0] ){
    case 0: //  T period length
      cout << "deri w.r.t. T: not implemented\n";
      break;
    case 1:
      switch( vx[0] ){
      case 0:
	out(0,0) = 0.0; 
	out(0,1) = 0.0;
	out(1,0) = g*d_scf_x0( x(0,0), x(0,1), x(0,2) );
	out(1,1) = 0.0;
	break;
      case 1:
	out(0,0) = 0.0; 
	out(0,1) = 0.0;
	out(1,0) = g*d_scf_x1( x(0,0), x(0,1), x(0,2) );
	out(1,1) = 0.0;
	break;
      case 2:
	out(0,0) = 0.0; 
	out(0,1) = 0.0;
	out(1,0) = g*d_scf_x2( x(0,0), x(0,1), x(0,2) );
	out(1,1) = 0.0;
	break;
      default:
	cout << "deri: not implemented\n";
	break;
      }		
      break;
    default:
      cout << "deri: not implemented\n";
      break;
    }
  }
}