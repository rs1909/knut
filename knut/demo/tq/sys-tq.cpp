/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <cmath>
#include <iostream>
#include "knutsys.h"

using namespace std;

// system definition
// par(0) = T the period length (=tau)
// par(1) = k1 cutting coefficient
// par(3) = which period are we in
// ...
// x[0] = xdot(t)
// x[1] = x(t)
// x[2] = x(t-tau_1)
// ...

// for the pd-analytic
// #define ZETA (0.0331166)
// #define H0 (1.0)
// #define TAU2 (0.0631373)

// #define cf_a (2.13727)
// #define cf_b (-1.28236/H0)
// #define cf_c (0.467529/H0/H0)

#define ZETA (0.015)
#define H0 (1.0)
#define TAU2 (0.1)

#define cf_a (120.0/59.0)
#define cf_b (-48.0/59.0)
#define cf_c (35.0/177.0)

static inline double cf_p( double x )
{
// 	if( x > 0 ) return pow(x,0.75)/0.75;
// 	else return 0.0;
  return x*( cf_a + cf_b*x + cf_c*x*x );
}

static inline double d_cf_p( double x )
{
// 	if( x > 0 ) return pow(x,-0.25);
// 	else return 0.0;
  return cf_a + 2*cf_b*x + 3*cf_c*x*x;
}

static inline double dd_cf_p( double x )
{
// 	if( x > 0 ) return -0.25*pow(x,-1.25);
// 	else return 0.0;
  return 2*cf_b + 6*cf_c*x;
}

// SMOOTHING THE PER 2 THINGS
#define CC (20.0)

static inline double cg2( double x )
{
// 	if( x > 0 ) return 1.0;
// 	else return 0.0;
  return (1.0+tanh(CC*x))/2.0;
}

static inline double d_cg2( double x )
{
// 	return 0.0;
  return CC*(1.0/pow(cosh(CC*x),2.0))/2.0;
}

static inline double dd_cg2( double x )
{
// 	return 0.0;
  return -CC*CC*(1.0/pow(cosh(CC*x),2.0))*tanh(CC*x);
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

extern "C"
{

int sys_ndim(){ return 2; }
int sys_npar(){ return 2; }
int sys_ntau(){ return 3; }
int sys_nderi() { return 2; }

void sys_p_tau( Array2D<double>& vout, const Array1D<double>& time, const Array1D<double>& par )
{
#define out(i) vout(i,idx)
  for (int idx=0; idx < time.Size(); ++idx)
  {
    out(0) = 0.0;
    out(1) = 1.0*par(0)/par(3);
    out(2) = 2.0*par(0)/par(3);
  }
#undef out
}

void sys_p_dtau( Array2D<double>& vout, const Array1D<double>& time, const Array1D<double>& par, int vp )
{
#define out(i) vout(i,idx)
  for (int idx=0; idx < time.Size(); ++idx)
  {
    switch( vp )
    {
      case 0:
        out(0) = 0.0;
        out(1) = 1.0/par(3);
        out(2) = 2.0/par(3);
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
#undef out
}

void sys_p_rhs( Array2D<double>& vout, const Array1D<double>& time, const Array3D<double>& yy, const Array1D<double>& par, int sel )
{
#define out(i) vout(i,idx)
#define xx(i,j) yy(i,j,idx)
  for (int idx=0; idx < time.Size(); ++idx)
  {
    const double t = time(idx);
    double g;
    if( (t < 0)||(t > 1) ) cout << "rhs: t is not element of the interval\n";
    double tt = par(3)*t - floor(par(3)*t);
    if( tt <= TAU2 )
    {
      g = 1.0;
    }else
    {
      g = 0.0;
    }

    out(0) = xx(1,0);
    out(1) = -xx(0,0) - 2*ZETA*xx(1,0) + g*par(1)*scf( xx(0,0), xx(0,1), xx(0,2) );
  }
#undef out
#undef xx
}

void sys_p_deri( Array3D<double>& mout, const Array1D<double>& time, const Array3D<double>& yy, const Array1D<double>& par, int sel, int nx, const int* vx, int np, const int* vp, const Array3D<double>& ww )
{
#define out(i,j) mout(i,j,idx)
#define xx(i,j) yy(i,j,idx)
#define vv(i,j) ww(i,j,idx)
  for (int idx=0; idx < time.Size(); ++idx)
  {
    const double t = time(idx);
    double g;
    if( (t < 0)||(t > 1) ) cout << "deri: t is not element of the interval\n";
    double tt = par(3)*t - floor(par(3)*t);
    if( tt <= TAU2 ){
      g = 1.0;
    }else{
      g = 0.0;
    }

    // derivatives w.r.t. the dependent variables: xx(t), xx(t-tau1), etc.
    if( (nx == 1) && (np == 0) )
    {
      switch( vx[0] )
      {
        case 0:
          out(0,0) = 0.0; 
          out(0,1) = 1.0;
          out(1,0) = -1.0 + g*par(1)*d_scf_x0( xx(0,0), xx(0,1), xx(0,2) );
          out(1,1) = -2.0*ZETA;
          break;
        case 1:
          out(0,0) = 0.0; 
          out(0,1) = 0.0;
          out(1,0) = g*par(1)*d_scf_x1( xx(0,0), xx(0,1), xx(0,2) );
          out(1,1) = 0.0;
          break;
        case 2:
          out(0,0) = 0.0; 
          out(0,1) = 0.0;
          out(1,0) = g*par(1)*d_scf_x2( xx(0,0), xx(0,1), xx(0,2) );
          out(1,1) = 0.0;
          break;
        default:
          std::cout << "deri: not implemented\n";
          break;
      }
    }
    // derivatives w.r.t. the parameters, purely, so this results a vector
    if( (nx == 0) && (np == 1) )
    {
      switch( vp[0] )
      {
        case 0: //  T period length
          std::cout << "deri w.r.t. T: not implemented\n";
          break;
        case 1:
          out(0,0) = 0.0;
          out(1,0) = g*scf( xx(0,0), xx(0,1), xx(0,2) );
          break;
        default:
          cout << "deri: not implemented\n";
          break;
      }
    }

    // second derivatives wrt. xx
    if( (nx == 2) && (np == 0) )
    {
      switch( vx[0] )
      { // phi(.,vx[0])
        case 0:
          switch( vx[1] )
          {
            case 0:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x0x0( xx(0,0), xx(0,1), xx(0,2) )*vv(0,0);
              out(1,1) = 0.0;
              break;
            case 1:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x0x1( xx(0,0), xx(0,1), xx(0,2) )*vv(0,0);
              out(1,1) = 0.0;
              break;
            case 2:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x0x2( xx(0,0), xx(0,1), xx(0,2) )*vv(0,0);
              out(1,1) = 0.0;
              break;
            default:
              cout << "deri: not implemented\n";
              break;
          }
          break;
        case 1:
          switch( vx[1] )
          {
            case 0:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x0x1( xx(0,0), xx(0,1), xx(0,2) )*vv(0,1);
              out(1,1) = 0.0;
              break;
            case 1:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x1x1( xx(0,0), xx(0,1), xx(0,2) )*vv(0,1);
              out(1,1) = 0.0;
              break;
            case 2:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x1x2( xx(0,0), xx(0,1), xx(0,2) )*vv(0,1);
              out(1,1) = 0.0;
              break;
            default:
              cout << "deri: not implemented\n";
              break;
          }
          break;
        case 2:
          switch( vx[1] )
          {
            case 0:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x0x2( xx(0,0), xx(0,1), xx(0,2) )*vv(0,2);
              out(1,1) = 0.0;
              break;
            case 1:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x1x2( xx(0,0), xx(0,1), xx(0,2) )*vv(0,2);
              out(1,1) = 0.0;
              break;
            case 2:
              out(0,0) = 0.0;
              out(0,1) = 0.0;
              out(1,0) = g*par(1)*d_scf_x2x2( xx(0,0), xx(0,1), xx(0,2) )*vv(0,2);
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
    if( (nx == 1) && (np == 1) )
    {
      switch( vp[0] )
      {
        case 0: //  T period length
          cout << "deri w.r.t. T: not implemented\n";
          break;
        case 1:
          switch( vx[0] )
          {
            case 0:
              out(0,0) = 0.0; 
              out(0,1) = 0.0;
              out(1,0) = g*d_scf_x0( xx(0,0), xx(0,1), xx(0,2) );
              out(1,1) = 0.0;
              break;
            case 1:
              out(0,0) = 0.0; 
              out(0,1) = 0.0;
              out(1,0) = g*d_scf_x1( xx(0,0), xx(0,1), xx(0,2) );
              out(1,1) = 0.0;
              break;
            case 2:
              out(0,0) = 0.0; 
              out(0,1) = 0.0;
              out(1,0) = g*d_scf_x2( xx(0,0), xx(0,1), xx(0,2) );
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
#undef out
#undef xx
#undef ww
}

void sys_stpar( Vector& par )
{
  par(0) = 14.75;
  par(1) = 0.00;
}

void sys_stsol( Vector& out, double t )
{
  out(0) = 0.0;
  out(1) = 0.0;
}

} // extern "C"
