#include <cmath>
#include "knutsys.h"

extern "C"
{

size_t sys_ndim(){ return 1; }
size_t sys_npar(){ return 4; }
size_t sys_ntau(){ return 2; }
size_t sys_nderi(){ return 2; }

void sys_p_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par )
{
  for (size_t idx=0; idx < time.size(); ++idx)
  {
    out(0,idx) = 0.0; out(1,idx) = par(3);
  }
}

void sys_p_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par, size_t vp )
{
  for (size_t idx=0; idx < time.size(); ++idx)
  {
    switch( vp ) {
      case 0:
        out(0,idx) = 0.0; out(1,idx) = 0.0; break;
      case 1:
        out(0,idx) = 0.0; out(1,idx) = 0.0; break;
      case 2:
        out(0,idx) = 0.0; out(1,idx) = 0.0; break;
      case 3:
        out(0,idx) = 0.0; out(1,idx) = 1.0; break;
    }
  }
}

void sys_p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel)
{
  for (size_t idx=0; idx < time.size(); ++idx)
  {
    out(0,idx) = par(1)*x(0,0,idx) + par(2)*x(0,1,idx)/(1+pow(x(0,1,idx), 10));
  }
}

void sys_p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel,
         size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv )
{
  for (size_t idx=0; idx < time.size(); ++idx)
  {
    // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
    if( (nx == 1) && (np == 0) ) {
      switch( vx[0] ) {
        case 0:
          out(0,0,idx) = par(1);
          break;
        case 1:
          out(0,0,idx) = par(2)*(1.0/(1+pow(x(0,1,idx), 10)) - 10.0*pow(x(0,1,idx), 10)/( (1+pow(x(0,1,idx), 10))*(1+pow(x(0,1,idx), 10)) ) );
          break;
      }
    }
    // derivatives w.r.t. the parameters, purely, so this results a vector
    if( (nx == 0) && (np == 1) )
    {
      switch( vp[0] )
      {
        case 1:
          out(0,0,idx) = x(0,0,idx); break;
        case 2:
          out(0,0,idx) = x(0,1,idx)/(1+pow(x(0,1,idx), 10)); break;
        case 3:
          out(0,0,idx) = 0.0; break;
      }
    }
    // second derivatives w.r.t. x
    if( (nx == 2) && (np == 0) ) {
      switch( vx[0] ) {
        case 0:
          switch( vx[1] ) {
            case 0:
              out(0,0,idx) = 0.0; break;
            case 1:
              out(0,0,idx) = 0.0; break;
          } break;
        case 1:
          switch( vx[1] ) {
            case 0:
              out(0,0,idx) = 0.0; break;
            case 1:
              out(0,0,idx) = 10.0*par(2)*pow(x(0,1,idx), 9)*(9.0*pow(x(0,1,idx), 10) - 11.0)/
                                    pow((1.0+pow(x(0,1,idx), 10)), 3) * vv(0,1,idx); break;
          } break;
      }
    }
    // mixed derivative w.r.t. x and the parameters
    if( (nx == 1) && (np == 1) ) {
      switch( vp[0] ) {
        case 1:
          switch( vx[0] ) {
            case 0:
              out(0,0,idx) = 1.0; break;
            case 1:
              out(0,0,idx) = 0.0; break;
          } break;
        case 2:
          switch( vx[0] ) {
            case 0:
              out(0,0,idx) = 0.0; break;
            case 1:
              out(0,0,idx) = (1.0 - 9.0*pow(x(0,1,idx), 10.0) )/
                        ((1.0 + pow(x(0,1,idx), 10.0))*(1.0+pow(x(0,1,idx), 10.0)) ); break;
          } break;
        case 3:
          switch( vx[0] ) {
            case 0:
              out(0,0,idx) = 0.0; break;
            case 1:
              out(0,0,idx) = 0.0; break;
          } break;
      }
    }
  }
}

void sys_stpar( KNVector& par )
{
  par(0) = 2.0; // 5.3918; 
  par(1) = -1.0; par(2) = 1.5; par(3) = 2.0;
}

void sys_stsol( KNVector& out, double t )
{
  out(0) = pow((1.0-1.5)/(-1.0), 1.0/10.0);
}

} // extern "C"
