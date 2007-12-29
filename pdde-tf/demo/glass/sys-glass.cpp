#include <cmath>
#include "pddesys.h"

extern "C"
int sys_ndim(){ return 1; }

extern "C"
{

int sys_npar(){ return 4; }
int sys_ntau(){ return 2; }
int sys_nderi(){ return 2; }

void sys_p_tau( Array2D<double>& out, Array1D<double>& time, const Array1D<double>& par )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    out(0,idx) = 0.0; out(1,idx) = par(3);
  }
}

void sys_p_dtau( Array2D<double>& out, Array1D<double>& time, const Array1D<double>& par, int vp )
{
  for (int idx=0; idx < time.Size(); ++idx)
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

void sys_p_rhs( Array2D<double>& out, Array1D<double>& time, const Array3D<double>& x, const Array1D<double>& par )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    out(0,idx) = par(1)*x(0,0,idx) + par(2)*x(0,1,idx)/(1+pow(x(0,1,idx), 10.0));
  }
}

void sys_p_deri( Array3D<double>& out, Array1D<double>& time, const Array3D<double>& x, const Array1D<double>& par,
         int nx, const int* vx, int np, const int* vp, const Array3D<double>& vv )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
    if( (nx == 1) && (np == 0) ) {
      switch( vx[0] ) {
        case 0:
          out(0,0,idx) = par(1);
          break;
        case 1:
          out(0,0,idx) = par(2)*(1.0 - 9.0*pow(x(0,1,idx), 10.0) )/( (1+pow(x(0,1,idx), 10.0))*(1+pow(x(0,1,idx), 10.0)) );
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
              out(0,0,idx) = 10.0*par(2)*pow(x(0,1,idx), 9.0)*(9.0*pow(x(0,1,idx), 10.0) - 11.0)/
                                    pow((1.0+pow(x(0,1,idx), 10.0)), 3.0) * vv(0,1,idx); break;
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

void sys_stpar( Vector& par )
{
  par(0) = 2.0; // 5.3918; 
  par(1) = -1.0; par(2) = 1.5; par(3) = 2.0;
}

void sys_stsol( Vector& out, double t )
{
  out(0) = pow((1.0-1.5)/(-1.0), 1.0/10.0);
}

} // extern "C"
