#include <cmath>
#include "pddesys.h"

// xx:  Ax Ay    N
//      0  1     2  3     4  5  6   7
// par: T  alpha TT kappa PP bb tau omega_0

extern "C"
{

int sys_ndim(){ return 3; }
int sys_npar(){ return 8; }
int sys_ntau(){ return 2; }
int sys_nderi(){ return 0; }


void sys_p_tau( Array2D<double>& out, Array1D<double>& time, const Array1D<double>& par )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    out(0,idx) = 0.0;
    out(1,idx) = par(6);
  }
}

void sys_p_dtau( Array2D<double>& out, Array1D<double>& time, const Array1D<double>& par, int vp )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    // it is a constant delay
    out(0,idx) = 0.0;
    if( vp == 6 ) out(1,idx) = 1.0;
    else out(1,idx) = 0.0;
  }
}

void sys_p_rhs( Array2D<double>& out, Array1D<double>& time, const Array3D<double>& xx, const Array1D<double>& par )
{
  // parameters
  const double alpha = par(1);
  const double TT = par(2);
  const double kappa = par(3);
  const double PP = par(4);
  const double bb = par(5);
  const double tau = par(6);
  const double omega_0 = par(7);

  for (int idx=0; idx < time.Size(); ++idx)
  {
    // state variables
    const double Are = xx(0,0,idx);
    const double Aim = xx(1,0,idx);
    const double NN = xx(2,0,idx);
    const double AreD = xx(0,1,idx);
    const double AimD = xx(1,1,idx);

    out(0,idx) = Aim*bb - Aim*alpha*NN + Are*NN + AreD*kappa*cos((bb + omega_0)*tau) + AimD*kappa*sin((bb + omega_0)*tau);
    out(1,idx) = -(Are*bb) + Aim*NN + alpha*Are*NN + AimD*kappa*cos((bb + omega_0)*tau) - AreD*kappa*sin((bb + omega_0)*tau);
    out(2,idx) = (-NN - (Aim*Aim + Are*Are)*(1 + 2*NN) + PP)/TT;
  }
}

void sys_p_deri( Array3D<double>& out, Array1D<double>& time, const Array3D<double>& xx, const Array1D<double>& par, int nx, const int* vx, int np, const int* vp, const Array3D<double>& vv )
{

}

void sys_stpar( Vector& par )
{
  // bb controls the solution
  const double bb      = -0.00483423*0.98;
  const double period  = 501.0;
  const double alpha   = 4.0;
  const double TT      = 1000.0;
  const double PP      = 0.001;
  const double tau     = 1000.0;
  const double omega_0 = -1.0/tau;
  const double kappa = -bb/(alpha*cos((omega_0+bb)*tau)+sin((omega_0+bb)*tau));

  par(0) = period;
  par(1) = alpha;
  par(2) = TT;
  par(3) = kappa;
  par(4) = PP;
  par(5) = bb;
  par(6) = tau;
  par(7) = omega_0;

//   std::cout<<"kappa "<< kappa << "\n";
}

void sys_stsol( Vector& out, double t )
{
  Vector par(8);
  sys_stpar(par);
  // parameters
  const double alpha = par(1);
  const double TT = par(2);
  const double kappa = par(3);
  const double PP = par(4);
  const double bb = par(5);
  const double tau = par(6);
  const double omega_0 = par(7);
  const double cs = cos(omega_0*tau+bb);
  const double sn = sin(omega_0*tau+bb);

  out(0) = sqrt( ((PP*alpha-bb)*cs + PP*sn ) / ((2*bb+alpha)*cs + sn) );
  out(1) = 0.0;
  out(2) = bb/(alpha+tan(omega_0*tau+bb));
//   std::cout<<"Are "<< out(0) << " Aim "<< out(1) << " N " << out(2) <<"\n";
}

} // extern "C"
