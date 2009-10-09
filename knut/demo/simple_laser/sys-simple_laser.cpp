#include <cmath>
#include "knutsys.h"

// xx:  Ax Ay    N
//      0  1     2  3     4  5  6   7
// par: T  alpha TT kappa PP bb tau omega_0

extern "C"
{

int sys_ndim(){ return 3; }
int sys_npar(){ return 8; }
int sys_ntau(){ return 2; }
int sys_nderi(){ return 2; }

// VECTORIZED DEF.

void sys_p_tau( Array2D<double>& out, const Array1D<double>& time, const Array1D<double>& par )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    out(0,idx) = 0.0;
    out(1,idx) = par(6);
  }
}

void sys_p_dtau( Array2D<double>& out, const Array1D<double>& time, const Array1D<double>& par,                       const int vp )
{
  for (int idx=0; idx < time.Size(); ++idx)
  {
    // it is a constant delay
    out(0,idx) = 0.0;
    if( vp == 6 ) out(1,idx) = 1.0;
    else out(1,idx) = 0.0;
  }
}

void sys_p_rhs( Array2D<double>& out, const Array1D<double>& time,
                const Array3D<double>& xx, const Array1D<double>& par, int sel )
{
  // parameters
  const double alpha = par(1);
  const double TT = par(2);
  const double kappa = par(3);
  const double PP = par(4);
  const double bb = par(5);
  const double tau = par(6);
  const double omega_0 = par(7);
  const double SN = sin((bb + omega_0)*tau);
  const double CS = cos((bb + omega_0)*tau);

  for (int idx=0; idx < time.Size(); ++idx)
  {
    // state variables
    const double Are = xx(0,0,idx);
    const double Aim = xx(1,0,idx);
    const double NN = xx(2,0,idx);
    const double AreD = xx(0,1,idx);
    const double AimD = xx(1,1,idx);

    out(0,idx) = Aim*bb + AreD*CS*kappa - Aim*alpha*NN + Are*NN + AimD*kappa*SN;
    out(1,idx) = -(Are*bb) + AimD*CS*kappa + Aim*NN + alpha*Are*NN - AreD*kappa*SN;
    out(2,idx) = (-NN - (Aim*Aim + Are*Are)*(1 + 2*NN) + PP)/TT;
  }
}

void sys_p_deri( Array3D<double>& mout, const Array1D<double>& time,
                 const Array3D<double>& xx, const Array1D<double>& par, int sel,
                 const int nx, const int* vx, const int np, const int* vp,
                 const Array3D<double>& vv )
{
#define out(i,j) mout(i,j,idx)
const double alpha = par(1);
const double TT = par(2);
const double kappa = par(3);
const double PP = par(4);
const double bb = par(5);
const double tau = par(6);
const double omega_0 = par(7);
const double SN = sin((bb + omega_0)*tau);
const double CS = cos((bb + omega_0)*tau);

for (int idx=0; idx<time.Size(); ++idx)
{
  const double Are = xx(0,0,idx);
  const double Aim = xx(1,0,idx);
  const double NN = xx(2,0,idx);
  const double AreD = xx(0,1,idx);
  const double AimD = xx(1,1,idx);

  if ((nx==1)&&(np==0))
  {
    switch (vx[0])
    {
      case 0:
        //   w.r.t. x(t)
        out(0,0) = NN;
        out(0,1) = bb - NN*alpha;
        out(0,2) = Are - Aim*alpha;
        //
        out(1,0) = NN*alpha - bb;
        out(1,1) = NN;
        out(1,2) = Aim + Are*alpha;
        //
        out(2,0) = -2*Are*(1 + 2*NN)/TT;
        out(2,1) = -2*Aim*(1 + 2*NN)/TT;
        out(2,2) = (-1 - 2*(pow(Aim, 2) + pow(Are, 2)))/TT;
        break;
      case 1:
        //   w.r.t. x(t-tau)
        out(0,0) = kappa*CS;
        out(0,1) = kappa*SN;
        out(0,2) = 0;
        //
        out(1,0) = -kappa*SN;
        out(1,1) = kappa*CS;
        out(1,2) = 0;
        //
        out(2,0) = 0;
        out(2,1) = 0;
        out(2,2) = 0;
        break;
      default:
        std::cout<< "not implemented nx=1, vx[0]="<<vx[0]<<"\n";
        break;
    }
  }
  if ((nx==0)&&(np==1))
  {
    switch (vp[0])
    {
      case 0:
        // w.r.t period = p1
        out(0,0) = 0;
        out(1,0) = 0;
        out(2,0) = 0;
        break;
      case 1:
        // w.r.t alpha = p1
        out(0,0) = -Aim*NN;
        out(1,0) = Are*NN;
        out(2,0) = 0;
        break;
      case 2:
        // w.r.t TT = p2
        out(0,0) = 0;
        out(1,0) = 0;
        out(2,0) = (Aim*Aim+Are*Are+NN+2*(Aim*Aim+Are*Are)*NN -PP)/(TT*TT);
        break;
      case 3:
        // w.r.t kappa = p3
        out(0,0) = AreD*CS + AimD*SN;
        out(1,0) = AimD*CS - AreD*SN;
        out(2,0) = 0;
        break;
      case 4:
        // w.r.t PP = p4
        out(0,0) = 0;
        out(1,0) = 0;
        out(2,0) = 1.0/TT;
        break;
      case 5:
        // w.r.t bb = p5
        out(0,0) = Aim + kappa*(AimD*CS - AreD*SN)*tau;
        out(1,0) = -Are - AreD*CS*kappa*tau - AimD*kappa*SN*tau;
        out(2,0) = 0;
        break;
      case 6:
        // w.r.t tau = p6
        out(0,0) = kappa*(bb + omega_0)*(AimD*CS - AreD*SN);
        out(1,0) = -(kappa*(bb + omega_0)*(AreD*CS + AimD*SN));
        out(2,0) = 0;
        break;
      case 7:
        // w.r.t omega_0 = p7
        out(0,0) = kappa*(AimD*CS - AreD*SN)*tau;
        out(1,0) = -(kappa*(AreD*CS + AimD*SN)*tau);
        out(2,0) = 0;
        break;
      default:
        std::cout<< "not implemented np=1, vp[0]="<<vp[0]<<"\n";
        break;
    }
  }
  if ((nx==2)&&(np==0))
  {
    const double vv0 = vv(0,vx[0],idx);
    const double vv1 = vv(1,vx[0],idx);
    const double vv2 = vv(2,vx[0],idx);
    switch (vx[0])
    {
      case 0:
        switch (vx[1])
        {
          case 0:
            //   x(t), x(t)
            out(0,0) = vv2;
            out(0,1) = -(alpha*vv2);
            out(0,2) = vv0 - alpha*vv1;
            //
            out(1,0) = alpha*vv2;
            out(1,1) = vv2;
            out(1,2) = alpha*vv0 + vv1;
            //
            out(2,0) = (-2*(vv0 + 2*NN*vv0 + 2*Are*vv2))/TT;
            out(2,1) = (-2*(vv1 + 2*NN*vv1 + 2*Aim*vv2))/TT;
            out(2,2) = (-4*(Are*vv0 + Aim*vv1))/TT;
            break;
          case 1:
            //   x(t), x(t-tau)
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          default:
            std::cout<< "not implemented nx=2, vx[0]=0, vx[1]="<<vx[1]<<"\n";
            throw(-1);
            break;
        }
        break;
      case 1:
        switch (vx[1])
        {
          case 0:
            //   x(t-tau), x(t)
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 1:
            //   x(t-tau), x(t-tau)
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          default:
            std::cout<< "not implemented nx=2, vx[0]=1, vx[1]="<<vx[1]<<"\n";
            break;
        }
        break;
      default:
        std::cout<< "not implemented nx=2, vx[0]="<<vx[1]<<"\n";
        break;
    }
  }
  if ((nx==1)&&(np==1))
  {
    switch (vx[0])
    {
      case 0:
        switch (vp[0])
        {
          case 0:
            // x(t), p1 = alpha
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;

          case 1:
            // x(t) , p1 = alpha
            out(0,0) = 0;
            out(0,1) = -NN;
            out(0,2) = -Aim;
            //
            out(1,0) = NN;
            out(1,1) = 0;
            out(1,2) = Are;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 2:
            // x(t) , p2 = TT
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = (2*Are*(1 + 2*NN))/pow(TT, 2);
            out(2,1) = (2*Aim*(1 + 2*NN))/pow(TT, 2);
            out(2,2) = (1 + 2*pow(Aim, 2) + 2*pow(Are, 2))/pow(TT, 2);
            break;
          case 3:
            // x(t) , p3 = kappa
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 4:
            // x(t) , p4 = PP
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 5:
            // x(t) , p5 = bb
            out(0,0) = 0;
            out(0,1) = 1;
            out(0,2) = 0;
            //
            out(1,0) = -1;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 6:
            // x(t) , p6 = tau
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 7:
            // x(t) , p7 = omega_0
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          default:
            std::cout<< "not implemented\n";
            break;
        }
        break;
      case 1:
        switch (vp[0])
        {
          case 0:
            // x(t-tau), p1 = period
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 1:
            // x(t-tau), p1 = alpha
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 2:
            // x(t-tau) , p2 = TT
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 3:
            // x(t-tau) , p3 = kappa
            out(0,0) = CS;
            out(0,1) = SN;
            out(0,2) = 0;
            //
            out(1,0) = -SN;
            out(1,1) = CS;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 4:
            // x(t-tau) , p4 = PP
            out(0,0) = 0;
            out(0,1) = 0;
            out(0,2) = 0;
            //
            out(1,0) = 0;
            out(1,1) = 0;
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 5:
           // x(t-tau) , p5 = bb
            out(0,0) = -(kappa*SN*tau);
            out(0,1) = CS*kappa*tau;
            out(0,2) = 0;
            //
            out(1,0) = -(CS*kappa*tau);
            out(1,1) = -(kappa*SN*tau);
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 6:
            // x(t-tau) , p6 = tau
            out(0,0) = -(kappa*(bb + omega_0)*SN);
            out(0,1) = CS*kappa*(bb + omega_0);
            out(0,2) = 0;
            //
            out(1,0) = -(CS*kappa*(bb + omega_0));
            out(1,1) = -(kappa*(bb + omega_0)*SN);
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          case 7:
            // x(t) , p7 = omega_0
            out(0,0) = -(kappa*SN*tau);
            out(0,1) = CS*kappa*tau;
            out(0,2) = 0;
            //
            out(1,0) = -(CS*kappa*tau);
            out(1,1) = -(kappa*SN*tau);
            out(1,2) = 0;
            //
            out(2,0) = 0;
            out(2,1) = 0;
            out(2,2) = 0;
            break;
          default:
            std::cout<< "not implemented\n";
            break;
        }
        break;
      default:
        std::cout<< "not implemented\n";
        break;
    }
  }
}
#undef out
}

void sys_stpar( Vector& par )
{
  // bb controls the solution
  const double bb      = -0.0002;//;-0.00483423*0.98;
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

// 	std::cout<<"kappa "<< kappa << "\n";
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
  const double SN = sin((bb + omega_0)*tau);
  const double CS = cos((bb + omega_0)*tau);

  out(0) = -sqrt( ((PP*alpha-bb)*CS + PP*SN ) / ((2*bb+alpha)*CS + SN) );
  out(1) = 0.0;
  out(2) = bb*CS/(alpha*CS+SN);
// 	std::cout<<"Are "<< out(0) << " Aim "<< out(1) << " N " << out(2) <<"\n";
}

} // extern "C"
