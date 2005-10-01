#include <cmath>
#include "pddesys.h"

// xx:  Ax Ay    ZZ  Xx Xy
//      0  1     2   3  4  5      6     7       8          9      10     11 12  13
// par: T  alpha eta TT Pc lambda sigma omega_m thetascale Escale Nscale bb tau omega_0

extern "C"
{

int sys_ndim(){ return 5; }
int sys_npar(){ return 14; }
int sys_ntau(){ return 2; }
int sys_nderi(){ return 0; }
 
void sys_tau( Vector& out, double t, const Vector& par )
{
	out(0) = 0.0;
	out(1) = par(12);
}

void sys_dtau( Vector& out, double t, const Vector& par, int vp )
{
	// it is a constant delay
	out(0) = 0.0;
	if( vp == 12 ) out(1) = 1.0;
	else out(1) = 0.0;	
}

void sys_rhs( Vector& out, double t, const Matrix& yy, const Vector& par )
{

#define xx(a,b) yy(a-1,b-1)

	// some calculations
const double phi=par(6) + par(8)*par(11)*par(12) + par(13)*par(12);


// the equations
#define tmp par(8)

out(0)=(tmp*(par(9)*par(11)*xx(2,1)
                +par(9)*par(10)*xx(1,1)*xx(3,1)
                -par(1)*par(9)*par(10)*xx(2,1)*xx(3,1)
               +par(2)*par(9)*xx(4,1)))/par(9);

out(1)=(tmp*(-(par(9)*par(11)*xx(1,1))
                +par(1)*par(9)*par(10)*xx(1,1)*xx(3,1)
                +par(9)*par(10)*xx(2,1)*xx(3,1)
                +par(2)*par(9)*xx(5,1)))/par(9);

out(2)=(tmp*(par(4)
                -par(10)*xx(3,1)
                -(pow(par(9)*xx(1,1),2.0) + pow(par(9)*xx(2,1),2.0))*(1 + 2*par(10)*xx(3,1))))/(par(3)*par(10));

out(3)=(tmp*(cos(phi)*par(5)*par(9)*xx(1,2)
                +par(5)*par(9)*sin(phi)*xx(2,2)
                -par(5)*par(9)*xx(4,1)
                -(par(7)-par(13))*par(9)*xx(5,1)
                +par(9)*par(11)*xx(5,1)))/par(9);

out(4)=(tmp*(-(par(5)*par(9)*sin(phi)*xx(1,2))
                +cos(phi)*par(5)*par(9)*xx(2,2)
                +(par(7)-par(13))*par(9)*xx(4,1)
                -par(9)*par(11)*xx(4,1)
                -par(5)*par(9)*xx(5,1)))/par(9);
}

void sys_deri( Matrix &out, double t, const Matrix& xx, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{

}

void sys_stpar( Vector& par )
{
	// par(0) = 16.0;
}

void sys_stsol( Vector& out, double t )
{
	// out(0)   = 0.9+0.3*sin(2*M_PI*t);
}

} // extern "C"
