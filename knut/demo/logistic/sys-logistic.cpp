/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <cmath>
#include "knutsys.h"
extern "C"
{

size_t sys_ndim(){ return 1; }
size_t sys_npar(){ return 2; }
size_t sys_ntau(){ return 2; }

void sys_tau( KNVector& out, double t, const KNVector& par )
{
  out(0) = 0.0;
  out(1) = 1.0;
}

void sys_dtau( KNVector& out, double t, const KNVector& par, size_t vp )
{
	switch( vp )
	{
		case 0:
			out(0) = 0.0;
			out(1) = 0.0;
			break;
		case 1:
			out(0) = 0.0;
			out(1) = 0.0;
			break;
		default:
			std::cout << "dtau: not implemented\n";
			break;
	}
}

void sys_rhs( KNVector& out, double t, const KNMatrix& x, const KNVector& par )
{
  out(0) = ( par(1)-x(0,1) )*x(0,0);
}

void sys_deri( KNMatrix &out, double t, const KNMatrix& x, const KNVector& par, 
	       size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNMatrix& vv )
{
  if( (t < 0)||(t > 1) ){ std::cout << "deri: t is not element of the interval\n"; }
  
  // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
  if( (nx == 1) && (np == 0) ){
    switch( vx[0] ){
    case 0:
      out(0,0) = par(1)-x(0,1);
      break;
    case 1:
      out(0,0) = - x(0,0);
      break;
    default:
      std::cout << "deri: not implemented\n";
      break;
    }
  }     
  // derivatives w.r.t. the parameters, purely, so this results a vector
  if( (nx == 0) && (np == 1) ){
    switch( vp[0] ){
    case 0: //  T period length
      std::cout << "deri w.r.t. T: not implemented\n";
      break;
    case 1:
      out(0) = x(0,0);
      break;
    default:
      std::cout << "deri: not implemented\n";
      break;
    }
  }
  // second derivatives wrt. xx
	if( (nx == 2) && (np == 0) )
	{
		switch( vx[0] )
		{
			case 0:
				switch( vx[1] )
				{
					case 0:
						out(0,0) = 0.0;
						break;
					case 1:
						out(0,0) = -1.0*vv(0,0);
						break;
					default:
						std::cout << "deri: not implemented\n";
						break;
				}
				break;
			case 1:
				switch( vx[1] )
				{
					case 0:
						out(0,0) = -1.0*vv(0,1);
						break;
					case 1:
						out(0,0) = 0.0;
						break;
					default:
						std::cout << "deri: not implemented\n";
						break;
				}
				break;
			default:
				std::cout << "deri: not implemented\n";
				break;
		}
	}
  // mixed derivative wrt to xx and par
  if( (nx == 1) && (np == 1) )
  {
		switch( vp[0] )
		{
			case 0: //  T period length
				std::cout << "deri w.r.t. T: not implemented\n";
				break;
			case 1:
				switch( vx[0] )
				{
					case 0:
						out(0,0) = 1.0;
						break;
					case 1:
						out(0,0) = 0.0;
						break;
					default:
						std::cout << "deri: not implemented\n";
						break;
				}	
				break;
			default:
				std::cout << "deri: not implemented\n";
				break;
		}
	}
}

} // extern "C"
