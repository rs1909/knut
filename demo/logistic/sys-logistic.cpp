/***************************************************************************
 *   Copyright (C) 2004 by Robert Szalai                                   *
 *   szalai@localhost.localdomain                                          *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <cmath>
#include "pddesys.h"

int Sys::ndim(){ return 1; }
int Sys::npar(){ return 2; }
int Sys::ntau(){ return 2; }

void Sys::tau( Vector& out, double t, const Vector& par )
{
  out(0) = 0.0;
  out(1) = 1.0;
}

void Sys::dtau( Vector& out, double t, const Vector& par, int vp )
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
			cout << "dtau: not implemented\n";
			break;
	}
}

void Sys::rhs( Vector& out, double t, const Matrix& x, const Vector& par )
{
  out(0) = ( par(1)-x(0,1) )*x(0,0);
}

void Sys::deri( Matrix &out, double t, const Matrix& x, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{
  if( (t < 0)||(t > 1) ){ cout << "deri: t is not element of the interval\n"; }
  
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
      out(0) = x(0,0);
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
						cout << "deri: not implemented\n";
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
						out(0,0) = 1.0;
						break;
					case 1:
						out(0,0) = 0.0;
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
