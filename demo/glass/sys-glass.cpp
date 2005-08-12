#include <cmath>
#include "pddesys.h"

int Sys::ndim(){ return 1; }
int Sys::npar(){ return 4; }
int Sys::ntau(){ return 2; }
int Sys::nderi(){ return 0; }
 
void Sys::tau( Vector& out, double t, const Vector& par )
{
	out(0) = 0.0; out(1) = par(3);
}

void Sys::dtau( Vector& out, double t, const Vector& par, int vp )
{
	switch( vp ) {
		case 0:
			out(0) = 0.0; out(1) = 0.0; break;
		case 1:
			out(0) = 0.0; out(1) = 0.0; break;
		case 2:
			out(0) = 0.0; out(1) = 0.0; break;
		case 3:
			out(0) = 0.0; out(1) = 1.0; break;
	}
}

void Sys::rhs( Vector& out, double t, const Matrix& x, const Vector& par )
{
	out(0) = par(1)*x(0,0) + par(2)*x(0,1)/(1+pow(x(0,1), 10.0));
}

void Sys::deri( Matrix &out, double t, const Matrix& x, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{
	// derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
	if( (nx == 1) && (np == 0) ) {
		switch( vx[0] ) {
			case 0:
				out(0,0) = par(1);
				break;
			case 1:
				out(0,0) = par(2)*(1.0 - 9.0*pow(x(0,1), 10.0) )/( (1+pow(x(0,1), 10.0))*(1+pow(x(0,1), 10.0)) );
				break;
		}
	}
	// derivatives w.r.t. the parameters, purely, so this results a vector
	if( (nx == 0) && (np == 1) )
	{
		switch( vp[0] )
		{
			case 1:
				out(0) = x(0,0); break;
			case 2:
				out(0) = x(0,1)/(1+pow(x(0,1), 10)); break;
			case 3:
				out(0) = 0.0; break;
		}
	}
	// second derivatives w.r.t. x
	if( (nx == 2) && (np == 0) ) {
		switch( vx[0] ) {
			case 0:
				switch( vx[1] ) {
					case 0:
						out(0,0) = 0.0; break;
					case 1:
						out(0,0) = 0.0; break;
				} break;
			case 1:
				switch( vx[1] ) {
					case 0:
						out(0,0) = 0.0; break;
					case 1:
						out(0,0) = 10.0*par(2)*pow(x(0,1), 9.0)*(9.0*pow(x(0,1), 10.0) - 11.0)/
						                       pow((1.0+pow(x(0,1), 10.0)), 3.0) * vv(0,1); break;
				} break;
		}
	}
	// mixed derivative w.r.t. x and the parameters
	if( (nx == 1) && (np == 1) ) {
		switch( vp[0] ) {
			case 1:
				switch( vx[0] ) {
					case 0:
						out(0,0) = 1.0; break;
					case 1:
						out(0,0) = 0.0; break;
				} break;
			case 2:
				switch( vx[0] ) {
					case 0:
						out(0,0) = 0.0; break;
					case 1:
						out(0,0) = (1.0 - 9.0*pow(x(0,1), 10.0) )/
						           ((1.0 + pow(x(0,1), 10.0))*(1.0+pow(x(0,1), 10.0)) ); break;
				} break;
			case 3:
				switch( vx[0] ) {
					case 0:
						out(0,0) = 0.0; break;
					case 1:
						out(0,0) = 0.0; break;
				} break;
		}
	}
}

void Sys::stpar( Vector& par )
{
	par(0) = 5.3918; par(1) = -1.0; par(2) = 1.74; par(3) = 2.0;
}

void Sys::stsol( Vector& out, double t )
{
	out(0)   = 0.9+0.3*sin(2*M_PI*t);
}
