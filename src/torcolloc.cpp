// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "system.h"
#include "matrix.h"
#include "spmatrix.h"
#include "hypermatrix.h"
#include "plot.h"
#include <cmath>

// support routines

// it is used inly by poly_coeff_lgr
static inline void poly_linmul( Array1D<double>& pp, double aa, double bb )
{
	// if pp(pp.Size()-1) != 0.0 ) std::cout<<"Error in \"poly_linmul\"\n";
	// pp * ( aa + bb*x )
	if( pp(pp.Size()-1) != 0.0 ) std::cout<<"Warning: \"poly_linmul\" is truncating the highest order term!\n";
	for( int i = pp.Size()-1; i > 0; --i )
	{
		pp(i) = aa*pp(i) + bb*pp(i-1);
	}
	pp(0) = aa*pp(0);
}

static void poly_coeff_lgr( Array1D<double>& out, Array1D<double>& t, int i )
{
	out.Clear();
	out(0) = 1.0;
	for( int j = 0; j < t.Size(); j++ )
	{
		if( j != i )
		{
			poly_linmul( out, -t(j)/(t(i)-t(j)), 1.0/(t(i)-t(j)) );
		}
	}
}

static void poly_coeff_mul( Array1D<double>& out, Array1D<double>& in1, Array1D<double>& in2 )
{
	if( out.Size() < (in1.Size()-1)+(in2.Size()-1)+1 ) std::cout<<"Error in \"poly_mul\"\n";
	out.Clear();
	for( int i=0; i<in1.Size(); i++ )
	{
		for( int j=0; j<in2.Size(); j++ )
		{
			out(i+j) = in1(i)*in2(j);
		}
	}
}

static void poly_coeff_int( Array1D<double>& out, Array1D<double>& in )
{
	out.Clear();
	out(0) = 0.0;
	for( int i=0; i<in.Size()-1; i++ )
	{
		out(i+1) += in(i)/(i+1);
	}
}

static void poly_coeff_diff( Array1D<double>& out, Array1D<double>& in )
{
	out.Clear();
	out(in.Size()-1) = 0.0;
	for( int i=1; i<in.Size(); i++ )
	{
		out(i-1) += in(i)*i;
	}
}

static double poly_eval( Array1D<double>& in, double c )
{
	double tmp = in(in.Size()-1);
	for( int i = in.Size()-1; i > 0; --i )
	{
		tmp = in(i-1) + c*tmp;
	}
	return tmp;
}

static void lobatto( Array1D<double>& mesh )
{
	const int N = mesh.Size()-1;
	if( mesh.Size() > 1 ) for( int i=0; i<mesh.Size(); i++ ) mesh(i) = 1.0/2.0*(1-cos((i*M_PI)/N));
}

static void gauss( Array1D<double>& mesh )
{
	const int N = mesh.Size();
	if( mesh.Size() > 1 ) for( int i=0; i<mesh.Size(); i++ ) mesh(i) = 1.0/2.0*(1-cos((2.0*i+1.0)/(2.0*N+2)*M_PI));
}

#include <functional>
#include <algorithm>

struct comp : public std::binary_function<int, int, bool> 
{
	const int* V;
	comp(int *V_) : V(V_) {} 
	bool operator()(int i, int j) { return V[i] < V[j]; }
};

// kk[.]     : i1,i2,j1,j2 -> idx2
// kk[ee[.]] :      s      -> idx2
// ee[.]     : i1,i2,j1,j2 ->  s
// rr[.]     : i1,i2,j1,j2 ->  index of Ai 

// inline int idxmap( int j1, int j2, int i1, int i2 )
// {
// 	if( (i2 >= NINT2)||(i1 >= NINT1) ) std::cout<<"IDER ";
// 	if( (j2 > NDEG2)||(j1 > NDEG1) ){ std::cout<<"IDR "; PDError(-1); }
// 	
// 	return j1 + i1*NDEG1 + (NDEG1*NINT1+1)*(j2 + i2*NDEG1);
// }

//-----------------------------------------------------------------------------------------------
//
//    CollocTR CLASS
//
//-----------------------------------------------------------------------------------------------

#define NTAU sys->ntau()
#define NDIM sys->ndim()
#define NPAR sys->npar()

#define RHO (NPAR+ParRot)

#include "torcolloc.h"
#include "pointtype.h"

CollocTR::CollocTR( System& sys_, int ndeg1_, int ndeg2_, int nint1_, int nint2_ ) :
	sys( &sys_ ),
	ndeg1(ndeg1_), ndeg2(ndeg2_), nint1(nint1_), nint2(nint2_),
	col1(ndeg1_),    col2(ndeg2_),
	mesh1(ndeg1_+1), mesh2(ndeg2_+1),
	lgr1( ndeg1_+1, ndeg1_+1 ), lgr2( ndeg2_+1, ndeg2_+1 ),
	dlg1( ndeg1_+1, ndeg1_+1 ), dlg2( ndeg2_+1, ndeg2_+1 ),
	I1( (ndeg1_+1), (ndeg1_+1) ),
	ID1( (ndeg1_+1), (ndeg1_+1) ),
	I2( (ndeg2_+1), (ndeg2_+1) ),
	ID2( (ndeg2_+1), (ndeg2_+1) ),
	mlg1( (ndeg1_+1)*(ndeg1_+1) ),
	mlg2( (ndeg2_+1)*(ndeg2_+1) ),
	mlgd1( (ndeg1_+1)*(ndeg1_+1) ),
	mlgd2( (ndeg2_+1)*(ndeg2_+1) ),
	ilg1( (ndeg1_+1)*(ndeg1_+1)+1 ),
	ilg2( (ndeg2_+1)*(ndeg2_+1)+1 ),
	ilgd1( (ndeg1_+1)*(ndeg1_+1)+1 ),
	ilgd2( (ndeg2_+1)*(ndeg2_+1)+1 )
{
	lobatto( mesh1 );
	lobatto( mesh2 );
	gauss( col1 );
	gauss( col2 );
// 	col1.Print();
// 	col2.Print();
	for( int i=0; i<mesh1.Size(); i++ ){ poly_coeff_lgr( lgr1(i), mesh1, i ); poly_coeff_diff( dlg1(i), lgr1(i) ); }
	for( int i=0; i<mesh2.Size(); i++ ){ poly_coeff_lgr( lgr2(i), mesh2, i ); poly_coeff_diff( dlg2(i), lgr2(i) ); }

	// here comes the phase condition for par(OMEGA0) and par(OMEGA1) 
	// the integration in the bottom border
	// construct the diffint matrix
	for( int i = 0; i < ndeg1+1; i++ )
	{
		for( int k = 0; k < ndeg1+1; k++ )
		{
			mlg1.Clear(); mlgd1.Clear(); ilg1.Clear(); ilgd1.Clear();
			poly_coeff_mul( mlg1,  lgr1(i), lgr1(k) ); 
			poly_coeff_mul( mlgd1, dlg1(i), lgr1(k) );
			poly_coeff_int( ilg1,  mlg1 );
			poly_coeff_int( ilgd1,  mlgd1 );
			I1( i, k ) = poly_eval( ilg1, 1.0 ) / nint1;
			ID1( i, k ) = poly_eval( ilgd1, 1.0 );
		}
	}
	for( int j = 0; j < ndeg2+1; j++ )
	{
		for( int l = 0; l < ndeg2+1; l++ )
		{
			mlg2.Clear(); mlgd2.Clear(); ilg2.Clear(); ilgd2.Clear();
			poly_coeff_mul( mlg2,  lgr2(j), lgr2(l) ); 
			poly_coeff_mul( mlgd2, dlg2(j), lgr2(l) );
			poly_coeff_int( ilg2,  mlg2 );
			poly_coeff_int( ilgd2,  mlgd2 );
			I2( j, l ) = poly_eval( ilg2, 1.0 ) / nint2;
			ID2( j, l ) = poly_eval( ilgd2, 1.0 );
		}
	}
}



void CollocTR::index( int* kk, int* ee, int* rr, Vector& par, Vector& tau, double* t1, double* t2 )
{
	sys->tau( tau, t1[0], par );
	
	for( int k=0; k<NTAU; k++ )
	{
		t1[1+k] = (t1[0] - tau(k)/par(0)) - floor(t1[0] - tau(k)/par(0));
		t2[1+k] = (t2[0] - par(RHO)*tau(k)/par(0)) - floor(t2[0] - par(RHO)*tau(k)/par(0));
	}
	
	for( int k=0; k<NTAU+1; k++ )
	{
		int i1 = static_cast<int>(floor(nint1*t1[k]));
		int i2 = static_cast<int>(floor(nint2*t2[k]));
		// some error checking
		if( (t1[k] > 1.0)||(t2[k] > 1.0)||
		    (t1[k] < 0.0)||(t2[k] < 0.0) ) std::cout<<"Er ";
		if( (i1>=nint1)||(i2>=nint2) ) std::cout<<"Ei ";
		// end error checking
		for( int j2=0; j2 < ndeg2+1; j2++ )
		{
			for( int j1=0; j1 < ndeg1+1; j1++ )
			{
				const int idx = idxkk( j1, j2, k );
				kk[idx] = idxmap( j1, j2, i1, i2 );
				ee[idx] = idx;
			}
		}
	}
	
	// sorting
	comp aa(kk);
	std::sort(ee, ee + (NTAU+1)*(ndeg1+1)*(ndeg2+1), aa );
	// with the quadratic algorithm
// 	for( int i=1; i<(NTAU+1)*(NDEG1+1)*(NDEG2+1); i++ )
// 	{
// 		for( int j=1; j<(NTAU+1)*(NDEG1+1)*(NDEG2+1); j++ )
// 		{
// 			if( kk[ee[j-1]] > kk[ee[j]] )
// 			{
// 				double tmp = ee[j-1];
// 				ee[j-1] = ee[j];
// 				ee[j] = tmp;
// 			}
// 		}
// 	}
	
	// filtering same indices
	rr[ee[0]] = 0;
	int nz = 0;
	for( int i=1; i < (NTAU+1)*(ndeg1+1)*(ndeg2+1); i++ )
	{
		if( kk[ee[i-1]] != kk[ee[i]] ) nz++;
		rr[ee[i]] = nz;
	}
// 	std::cout<<"BEGIN\n";
// 	for( int i=0; i < (NTAU+1)*(NDEG1+1)*(NDEG2+1); i++ )
// 	{
// 		std::cout<<" k "<<kk[ee[i]]<<" e "<<ee[ee[i]]<<" r "<<rr[ee[i]];
// 	}
// 	std::cout<<"\nEND\n";
}

void CollocTR::indexSep( int* kk, int* ee, int* rr, Vector& par, Vector& tau, double* t1, double* t2 )
{
	sys->tau( tau, t1[0], par );
	
	for( int k=0; k<NTAU; k++ )
	{
		t1[1+k] = (t1[0] - tau(k)/par(0));
		t2[1+k] = (t2[0] - par(RHO)*tau(k)/par(0)) - floor(t2[0] - par(RHO)*tau(k)/par(0));
	}
	
	for( int k=0; k<NTAU+1; k++ )
	{
		int i1 = static_cast<int>(floor(nint1*t1[k]));
		int i2 = static_cast<int>(floor(nint2*t2[k]));
		// some error checking
		if( (t1[k] > 1.0)||(t2[k] > 1.0)||
		    (t1[k] < -1.0)||(t2[k] < 0.0) ) std::cout<<"Er ";
		if( (i1>=nint1)||(i2>=nint2) ) std::cout<<"Ei ";
		// end error checking
		for( int j2=0; j2 < ndeg2+1; j2++ )
		{
			for( int j1=0; j1 < ndeg1+1; j1++ )
			{
				const int idx = idxkk( j1, j2, k );
				kk[idx] = idxmapSep( j1, j2, i1, i2 );
				ee[idx] = idx;
			}
		}
	}
	
	// sorting
	comp aa(kk);
	std::sort(ee, ee + (NTAU+1)*(ndeg1+1)*(ndeg2+1), aa );
	// with the quadratic algorithm
// 	for( int i=1; i<(NTAU+1)*(NDEG1+1)*(NDEG2+1); i++ )
// 	{
// 		for( int j=1; j<(NTAU+1)*(NDEG1+1)*(NDEG2+1); j++ )
// 		{
// 			if( kk[ee[j-1]] > kk[ee[j]] )
// 			{
// 				double tmp = ee[j-1];
// 				ee[j-1] = ee[j];
// 				ee[j] = tmp;
// 			}
// 		}
// 	}
	
	// filtering same indices
	rr[ee[0]] = 0;
	int nz = 0;
	for( int i=1; i < (NTAU+1)*(ndeg1+1)*(ndeg2+1); i++ )
	{
		if( kk[ee[i-1]] != kk[ee[i]] ) nz++;
		rr[ee[i]] = nz;
	}
// 	std::cout<<"BEGIN\n";
// 	for( int i=0; i < (NTAU+1)*(NDEG1+1)*(NDEG2+1); i++ )
// 	{
// 		std::cout<<" k "<<kk[ee[i]]<<" e "<<ee[ee[i]]<<" r "<<rr[ee[i]];
// 	}
// 	std::cout<<"\nEND\n";
}

//seems to be right
void CollocTR::interpolate( Matrix& xx, Array1D<double>& sol, int* kk, double* t1, double* t2 )
{
	xx.Clear();
	for( int k = 1; k < NTAU+1; k++ )
	{
		const double c1 = nint1*t1[k] - floor(nint1*t1[k]);
		const double c2 = nint2*t2[k] - floor(nint2*t2[k]);
		for( int j2 = 0; j2 < ndeg2+1; j2++ )
		{
			for( int j1 = 0; j1 < ndeg1+1; j1++ )
			{
				const int idxK = idxkk( j1, j2, k );
				const double cf = poly_eval( lgr1(j1), c1 ) * poly_eval( lgr2(j2), c2 );
				for( int p = 0; p < NDIM; p++ )
				{
					xx(p, k-1) += cf*sol(p+NDIM*kk[idxK]);
				}
			}
		}
	}
	// derivatives at all delays
	for( int k = 0; k < NTAU+1; k++ )
	{
		const double c1 = nint1*t1[k] - floor(nint1*t1[k]);
		const double c2 = nint2*t2[k] - floor(nint2*t2[k]);
		for( int j2 = 0; j2 < ndeg2+1; j2++ )
		{
			for( int j1 = 0; j1 < ndeg1+1; j1++ )
			{
				const int idx = idxkk( j1, j2, k );
				const double cf1 = poly_eval( dlg1(j1), c1 ) * nint1 * poly_eval( lgr2(j2), c2 );
				const double cf2 = poly_eval( lgr1(j1), c1 ) * poly_eval( dlg2(j2), c2 ) * nint2;
				for( int p = 0; p < NDIM; p++ )
				{
					xx(p, NTAU + 2*k )   += cf1 * sol(p+NDIM*kk[idx]);
					xx(p, NTAU + 2*k+1 ) += cf2 * sol(p+NDIM*kk[idx]);
				}
			}
		}
	}
}

// this is the main attraction

//***************************************************************************
//
// NOT Just for periodic equations
//
//  D[u[x,y],x] + \rho*D[u[x,y],y] = T * f[t, u(x-tau(T)/T, y-rho*tau(T)/T), ... ]
//
//  ParRot is RHO
//
//***************************************************************************

void CollocTR::Jacobian( SpMatrix& A, Array1D< Vector* > Avar, Vector& rhs, Vector& par, Vector& sol, Array1D<int>& var )
{
	Vector tau( NTAU );
	Vector dtau( NTAU );
	Matrix dfx(NDIM,NDIM);
	Matrix dfp(NDIM,1);
	Vector fx(NDIM);
	Matrix dummy(0,0);
	Matrix xx(NDIM,NTAU+2*(NTAU+1));
	
	double* t1 = new double[NTAU+1];
	double* t2 = new double[NTAU+1];
	int*    kk = new int[(NTAU+1)*(ndeg1+1)*(ndeg2+1)];
	int*    ee = new int[(NTAU+1)*(ndeg1+1)*(ndeg2+1)];
	int*    rr = new int[(NTAU+1)*(ndeg1+1)*(ndeg2+1)];
	
	A.Clear();
	for( int r=0; r<var.Size(); r++ ) Avar(r)->Clear();
	// rhs doesn't need to be cleared

	// fill the matrix
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1; j1++ )
				{
					const int idx1 = idxmap( j1, j2, i1, i2 );
					if( idx1 >= (ndeg1*nint1)*(ndeg2*nint2) ) std::cout<<"IDXOVERFLOW\n";
					
					t1[0] = ((double)i1 + col1(j1))/nint1;
					t2[0] = ((double)i2 + col2(j2))/nint2;
					
					index( kk, ee, rr, par, tau, t1, t2 );
					const int lend = NDIM*(rr[ee[(NTAU+1)*(ndeg1+1)*(ndeg2+1)-1]] + 1);

					for( int p = 0; p < NDIM; p++ ) A.NewL( lend );
					
					interpolate( xx, sol, kk, t1, t2 );
					// righ-hand side
					sys->rhs( fx, t1[0], xx, par );
					for( int p=0; p < NDIM; p++ ) rhs( p+NDIM*idx1 ) = par(0)*fx(p) - xx(p,NTAU) - par(RHO)*xx(p,NTAU+1);
					
					// derivatives w.r.t the parameters
					for( int r=0; r<var.Size(); r++ )
					{
						Vector& deri = *Avar(r);
//!!!!!!!!!!!!!!!!!
//!BEGIN w.r.t the period
//!!!!!!!!!!!!!!!!!
						if( var(r) == 0 )
						{
							for( int p=0; p < NDIM; p++ ) deri( p+NDIM*idx1 ) = -fx(p);
							sys->dtau( dtau, t1[0], par, 0 );
							for( int k = 0; k < NTAU; k++ )
							{
								const double d = tau(k)/par(0) - dtau(k);
								if( d != 0 )
								{
									dfx.Clear();
									int nx = 1, vx = k, np = 0, vp = 0;
									sys->deri( dfx, t1[0], xx, par, nx, &vx, np, &vp, dummy );
									for( int p=0; p < NDIM; p++ )
									{
										for( int q=0; q < NDIM; q++ )
										{
											deri( p+NDIM*idx1 ) -= dfx(p,q) * d *(xx(q,NTAU+2*(k+1)) + par(RHO)*xx(q,NTAU+2*(k+1)+1) );
										}
									}
								}
							}
						}else
//!!!!!!!!!!!!!!!!!
//!END w.r.t the period
//!!!!!!!!!!!!!!!!!
//
//!!!!!!!!!!!!!!!!!
//!BEGIN w.r.t the ordinary parameters
//!!!!!!!!!!!!!!!!!
						// derivatives w.r.t. the real parameters
						// assuming that delays do not depend on parameters
						if( var(r) < NPAR )
						{
							dfp.Clear();
							int nx = 0, vx = 0, np = 1, vp = var(r);
							sys->deri( dfp, t1[0], xx, par, nx, &vx, np, &vp, dummy );
							for( int p=0; p < NDIM; p++ ) deri( p+NDIM*idx1 ) = - par(0) * dfp(p);
						}else
//!!!!!!!!!!!!!!!!!
//!END w.r.t the ordinary parameters
//!!!!!!!!!!!!!!!!!
//
//!!!!!!!!!!!!!!!!!
//!BEGIN w.r.t the rotation number
//!!!!!!!!!!!!!!!!!
						if( var(r) == RHO )
						{
							// the derivatives for the phase condition
							// it is much more difficult
							// the first frequency, which is the inverse of the period
							// std::cout<<"JAC:OM0\n";
							for( int p=0; p < NDIM; p++ ) deri(p + NDIM*idx1) = xx(p,NTAU+1);
							sys->dtau( dtau, t1[0], par, 0 );
							for( int k = 0; k < NTAU; k++ )
							{
								const double d = -tau(k);
								
								dfx.Clear();
								const int nx = 1, vx = k, np = 0, vp = 0;
								sys->deri( dfx, t1[0], xx, par, nx, &vx, np, &vp, dummy );
								for( int p=0; p < NDIM; p++ )
								{
									for( int q=0; q < NDIM; q++ )
									{
										deri(p + NDIM*idx1) -= dfx(p,q) * d * xx(q,NTAU+2*(k+1)+1);
									}
								}
							}
						}
//!!!!!!!!!!!!!!!!!
//!END w.r.t the rotation number
//!!!!!!!!!!!!!!!!!
						else std::cout<<"Jac:NNN\n";
					}
					// the matrix
					for( int k = 0; k < NTAU+1; k++ )
					{
						if( k != 0 )
						{
							// evaluate the solution
							dfx.Clear();
							const int nx = 1, vx = k-1, np = 0, vp = 0;
							sys->deri( dfx, t1[0], xx, par, nx, &vx, np, &vp, dummy );
						}
						
						const double c1 = nint1*t1[k] - floor(nint1*t1[k]);
						const double c2 = nint2*t2[k] - floor(nint2*t2[k]);
						
						// the real jacobian
						for( int l2 = 0; l2 < ndeg2+1; l2++ )
						{
							for( int l1 = 0; l1 < ndeg1+1; l1++ )
							{
								const int idxK = idxkk( l1, l2, k );
								if( k != 0 )
								{
									const double cf = poly_eval( lgr1(l1), c1 ) * poly_eval( lgr2(l2), c2 );
									for( int p = 0; p < NDIM; p++ )
									{
										for( int q = 0; q < NDIM; q++ )
										{
											// jacobian of rhs
											A.WrLi( p + NDIM*idx1, q+NDIM*rr[idxK] ) = q+NDIM*kk[idxK];
											A.WrLx( p + NDIM*idx1, q+NDIM*rr[idxK] ) -= cf*par(0)*dfx(p,q);
										}
									}
								}else{
									const double cf1 = poly_eval( dlg1(l1), c1 ) * nint1 * poly_eval( lgr2(l2), c2 );
									const double cf2 = poly_eval( lgr1(l1), c1 ) * poly_eval( dlg2(l2), c2 ) * nint2;
									for( int p = 0; p < NDIM; p++ )
									{
										for( int q = 0; q < NDIM; q++ )
										{
											// derivative part of the jacobian
											A.WrLi( p + NDIM*idx1, q+NDIM*rr[idxK] ) = q+NDIM*kk[idxK];
											if( p == q ) 
												A.WrLx( p + NDIM*idx1, q+NDIM*rr[idxK] ) += cf1 + par(RHO)*cf2;
										}
									}
								}
								// error check
								if( kk[idxK] > (ndeg1*nint1)*(ndeg2*nint2) ) std::cout<<"D"<<kk[idxK];
							}
						}
					}
				}
			} 
		}
	}
	delete[] rr;
	delete[] ee;
	delete[] kk;
	delete[] t2;
	delete[] t1;
}

//***************************************************************************************
//
// End of the periodic case
//
//***************************************************************************************


void CollocTR::PhaseONE( Vector& ph, Vector& presol )
{
	ph.Clear();
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2+1; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1+1; j1++ )
				{
					int idx1 = idxmap( j1, j2, i1, i2 );
					// matrix multiplication 
					for( int l2 = 0; l2 < ndeg2+1; l2++ )
					{
						for( int l1 = 0; l1 < ndeg1+1; l1++ )
						{
							const int idx2 = idxmap( l1, l2, i1, i2 );
							for( int p = 0; p < NDIM; p++ )
							{
								ph( p+NDIM*idx1 ) += presol( p+NDIM*idx2 ) * I1( l1, j1 ) * ID2( l2, j2 );
							}
						}
					}
				}
			}
		}
	}
}


void CollocTR::PhaseBOTH( Vector& ph0, Vector& ph1, Vector& presol )
{
	ph0.Clear();
	ph1.Clear();
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2+1; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1+1; j1++ )
				{
					int idx1 = idxmap( j1, j2, i1, i2 );
					// matrix multiplication 
					for( int l2 = 0; l2 < ndeg2+1; l2++ )
					{
						for( int l1 = 0; l1 < ndeg1+1; l1++ )
						{
							const int idx2 = idxmap( l1, l2, i1, i2 );
							for( int p = 0; p < NDIM; p++ )
							{
								ph0( p+NDIM*idx1 ) += presol( p+NDIM*idx2 ) * ID1( l1, j1 ) * I2( l2, j2 );
								ph1( p+NDIM*idx1 ) += presol( p+NDIM*idx2 ) * I1( l1, j1 ) * ID2( l2, j2 );
							}
						}
					}
				}
			}
		}
	}
}


double CollocTR::Integrate( Vector& ph1, Vector& ph2 )
{
	double res = 0.0;
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2+1; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1+1; j1++ )
				{
					int idx1 = idxmap( j1, j2, i1, i2 );
					// matrix multiplication 
					for( int l2 = 0; l2 < ndeg2+1; l2++ )
					{
						for( int l1 = 0; l1 < ndeg1+1; l1++ )
						{
							const int idx2 = idxmap( l1, l2, i1, i2 );
							for( int p = 0; p < NDIM; p++ )
							{
								res += ph1( p+NDIM*idx1 ) * I1( j1, l1 ) * I2( j2, l2 ) * ph2( p+NDIM*idx2 );
							}
						}
					}
				}
			}
		}
	}
	return res;
}

void CollocTR::Star( Vector& ph1, Vector& ph2 )
{
	ph1.Clear();
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2+1; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1+1; j1++ )
				{
					int idx1 = idxmap( j1, j2, i1, i2 );
					// matrix multiplication 
					for( int l2 = 0; l2 < ndeg2+1; l2++ )
					{
						for( int l1 = 0; l1 < ndeg1+1; l1++ )
						{
							const int idx2 = idxmap( l1, l2, i1, i2 );
							for( int p = 0; p < NDIM; p++ )
							{
								ph1( p+NDIM*idx2 ) += I1( l1, j1 ) * I2( l2, j2 ) * ph2( p+NDIM*idx1 );
							}
						}
					}
				}
			}
		}
	}
}


double CollocTR::IntegrateDIFF( Vector& ph1, Vector& ph2, Vector& ph3 )
{
	double res = 0.0;
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2+1; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1+1; j1++ )
				{
					int idx1 = idxmap( j1, j2, i1, i2 );
					// matrix multiplication 
					for( int l2 = 0; l2 < ndeg2+1; l2++ )
					{
						for( int l1 = 0; l1 < ndeg1+1; l1++ )
						{
							const int idx2 = idxmap( l1, l2, i1, i2 );
							for( int p = 0; p < NDIM; p++ )
							{
								res += (ph1( p+NDIM*idx1 )-ph2( p+NDIM*idx1 )) * I1( j1, l1 ) * I2( j2, l2 ) * ph3( p+NDIM*idx2 );
							}
						}
					}
				}
			}
		}
	}
	return res;
}

void CollocTR::ImportSol( Vector& out, Vector& in )
{
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1; j1++ )
				{
					const int idx1 = idxmap( j1, j2, i1, i2 );
					for( int p=0; p < NDIM; p++ )
					{
						out( p + NDIM*idx1 ) = in( p + NDIM*(j1 + i1*ndeg1) );
					}
				}
			}
		}
	}
}

void CollocTR::ImportTan( Vector& out, Vector& Re, Vector& Im, double alpha )
{
	// alpha /= 2.0;
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int i1 = 0; i1 < nint1; i1++ )
		{
			for( int j2 = 0; j2 < ndeg2; j2++ )
			{
				for( int j1 = 0; j1 < ndeg1; j1++ )
				{
					const int idx1 = idxmap( j1, j2, i1, i2 );
					const double t1 = (mesh1(j1) + i1)/((double)nint1);
					const double t2 = (mesh2(j2) + i2)/((double)nint2);
					for( int p=0; p < NDIM; p++ )
					{
						out( p + NDIM*idx1 ) = cos(-alpha*t1+2*M_PI*t2)*Re( p + NDIM*(j1 + i1*ndeg1) ) + sin(-alpha*t1+2*M_PI*t2)*Im( p + NDIM*(j1 + i1*ndeg1) );
					}
				}
			}
		}
	}
// 	std::cout<<"ImportTan : NORMing\n";
	double norm = sqrt( Integrate( out, out ) );
// 	std::cout<<"NORMing2\n";
	out /= norm;
}

#include <fstream>

void CollocTR::Save( char* dat, char* idx, const Vector& in )
{
	// writing to file for plotting
	std::ofstream ff(dat);
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int j2 = 0; j2 < ndeg2; j2++ )
		{
			for( int i1 = 0; i1 < nint1; i1++ )
			{
				for( int j1 = 0; j1 < ndeg1; j1++ )
				{
					const int idx1 = idxmap(j1, j2, i1, i2);
					for( int p = 0; p < NDIM; p++ )
					{
						ff<<in( p + NDIM*idx1 )<<"\t";
					}
				}
			}
			ff<<"\n";
		}
	}

	std::ofstream gg(idx);
	for( int i1 = 0; i1 < nint1; i1++ )
	{
		for( int j1 = 0; j1 < ndeg1; j1++ )
		{
			const double t1 = (mesh1(j1) + i1)/((double)nint1);
			gg<<t1<<"\t";
		}
	}
	gg<<"\n";
	for( int i2 = 0; i2 < nint2; i2++ )
	{
		for( int j2 = 0; j2 < ndeg2; j2++ )
		{
			const double t2 = (mesh2(j2) + i2)/((double)nint2);
			gg<<t2<<"\t";
		}
	}
	gg<<"\n";
}

#undef NTAU
#undef NDIM

