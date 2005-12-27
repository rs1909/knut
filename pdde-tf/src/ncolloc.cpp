// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "pderror.h"
#include "ncolloc.h"
#include "system.h"

#include "matrix.h"
#include <cmath>

// #define MADD

// support routines

static void equidist( Vector& mesh )
{
	const int m = mesh.Size();
	for( int i = 0; i < m; i++ ) mesh(i) = (double)i/((double)m-1);
}

static void poly_gau( Vector& roots )
{
	const int m=roots.Size();

	Matrix a( m, m );
	/* construct the matrix */
	switch( m )
	{
		case 1:
			a(0,0) = 0.0;
			break;
		case 2:
			a(0,0) = -1.0/3; a(0,1) = 0.0;
			break;
		case 3:
			a(0,0) = 0.0; a(0,1) = -3.0/5; a(0,2) = 0.0;
			break;
		case 4:
			a(0,0) = 3.0/35; a(0,1) = 0.0; a(0,2) = -6.0/7; a(0,3) = 0.0;
			break;
		case 5:
			a(0,0) = 0.0; a(0,1) = 5.0/21; a(0,2) = 0.0; a(0,3) = -10.0/9; a(0,4) = 0.0;
			break;
		case 6:
			a(0,0) = -5.0/231; a(0,1) = 0.0; a(0,2) = 5.0/11; a(0,3) = 0.0; 
			a(0,4) = -15.0/11; a(0,5) = 0.0;
			break;
		case 7:
			a(0,0) = 0.0; a(0,1) = -35.0/429; a(0,2) = 0.0; a(0,3) = 105.0/143; 
			a(0,4) = 0.0; a(0,5) = -21.0/13; a(0,6) = 0;
			break;
		default:
			std::cout << "Something wrong! \n";
			PDError(-1);
			return;
			break;
	}
	
	for(int i = 0; i < m; i++) a(m-1,i) = -a(0,i);
	for(int i = 0; i < m-1; i++)
	{
		for( int j = 0; j < m; j++)
		{
			if( i+1 == j ) a(i,j) = 1.0;
			else a(i,j) = 0.0;
		}
	}
	
	Vector wi(m);
	a.Eigval( roots, wi );
	
	/* scaling */
	for(int i = 0; i < m; i++) roots(i) = (roots(i)+1.0)/2.0;
	
	/* sorting */
	double tmp;
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < (m-1); j++)
		{
			if( roots(j) > roots(j+1) )
			{
				tmp = roots(j); roots(j) = roots(j+1); roots(j+1) = tmp;
			}
		}
	}
}

static void lobatto( Array1D<double>& mesh )
{
	const int N = mesh.Size()-1;
	if( mesh.Size() > 1 ) for( int i=0; i<mesh.Size(); i++ ) mesh(i) = 1.0/2.0*(1-cos((i*M_PI)/N));
	for( int i=0; i<mesh.Size(); i++ ) std::cout<<mesh(i)<<" ";
	std::cout<<"LOB\n";
}

static void gauss( Array1D<double>& mesh )
{
	const int N = mesh.Size();
	if( mesh.Size() > 1 ) for( int i=0; i<mesh.Size(); i++ ) mesh(i) = 1.0/2.0*(1-cos((2.0*i+1.0)/(2.0*N+2)*M_PI));
	for( int i=0; i<mesh.Size(); i++ ) std::cout<<mesh(i)<<" ";
	std::cout<<"GAU\n";
}

inline static void repr_mesh( Vector& V )
{
	equidist( V );
}

inline static void col_mesh( Vector& V )
{
	poly_gau( V );
}

static void poly_lgr( Vector& t, Vector &out, double c )
{

	if( t.Size() != out.Size() )
	{
		std::cout << "poly_lgr: wrong dimensions" << '\n';
		PDError(-1);
	}
	for( int i = 0; i < t.Size(); i++)
	{
		out(i) = 1.0;
		for( int j = 0; j < t.Size(); j++)
		{
			if( i != j )
			{
				out(i) *= ( c - t(j) )/( t(i) - t(j) );
			}
		}
	}
}

static void poly_dlg( Vector& t, Vector& out, double c)
{
	int j,k,l;
	double f;

	if( t.Size() != out.Size() )
	{
		std::cout << "poly_dlg: wrong dimensions" << '\n';
		PDError(-1);
	}

	for(j = 0; j < t.Size(); j++ )
	{
		out(j) = 0.0;
		for(k = 0; k < t.Size(); k++ )
		{
			if( k != j )
			{
				f = 1.0;
				for(l = 0; l < t.Size(); l++ )
				{
					if( (l != k)&&(l != j) ) f *= ( c - t(l) )/( t(j) - t(l) );
				}
				out(j) += f/( t(j) - t(k) );
			}
		}
	}
}

static inline void poly_mul( Vector& pp, double bb, double aa )
{
	// pp * ( aa + bb*x )
	Vector tmp(pp.Size());
	for( int i = 0; i < pp.Size(); i++ )
	{
		tmp(i) = pp(i)*bb;
		pp(i) *= aa;
	}
	for( int i = 1; i < pp.Size(); i++ ) pp(i) += tmp(i-1);
}

static void poly_int( Matrix& out, Vector& t )
{
	int i,j,k;
	Vector poly(2*t.Size());
	
	for( i = 0; i < t.Size(); i++ )
	{
		for( j = 0; j < t.Size(); j++ )
		{
			poly(0) = 1.0;
			for( k = 1; k < poly.Size(); k++ ) poly(k) = 0.0;
			//      poly.Print();
			// i,j az out matrix indexe
			//      cout<<"in:poly_mul\n";
			for( k = 0; k < t.Size(); k++ )
			{
				if( k != i ) 
					poly_mul( poly, 1.0/(t(i)-t(k)), -t(k)/(t(i)-t(k)) );
				if( k != j ) 
					poly_mul( poly, 1.0/(t(j)-t(k)), -t(k)/(t(j)-t(k)) );
			}
			//      cout<<"out:poly_mul\n";
			//      t.Print();
			//      poly.Print();
			// integrate
			for( k = 0; k < poly.Size(); k++ ) poly(k) /= k+1.0;
			out(i,j) = 0.0;
			// evaluate at x = 0..1
			for( k = 0; k < poly.Size(); k++ ) out(i,j) += poly(k);
		}
	}
} 

static void poly_diff_int( Matrix& out, Vector& t )
{
	Vector poly(2*t.Size());
	Vector poly_fin(2*t.Size());
	
	for( int i = 0; i < t.Size(); i++ )
	{
		for( int j = 0; j < t.Size(); j++ )
		{
			poly_fin.Clear();
			for( int s = 0; s < t.Size(); s++ )
			{
				if( s != i )
				{
					poly(0) = 1.0; 
					for( int r = 1; r < poly.Size(); r++ ) poly(r) = 0.0;
					
					for( int r = 0; r < t.Size(); r++ )
					{
						if( (i != r)&&(i != s) ) poly_mul( poly, 1.0/(t(i)-t(r)), -t(r)/(t(i)-t(r)) );
						if( j != r )             poly_mul( poly, 1.0/(t(j)-t(r)), -t(r)/(t(j)-t(r)) );
					}
					// adding
					for( int r = 0; r < poly.Size(); r++ ) poly_fin(r) += poly(r)/(t(i)-t(s));
				}
			}
			// integrate
			for( int k = 0; k < poly_fin.Size(); k++ ) poly_fin(k) /= k+1.0;
			out(i,j) = 0.0;
			// evaluate at x = 0..1
			for( int k = 0; k < poly_fin.Size(); k++ ) out(i,j) += poly_fin(k);
		}
	}
}

//
// the NColloc class
//

static inline int meshlookup( const Vector& mesh, double t )
{
	// binary search for in which interval is t-tau(k)
	int mid, low=0, up=mesh.Size()-1;
	while( up - low > 1 )
	{
		mid = low + (up - low)/2;
		if( (mesh(low) <= t) && (mesh(mid) > t) ) up = mid;
		else low = mid;
	}
	return low;
}

#define NDIM ndim
#define NPAR npar
#define NTAU ntau
#define NINT nint
#define NDEG ndeg
#define NMAT nmat

NColloc::NColloc( System& _sys, const int _nint, const int _ndeg, int _nmat )
	:
		ndim(_sys.ndim()), npar(_sys.npar()), ntau(_sys.ntau()),
		nint(_nint), ndeg(_ndeg), nmat(_nmat),
		mesh(nint+1), time(nint*ndeg),
		
		kk(ntau+1, nint*ndeg), ee(ntau+1, nint*ndeg),
		dd(ntau+1, nint*ndeg), rr(ntau+1, nint*ndeg),
		
		kkS(ntau+1, nint*ndeg), eeS(ntau+1, nint*ndeg),
		rrS(ntau+1, nint*ndeg), ddS(ntau+1, nint*ndeg),
		mmS(ntau+1, nint*ndeg), szS(2, nint*ndeg),
		
		kkI(ntau+1, nint*ndeg), eeI(ntau+1, nint*ndeg),
		rrI(ntau+1, nint*ndeg), ddI(ntau+1, nint*ndeg),
		mmI(ntau+1, nint*ndeg), szI(nmat+1, nint*ndeg),
		
		tt(2*ntau+1, ndeg+1, ndeg*nint),
		ttMSH( 2*ntau, ndeg+1, ndeg*nint + 1 ),
		kkMSH(ntau, nint*ndeg+1),
		metric(ndeg+1, ndeg+1),
		metricPhase(ndeg+1, ndeg+1),
		col(ndeg),
		out(ndeg+1),
		meshINT(ndeg+1)
{
	sys = &_sys;
	for( int i = 0; i < NINT+1; i++ ) mesh(i) = i*1.0/NINT;
	repr_mesh( meshINT );

	col_mesh( col );
	poly_int( metric, meshINT ); // works for any meshes
	poly_diff_int( metricPhase, meshINT );
}

void NColloc::Init( const Vector& par, const Vector& /*sol*/ )
{
	double *t = new double[NTAU+1];
	double *tMSH = new double[NTAU];
	Vector ttau(NTAU);

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			int idx = j+i*NDEG;
			const double h0 = mesh(i+1) - mesh(i);

			/* first, the derivative: x'(t)*/
			t[0] = mesh(i) + h0*col(j); // xdot(c_i,j)
			time(idx) = t[0];
			kk(0, idx) = i;
			kkS(0, idx) = i;
			kkI(0, idx) = i;
			poly_dlg( meshINT, out, col(j) );
			// out.Print(); std::cout<<"SS\n";
			for( int k = 0; k < NDEG+1; k++ )
			{
				tt( 0, k, idx ) = out(k)/h0;
			}

			/* second, the right-hand side */
			sys->tau( ttau, t[0], par );
			for( int k = 0; k < NTAU; k++ )
			{
				ttau(k) /= par(0);
				if( ttau(k) > NMAT ){ std::cout<<"\n NColloc::Init: DELAY > NMAT*PERIOD  "<<k<<"\n\n"; PDError(12); }
				if( ttau(k) < 0.0 ){ std::cout<<"\n NColloc::Init: Negative DELAY "<<k<<"\n\n"; PDError(12); }
				t[1+k] = (t[0] - ttau(k)) - floor(t[0] - ttau(k));  // nem szetvalasztott
				
				// binary search for in which interval is t-tau(k)
				const int low = meshlookup( mesh, t[1+k] );
				kk(k+1, idx) = low;

				if( t[0] - ttau(k) >= 0 ) kkS(k+1, idx) = low;
				else kkS(k+1, idx) = low - NINT;
				
				kkI(k+1, idx ) = low + NINT*static_cast<int>(floor(t[0] - ttau(k)));
				
				const double hk = mesh(low+1) - mesh(low);
				// x(t-\tau_i)
				poly_lgr( meshINT, out, (t[1+k] - mesh(low))/hk );
				for( int l = 0; l < NDEG+1; l++ ) 
				{
					tt( 1+k, l, idx ) = out(l);
				}
				// x'(t-\tau_i)
				poly_dlg( meshINT, out, (t[1+k] - mesh(low))/hk );
				for( int l = 0; l < NDEG+1; l++ ) 
				{
					tt( NTAU+1+k, l, idx ) = out(l)/hk;
				}
				// creating interpolation at the representation points
				tMSH[k] = mesh(i) + j*h0/(double)NDEG - ttau(k) - floor(mesh(i) + j*h0/(double)NDEG - ttau(k));
				const int lowMSH = meshlookup( mesh, tMSH[k] );
				kkMSH(k, idx) = lowMSH;
				const double hkMSH = mesh(lowMSH+1) - mesh(lowMSH);
				poly_lgr( meshINT, out, (tMSH[k] - mesh(lowMSH))/hkMSH );
				for( int l = 0; l < NDEG+1; l++ ) 
				{
					ttMSH( k, l, idx ) = out(l);
				}
				// x'(t-\tau_i)
				poly_dlg( meshINT, out, (tMSH[k] - mesh(lowMSH))/hkMSH );
				for( int l = 0; l < NDEG+1; l++ ) 
				{
					ttMSH( NTAU+k, l, idx ) = out(l)/hkMSH;
				}
				// REPR ENDS
			}

			// sorting kk

			for( int r=0; r < NTAU+1; r++ )
			{
				ee(r,idx) = r;
				eeS(r,idx) = r;
				eeI(r,idx) = r;
			}
			
			int tmp;
			// probably there is a better algorithm. It is n^2, but NTAU is generally small
			for( int r=1; r < NTAU+1; r++ )
			{
				for( int s=r; s < NTAU+1; s++ )
				{
					if( kk(ee(r-1,idx),idx) > kk(ee(s,idx),idx) )
					{
						tmp = ee(r-1,idx); ee(r-1,idx) = ee(s,idx); ee(s,idx) = tmp;
					}
					if( kkS(eeS(r-1,idx),idx) > kkS(eeS(s,idx),idx) )
					{
						tmp = eeS(r-1,idx); eeS(r-1,idx) = eeS(s,idx); eeS(s,idx) = tmp;
					}
					if( kkI(eeI(r-1,idx),idx) > kkI(eeI(s,idx),idx) )
					{
						tmp = eeI(r-1,idx); eeI(r-1,idx) = eeI(s,idx); eeI(s,idx) = tmp;
					}
				}
			}

			// filter out elements at the same position and find adjacent elements
			// for one matrix
			int idp = 0, del = 0;
			rr(0,idx) = 0; dd(0,idx) = 0;
			
			// for two separate matrices
			for( int r = 0; r < 2; r++ ) szS(r,idx) = 0;
			for( int r = 0; r < NTAU+1; r++ ) mmS(r,idx) = (-kkS(eeS(r,idx),idx) + NINT-1)/NINT;
			int idpS = 0, delS = 0;
			rrS(0,idx) = 0; ddS(0,idx) = 0;
			szS(mmS(0,idx),idx) = NDEG+1;
			
			// for many separate matrices
			for( int r = 0; r < NMAT; r++ ) szI(r,idx) = 0;
			for( int r = 0; r < NTAU+1; r++ ) mmI(r,idx) = (-kkI(eeI(r,idx),idx) + NINT-1)/NINT;
			int idpI = 0, delI = 0;
			rrI(0,idx) = 0; ddI(0,idx) = 0;
			szI(mmI(0,idx),idx) = NDEG+1;
			
			for( int r = 1; r < NTAU+1; r++ )
			{
				// for one matrix
				if( kk(ee(r-1,idx),idx) != kk(ee(r,idx),idx) ) idp++;
				rr(r,idx) = idp;
				if( kk(ee(r,idx),idx) - kk(ee(r-1,idx),idx) == 1 ) del++;
				dd(r,idx) = del;
				
				// for two separate matrices
				if( mmS(r-1,idx) != mmS(r,idx) )
				{
					idpS = 0, delS = 0;
					szS(mmS(r,idx),idx) = NDEG+1;
				}else
				{
					if( kkS(eeS(r-1,idx),idx) != kkS(eeS(r,idx),idx) )
					{
						idpS++;
						szS(mmS(r,idx),idx) += NDEG+1;
					}
					if( kkS(eeS(r,idx),idx) - kkS(eeS(r-1,idx),idx) == 1 )
					{
						delS++;
						szS(mmS(r,idx),idx) -= 1;
					}
				}
				rrS(r,idx) = idpS;
				ddS(r,idx) = delS;
				
				// for many separate matrices
				if( mmI(r-1,idx) != mmI(r,idx) )
				{
					idpI = 0, delI = 0;
					szI(mmI(r,idx),idx) = NDEG+1;
				}else
				{
					if( kkI(eeI(r-1,idx),idx) != kkI(eeI(r,idx),idx) )
					{
						idpI++;
						szI(mmI(r,idx),idx) += NDEG+1;
					}
					if( kkI(eeI(r,idx),idx) - kkI(eeI(r-1,idx),idx) == 1 )
					{
						delI++;
						szI(mmI(r,idx),idx) -= 1;
					}
				}
				rrI(r,idx) = idpI;
				ddI(r,idx) = delI;
			}
		}
	}
	// making periodic
	for( int k = 0; k < NTAU; k++ )
	{
		kkMSH(k, NDEG*NINT) = kkMSH(k, 0);
		for( int l = 0; l < NDEG+1; l++ ) 
		{
			ttMSH( k, l, NDEG*NINT ) = ttMSH( k, l, 0 );
			ttMSH( NTAU+k, l, NDEG*NINT ) = ttMSH( NTAU+k, l, 0 );
		}
	}
	delete[] tMSH;
	delete[] t;
}

void NColloc::getMesh( Vector& msh )
{
	if( msh.Size() != NDEG*NINT+1 ) { std::cout<<"Error in NColloc::getMesh : Bad dimensions\n"; PDError(-1); }
	for( int i = 0; i < NINT; i++ )
		for( int j = 0; j < NDEG; j++ )
			msh( j + i*NDEG ) = mesh(i) + meshINT(j)*(mesh(i+1)-mesh(i));
	msh( NDEG*NINT ) = 1.0;
}

void NColloc::setMesh( const Vector& msh )
{
	if( msh.Size() != NDEG*NINT+1 ) { std::cout<<"Error in NColloc::setMesh : Bad dimensions\n"; PDError(-1); }
	for( int i = 0; i < NINT; i++ ) mesh(i) = msh( i*NDEG );
	mesh( NINT ) = 1.0;
}

void NColloc::Interpolate( JagMatrix3D& jag, const Vector& sol )
{
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int p = 0; p < NTAU; p++ )
			{
				for( int k = 0; k < NDIM; k++ )
				{
					jag(idx)( k, p ) = 0.0;
					for( int l = 0; l < NDEG+1; l++ )
					{
						jag(idx)( k, p ) += sol( k + NDIM*(l+kk(p+1,idx)*NDEG) )*tt( p+1, l, idx );
					}
				}
			}
			// x'(t)
			for( int k = 0; k < NDIM; k++ )
			{
				jag(idx)( k, NTAU ) = 0.0;
				for( int l = 0; l < NDEG+1; l++ )
				{
					jag(idx)( k, NTAU ) += sol( k + NDIM*(l+kk(0,idx)*NDEG) )*tt( 0, l, idx );
				}
			}
			// x'(t-\tau_i)
			for( int p = 0; p < NTAU; p++ )
			{
				for( int k = 0; k < NDIM; k++ )
				{
					jag(idx)( k, NTAU+1+p ) = 0.0;
					for( int l = 0; l < NDEG+1; l++ )
					{
						jag(idx)( k, NTAU+1+p ) += sol( k + NDIM*(l+kk(p+1,idx)*NDEG) )*tt( NTAU+1+p, l, idx );
					}
				}
			}
		}
	}
}

void NColloc::InterpolateREAL( JagMatrix3D& jag, const Vector& sol )
{
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int p = 0; p < NTAU; p++ )
			{
				for( int k = 0; k < NDIM; k++ )
				{
					// std::cout<<"InterR: "<<idx<<", "<<out(idx).Row()<<", "<<out(idx).Col()<<", "<<k<<", "<<p<<"\n";
					jag(idx)( k, p ) = 0.0;
					for( int l = 0; l < NDEG+1; l++ )
					{
						jag(idx)( k, p ) += sol( k + NDIM*(l+kk(p+1,idx)*NDEG) )*tt( p+1, l, idx );
					}
				}
			}
			// x'(t)
			for( int k = 0; k < NDIM; k++ )
			{
				jag(idx)( k, NTAU ) = 0.0;
				for( int l = 0; l < NDEG+1; l++ )
				{
					jag(idx)( k, NTAU ) += sol( k + NDIM*(l+kk(0,idx)*NDEG) )*tt( 0, l, idx );
				}
			}
		}
	}
}

void NColloc::InterpolateMSH( JagMatrix3D& jag, const Vector& sol )
{
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			for( int p = 0; p < NTAU; p++ )
			{
				for( int k = 0; k < NDIM; k++ )
				{
					// std::cout<<"InterR: "<<idx<<", "<<out(idx).Row()<<", "<<out(idx).Col()<<", "<<k<<", "<<p<<"\n";
					jag(idx)( k, p )      = 0.0;
					jag(idx)( k, NTAU+p ) = 0.0;
					for( int l = 0; l < NDEG+1; l++ )
					{
						jag(idx)( k, p )      += sol( k + NDIM*(l+kkMSH(p,idx)*NDEG) )*ttMSH( p, l, idx );
						jag(idx)( k, NTAU+p ) += sol( k + NDIM*(l+kkMSH(p,idx)*NDEG) )*ttMSH( NTAU+p, l, idx );
					}
				}
			}
		}
	}
	for( int p = 0; p < NTAU; p++ )
	{
		for( int k = 0; k < NDIM; k++ )
		{
			// std::cout<<"InterR: "<<idx<<", "<<out(idx).Row()<<", "<<out(idx).Col()<<", "<<k<<", "<<p<<"\n";
			jag(NDEG*NINT)( k, p )        = 0.0;
			jag(NDEG*NINT)( k, NTAU + p ) = 0.0;
			for( int l = 0; l < NDEG+1; l++ )
			{
				jag(NDEG*NINT)( k, p )      += sol( k + NDIM*(l+kkMSH(p,NDEG*NINT)*NDEG) )*ttMSH( p, l, NDEG*NINT );
				jag(NDEG*NINT)( k, NTAU+p ) += sol( k + NDIM*(l+kkMSH(p,NDEG*NINT)*NDEG) )*ttMSH( NTAU+p, l, NDEG*NINT );
			}
		}
	}
}

// in complex form on 2*i th places are the reals and on 2*i+1 th places the imaginary parts

void NColloc::InterpolateCPLX( JagMatrix3D& jag, const Vector& sol ) 
{
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			const int idxRe = 2*(j+i*NDEG);
			const int idxIm = 2*(j+i*NDEG) + 1;

			for( int p = 0; p < NTAU; p++ )
			{
				for( int k = 0; k < NDIM; k++ )
				{
					jag(idxRe)( k, p ) = 0.0;
					jag(idxIm)( k, p ) = 0.0;
					for( int l = 0; l < NDEG+1; l++ )
					{
						jag(idxRe)( k, p ) += sol( 2*(k + NDIM*(l+kk(p+1,idx)*NDEG)) )*tt( p+1, l, idx );
						jag(idxIm)( k, p ) += sol( 2*(k + NDIM*(l+kk(p+1,idx)*NDEG)) + 1 )*tt( p+1, l, idx );
					}
				}
			}
			for( int k = 0; k < NDIM; k++ )
			{
				jag(idxRe)( k, NTAU ) = 0.0;
				jag(idxIm)( k, NTAU ) = 0.0;
				for( int l = 0; l < NDEG+1; l++ )
				{
					jag(idxRe)( k, NTAU ) += sol( 2*(k + NDIM*(l+kk(0,idx)*NDEG)) )*tt( 0, l, idx );
					jag(idxIm)( k, NTAU ) += sol( 2*(k + NDIM*(l+kk(0,idx)*NDEG)) + 1 )*tt( 0, l, idx );
				}
			}
		}
	}
}

void NColloc::RHS( Vector& rhs, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	Vector fx(NDIM); // it should be a variable in the class itself

	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		rhs(r) = sol(r+NDIM*NDEG*NINT) - sol(r);;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			sys->rhs( fx, time(idx), solData(idx), par );
			for( int k = 0; k < NDIM; k++ )
				rhs( NDIM + k + NDIM*(j+NDEG*i) ) = par(0)*fx(k) - solData(idx)(k, NTAU);
		}
	}
}

void NColloc::RHS_p( Vector& rhs, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData, int alpha )
{

	Vector tau(NTAU);
	Vector dtau(NTAU);
	Vector fx(NDIM);
	Matrix dfx(NDIM,NDIM);
	Matrix dfx2(NDIM,1);
	Matrix dummy(0,0);

	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		rhs(r) = 0.0;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			sys->tau( tau, time(idx), par );
			sys->dtau( dtau, time(idx), par, alpha );

			int nx,vx,np,vp;

			if( alpha == 0 )
			{
				sys->rhs( fx, time(idx), solData(idx), par );
				for( int p = 0; p < NDIM; p++ )
				{
					rhs( NDIM + p + NDIM*idx ) = -fx(p); 
					// if( fabs(fx(p)) >= 1e-4 ) std::cout<<"Bb";
				}
				
				nx = 1; np = 0;
				for( int r=0; r<NTAU; r++ )
				{
					double d = (dtau(r)-tau(r)/par(0));
					if( d != 0.0 )
					{
// 						std::cout<<"P0"; // working for the glass an logistic eqns
						vx = r;
						sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
						for( int p=0; p<NDIM; p++ )
						{
							for( int q=0; q<NDIM; q++ )
							{
								rhs( NDIM + p + NDIM*idx ) += d*dfx(p,q)*solData(idx)(q,NTAU+1+r);
							}
						}
					}
				}
			}else
			{
				nx = 0, np = 1; vp = alpha;
				sys->deri( dfx2, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
				for( int k = 0; k < NDIM; k++ ) rhs( NDIM + k + NDIM*idx ) = -par(0)*dfx2(k);
				
				nx = 1, np = 0;
				for( int r=0; r<NTAU; r++ )
				{
					if( dtau(r) != 0.0 )
					{
// 						std::cout<<"P"<<alpha; // it is workng for the Glass and logistic eqns
						vx = r;
						sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
						for( int p=0; p<NDIM; p++ )
						{
							for( int q=0; q<NDIM; q++ )
							{
								rhs( NDIM + p + NDIM*idx ) += dtau(r)*dfx(p,q)*solData(idx)(q,NTAU+1+r);
							}
						}
					}
				}
			}
		}
	}
}


void NColloc::RHS_x( SpMatrix& A, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{

	Matrix dfx(NDIM,NDIM);
	Matrix dummy(0,0);
	
	A.Clear('R');
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		A.NewL( 2 );
		A.WrLi(r,0) = r;   A.WrLi(r,1) = r+NDIM*NDEG*NINT;
		A.WrLx(r,0) = 1.0; A.WrLx(r,1) = -1.0;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			for( int r = 0; r < NDIM; r++ ){
				A.NewL( NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // check the line
			}
			
			int nx=1, vx, np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				if( ee(k,idx) != 0 )
				{
					vx = ee(k,idx)-1;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
				}
				for( int l = 0; l < NDEG+1; l++) // degree
				{
					for( int p = 0; p < NDIM; p++ )  // row
					{
						for( int q = 0; q < NDIM; q++ )  //column
						{
							WRIDX(A,idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(ee(k,idx),idx));
							if( ee(k,idx) == 0 )
							{
								if( p == q ) WRDAT(A,idx, i,j,p, k,l,q) += tt(0,l,idx);
							}else
							{
								WRDAT(A, idx, i,j,p, k,l,q) -= par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
							}
						}
					}
				}
			}
		}
	}

}


//! its very differnt from all of them
void NColloc::StabJac( StabMatrix& AB, const Vector& par, const JagMatrix3D& solData )
{
	Matrix dfx(NDIM,NDIM);
	Matrix dummy(0,0);

	AB.getA0().Clear('R');
	for( int s = 0; s < NMAT; s++ ) AB.getAI(s).Clear('R');
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		AB.getA0().NewL( 1 );						// L
		AB.getA0().WrLi(r,0) = r;
		AB.getA0().WrLx(r,0) = 1.0;

		for( int s = 0; s < NMAT; s++ )
		{
			if( s == 0 )
			{
				AB.getAI(s).NewL( 1 );						// M
				AB.getAI(s).WrLi(r,0) = r+NDIM*NDEG*NINT;
				AB.getAI(s).WrLx(r,0) = 1.0;
			}else
			{
				AB.getAI(s).NewL( 0 );
			}
		}
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int r = 0; r < NDIM; r++ )
			{
				AB.getA0().NewL( NDIM*szI(0, idx) );
			}
			for( int s = 0; s < NMAT; s++ )
			{
				for( int r = 0; r < NDIM; r++ )
				{
					AB.getAI(s).NewL( NDIM*szI(s+1, idx) );
				}
			}
			
			int nx=1, vx, np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				if( eeI(k,idx) != 0 )
				{
					vx = eeI(k,idx)-1;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
				}
				for( int l = 0; l < NDEG+1; l++)
				{
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							if( mmI(k,idx) == 0  ) // A matrix -- kkS(eeS(k,idx),idx) >= 0 
							{
								WRIDXI(AB.getA0(), idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(eeI(k,idx),idx));
								if( eeI(k,idx) == 0 )
								{
									if( p == q ) WRDATI(AB.getA0(), idx, i,j,p, k,l,q) += tt(0,l,idx);
								}
								else
								{
									WRDATI(AB.getA0(), idx, i,j,p, k,l,q) -= par(0)*dfx(p,q)*tt(eeI(k,idx),l,idx);
								}
							}
							else // B matrices
							{
								WRIDXI(AB.getAI(mmI(k,idx)-1), idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(eeI(k,idx),idx));
								if( eeS(k,idx) == 0 )
								{
// 									if( p == q ) WRDATSB(B,idx, i,j,p, k,l,q) -= tt(0,l,idx);
								}
								else
								{
									WRDATI(AB.getAI(mmI(k,idx)-1), idx, i,j,p, k,l,q) += par(0)*dfx(p,q)*tt(eeI(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}

// similar to RHS_x but with one boundary condition only and multiplied by Z
void NColloc::CharJac_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, double Z )
{

	Matrix dfx(NDIM,NDIM);
	Matrix dummy(0,0);

	A.Clear('R');
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		A.NewL( 2 );
		A.WrLi(r,0) = r;
		A.WrLx(r,0) = 1.0;
		// just for the test function
		A.WrLi(r,1) = NDIM*NINT*NDEG + r;
		A.WrLx(r,1) = - 1.0 * Z;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int r = 0; r < NDIM; r++ ){
				A.NewL( NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // check the line
			}
  
			int nx=1, vx, np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				const int zpow = (-kkI(ee(k,idx),idx) + NINT-1)/NINT;
				const double ZP = pow(Z,zpow);
// 				std::cout<<kkI(ee(k,idx),idx)<<" Z "<<Z<<" zpow "<<zpow<<" ZP "<<ZP<<"\n";
				if( ee(k,idx) != 0 )
				{
					vx = ee(k,idx)-1;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
				}
				for( int l = 0; l < NDEG+1; l++) // degree
				{
					for( int p = 0; p < NDIM; p++ )  // row
					{
						for( int q = 0; q < NDIM; q++ )  //column
						{
							WRIDX(A,idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(ee(k,idx),idx));
							if( ee(k,idx) == 0 )             // A
							{
								if( p == q ) WRDAT(A,idx, i,j,p, k,l,q) += tt(0,l,idx);
							}
							else
							{
								if( zpow > 0 )  // -zB0 - z^2*B1 ... - z^n*BN
								{
									WRDAT(A, idx, i,j,p, k,l,q) -= ZP*par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
								}
								else                         // A
								{
									WRDAT(A, idx, i,j,p, k,l,q) -= par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}


// this has to be changed only to packed complex.
void NColloc::CharJac_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, double Re, double Im )
{
	Matrix dfx(NDIM,NDIM);
	Matrix dummy(0,0);

	A.Clear('R');
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		A.NewL( 3 );    // Re
		// L
		A.WrLi(2*r,0) = 2*r;
		A.WrLx(2*r,0) = 1.0;
		// -z M
		A.WrLi(2*r,1) = 2*NDIM*NINT*NDEG + 2*r;
		A.WrLi(2*r,2) = 2*NDIM*NINT*NDEG + 2*r + 1;
		A.WrLx(2*r,1) = -1.0 * Re;
		A.WrLx(2*r,2) = 1.0 * Im;
		A.NewL( 3 );    // Im
		// L
		A.WrLi(2*r+1,0) = 2*r+1;
		A.WrLx(2*r+1,0) = 1.0;
		// -z M
		A.WrLi(2*r+1,1) = 2*NDIM*NINT*NDEG + 2*r;
		A.WrLi(2*r+1,2) = 2*NDIM*NINT*NDEG + 2*r + 1;
		A.WrLx(2*r+1,1) = -1.0 * Im;
		A.WrLx(2*r+1,2) = -1.0 * Re;
	}
	
	// computing the powers of the multiplier
	Vector ZReP(NMAT+1), ZImP(NMAT+1);
	ZReP(0) = 1.0; ZImP(0) = 0.0;
	ZReP(1) = Re;  ZImP(1) = Im;
	for( int r = 2; r < NMAT+1; r++ )
	{
		ZReP(r) = Re*ZReP(r-1) - Im*ZImP(r-1);
		ZImP(r) = Re*ZImP(r-1) + Im*ZReP(r-1);
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int r = 0; r < NDIM; r++ ){
				A.NewL( 2*NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // Re
				A.NewL( 2*NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // Im
			}
			
			int nx=1, vx, np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				const int zpow = (-kkI(ee(k,idx),idx) + NINT-1)/NINT;
				if( ee(k,idx) != 0 )
				{
					vx = ee(k,idx)-1;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
				}
				for( int l = 0; l < NDEG+1; l++) // degree
				{
					for( int p = 0; p < NDIM; p++ )  // row
					{
						for( int q = 0; q < NDIM; q++ )  //column
						{
							WRIDXCPLX(A,idx, i,j,p,0, k,l,q,0) = 2*(q+NDIM*(l+NDEG*kk(ee(k,idx),idx)));
							WRIDXCPLX(A,idx, i,j,p,0, k,l,q,1) = 2*(q+NDIM*(l+NDEG*kk(ee(k,idx),idx))) + 1;
							WRIDXCPLX(A,idx, i,j,p,1, k,l,q,0) = 2*(q+NDIM*(l+NDEG*kk(ee(k,idx),idx)));
							WRIDXCPLX(A,idx, i,j,p,1, k,l,q,1) = 2*(q+NDIM*(l+NDEG*kk(ee(k,idx),idx))) + 1;
							if( ee(k,idx) == 0 )             // A
							{
								if( p == q )
								{
									WRDATCPLX(A,idx, i,j,p,0, k,l,q,0) += tt(0,l,idx);
									WRDATCPLX(A,idx, i,j,p,1, k,l,q,1) += tt(0,l,idx);
								}
							}
							else
							{ //kkS(ee(k,idx),idx) < 0
								if( zpow > 0 )  // -zB0 - z^2*B1 ... - z^n*BN
								{
									WRDATCPLX(A, idx, i,j,p,0, k,l,q,0) -= ZReP(zpow)*par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
									WRDATCPLX(A, idx, i,j,p,0, k,l,q,1) += ZImP(zpow)*par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
									WRDATCPLX(A, idx, i,j,p,1, k,l,q,0) -= ZImP(zpow)*par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
									WRDATCPLX(A, idx, i,j,p,1, k,l,q,1) -= ZReP(zpow)*par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
								}
								else                         // A
								{
									WRDATCPLX(A, idx, i,j,p,0, k,l,q,0) -= par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
									WRDATCPLX(A, idx, i,j,p,1, k,l,q,1) -= par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}

// REAL

void NColloc::CharJac_x_p( Vector& V, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z, int alpha )
{
	Vector tau(NTAU);
	Vector dtau(NTAU);
	Vector fx(NDIM);
	Matrix dfx(NDIM,NDIM);
	Matrix t_dfx(NDIM,NDIM);
	Matrix dummy(0,0);
	
	V.Clear();
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		V(r) = 0.0;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			sys->tau( tau, time(idx), par );
			sys->dtau( dtau, time(idx), par, alpha );
			
			if( alpha == 0 )  // a periodusido szerint deriv...
			{
				int nx=1, vx[2], np=1, vp=alpha;
				for( int k = 0; k < NTAU; k++ )
				{
					const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
					const double ZP = pow(Z,zpow);
					nx = 1; np = 0;
					vx[0] = k;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, dummy );
					fx.Clear(); //OK this is cleared, but why here?
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							fx(p) -= dfx(p,q)*phiData(idx)(q,k);
						}
					}
					
					nx = 2;
					np = 0;
					for( int r = 0; r < NTAU; r++ )
					{
						double d = (dtau(r)-tau(r)/par(0));
						if( d != 0.0 )
						{
// 							std::cout<<"dP0"<<d;  // working for the glass and logistic eqn
							vx[1] = r; vx[0] = k; // CHANGE THIS 0 1
							sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idx) );
							for( int p = 0; p < NDIM; p++ )
							{
								for( int q = 0; q < NDIM; q++ )
								{
									fx(p) += d*dfx(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
								}
							}
						}
					}
					
					for( int p = 0; p < NDIM; p++ )
					{
						if( zpow > 0 )
						{
							V( NDIM +p+NDIM*(j+NDEG*i) ) += ZP*fx(p);
						}
						else
						{
							V( NDIM +p+NDIM*(j+NDEG*i) ) += fx(p);
						}
					}
				}
			}
			else
			{
				int nx=1, vx[2], np=1, vp=alpha;
				for( int k = 0; k < NTAU; k++ )
				{
					const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
					const double ZP = pow(Z,zpow);
					nx = 1; np = 1;
					vx[0] = k; vp=alpha; // though, this last is never changed...
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, dummy );
					fx.Clear();
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							fx(p) -= par(0)*dfx(p,q)*phiData(idx)(q,k);
						}
					}
					
					nx = 2; np = 0;
					for( int r = 0; r < NTAU; r++ )
					{
						if( dtau(r) != 0.0 )
						{
// 							std::cout<<"dP"<<alpha; // NOT TESTED !!! it is working for the glass eq
							vx[1] = r; vx[0] = k; // CHANGE THIS to 0, 1
							sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idx) );
							for( int p = 0; p < NDIM; p++ )
							{
								for( int q = 0; q < NDIM; q++ )
								{
									fx(p) += dtau(r)*dfx(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
								}
							}
						}
					}
					
					for( int p = 0; p < NDIM; p++ )
					{
						if( zpow > 0 )
						{
							V( NDIM +p+NDIM*(j+NDEG*i) ) += ZP*fx(p);
						}
						else
						{
							V( NDIM +p+NDIM*(j+NDEG*i) ) += fx(p);
						}
					}
				}
			}
		}
	}
}

// COMPLEX

void NColloc::CharJac_x_p( Vector& V, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Re, double Im, int alpha )
{
	Vector tau(NTAU);
	Vector dtau(NTAU);
	Matrix dfx(NDIM,NDIM);
	Matrix dfxRe(NDIM,NDIM);
	Matrix dfxIm(NDIM,NDIM);
	Vector fxRe(NDIM);
	Vector fxIm(NDIM);
	Matrix dummy(0,0);
	
	V.Clear();
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		V(2*r) = 0.0;
		V(2*r+1) = 0.0;
	}

	// computing the powers of the multiplier
	Vector ZReP(NMAT+1), ZImP(NMAT+1);
	ZReP(0) = 1.0; ZImP(0) = 0.0;
	ZReP(1) = Re;  ZImP(1) = Im;
	for( int r = 2; r < NMAT+1; r++ )
	{
		ZReP(r) = Re*ZReP(r-1) - Im*ZImP(r-1);
		ZImP(r) = Re*ZImP(r-1) + Im*ZReP(r-1);
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			const int idxRe = 2*(j+i*NDEG);
			const int idxIm = 2*(j+i*NDEG)+1;
			
			sys->tau( tau, time(idx), par );
			sys->dtau( dtau, time(idx), par, alpha );

			if( alpha == 0 )  // a periodusido szerint deriv...
			{
				int nx=1, vx[2], np=1, vp=alpha;
				for( int k = 0; k < NTAU; k++ )
				{
					const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
					nx = 1; np = 0;
					vx[0] = k;
					fxRe.Clear();
					fxIm.Clear();
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, dummy );
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							fxRe(p) -= dfx(p,q)*phiData(idxRe)(q,k);
							fxIm(p) -= dfx(p,q)*phiData(idxIm)(q,k);
						}
					}
					
					nx = 2; 
					np = 0;
					for( int r = 0; r < NTAU; r++ )
					{
						double d = (dtau(r)-tau(r)/par(0));
						if( d != 0.0 )
						{
							// std::cout<<"cP0";  // NOT TESTED !!!
							vx[1] = r; vx[0] = k; // CHANGE THIS to 0, 1
							sys->deri( dfxRe, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idxRe) );
							sys->deri( dfxIm, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idxIm) );
							for( int p = 0; p < NDIM; p++ )
							{
								for( int q = 0; q < NDIM; q++ )
								{
									fxRe(p) += d*dfxRe(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
									fxIm(p) += d*dfxIm(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
								}
							}
						}
					}
					
					for( int p = 0; p < NDIM; p++ )
					{
						if( zpow > 0 )
						{
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) )     += ZReP(zpow)*fxRe(p) - ZImP(zpow)*fxIm(p);
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) + 1 ) += ZImP(zpow)*fxRe(p) + ZReP(zpow)*fxIm(p);
						}
						else
						{
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) )     += fxRe(p);
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) + 1 ) += fxIm(p);
						}
					}
				}
			}
			else
			{
				int nx=1, vx[2], np=1, vp=alpha;
				for( int k = 0; k < NTAU; k++ )
				{
					const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
					nx = 1; np = 1;
					vx[0] = k; vp=alpha;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, dummy );
					fxRe.Clear();
					fxIm.Clear();
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							fxRe(p) -= par(0)*dfx(p,q)*phiData(idxRe)(q,k);
							fxIm(p) -= par(0)*dfx(p,q)*phiData(idxIm)(q,k);
						}
					}
					
					nx = 2;
					np = 0;
					for( int r = 0; r < NTAU; r++ )
					{
						if( dtau(r) != 0.0 )
						{
							// std::cout<<"cP"<<alpha; // NOT TESTED !!!
							vx[1] = r; vx[0] = k; // CHANGE THIS to 0, 1
							sys->deri( dfxRe, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idxRe) );
							sys->deri( dfxIm, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idxIm) );
							for( int p = 0; p < NDIM; p++ )
							{
								for( int q = 0; q < NDIM; q++ )
								{
									fxRe(p) += dtau(r)*dfxRe(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
									fxIm(p) += dtau(r)*dfxIm(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
								}
							}
						}
					}
					for( int p = 0; p < NDIM; p++ )
					{
						if( zpow > 0 )
						{
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) )     += ZReP(zpow)*fxRe(p) - ZImP(zpow)*fxIm(p);
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) + 1 ) += ZImP(zpow)*fxRe(p) + ZReP(zpow)*fxIm(p);
						}
						else
						{
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) )     += fxRe(p);
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) + 1 ) += fxIm(p);
						}
					}
				}
			}
		}
	}
}

void NColloc::CharJac_x_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z )
{
	Matrix dfx(NDIM,NDIM);
	Matrix t_dfx(NDIM,NDIM);
	Matrix dummy(0,0);
	
	A.Clear('R');

	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		A.NewL( 1 );
		A.WrLi(r,0) = r;
		A.WrLx(r,0) = 0.0;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int r = 0; r < NDIM; r++ ){
				A.NewL( NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // check the line
			}
  
			int nx=2, vx[2], np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				const int zpow = (-kkI(ee(k,idx),idx) + NINT-1)/NINT;
				const double ZP = pow(Z,zpow);
				if( ee(k,idx) != 0 )
				{
					vx[1] = ee(k,idx)-1;  // CHANGE THIS to 1
					
					dfx.Clear();
					t_dfx.Clear();
					for( int r = 0; r < NTAU; r++ )
					{
						vx[0] = r;        // CHANGE THIS to 0
						sys->deri( t_dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idx) );
						//std::cout<<"t_dfx "; t_dfx.Print();
						for( int ra = 0; ra < NDIM; ra++ )
						{
							for( int rb = 0; rb < NDIM; rb++ )
							{
								dfx(ra,rb) += t_dfx(ra,rb);
							}
						}
					}
					//std::cout<<"dfx "; dfx.Print();
				}
				for( int l = 0; l < NDEG+1; l++) // degree
				{
					for( int p = 0; p < NDIM; p++ )  // row
					{
						for( int q = 0; q < NDIM; q++ )  //column
						{
							WRIDX(A,idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(ee(k,idx),idx));
							if( ee(k,idx) != 0 )
							{
								if( zpow > 0 )  // -zB
								{
									WRDAT(A, idx, i,j,p, k,l,q) -= ZP*par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
								}
								else                         // A
								{
									WRDAT(A, idx, i,j,p, k,l,q) -= par(0)*dfx(p,q)*tt(ee(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}

void NColloc::CharJac_x_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Re, double Im )
{
	Matrix dfxRe(NDIM,NDIM);
	Matrix dfxIm(NDIM,NDIM);
	Matrix t_dfxRe(NDIM,NDIM);
	Matrix t_dfxIm(NDIM,NDIM);
	
	A.Clear('R');

	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		A.NewL( 1 );    // Re
		A.WrLi(2*r,0) = 2*r;
		A.WrLx(2*r,0) = 0.0;
		A.NewL( 1 );    // Im
		A.WrLi(2*r+1,0) = 2*r+1;
		A.WrLx(2*r+1,0) = 0.0;
	}
	
	// computing the powers of the multiplier
	Vector ZReP(NMAT+1), ZImP(NMAT+1);
	ZReP(0) = 1.0; ZImP(0) = 0.0;
	ZReP(1) = Re;  ZImP(1) = Im;
	for( int r = 2; r < NMAT+1; r++ )
	{
		ZReP(r) = Re*ZReP(r-1) - Im*ZImP(r-1);
		ZImP(r) = Re*ZImP(r-1) + Im*ZReP(r-1);
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			const int idxRe = 2*(j+i*NDEG);
			const int idxIm = 2*(j+i*NDEG)+1;

			for( int r = 0; r < NDIM; r++ )
			{
				A.NewL( NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // Re
				A.NewL( NDIM*( (NDEG+1)*(rr(NTAU,idx)+1) - dd(NTAU,idx) ) ); // Im
			}
  
			int nx=2, vx[2], np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				const int zpow = (-kkI(ee(k,idx),idx) + NINT-1)/NINT;
				if( ee(k,idx) != 0 )
				{
					vx[1] = ee(k,idx)-1; // CHANGE THIS to 1
					
					dfxRe.Clear();
					dfxIm.Clear();
					for( int r = 0; r < NTAU; r++ )
					{
						vx[0] = r;       // CHANGE THIS to 0
						sys->deri( t_dfxRe, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idxRe) );
						sys->deri( t_dfxIm, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idxIm) );
						//t_dfx.Print();
						for( int ra = 0; ra < NDIM; ra++ )
						{
							for( int rb = 0; rb < NDIM; rb++ )
							{
								dfxRe(ra,rb) += t_dfxRe(ra,rb);
								dfxIm(ra,rb) += t_dfxIm(ra,rb);
							}
						}
					}
					//std::cout<<"dfx "; dfx.Print();
				}
				for( int l = 0; l < NDEG+1; l++) // degree
				{
					for( int p = 0; p < NDIM; p++ )  // row
					{
						for( int q = 0; q < NDIM; q++ )  //column
						{
							WRIDXCPLXM(A,idx, i,j,p,0, k,l,q) = q+NDIM*(l+NDEG*kk(ee(k,idx),idx));
							WRIDXCPLXM(A,idx, i,j,p,1, k,l,q) = q+NDIM*(l+NDEG*kk(ee(k,idx),idx));
							if( ee(k,idx) == 0 )             // A
							{
								// WRDATCPLX(A,idx, i,j,p,0, k,l,q,0) += tt(0,l,idx); // no derivative !!!!!!
								// WRDATCPLX(A,idx, i,j,p,1, k,l,q,1) += tt(0,l,idx);
							}
							else
							{
								if( zpow > 0 )  // -zB
								{
									WRDATCPLXM(A, idx, i,j,p,0, k,l,q) -= 
										ZReP(zpow)*par(0)*dfxRe(p,q)*tt(ee(k,idx),l,idx) - ZImP(zpow)*par(0)*dfxIm(p,q)*tt(ee(k,idx),l,idx);
									WRDATCPLXM(A, idx, i,j,p,1, k,l,q) -= 
										ZImP(zpow)*par(0)*dfxRe(p,q)*tt(ee(k,idx),l,idx) + ZReP(zpow)*par(0)*dfxIm(p,q)*tt(ee(k,idx),l,idx);
								}
								else                         // A
								{
									WRDATCPLXM(A, idx, i,j,p,0, k,l,q) -= par(0)*dfxRe(p,q)*tt(ee(k,idx),l,idx);
									WRDATCPLXM(A, idx, i,j,p,1, k,l,q) -= par(0)*dfxIm(p,q)*tt(ee(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}

// this is for CharmatCPLX

void NColloc::CharJac_x_z( Vector& V, const Vector& par, const JagMatrix3D& solData, const Vector& phi, const JagMatrix3D& phiData, double Re, double Im )
{
	Matrix dfx(NDIM,NDIM);
	Matrix dummy(0,0);
	
	V.Clear();
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		V(2*r)   = -phi( 2*(r + NDIM*NDEG*NINT) );
		V(2*r+1) = -phi( 2*(r + NDIM*NDEG*NINT)+1 );
	}
	
	// computing the powers of the multiplier
	Vector ZReP(NMAT+1), ZImP(NMAT+1);
	ZReP(0) = 1.0; ZImP(0) = 0.0;
	ZReP(1) = Re;  ZImP(1) = Im;
	for( int r = 2; r < NMAT+1; r++ )
	{
		ZReP(r) = Re*ZReP(r-1) - Im*ZImP(r-1);
		ZImP(r) = Re*ZImP(r-1) + Im*ZReP(r-1);
	}
	
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			const int idxRe = 2*(j+i*NDEG);
			const int idxIm = 2*(j+i*NDEG)+1;

			int nx=1, vx, np=0, vp;
			for( int k = 0; k < NTAU; k++ )
			{
				const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
				vx = k; 
				sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
				for( int p = 0; p < NDIM; p++ )
				{
					for( int q = 0; q < NDIM; q++ )
					{
						if( zpow > 0 )
						{
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) )     -= zpow*(ZReP(zpow-1)*par(0)*dfx(p,q)*phiData(idxRe)(q,k) -
							                                              ZImP(zpow-1)*par(0)*dfx(p,q)*phiData(idxIm)(q,k));
							V( 2*(NDIM +p+NDIM*(j+NDEG*i)) + 1 ) -= zpow*(ZImP(zpow-1)*par(0)*dfx(p,q)*phiData(idxRe)(q,k) +
							                                              ZReP(zpow-1)*par(0)*dfx(p,q)*phiData(idxIm)(q,k));
						}
					}
				}
			}
		}
	}
}

//!
//! from now CharmatLPAUT
//!

void NColloc::CharJac_mB( SpMatrix& B, const Vector& par, const JagMatrix3D& solData, double Z )
{

	Matrix dfx(NDIM,NDIM);
	Matrix dummy(0,0);

	B.Clear('R');
	
	Vector ZP(NMAT+1);
	ZP(0) = 1.0;
	for( int r = 1; r < NMAT+1; r++ )
	{
		ZP(r) = Z*ZP(r-1);
	}
	
	// boundary conditions: no boundary condition
	for( int r = 0; r < NDIM; r++ )
	{
		B.NewL( 1 );						// M
		B.WrLi(r,0) = r+NDIM*NDEG*NINT;
		B.WrLx(r,0) = -1.0;
	}

	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int r = 0; r < NDIM; r++ )
			{
				B.NewL( NDIM*szS(1,idx) );
			}
			
			int nx=1, vx, np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				const int zpow = (-kkI(eeS(k,idx),idx) + NINT-1)/NINT;
				if( eeS(k,idx) != 0 )
				{
					vx = eeS(k,idx)-1;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx, np, &vp, dummy );
					
					for( int l = 0; l < NDEG+1; l++)
					{
						for( int p = 0; p < NDIM; p++ )
						{
							for( int q = 0; q < NDIM; q++ )
							{
								if( zpow > 0 ) // B matrix -- kkS(eeS(k,idx),idx) < 0 
								{
									WRIDXS(B, idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(eeS(k,idx),idx));
									WRDATS(B, idx, i,j,p, k,l,q) -= zpow*ZP(zpow-1)*par(0)*dfx(p,q)*tt(eeS(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}

// same as CharJac_x_p, but only writes the B part

void NColloc::CharJac_mB_p( Vector& V, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z, int alpha )
{
	Vector tau(NTAU);
	Vector dtau(NTAU);
	Vector fx(NDIM);
	Matrix dfx(NDIM,NDIM);
	Matrix t_dfx(NDIM,NDIM);
	Matrix dummy(0,0);
	
	V.Clear();
	
	// boundary conditions
	for( int r = 0; r < NDIM; r++ )
	{
		V(r) = 0.0;
	}
	
	Vector ZP(NMAT+1);
	ZP(0) = 1.0;
	for( int r = 1; r < NMAT+1; r++ )
	{
		ZP(r) = Z*ZP(r-1);
	}
	
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			sys->tau( tau, time(idx), par );
			sys->dtau( dtau, time(idx), par, alpha );
			
			if( alpha == 0 )  // a periodusido szerint deriv...
			{
				int nx=1, vx[2], np=1, vp=alpha;
				for( int k = 0; k < NTAU; k++ )
				{
					const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
					nx = 1; np = 0;
					vx[0] = k;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, dummy );
					fx.Clear();
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							fx(p) -= dfx(p,q)*phiData(idx)(q,k);
						}
					}
					
					nx = 2;
					np = 0;
					for( int r = 0; r < NTAU; r++ )
					{
						double d = (dtau(r)-tau(r)/par(0));
						if( d != 0.0 )
						{
// 							std::cout<<"mP0";  // working for the glass and logistic eqn
							vx[1] = r; vx[0] = k;
							sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idx) );
							for( int p = 0; p < NDIM; p++ )
							{
								for( int q = 0; q < NDIM; q++ )
								{
									fx(p) += d*dfx(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
								}
							}
						}
					}
					
					for( int p = 0; p < NDIM; p++ )
					{
						if( zpow > 0 )
						{
							V( NDIM+p+NDIM*(j+NDEG*i) ) += zpow*ZP(zpow-1)*fx(p);
						}
					}
				}
			}
			else
			{
				int nx=1, vx[2], np=1, vp=alpha;
				for( int k = 0; k < NTAU; k++ )
				{
					const int zpow = (-kkI(k+1,idx) + NINT-1)/NINT;
					nx = 1; np = 1;
					vx[0] = k;
					sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, dummy );
					fx.Clear();
					for( int p = 0; p < NDIM; p++ )
					{
						for( int q = 0; q < NDIM; q++ )
						{
							fx(p) -= par(0)*dfx(p,q)*phiData(idx)(q,k);
						}
					}
					
					nx = 2;
					np = 0;
					for( int r = 0; r < NTAU; r++ )
					{
						if( dtau(r) != 0.0 )
						{
// 							std::cout<<"mP"<<alpha; // NOT TESTED !!! it is working for the glass eq
							vx[1] = r; vx[0] = k;
							sys->deri( dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idx) );
							for( int p = 0; p < NDIM; p++ )
							{
								for( int q = 0; q < NDIM; q++ )
								{
									fx(p) += dtau(r)*dfx(p,q)*solData(idx)(q,NTAU+1+r); // this MUST be k, was r
								}
							}
						}
					}
					
					for( int p = 0; p < NDIM; p++ )
					{
						if( zpow > 0 )
						{
							V( NDIM +p+NDIM*(j+NDEG*i) ) += zpow*ZP(zpow-1)*fx(p);
						}
					}
				}
			}
		}
	}
}

// like x_x but write bpart only

void NColloc::CharJac_mB_x( SpMatrix& B, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z )
{
	Matrix dfx(NDIM,NDIM);
	Matrix t_dfx(NDIM,NDIM);
	Matrix dummy(0,0);

	B.Clear('R');
	
	// boundary conditions: no boundary condition
	for( int r = 0; r < NDIM; r++ )
	{
		B.NewL( 0 );						// M
// 		B.WrLi(r,0) = r+NDIM*NDEG*NINT;
// 		B.WrLx(r,0) = 0.0;
	}
	
	Vector ZP(NMAT+1);
	ZP(0) = 1.0;
	for( int r = 1; r < NMAT+1; r++ )
	{
		ZP(r) = Z*ZP(r-1);
	}
	
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;

			for( int r = 0; r < NDIM; r++ )
			{
				B.NewL( NDIM*szS(1,idx) );
			}

			int nx=2, vx[2], np=0, vp;
			for( int k = 0; k < NTAU+1; k++ )
			{
				const int zpow = (-kkI(eeS(k,idx),idx) + NINT-1)/NINT;
				if( eeS(k,idx) != 0 )
				{
					vx[1] = eeS(k,idx)-1; // CHANGE THIS to 1
					
					dfx.Clear();
					t_dfx.Clear();
					for( int r = 0; r < NTAU; r++ )
					{
						vx[0] = r; // CHANGE THIS to 0
						sys->deri( t_dfx, time(idx), solData(idx), par, nx, &vx[0], np, &vp, phiData(idx) );
						//std::cout<<"t_dfx "; t_dfx.Print();
						for( int ra = 0; ra < NDIM; ra++ )
						{
							for( int rb = 0; rb < NDIM; rb++ )
							{
								dfx(ra,rb) += t_dfx(ra,rb);
							}
						}
					}
					//std::cout<<"dfx "; dfx.Print();
					for( int l = 0; l < NDEG+1; l++) // degree
					{
						for( int p = 0; p < NDIM; p++ )  // row
						{
							for( int q = 0; q < NDIM; q++ )  //column
							{
								if( zpow > 0 ) // A matrix -- kkS(eeS(k,idx),idx) >= 0 
								{
									WRIDXS(B, idx, i,j,p, k,l,q) = q+NDIM*(l+NDEG*kk(eeS(k,idx),idx));
									WRDATS(B, idx, i,j,p, k,l,q) -= zpow*ZP(zpow-1)*par(0)*dfx(p,q)*tt(eeS(k,idx),l,idx);
								}
							}
						}
					}
				}
			}
		}
	}
}

//!
//! until this CharmatLPAUT
//!

//!
//! this is also for CharmatLPAUT: for computing q_0 and its derivatives: q_0, D_x q_0, D_p q_0
//!

void NColloc::CharJac_MSHphi( Vector& V, const Vector& par, const JagMatrix3D& solData )
{
	Vector fx(NDIM); // it should be a variable in the class itself
	
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		const double h = mesh(i+1) - mesh(i);
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			sys->rhs( fx, mesh(i) + j*h/(double)NDEG, solData(idx), par );
			for( int k = 0; k < NDIM; k++ ) V( k + NDIM*(j+NDEG*i) ) = - fx(k);
		}
	}
	// make it periodic
	sys->rhs( fx, 1.0, solData(NDEG*NINT), par );
	for( int r = 0; r < NDIM; r++ )
	{
		V(NDIM*NDEG*NINT + r) = - fx(r);
	}
}

void NColloc::CharJac_MSHphi_p( Vector& V, const Vector& par, const JagMatrix3D& solData, int alpha )
{
	Vector tau(NTAU);
	Vector dtau(NTAU);
	Vector fx(NDIM);
	Matrix dfx(NDIM,NDIM);
	Matrix dfx2(NDIM,1);
	Matrix dummy(0,0);
	
	V.Clear(); /// it is not cleared otherwise!!!!
	
	// boundary conditions
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		const double h = mesh(i+1) - mesh(i);
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			
			sys->tau( tau, time(idx), par );
			sys->dtau( dtau, time(idx), par, alpha );

			int nx,vx,np,vp;

			if( alpha == 0 )
			{
// 				sys->rhs( fx, time(idx), solData(idx), par );
// 				for( int p = 0; p < NDIM; p++ )
// 				{
// 					V( p + NDIM*idx ) = -fx(p); 
// 					// if( fabs(fx(p)) >= 1e-4 ) std::cout<<"Bb";
// 				}
				nx = 1; np = 0;
				for( int r=0; r<NTAU; r++ )
				{
					const double d = (dtau(r)-tau(r)/par(0))/par(0);
					if( d != 0.0 )
					{
// 						std::cout<<"P0"; // working for the glass an logistic eqns
						vx = r;
						sys->deri( dfx, mesh(i) + j*h/(double)NDEG, solData(idx), par, nx, &vx, np, &vp, dummy );
						for( int p=0; p<NDIM; p++ )
						{
							for( int q=0; q<NDIM; q++ )
							{
								V( p + NDIM*idx ) += d*dfx(p,q)*solData(idx)(q,NTAU+r);
							}
						}
					}
				}
			}else
			{
				nx = 0, np = 1; vp = alpha;
				sys->deri( dfx2, mesh(i) + j*h/(double)NDEG, solData(idx), par, nx, &vx, np, &vp, dummy );
				for( int k = 0; k < NDIM; k++ ) V( k + NDIM*idx ) = - dfx2(k);
				
				nx = 1, np = 0;
				for( int r=0; r<NTAU; r++ )
				{
					const double d = dtau(r)/par(0);
					if( d != 0.0 )
					{
// 						std::cout<<"P"<<alpha; // it is workng for the Glass and logistic eqns
						vx = r;
						sys->deri( dfx, mesh(i) + j*h/(double)NDEG, solData(idx), par, nx, &vx, np, &vp, dummy );
						for( int p=0; p<NDIM; p++ )
						{
							for( int q=0; q<NDIM; q++ )
							{
								V( p + NDIM*idx ) += d*dfx(p,q)*solData(idx)(q,NTAU+r);
							}
						}
					}
				}
			}
		}
	}
	for( int r = 0; r < NDIM; r++ )
	{
		V(r+NDIM*NDEG*NINT) = V(r);
	}
}


//!
//! until this CharmatLPAUT
//!

//----------------------------------------------------------------------------
//
// here are the integration routines
//
//----------------------------------------------------------------------------

void NColloc::Star( Vector& V1, Vector& V2 )
{
	V1.Clear();
	for( int i=0; i < NINT; i++ )
	{
		const double dx = mesh(i+1) - mesh(i);
		for( int j=0; j < NDIM; j++ )
		{
			// ez itt a matrixszorzas
			for( int k=0; k < NDEG+1; k++ )
			{
				for( int l=0; l < NDEG+1; l++ )
				{
					V1( j+NDIM*(k+i*NDEG) ) += dx*metric(k,l)*V2( j+NDIM*(l+i*NDEG) );
				}
			}
		}
	}

#ifdef MADD // whether we need to add the headpoint ?

	// Now we add (M udot)^* M\phi + <udot, \phi > ...
	for( int j=0; j < NDIM; j++ ){
		V1( NINT*NDEG*NDIM + j ) += V2( NINT*NDEG*NDIM + j );
	}

#endif
}

double NColloc::Integrate( Vector& V1, Vector& V2 )
{
	double res=0.0,head=0.0;
	for( int i=0; i < NINT; i++ )
	{
		const double dx = mesh(i+1) - mesh(i);
		for( int j=0; j < NDIM; j++ )
		{
			// ez itt a matrixszorzas
			for( int k=0; k < NDEG+1; k++ )
			{
				for( int l=0; l < NDEG+1; l++ )
				{
					res += dx*V1( j+NDIM*(k+i*NDEG) )*metric(k,l)*V2( j+NDIM*(l+i*NDEG) );
				}
			}
		}
	}
	
#ifdef MADD
	for( int j=0; j < NDIM; j++ ){
		head += V1( NINT*NDEG*NDIM + j )*V2( NINT*NDEG*NDIM + j );
	}
#endif
	
	return res + head;
}

double NColloc::IntegrateCont( Vector& V1, Vector& V2, Vector& V3 )
{
	double res=0.0,head=0.0;
	for( int i=0; i < NINT; i++ )
	{
		const double dx = mesh(i+1) - mesh(i);
		for( int j=0; j < NDIM; j++ ){
		// ez itt a matrixszorzas
			for( int k=0; k < NDEG+1; k++ ){
				for( int l=0; l < NDEG+1; l++ ){
					res += dx*V1( j+NDIM*(k+i*NDEG) )*metric(k,l)*( V2( j+NDIM*(l+i*NDEG) ) - V3( j+NDIM*(l+i*NDEG) ) );
				}
			}
		}
	}
#ifdef MADD
	for( int j=0; j < NDIM; j++ ){
		head += V1( NINT*NDEG*NDIM + j )*( V2( NINT*NDEG*NDIM + j ) - V3( NINT*NDEG*NDIM + j ) );
	}
#endif
	return res + head;
}

void NColloc::PhaseStar( Vector& V1, Vector& V2 )
{
	V1.Clear();
	for( int i=0; i < NINT; i++ )
	{
		for( int j=0; j < NDIM; j++ )
		{
			// ez itt a matrixszorzas
			for( int k=0; k < NDEG+1; k++ )
			{
				for( int l=0; l < NDEG+1; l++ )
				{
					V1( j+NDIM*(k+i*NDEG) ) += metricPhase(k,l)*V2( j+NDIM*(l+i*NDEG) );
				}
			}
		}
	}
}

void NColloc::PhaseRotStar( Vector& V1, Vector& V2, Array1D<int>& Re, Array1D<int>& Im )
{
	V1.Clear();
	for( int i=0; i < NINT; i++ )
	{
		const double dx = mesh(i+1) - mesh(i);
		for( int j=0; j < Re.Size(); j++ )
		{
			// ez itt a matrixszorzas
			for( int k=0; k < NDEG+1; k++ )
			{
				for( int l=0; l < NDEG+1; l++ )
				{
					V1( Re(j)+NDIM*(k+i*NDEG) ) -= dx*metric(k,l)*V2( Im(j)+NDIM*(l+i*NDEG) );
					V1( Im(j)+NDIM*(k+i*NDEG) ) += dx*metric(k,l)*V2( Re(j)+NDIM*(l+i*NDEG) );
				}
			}
		}
	}
}

void NColloc::Import( Vector& outs, const Vector& in, const Vector& msh_, int deg_ )
{
	if( (msh_.Size()-1) % deg_ != 0 ) { std::cout<<"NColloc::Import: bad mesh"; PDError(-1); }
	int int_ = (msh_.Size()-1)/deg_;
	Vector msh( int_+1 );
	Vector in_mesh(deg_+1);
	Vector in_lgr(deg_+1);
	
	std::cout<<"Import\n";
	repr_mesh( in_mesh ); /// now we can use chebyshev
	for( int i = 0; i < int_; i++ ) msh(i) = msh_( deg_*i );
	msh( int_ ) = 1.0;
	
	for( int i = 0; i < NINT; i++ )
	{
		for( int j = 0; j < NDEG; j++ )
		{
			double t = mesh(i) + meshINT(j)*(mesh(i+1)-mesh(i));
			int k = meshlookup( msh, t );
			// std::cout<<"int "<<i<<" "<<k<<"\n";
			double c = (t-msh(k))/(msh(k+1)-msh(k));
			if( c < 0.0 || c > 1.0 ) { std::cout<<"NColloc::Import: FATAL ERROR"; PDError(-1); }
			
			poly_lgr( in_mesh, in_lgr, c );
			// in_lgr.Print();
			for( int p = 0; p < NDIM; p++ )
			{
				outs( p + NDIM*(NDEG*i + j) ) = 0.0;
				for( int r = 0; r < deg_+1; r++ )
				{
					outs( p + NDIM*(NDEG*i + j) ) += in( p + NDIM*(deg_*k + r) ) * in_lgr(r);
				}
			}
		}
	}
	for( int p = 0; p < NDIM; p++ )
	{
		outs( p + NDIM*(NINT*NDEG) ) = in( p + NDIM*(int_*deg_) );
	}
}


// it exports for CollocTR and PointTR, so no last value is necessary
void   NColloc::Export( Vector& outs, const Vector& mshint, const Vector& mshdeg, const Vector& in )
{
	int nint_ = mshint.Size()-1;
	int ndeg_ = mshdeg.Size()-1;
	Vector in_mesh(NDEG+1);
	Vector in_lgr(NDEG+1);
	
	for( int i = 0; i < NDEG+1; i++ ) in_mesh(i) = i*1.0/NDEG;
	
	for( int i = 0; i < nint_; i++ )
	{
		for( int j = 0; j < ndeg_; j++ )
		{
			double t = mshint(i) + mshdeg(j)/nint_;
			int k = meshlookup( mesh, t );
			// std::cout<<"int "<<i<<" "<<k<<"\n";
			double c = (t - mesh(k))/(mesh(k+1)-mesh(k));  // mesh is the interval mesh in the class
			
			poly_lgr( in_mesh, in_lgr, c );
			// in_lgr.Print();
			for( int p = 0; p < NDIM; p++ )
			{
				outs( p + NDIM*(j + i*ndeg_) ) = 0.0;
				for( int r = 0; r < NDEG+1; r++ )
				{
					outs( p + NDIM*(j + i*ndeg_) ) += in( p + NDIM*(r + k*NDEG) ) * in_lgr(r);
				}
			}
		}
	}
}
