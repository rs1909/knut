// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "config.h"

#include "testfunct.h"
#include "matrix.h"
#include "spmatrix.h"
#include "hypermatrix.h"
#include "ncolloc.h"

#define NKERNITER 20
#define KERNEPS 1.0e-12

#define NDIM (col.Ndim())
#define NDEG (col.Ndeg())
#define NINT (col.Nint())
#define NTAU (col.Ntau())
#define NPAR (col.Npar())

template<bool trans> inline void rotbord( Vector& V, NColloc& col, const Vector& IN, Array1D<int>& Re, Array1D<int>& Im )
{
	V.Clear();
	for( int idx = 0; idx < NINT*NDEG+1; idx++ )
	{
		for( int k = 0; k < Re.Size(); k++ )
		{
			if(trans)
			{
				V( Re(k) + NDIM*idx ) = - IN( Im(k) + NDIM*idx );
				V( Im(k) + NDIM*idx ) = IN( Re(k) + NDIM*idx );
			}else
			{
				V( Re(k) + NDIM*idx ) = IN( Im(k) + NDIM*idx );
				V( Im(k) + NDIM*idx ) = - IN( Re(k) + NDIM*idx );
			}
		}
	}
}

inline void conjugate( Vector& out, const Vector& inp )
{
	for( int i = 0; i < out.Size()/2; ++i )
	{
		out( 2*i )   = -inp( 2*i+1 );
		out( 2*i+1 ) =  inp( 2*i );
    }
}

TestFunct::TestFunct( NColloc& col, double Z ) :
	ZZ(Z),
	AHAT( NDIM*(NDEG*NINT+1), 0, 1, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	rhs( NDIM*(NDEG*NINT+1) ),
	uu( NDIM*(NDEG*NINT+1) ), vv( NDIM*(NDEG*NINT+1) ),
	uudiff( NDIM*(NDEG*NINT+1) ), vvdiff( NDIM*(NDEG*NINT+1) ),
	vvData( NDEG*NINT, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunct::~TestFunct()
{
}

void TestFunct::Init( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	AHAT.getA13(0).Rand();
	AHAT.getA31(0).Rand();
	AHAT.getA33()(0,0) = 0.0;
	// norming the borders
	AHAT.getA31(0) /= sqrt( AHAT.getA31(0) * AHAT.getA31(0) );
	AHAT.getA13(0) /= sqrt( AHAT.getA13(0) * AHAT.getA13(0) );
	vv = AHAT.getA31(0);
	uu = AHAT.getA13(0);
	double one = 1.0;
	double gg = 0.0;
	double hh = 0.0;
	
	rhs.Clear();
	
	double ggdiff, hhdiff;
	double diffnorm = 1.0;
	int it = 0;
	do {
		AHAT.Multiply<false>( rhs, one, vv, gg );
		one -= 1.0;
		AHAT.Solve  ( vvdiff, ggdiff, rhs, one );
		AHAT.Multiply<true>( rhs, one, uu, hh );
		one -= 1.0;
		AHAT.SolveTR( uudiff, hhdiff, rhs, one );
		vv -= vvdiff;
		gg -= ggdiff;
		uu -= uudiff;
		hh -= hhdiff;
		vv /= sqrt( vv * vv );
		uu /= sqrt( uu * uu );
		AHAT.getA31(0) = vv;
		AHAT.getA13(0) = uu;
		diffnorm = std::max<double>(sqrt(vvdiff*vvdiff+ggdiff*ggdiff),sqrt(uudiff*uudiff+hhdiff*hhdiff));
		std::cout<<"dnor "<<diffnorm<<"\n";
		std::cout<<"vv "<<sqrt(AHAT.getA31(0)*AHAT.getA31(0))<<", uu "<<sqrt(AHAT.getA13(0)*AHAT.getA13(0))<<"\n";
	} while( (++it < NKERNITER)&&(diffnorm > KERNEPS) );
	if( diffnorm > KERNEPS ) std::cout<<"TestFunct::Init: warning: No convergence in finding the singular vector. Residual = "<<diffnorm<<"\n";
// 	std::cout<<"TF: "<<gg<<", "<<hh<<"\n";
}

double TestFunct::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	double one = 1.0;
	double gg;
	double hh;
	double ggdiff, hhdiff;
	
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	
	AHAT.Multiply<false>( rhs, one, vv, gg );
	one -= 1.0;
	AHAT.Solve  ( vvdiff, ggdiff, rhs, one );
	AHAT.Multiply<true>( rhs, one, uu, hh );
	one -= 1.0;
	AHAT.SolveTR( uudiff, hhdiff, rhs, one );
	vv -= vvdiff;
	gg -= ggdiff;
	uu -= uudiff;
	hh -= hhdiff;
	AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
	AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
// 	std::cout<<"TF: "<<gg<<", "<<hh;
// 	if( gg > 0.0 ) std::cout<<"\t+++\n";
// 	else           std::cout<<"\t---\n";
	return gg;
}

double TestFunct::Funct_p( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData, int alpha )
{
	col.Interpolate( vvData, vv );
	col.CharJac_x_p( A_p, par, solData, vvData, ZZ, alpha );
	const double gg_p = (uu * A_p);
// 	std::cout<<"GP "<<alpha<<":"<<gg_p;
	return gg_p; // positive, because the test function is not negated in the Jacobian
}

void   TestFunct::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	col.Interpolate( vvData, vv );
	col.CharJac_x_x( A_x, par, solData, vvData, ZZ );
	func = !A_x * uu; // positive, because the test function is not negated in the Jacobian
// 	func *= 100.0;
}

void   TestFunct::Switch( Vector& phi )
{
	phi = vv;
}

/// ---------------------------------------------------------
/// test function for TORUS BIFURCATIONS
/// ---------------------------------------------------------

TestFunctCPLX::TestFunctCPLX( NColloc& col ) :
	first(true),
	AHAT( 2*NDIM*(NDEG*NINT+1), 0, 2, 4*NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( 2*NDIM*(NDEG*NINT+1) ),
	A_x( 'R', 2*NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1), 2*NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	rhs( 2*NDIM*(NDEG*NINT+1) ),
	one( 2 ),
	uu( 2*NDIM*(NDEG*NINT+1) ),
	vv( 2*NDIM*(NDEG*NINT+1) ),
	gg( 2 ),
	hh( 2 ),
	vvData( 2*NDEG*NINT, NDIM, NTAU+1 )
{
}

TestFunctCPLX::~TestFunctCPLX( )
{
}

void TestFunctCPLX::Init( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData,
									double Re, double Im )
{
	ZRe = Re; ZIm = Im;
	col.CharJac_x( AHAT.getA11(), par, solData, Re, Im );
	AHAT.getA13(0).Rand();
	AHAT.getA31(0).Rand();
	AHAT.getA33().Clear();
	// norming the borders
	AHAT.getA31(0) /= sqrt( AHAT.getA31(0) * AHAT.getA31(0) );
	AHAT.getA13(0) /= sqrt( AHAT.getA13(0) * AHAT.getA13(0) );
	// conjugate
	conjugate( AHAT.getA31(1), AHAT.getA31(0) );
	conjugate( AHAT.getA13(1), AHAT.getA13(0) );
	
	rhs.Clear();
	one(0) = 1.0; one(1) = 0.0;
	
	Vector vdiff( vv ), udiff( uu );
	double diffnorm = 1.0;
	int it = 0;
	do {
		AHAT.Solve  ( 2, vv, gg, rhs, one );
		AHAT.SolveTR( 2, uu, hh, rhs, one );
		udiff = AHAT.getA13(0);
		udiff -= uu;
		vdiff = AHAT.getA31(0);
		vdiff -= vv;
		diffnorm = std::max<double>(sqrt(udiff*udiff),sqrt(vdiff*vdiff));
// 		std::cout<<"std::max<double>(udiff, vdiff) "<<diffnorm<<"\n";
		AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
		AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
		
		conjugate( AHAT.getA13(1), AHAT.getA13(0) );
		conjugate( AHAT.getA31(1), AHAT.getA31(0) );
	} while( (++it < NKERNITER)&&(diffnorm > KERNEPS) );
	if( diffnorm > KERNEPS ) std::cout<<"TestFunctCPLX::Init: warning: No convergence in finding the singular vector. Residual = "<<diffnorm<<"\n";
// 	std::cout<<"TF: "<<gg<<", "<<hh<<"\n";
}

void TestFunctCPLX::Funct  ( double& f1, double& f2,
					NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, double Re, double Im )
{
	if( first ){ Init( col, par, sol, solData, Re, Im ); first = false; }
	
	ZRe = Re; ZIm = Im;
	col.CharJac_x( AHAT.getA11(), par, solData, Re, Im );
	AHAT.Solve  ( 2, vv, gg, rhs, one );
	AHAT.SolveTR( 2, uu, hh, rhs, one );
	AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
	AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
	conjugate( AHAT.getA13(1), AHAT.getA13(0) );
	conjugate( AHAT.getA31(1), AHAT.getA31(0) );
	// for later use
	col.InterpolateCPLX( vvData, AHAT.getA31(0)/*vv*/ );
	
	f1 = gg(0);
	f2 = gg(1);
// 	std::cout<<"gg: "<<gg(0)<<", "<<gg(1)<<", "<<hh(0)<<", "<<hh(1)<<"\n";
// 	if( gg(0) > 0.0 ) std::cout<<"\t+++\n";
// 	else               std::cout<<"\t---\n";
// 	if( gg(1) > 0.0 ) std::cout<<"\t+++\n";
// 	else               std::cout<<"\t---\n";
}

void TestFunctCPLX::Funct_p( double& f1, double& f2,
					NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData,
					int alpha )
{
	col.CharJac_x_p( A_p, par, solData, vvData, ZRe, ZIm, alpha );
	f1 = (AHAT.getA13(0)/*uu*/ * A_p);
	f2 = (AHAT.getA13(1)/*uu^conj*/ * A_p);
}

void TestFunctCPLX::Funct_z( double& f1, double& f2,
					NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	col.CharJac_x_z( A_p, par, solData, AHAT.getA31(0), vvData, ZRe, ZIm );
	const double dzre = (AHAT.getA13(0)/*uu*/ * A_p);
	const double dzim = (AHAT.getA13(1)/*uu^conj*/ * A_p);
	f1 = ( -dzim*ZRe - dzre*ZIm );
	f2 = (  dzre*ZRe - dzim*ZIm );
}

void TestFunctCPLX::Funct_x( Vector& func1, Vector& func2,
					NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	col.CharJac_x_x( A_x, par, solData, vvData, ZRe, ZIm );
	func1 = !A_x * AHAT.getA13(0); /*uu*/
	func2 = !A_x * AHAT.getA13(1); /*uu^conj*/
}

void TestFunctCPLX::Switch( Vector& Re, Vector& Im, double& alpha )
{
	P_ASSERT_X( (2*Re.Size() == vv.Size())&&(2*Im.Size() == vv.Size()), "TestFunctCPLX::Switch: Bad sizes\n" );
	std::cout<<"zRe="<<ZRe<<", zIm="<<ZIm<<"\n";
	for( int i=0; i<Re.Size(); i++ )
	{
		Re(i) = vv( 2*i );
		Im(i) = vv( 2*i + 1 );
	}
	if( ZRe > 0.0 )
	{
		alpha = atan( fabs(ZIm/ZRe) );
	}
	else
	{
		alpha = atan( fabs(ZRe/ZIm) ) + M_PI/2.0;
	}
}

#ifdef DEBUG
	#include <iomanip>
	#include <fstream>
#endif

void TestFunctCPLX::SwitchHB( Vector& Re, Vector& Im, NColloc& col, const Vector& par )
{
// probably it is not the best point in the profile...
#ifdef DEBUG
	std::ofstream file( "eigenvec" );
	file<<std::scientific;
	file.precision(12);

	for( int i = 0; i < NINT; i++ )
	{
		for( int j = 0; j < NDEG+1; j++ )
		{
			const double t = col.Profile( i, j );
			for( int p = 0; p < NDIM; p++ )
			{
				file<<vv( 2*(p + (j+i*NDEG)*NDIM) )<<"\t";
			}
			file<<par(0)*t<<"\n";
		}
	}
#endif

	for( int i = 0; i < NDIM; i++ )
	{
		Re(i) = vv( 2*i );
		Im(i) = vv( 2*i + 1 );
	}
}

/// ---------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems
/// ---------------------------------------------------------

TestFunctLPAUT::TestFunctLPAUT( NColloc& col, double Z ) :
	ZZ(Z),
	AHAT( NDIM*(NDEG*NINT+1), 0, 2, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB_p( NDIM*(NDEG*NINT+1) ),
	mB_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	rhs( NDIM*(NDEG*NINT+1) ),
	phi( NDIM*(NDEG*NINT+1) ),
	temp( NDIM*(NDEG*NINT+1) ),
	DpPhi( NDIM*(NDEG*NINT+1) ),
	DxPhi( NDIM*(NDEG*NINT+1) ),
	uu2( NDIM*(NDEG*NINT+1) ),
	vv2( NDIM*(NDEG*NINT+1) ),
	gg2(2),
	hh2(2),
	one2(2),
	phiData( NDEG*NINT, NDIM, 2*NTAU+1 ),
	vv2Data( NDEG*NINT, NDIM, 2*NTAU+1 ),
	solMSHData( NDEG*NINT+1, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunctLPAUT::~TestFunctLPAUT() 
{
}

void TestFunctLPAUT::Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	col.InterpolateMSH( solMSHData, sol );
	col.CharJac_MSHphi( phi, par, solMSHData );
	
	col.Star( AHAT.getA31(0), phi );
	AHAT.getA13(0) = mB * phi;
	
	AHAT.getA13(1).Rand();
	AHAT.getA31(1).Rand();
	AHAT.getA33().Clear(); // 2x2 matrix
	// norming the borders
	AHAT.getA31(1) /= sqrt( AHAT.getA31(1) * AHAT.getA31(1) );
	AHAT.getA13(1) /= sqrt( AHAT.getA13(1) * AHAT.getA13(1) );
	
	// generating the non-trivial kernel	
	one2(0) = 0.0;
	one2(1) = 1.0;
	
	Vector v2diff( vv2 ), u2diff( uu2 );
	double g2diff, h2diff;
	double diffnorm = 1.0;
	int it = 0;
	do {
		AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
		AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
		u2diff = AHAT.getA13(1);
		u2diff -= uu2;
		h2diff = AHAT.getA33(0,1) - hh2(0);
		v2diff = AHAT.getA31(1);
		v2diff -= vv2;
		g2diff = AHAT.getA33(1,0) - gg2(0);
		diffnorm = std::max<double>(sqrt(u2diff*u2diff+h2diff*h2diff),sqrt(v2diff*v2diff+g2diff*g2diff));
// 		std::cout<<"std::max(u2diff, v2diff) "<<diffnorm<<"\n";
		const double nrm_v = (1.0/sqrt(vv2*vv2+gg2(0)*gg2(0)));
		const double nrm_u = (1.0/sqrt(uu2*uu2+hh2(0)*hh2(0)));
		AHAT.getA31(1) = nrm_v * vv2;
		AHAT.getA13(1) = nrm_u * uu2;
		AHAT.getA33(1,0) = nrm_v * gg2(0);
		AHAT.getA33(0,1) = nrm_v * hh2(0);
	} while( (++it < NKERNITER)&&(diffnorm > KERNEPS) );
	if( diffnorm > KERNEPS ) std::cout<<"TestFunctLPAUT::Init: warning: No convergence in finding the singular vector. Residual = "<<diffnorm<<"\n";
}

double TestFunctLPAUT::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	col.InterpolateMSH( solMSHData, sol );
	col.CharJac_MSHphi( phi, par, solMSHData );
	
	col.Star( AHAT.getA31(0), phi );
	AHAT.getA13(0) = mB * phi;
	
	// continuing the non-trivial kernel
	AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
	AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
	
	const double nrm_v = (1.0/sqrt(vv2*vv2+gg2(0)*gg2(0)));
	const double nrm_u = (1.0/sqrt(uu2*uu2+hh2(0)*hh2(0)));
	AHAT.getA31(1) = nrm_v * vv2;
	AHAT.getA13(1) = nrm_u * uu2;
	AHAT.getA33(1,0) = nrm_v * gg2(0);
	AHAT.getA33(0,1) = nrm_u * hh2(0);
	
	// for subsequent use
	col.Interpolate( phiData, phi );
	col.Interpolate( vv2Data, vv2 );
	
// 	std::cout<<"TF: "<<gg2(1)<<", "<<hh2(1)<<":"<<gg2(0)<<", "<<hh2(0)<<"\n";
// 	if( gg2(1) > 0.0 ) std::cout<<"\t+++\n";
// 	else               std::cout<<"\t---\n";
	return gg2(1);
}

double TestFunctLPAUT::Funct_p( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData, int alpha )
{
	col.CharJac_x_p( A_p, par, solData, vv2Data, ZZ, alpha );
	col.CharJac_mB_p( mB_p, par, solData, phiData, ZZ, alpha );
	
	col.CharJac_MSHphi_p( DpPhi, par, solMSHData, alpha );
	// check
// 	col.CharJac_x_p( temp, par, solData, phiData, ZZ, alpha );
// 	std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
// 	temp = AHAT.getA11() * DpPhi;
// 	std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
// 	std::cout<<"d: "<<alpha<<", "<<DpPhi*DpPhi<<"\n";
	// end check
	temp = mB * DpPhi;
	temp += mB_p;
	
	const double gg_p = (uu2 * A_p) + gg2(0) * (uu2 * temp) + hh2(0) * col.Integrate( DpPhi, vv2 );
// 	std::cout<<"GP "<<alpha<<":"<<gg_p;
	return gg_p;
}

void   TestFunctLPAUT::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	col.CharJac_x_x( A_x, par, solData, vv2Data, ZZ );
	col.CharJac_mB_x( mB_x, par, solData, phiData, ZZ );
	func = !A_x * uu2;
	
	temp = !mB * uu2;
	col.CharJac_MSHphi_x<true>( DxPhi, par, solMSHData, temp );
	func += gg2(0) * DxPhi;
	temp = !mB_x * uu2;
	func += gg2(0) * temp;

	col.CharJac_MSHphi_x<true>( DxPhi, par, solMSHData, vv2 );
	col.Star( temp, DxPhi );
	func += hh2(0) * temp;
}

void   TestFunctLPAUT::Switch( Vector& phi )
{
	phi = vv2;
}

/// -----------------------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems with SIMMETRY
/// -----------------------------------------------------------------------

TestFunctLPAUTROT::TestFunctLPAUTROT( NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z ) :
	ZZ(Z),
	Re(CRe), Im(CIm),
	AHAT( NDIM*(NDEG*NINT+1), 0, 3, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB_p( NDIM*(NDEG*NINT+1) ),
	mB_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	phi( NDIM*(NDEG*NINT+1) ),
	DpPhi( NDIM*(NDEG*NINT+1) ),
	DxPhi( NDIM*(NDEG*NINT+1) ),
	LAM( NDIM*(NDEG*NINT+1) ),
	DxLAM( NDIM*(NDEG*NINT+1) ),
	uu3( NDIM*(NDEG*NINT+1) ),
	vv3( NDIM*(NDEG*NINT+1) ),
	gg3(3),
	hh3(3),
	rhs3( NDIM*(NDEG*NINT+1) ),
	one3(3),
	temp( NDIM*(NDEG*NINT+1) ),
	phiData( NDEG*NINT, NDIM, 2*NTAU+1 ),
	LAMData( NDEG*NINT, NDIM, 2*NTAU+1 ),
	vv3Data( NDEG*NINT, NDIM, 2*NTAU+1 ),
	solMSHData( NDEG*NINT+1, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunctLPAUTROT::~TestFunctLPAUTROT() 
{
}

void TestFunctLPAUTROT::Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	JagMatrix3D solMSHData( NDEG*NINT+1, NDIM, 2*NTAU );
	col.InterpolateMSH( solMSHData, sol );
	col.CharJac_MSHphi( phi, par, solMSHData );

	temp = AHAT.getA11() * phi;
	std::cout<<"temp: "<<temp*temp<<" phi: "<<phi*phi<<"\n";
	col.Star( AHAT.getA31(0), phi );
	AHAT.getA13(0) = mB * phi;

	rotbord<false>( LAM, col, sol, Re, Im );
	col.Star( AHAT.getA31(1), LAM );
	AHAT.getA13(1) = mB * LAM;
	
	AHAT.getA13(2).Rand();
	AHAT.getA31(2).Rand();
	// norming the borders
	AHAT.getA31(2) /= sqrt( AHAT.getA31(2) * AHAT.getA31(2) );
	AHAT.getA13(2) /= sqrt( AHAT.getA13(2) * AHAT.getA13(2) );
	AHAT.getA33().Clear(); // 2x2 matrix
	
	// generating the non-trivial kernel	
	one3(0) = 0.0;
	one3(1) = 0.0;
	one3(2) = 1.0;
	rhs3.Clear();
	
	Vector v3diff( vv3 ), u3diff( uu3 );
	double g30diff, h30diff, g31diff, h31diff;
	double diffnorm = 1.0;
	int it = 0;
	do {
		AHAT.Solve  ( 3, vv3, gg3, rhs3, one3 );
		AHAT.SolveTR( 3, uu3, hh3, rhs3, one3 );
		u3diff = AHAT.getA13(2);
		u3diff -= uu3;
		h30diff = AHAT.getA33(0,2) - hh3(0);
		h31diff = AHAT.getA33(1,2) - hh3(1);
		v3diff = AHAT.getA31(2);
		v3diff -= vv3;
		g30diff = AHAT.getA33(2,0) - gg3(0);
		g31diff = AHAT.getA33(2,1) - gg3(1);
		diffnorm = std::max<double>(sqrt(u3diff*u3diff+h30diff*h30diff+h31diff*h31diff),
		                    sqrt(v3diff*v3diff+g30diff*g30diff+g31diff*g31diff));
		std::cout<<"std::max(u3diff, v3diff) "<<diffnorm<<"\n";
		const double nrm_v = (1.0/sqrt(vv3*vv3+gg3(0)*gg3(0)+gg3(1)*gg3(1)));
		const double nrm_u = (1.0/sqrt(uu3*uu3+hh3(0)*hh3(0)+hh3(1)*hh3(1)));
		AHAT.getA31(2)   = nrm_v * vv3;
		AHAT.getA33(2,0) = nrm_v * gg3(0);
		AHAT.getA33(2,1) = nrm_v * gg3(1);
		AHAT.getA13(2)   = nrm_u * uu3;
		AHAT.getA33(0,2) = nrm_u * hh3(0);
		AHAT.getA33(1,2) = nrm_u * hh3(1);
	} while( (++it < NKERNITER)&&(diffnorm > KERNEPS) );
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";
}

double TestFunctLPAUTROT::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	col.InterpolateMSH( solMSHData, sol );
	col.CharJac_MSHphi( phi, par, solMSHData );
	
// 	temp = AHAT.getA11() * phi;
// 	std::cout<<"temp: "<<temp*temp<<" phi: "<<phi*phi<<"\n";
	
	col.Star( AHAT.getA31(0), phi );
	AHAT.getA13(0) = mB * phi;

	rotbord<false>( LAM, col, sol, Re, Im );
	col.Star( AHAT.getA31(1), LAM );
	AHAT.getA13(1) = mB * LAM;
	
// 	temp = AHAT.getA11() * LAM;
// 	std::cout<<"temp: "<<temp*temp<<" LAM: "<<LAM*LAM<<"\n";
	
	AHAT.Solve  ( 3, vv3, gg3, rhs3, one3 );
	AHAT.SolveTR( 3, uu3, hh3, rhs3, one3 );
	const double nrm_v = (1.0/sqrt(vv3*vv3+gg3(0)*gg3(0)+gg3(1)*gg3(1)));
	const double nrm_u = (1.0/sqrt(uu3*uu3+hh3(0)*hh3(0)+hh3(1)*hh3(1)));
	AHAT.getA31(2)   = nrm_v * vv3;
	AHAT.getA33(2,0) = nrm_v * gg3(0);
	AHAT.getA33(2,1) = nrm_v * gg3(1);
	AHAT.getA13(2)   = nrm_u * uu3;
	AHAT.getA33(0,2) = nrm_u * hh3(0);
	AHAT.getA33(1,2) = nrm_u * hh3(1);
	
	// for subsequent use
	col.Interpolate( phiData, phi );
	col.Interpolate( LAMData, LAM );
	col.Interpolate( vv3Data, vv3 );
// 	std::cout<<" TF1: "<<gg3(0)<<", "<<hh3(0)<<" TF2: "<<gg3(1)<<", "<<hh3(1)<<" TF3: "<<gg3(2)<<", "<<hh3(2)<<"\n";

	if( gg3(2) > 0.0 ) std::cout<<"\t+++\n";
	else               std::cout<<"\t---\n";
	return gg3(2);
}

double TestFunctLPAUTROT::Funct_p( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData, int alpha )
{
	col.CharJac_x_p( A_p, par, solData, vv3Data, ZZ, alpha );
	col.CharJac_MSHphi_p( DpPhi, par, solMSHData, alpha );
	double gg_p = (uu3 * A_p);
	// check
// 	col.CharJac_x_p( temp, par, solData, phiData, ZZ, alpha );
// 	std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
// 	temp = AHAT.getA11() * DpPhi;
// 	std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
// 	std::cout<<"d: "<<alpha<<", "<<DpPhi*DpPhi<<"\n";
	// end check
	// first phi
	temp = mB * DpPhi;
	col.CharJac_mB_p( mB_p, par, solData, phiData, ZZ, alpha );
	temp += mB_p;
	gg_p += gg3(0) * (uu3 * temp) + hh3(0) * col.Integrate( DpPhi, vv3 );
	// second LAM ( mB * LAM -> only mB_p is nonzero )
	col.CharJac_mB_p( mB_p, par, solData, LAMData, ZZ, alpha );
	gg_p += gg3(1) * (uu3 * mB_p);
	
// 	if( alpha==0 ) gg_p *= 22.0;
	std::cout<<"\nGP "<<alpha<<":"<<gg_p;
	return gg_p;
}

void   TestFunctLPAUTROT::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	col.Interpolate( vv3Data, vv3 );
	col.CharJac_x_x( A_x, par, solData, vv3Data, ZZ );
	func = !A_x * uu3;
	
	temp = !mB * uu3;
	col.CharJac_MSHphi_x<true>( DxPhi, par, solMSHData, temp );
	func += gg3(0) * DxPhi;
	col.CharJac_mB_x( mB_x, par, solData, phiData, ZZ );
	temp = !mB_x * uu3;
	func += gg3(0) * temp;
	
	col.CharJac_mB_x( mB_x, par, solData, LAMData, ZZ );
	temp = !mB_x * uu3;
	func += gg3(1) * temp;
	
	rotbord<true>( DxLAM, col, uu3, Re, Im );
	func += gg3(1) * DxLAM;
	
	col.CharJac_MSHphi_x<true>( DxPhi, par, solMSHData, vv3 );
	col.Star( temp, DxPhi );
	func += hh3(0) * temp;
	
	rotbord<true>( DxLAM, col, vv3, Re, Im );
	col.Star( temp, DxLAM );
	func += hh3(1) * temp;
}

void   TestFunctLPAUTROT::Switch( Vector& phi )
{
	phi = vv3;
}

/// -----------------------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems with SIMMETRY
/// EXTENDED
/// -----------------------------------------------------------------------

TestFunctLPAUTROT_X::TestFunctLPAUTROT_X( NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z ) :
	ZZ(Z),
	Re(CRe), Im(CIm),
	AHAT( NDIM*(NDEG*NINT+1), 0, 3, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB_p( NDIM*(NDEG*NINT+1) ),
	mB_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	phi( NDIM*(NDEG*NINT+1) ),
	DpPhi( NDIM*(NDEG*NINT+1) ),
	DxPhi( NDIM*(NDEG*NINT+1) ),
	LAM( NDIM*(NDEG*NINT+1) ),
	DxLAM( NDIM*(NDEG*NINT+1) ),
	uu3( NDIM*(NDEG*NINT+1) ),
	vv3( NDIM*(NDEG*NINT+1) ),
	gg3(3),
	hh3(3),
	one3(3),
	uu1( NDIM*(NDEG*NINT+1) ),
	vv1( NDIM*(NDEG*NINT+1) ),
	gg1(2),
	hh1(2),
	one1(2),
	uu2( NDIM*(NDEG*NINT+1) ),
	vv2( NDIM*(NDEG*NINT+1) ),
	gg2(2),
	hh2(2),
	one2(2),
	rhs( NDIM*(NDEG*NINT+1) ),
	temp( NDIM*(NDEG*NINT+1) ),
	vv1Data( NDEG*NINT, NDIM, 2*NTAU+1 ),
	vv2Data( NDEG*NINT, NDIM, 2*NTAU+1 ),
	vv3Data( NDEG*NINT, NDIM, 2*NTAU+1 ),
	solMSHData( NDEG*NINT+1, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunctLPAUTROT_X::~TestFunctLPAUTROT_X() 
{
}

void TestFunctLPAUTROT_X::Init( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	vv1.Rand();
	uu1.Rand();
	vv2.Rand();
	uu2.Rand();
	
	vv1 /= sqrt( vv1 * vv1 );
	uu1 /= sqrt( uu1 * uu1 );
	vv2 -= (vv1 * vv2) * vv1;
	uu2 -= (uu1 * uu2) * uu1;
	vv2 /= sqrt( vv2 * vv2 );
	uu2 /= sqrt( uu2 * uu2 );
	
	// looking for the singular vectors (geometric multiplicity)
	AHAT.getA13(0) = uu1;
	AHAT.getA13(1) = uu2;
	AHAT.getA31(0) = vv1;
	AHAT.getA31(1) = vv2;
	
	AHAT.getA33().Clear();
	rhs.Clear();
	one1(0) = 1.0; one1(1) = 0.0;
	one2(0) = 0.0; one2(1) = 1.0;
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 2, vv1, gg1, rhs, one1 );
		AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
		AHAT.SolveTR( 2, uu1, hh1, rhs, one1 );
		AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
		const double nrm_v1 = (1.0/sqrt(vv1*vv1));
		const double nrm_v2 = (1.0/sqrt(vv2*vv2));
		const double nrm_u1 = (1.0/sqrt(uu1*uu1));
		const double nrm_u2 = (1.0/sqrt(uu2*uu2));
		AHAT.getA13(0) = nrm_u1 * uu1;
		AHAT.getA13(1) = nrm_u2 * uu2;
		AHAT.getA31(0) = nrm_v1 * vv1;
		AHAT.getA31(1) = nrm_v2 * vv2;
	}

	AHAT.getA13(0) = mB * vv1;
	AHAT.getA13(1) = mB * vv2;
	
	vv3.Rand();
	uu3.Rand();
	vv3 /= sqrt( vv1 * vv1 );
	uu3 /= sqrt( uu1 * uu1 );
	AHAT.getA13(2) = uu3;
	AHAT.getA31(2) = vv3;
	// generating the non-trivial kernel	
	one3(0) = 0.0; one3(1) = 0.0; one3(2) = 1.0;
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 3, vv3, gg3, rhs, one3 );
		AHAT.SolveTR( 3, uu3, hh3, rhs, one3 );
		const double nrm_v = (1.0/sqrt(vv3*vv3+gg3(0)*gg3(0)+gg3(1)*gg3(1)));
		const double nrm_u = (1.0/sqrt(uu3*uu3+hh3(0)*hh3(0)+hh3(1)*hh3(1)));
		AHAT.getA31(2)   = nrm_v * vv3;
		AHAT.getA33(2,0) = nrm_v * gg3(0);
		AHAT.getA33(2,1) = nrm_v * gg3(1);
		AHAT.getA13(2)   = nrm_u * uu3;
		AHAT.getA33(0,2) = nrm_u * hh3(0);
		AHAT.getA33(1,2) = nrm_u * hh3(1);
	}
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";
}

double TestFunctLPAUTROT_X::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	AHAT.getA13(0) = uu1;
	AHAT.getA13(1) = uu2;
	
	AHAT.Solve  ( 2, vv1, gg1, rhs, one1 );
	AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
	AHAT.SolveTR( 2, uu1, hh1, rhs, one1 );
	AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
	const double nrm_v1 = (1.0/sqrt(vv1*vv1));
	const double nrm_v2 = (1.0/sqrt(vv2*vv2));
	const double nrm_u1 = (1.0/sqrt(uu1*uu1));
	const double nrm_u2 = (1.0/sqrt(uu2*uu2));
	AHAT.getA13(0) = nrm_u1 * uu1;
	AHAT.getA13(1) = nrm_u2 * uu2;
	AHAT.getA31(0) = nrm_v1 * vv1;
	AHAT.getA31(1) = nrm_v2 * vv2;

	AHAT.getA13(0) = mB * vv1;
	AHAT.getA13(1) = mB * vv2;
	
// 	temp = AHAT.getA11() * LAM;
// 	std::cout<<"temp: "<<temp*temp<<" LAM: "<<LAM*LAM<<"\n";
	
	AHAT.Solve  ( 3, vv3, gg3, rhs, one3 );
	AHAT.SolveTR( 3, uu3, hh3, rhs, one3 );
	const double nrm_v = (1.0/sqrt(vv3*vv3+gg3(0)*gg3(0)+gg3(1)*gg3(1)));
	const double nrm_u = (1.0/sqrt(uu3*uu3+hh3(0)*hh3(0)+hh3(1)*hh3(1)));
	AHAT.getA31(2)   = nrm_v * vv3;
	AHAT.getA33(2,0) = nrm_v * gg3(0);
	AHAT.getA33(2,1) = nrm_v * gg3(1);
	AHAT.getA13(2)   = nrm_u * uu3;
	AHAT.getA33(0,2) = nrm_u * hh3(0);
	AHAT.getA33(1,2) = nrm_u * hh3(1);
	
	// for subsequent use
	col.Interpolate( vv1Data, vv1 );
	col.Interpolate( vv2Data, vv2 );
	col.Interpolate( vv3Data, vv3 );

	if( gg3(2) > 0.0 ) std::cout<<"\t+++\n";
	else               std::cout<<"\t---\n";
	return gg3(2);
}

double TestFunctLPAUTROT_X::Funct_p( NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData, int alpha )
{
	col.CharJac_x_p( A_p, par, solData, vv3Data, ZZ, alpha );
	col.CharJac_mB_p( mB_p, par, solData, vv1Data, ZZ, alpha );
	temp = gg3(0) * mB_p;
	col.CharJac_mB_p( mB_p, par, solData, vv2Data, ZZ, alpha );
	temp += gg3(1) * mB_p;
	
	double gg_p = (uu3 * A_p) + (uu3 * temp);
	std::cout<<"\nGP "<<alpha<<":"<<gg_p;
	return gg_p;
}

void   TestFunctLPAUTROT_X::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& /*sol*/, const JagMatrix3D& solData )
{
	col.CharJac_x_x( A_x, par, solData, vv3Data, ZZ );
	func = !A_x * uu3;
	col.CharJac_mB_x( mB_x, par, solData, vv1Data, ZZ );
	temp = !mB_x * uu3;
	func += gg3(0) * temp;
	col.CharJac_mB_x( mB_x, par, solData, vv2Data, ZZ );
	temp = !mB_x * uu3;
	func += gg3(1) * temp;
}

void   TestFunctLPAUTROT_X::Switch( Vector& phi )
{
	phi = vv3;
}
