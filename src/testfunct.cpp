// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "testfunct.h"
#include "matrix.h"
#include "spmatrix.h"
#include "hypermatrix.h"
#include "ncolloc.h"

#define NKERNITER 10

#define NDIM (col.Ndim())
#define NDEG (col.Ndeg())
#define NINT (col.Nint())
#define NTAU (col.Ntau())
#define NPAR (col.Npar())

template<bool trans> void inline rotbord( Vector& V, NColloc& col, const Vector& IN, Array1D<int>& Re, Array1D<int>& Im )
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

void inline conjugate( Vector& V, NColloc& col, const Vector& IN )
{
	for( int idx = 0; idx < NINT*NDEG+1; idx++ )
	{
		for( int k = 0; k < NDIM; k++ )
		{
			V( 2*(k+NDIM*idx) )   = -IN( 2*(k+NDIM*idx)+1 );
			V( 2*(k+NDIM*idx)+1 ) =  IN( 2*(k+NDIM*idx) );
		}
	}
}

TestFunct::TestFunct( NColloc& col, double Z ) :
	ZZ(Z),
	AHAT( NDIM*(NDEG*NINT+1), 0, 1, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	rhs( NDIM*(NDEG*NINT+1) ),
	uu( NDIM*(NDEG*NINT+1) ),
	vv( NDIM*(NDEG*NINT+1) ),
	vvData( NDEG*NINT, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunct::~TestFunct()
{
}

void TestFunct::Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	AHAT.getA13(0).Rand();
	AHAT.getA31(0).Rand();
	AHAT.getA33()(0,0) = 0.0;
	// norming the borders
	AHAT.getA31(0) /= sqrt( AHAT.getA31(0) * AHAT.getA31(0) );
	AHAT.getA13(0) /= sqrt( AHAT.getA13(0) * AHAT.getA13(0) );
	
	double one = 1.0;
	double gg;
	double hh;
	
	rhs.Clear();
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( vv, gg, rhs, one );
		AHAT.SolveTR( uu, hh, rhs, one );
		AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
		AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
	}
// 	std::cout<<"TF: "<<gg<<", "<<hh<<"\n";
}

double TestFunct::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	double one = 1.0;
	double gg;
	double hh;
	
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	
	AHAT.Solve  ( vv, gg, rhs, one );
	AHAT.SolveTR( uu, hh, rhs, one );
	AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
	AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
// 	std::cout<<"TF: "<<gg<<", "<<hh;
// 	if( gg > 0.0 ) std::cout<<"\t+++\n";
// 	else           std::cout<<"\t---\n";
	return gg;
}

double TestFunct::Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha )
{
	col.Interpolate( vvData, vv );
	col.CharJac_x_p( A_p, par, solData, vvData, ZZ, alpha );
	const double gg_p = (uu * A_p);
// 	std::cout<<"GP "<<alpha<<":"<<gg_p;
	return gg_p; // positive, because the test function is not negated in the Jacobian
}

void   TestFunct::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	col.Interpolate( vvData, vv );
	col.CharJac_x_x( A_x, par, solData, vvData, ZZ );
	func = !A_x * uu; // positive, because the test function is not negated in the Jacobian
// 	func *= 100.0;
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

void TestFunctCPLX::Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData,
									double Re, double Im )
{
	ZRe = Re; ZIm = Im;
	col.CharJac_x( AHAT.getA11(), par, solData, Re, Im, true );
	AHAT.getA13(0).Rand();
	AHAT.getA31(0).Rand();
	AHAT.getA33().Clear();
	// norming the borders
	AHAT.getA31(0) /= sqrt( AHAT.getA31(0) * AHAT.getA31(0) );
	AHAT.getA13(0) /= sqrt( AHAT.getA13(0) * AHAT.getA13(0) );
	// conjugate
	conjugate( AHAT.getA31(1), col, AHAT.getA31(0) );
	conjugate( AHAT.getA13(1), col, AHAT.getA13(0) );
	
	rhs.Clear();
	one(0) = 1.0; one(1) = 0.0;
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 2, vv, gg, rhs, one );
		AHAT.SolveTR( 2, uu, hh, rhs, one );
		AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
		AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
		conjugate( AHAT.getA13(1), col, AHAT.getA13(0) );
		conjugate( AHAT.getA31(1), col, AHAT.getA31(0) );
	}
// 	std::cout<<"TF: "<<gg<<", "<<hh<<"\n";
}

void TestFunctCPLX::Funct  ( double& f1, double& f2,
					NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, double Re, double Im )
{
	if( first ){ Init( col, par, sol, solData, Re, Im ); first = false; }
	
	ZRe = Re; ZIm = Im;
	col.CharJac_x( AHAT.getA11(), par, solData, Re, Im, true );
	AHAT.Solve  ( 2, vv, gg, rhs, one );
	AHAT.SolveTR( 2, uu, hh, rhs, one );
	AHAT.getA13(0) = (1.0/sqrt(uu*uu))*uu;
	AHAT.getA31(0) = (1.0/sqrt(vv*vv))*vv;
	conjugate( AHAT.getA13(1), col, AHAT.getA13(0) );
	conjugate( AHAT.getA31(1), col, AHAT.getA31(0) );
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
					NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData,
					int alpha )
{
	col.CharJac_x_p( A_p, par, solData, vvData, ZRe, ZIm, alpha );
	f1 = (AHAT.getA13(0)/*uu*/ * A_p);
	f2 = (AHAT.getA13(1)/*uu^conj*/ * A_p);
}

void TestFunctCPLX::Funct_z( double& f1, double& f2,
					NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	col.CharJac_x_z( A_p, par, solData, AHAT.getA31(0), vvData, ZRe, ZIm, true );
	const double dzre = (AHAT.getA13(0)/*uu*/ * A_p);
	const double dzim = (AHAT.getA13(1)/*uu^conj*/ * A_p);
	f1 = (- dzre * ZIm - dzim * ZRe);
	f2 = (  dzre * ZRe - dzim * ZIm);
}

void TestFunctCPLX::Funct_x( Vector& func1, Vector& func2,
					NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	col.CharJac_x_x( A_x, par, solData, vvData, ZRe, ZIm );
	func1 = !A_x * AHAT.getA13(0); /*uu*/
	func2 = !A_x * AHAT.getA13(1); /*uu^conj*/
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
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
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
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
		AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
		const double nrm_v = (1.0/sqrt(vv2*vv2+gg2(0)*gg2(0)));
		const double nrm_u = (1.0/sqrt(uu2*uu2+hh2(0)*hh2(0)));
		AHAT.getA31(1) = nrm_v * vv2;
		AHAT.getA13(1) = nrm_u * uu2;
		AHAT.getA33(1,0) = nrm_v * gg2(0);
		AHAT.getA33(0,1) = nrm_v * hh2(0);
	}
}

double TestFunctLPAUT::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
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

double TestFunctLPAUT::Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha )
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

void   TestFunctLPAUT::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
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
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
	JagMatrix3D solMSHData( NDEG*NINT+1, NDIM, 2*NTAU );
	col.InterpolateMSH( solMSHData, sol );
	col.CharJac_MSHphi( phi, par, solMSHData );

	temp = AHAT.getA11() * phi;
	std::cout<<"temp: "<<temp*temp<<" phi: "<<phi*phi<<"\n";
	col.Star( AHAT.getA31(0), phi );
	AHAT.getA13(0) = mB * phi;

	rotbord<false>( LAM, col, sol, Re, Im );
	col.Star( AHAT.getA31(1), LAM );
	AHAT.getA13(1) = LAM;
	
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
	
	for( int i = 0; i < NKERNITER; i++ )
	{
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
	}
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";
}

double TestFunctLPAUTROT::Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, sol, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
	col.InterpolateMSH( solMSHData, sol );
	col.CharJac_MSHphi( phi, par, solMSHData );
	
// 	temp = AHAT.getA11() * phi;
// 	std::cout<<"temp: "<<temp*temp<<" phi: "<<phi*phi<<"\n";
	
	col.Star( AHAT.getA31(0), phi );
	AHAT.getA13(0) = mB * phi;

	rotbord<false>( LAM, col, sol, Re, Im );
	col.Star( AHAT.getA31(1), LAM );
	AHAT.getA13(1) = LAM;
	
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
	col.Interpolate( vv3Data, vv3 );
// 	std::cout<<" TF1: "<<gg3(0)<<", "<<hh3(0)<<" TF2: "<<gg3(1)<<", "<<hh3(1)<<" TF3: "<<gg3(2)<<", "<<hh3(2)<<"\n";

	if( gg3(2) > 0.0 ) std::cout<<"\t+++\n";
	else               std::cout<<"\t---\n";
	return gg3(2);
}

double TestFunctLPAUTROT::Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha )
{
	col.CharJac_x_p( A_p, par, solData, vv3Data, ZZ, alpha );
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
	
	double gg_p = (uu3 * A_p) + gg3(0) * (uu3 * temp) + hh3(0) * col.Integrate( DpPhi, vv3 );
	if( alpha==0 ) gg_p = -gg3(0);
	std::cout<<"\nGP "<<alpha<<":"<<gg_p;
	return gg_p;
}

void   TestFunctLPAUTROT::Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData )
{
	col.Interpolate( vv3Data, vv3 );
	col.CharJac_x_x( A_x, par, solData, vv3Data, ZZ );
	col.CharJac_mB_x( mB_x, par, solData, phiData, ZZ );
	func = !A_x * uu3;
	
	temp = !mB * uu3;
	col.CharJac_MSHphi_x<true>( DxPhi, par, solMSHData, temp );
	func += gg3(0) * DxPhi;
	temp = !mB_x * uu3;
	func += gg3(0) * temp;
	
	rotbord<true>( DxLAM, col, uu3, Re, Im );
	func += gg3(1) * DxLAM;
	
	col.CharJac_MSHphi_x<true>( DxPhi, par, solMSHData, vv3 );
	col.Star( temp, DxPhi );
	func += hh3(0) * temp;
	
	rotbord<true>( DxLAM, col, vv3, Re, Im );
	col.Star( temp, DxLAM );
	func += hh3(1) * temp;
}
