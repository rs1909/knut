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

void TestFunct::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData )
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

double TestFunct::Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	double one = 1.0;
	double gg;
	double hh;
	
	if( first ){ Init( col, par, solData ); first = false; }
	
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

double TestFunct::Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha )
{
	col.Interpolate( vvData, vv );
	col.CharJac_x_p( A_p, par, solData, vvData, ZZ, alpha );
	const double gg_p = (uu * A_p);
// 	std::cout<<"GP "<<alpha<<":"<<gg_p;
	return gg_p; // positive, because the test function is not negated in the Jacobian
}

void   TestFunct::Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	col.Interpolate( vvData, vv );
	col.CharJac_x_x( A_x, par, solData, vvData, ZZ );
	func = !A_x * uu; // positive, because the test function is not negated in the Jacobian
// 	func *= 100.0;
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
	rhs( NDIM*(NDEG*NINT+1) ),
	uu1( NDIM*(NDEG*NINT+1) ),
	vv1( NDIM*(NDEG*NINT+1) ),
	uu2( NDIM*(NDEG*NINT+1) ),
	vv2( NDIM*(NDEG*NINT+1) ),
	gg2(2),
	hh2(2),
	one2(2),
	vv2Data( NDEG*NINT, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunctLPAUT::~TestFunctLPAUT() 
{
}

void TestFunctLPAUT::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	AHAT.getA13(0).Rand();
	AHAT.getA31(0).Rand();
	AHAT.getA13(1).Rand();
	AHAT.getA31(1).Rand();
	AHAT.getA33().Clear(); // 2x2 matrix
	// norming the borders
	AHAT.getA31(0) /= sqrt( AHAT.getA31(0) * AHAT.getA31(0) );
	AHAT.getA13(0) /= sqrt( AHAT.getA13(0) * AHAT.getA13(0) );
	AHAT.getA31(1) /= sqrt( AHAT.getA31(1) * AHAT.getA31(1) );
	AHAT.getA13(1) /= sqrt( AHAT.getA13(1) * AHAT.getA13(1) );

	// generating the trivial kernel
	double one1 = 1.0;
	double gg1;
	double hh1;
	
	rhs.Clear();
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( vv1, gg1, rhs, one1 );
		AHAT.SolveTR( uu1, hh1, rhs, one1 );
		AHAT.getA13(0) = (1.0/sqrt(uu1*uu1))*uu1;
		AHAT.getA31(0) = (1.0/sqrt(vv1*vv1))*vv1;
	}
// 	std::cout<<"T0: "<<gg1<<", "<<hh1<<"\n ";
	
	// generating the non-trivial kernel	
	one2(0) = 0.0;
	one2(1) = 1.0;
	
	AHAT.getA13(0) = mB * vv1;
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
		AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
		const double nrm_v = (1.0/sqrt(vv2*vv2));
		const double nrm_u = (1.0/sqrt(uu2*uu2));
		AHAT.getA31(1) = nrm_v * vv2;
		AHAT.getA13(1) = nrm_u * uu2;
	}
// 	std::cout<<"TF: "<<gg2(1)<<", "<<hh2(1)<<"\n";
}

double TestFunctLPAUT::Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
	// continuing the trivial kernel
	double one1 = 1.0;
	double gg1;
	double hh1;
	
	AHAT.getA13(0) = vv1;
	AHAT.Solve  ( vv1, gg1, rhs, one1 );
	AHAT.SolveTR( uu1, hh1, rhs, one1 );
	AHAT.getA13(0) = (1.0/sqrt(uu1*uu1))*uu1;
	AHAT.getA31(0) = (1.0/sqrt(vv1*vv1))*vv1;
// 	std::cout<<"T0: "<<gg1<<", "<<hh1;
	
	// continuing the non-trivial kernel
	AHAT.getA13(0) = mB * vv1;
	AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
	AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
	
	const double nrm_v = (1.0/sqrt(vv2*vv2));
	const double nrm_u = (1.0/sqrt(uu2*uu2));
	AHAT.getA31(1) = nrm_v * vv2;
	AHAT.getA13(1) = nrm_u * uu2;
	
// 	std::cout<<"TF: "<<gg2(1)<<", "<<hh2(1)<<"\n";
// 	if( gg2(1) > 0.0 ) std::cout<<"\t+++\n";
// 	else               std::cout<<"\t---\n";
	return gg2(1);
}

double TestFunctLPAUT::Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha )
{
// 	JagMatrix3D vv1Data( NDEG*NINT, NDIM, 2*NTAU+1 );
// 	col.Interpolate( vv1Data, vv1 );
	col.Interpolate( vv2Data, vv2 );
	Vector mB_p( NDIM*(NDEG*NINT+1) );
	col.CharJac_x_p( A_p, par, solData, vv2Data, ZZ, alpha );
// 	col.CharJac_mB_p( mB_p, par, solData, vv1Data, ZZ, alpha );
	const double gg_p = (uu2 * A_p) /*+ gg2(0) * (uu2 * mB_p)*/;
// 	std::cout<<"GP "<<alpha<<":"<<gg_p;
	return gg_p; // positive, because the test function is not negated in the Jacobian
}

void   TestFunctLPAUT::Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
// 	JagMatrix3D vv1Data( NDEG*NINT, NDIM, 2*NTAU+1 );
// 	col.Interpolate( vv1Data, vv1 );
	col.Interpolate( vv2Data, vv2 );
	col.CharJac_x_x( A_x, par, solData, vv2Data, ZZ );
// 	SpMatrix mB_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) );
// 	col.CharJac_mB_x( mB_x, par, solData, vv1Data, ZZ );
	func = !A_x * uu2;
// 	Vector temp( NDIM*(NDEG*NINT+1) );
// 	temp = !mB_x * uu2;
// 	func += gg2(0) * temp;
}

/// -----------------------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems with SIMMETRY
/// -----------------------------------------------------------------------

TestFunctLPAUTROT::TestFunctLPAUTROT( NColloc& col, double Z ) :
	ZZ(Z),
	AHAT( NDIM*(NDEG*NINT+1), 0, 3, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	rhs( NDIM*(NDEG*NINT+1) ),
	uu1( NDIM*(NDEG*NINT+1) ),
	vv1( NDIM*(NDEG*NINT+1) ),
	uu2( NDIM*(NDEG*NINT+1) ),
	vv2( NDIM*(NDEG*NINT+1) ),
	gg2(2),
	hh2(2),
	one2(2),
	uu3( NDIM*(NDEG*NINT+1) ),
	vv3( NDIM*(NDEG*NINT+1) ),
	gg3(3),
	hh3(3),
	one3(3),
	vv3Data( NDEG*NINT, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunctLPAUTROT::~TestFunctLPAUTROT() 
{
}

void TestFunctLPAUTROT::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	AHAT.getA13(0).Rand();
	AHAT.getA31(0).Rand();
	AHAT.getA13(1).Rand();
	AHAT.getA31(1).Rand();
	AHAT.getA13(2).Rand();
	AHAT.getA31(2).Rand();
	AHAT.getA33().Clear(); // 2x2 matrix
	// norming the borders
	AHAT.getA31(0) /= sqrt( AHAT.getA31(0) * AHAT.getA31(0) );
	AHAT.getA13(0) /= sqrt( AHAT.getA13(0) * AHAT.getA13(0) );
	AHAT.getA31(1) /= sqrt( AHAT.getA31(1) * AHAT.getA31(1) );
	AHAT.getA13(1) /= sqrt( AHAT.getA13(1) * AHAT.getA13(1) );
	AHAT.getA31(2) /= sqrt( AHAT.getA31(2) * AHAT.getA31(2) );
	AHAT.getA13(2) /= sqrt( AHAT.getA13(2) * AHAT.getA13(2) );

	// generating the trivial kernel
	double one1 = 1.0;
	double gg1;
	double hh1;
	
	rhs.Clear();
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( vv1, gg1, rhs, one1 );
		AHAT.SolveTR( uu1, hh1, rhs, one1 );
		AHAT.getA13(0) = (1.0/sqrt(uu1*uu1))*uu1;
		AHAT.getA31(0) = (1.0/sqrt(vv1*vv1))*vv1;
	}
	std::cout<<"T0: "<<gg1<<", "<<hh1<<"\n ";

	// generating the second trivial kernel
	one2(0) = 0.0;
	one2(1) = 1.0;
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
		AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
		const double nrm_v = (1.0/sqrt(vv2*vv2));
		const double nrm_u = (1.0/sqrt(uu2*uu2));
		AHAT.getA31(1) = nrm_v * vv2;
		AHAT.getA13(1) = nrm_u * uu2;
	}
	std::cout<<"T1: "<<gg2(1)<<", "<<hh2(1)<<"\n ";
	
	// generating the non-trivial kernel	
	one3(0) = 0.0;
	one3(1) = 0.0;
	one3(2) = 1.0;
	
	AHAT.getA13(0) = mB * vv1;
	AHAT.getA13(1) = mB * vv2;
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 3, vv3, gg3, rhs, one3 );
		AHAT.SolveTR( 3, uu3, hh3, rhs, one3 );
		const double nrm_v = (1.0/sqrt(vv3*vv3));
		const double nrm_u = (1.0/sqrt(uu3*uu3));
		AHAT.getA31(2) = nrm_v * vv3;
		AHAT.getA13(2) = nrm_u * uu3;
	}
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";
}

double TestFunctLPAUTROT::Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
	// continuing the trivial kernel
	double one1 = 1.0;
	double gg1;
	double hh1;
	
	AHAT.getA13(0) = vv1;
	AHAT.getA13(1) = vv2;
	AHAT.Solve  ( vv1, gg1, rhs, one1 );
	AHAT.SolveTR( uu1, hh1, rhs, one1 );
	AHAT.getA13(0) = (1.0/sqrt(uu1*uu1))*uu1;
	AHAT.getA31(0) = (1.0/sqrt(vv1*vv1))*vv1;
	std::cout<<"T0: "<<gg1<<", "<<hh1;
	
	// continuing the second trivial kernel
	AHAT.Solve  ( 2, vv2, gg2, rhs, one2 );
	AHAT.SolveTR( 2, uu2, hh2, rhs, one2 );
	const double nrm_v = (1.0/sqrt(vv2*vv2));
	const double nrm_u = (1.0/sqrt(uu2*uu2));
	AHAT.getA31(1) = nrm_v * vv2;
	AHAT.getA13(1) = nrm_u * uu2;
	std::cout<<"TF: "<<gg2(1)<<", "<<hh2(1)<<"\n";

	AHAT.getA13(0) = mB * vv1;
	AHAT.getA13(1) = mB * vv2;
	AHAT.Solve  ( 3, vv3, gg3, rhs, one3 );
	AHAT.SolveTR( 3, uu3, hh3, rhs, one3 );
	const double nrm3_v = (1.0/sqrt(vv3*vv3));
	const double nrm3_u = (1.0/sqrt(uu3*uu3));
	AHAT.getA31(2) = nrm3_v * vv3;
	AHAT.getA13(2) = nrm3_u * uu3;
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";

	if( gg3(2) > 0.0 ) std::cout<<"\t+++\n";
	else               std::cout<<"\t---\n";
	return gg3(2);
}

double TestFunctLPAUTROT::Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha )
{
	col.Interpolate( vv3Data, vv3 );
	col.CharJac_x_p( A_p, par, solData, vv3Data, ZZ, alpha );
	const double gg_p = (uu3 * A_p);
	std::cout<<"GP "<<alpha<<":"<<gg_p;
	if( alpha==0 ) return 22*gg_p;
	else return gg_p;
}

void   TestFunctLPAUTROT::Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	col.Interpolate( vv3Data, vv3 );
	col.CharJac_x_x( A_x, par, solData, vv3Data, ZZ );
	func = (1.0)*!A_x * uu3; // positive, because the test function is not negated in the Jacobian
}

/// -----------------------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems with SIMMETRY  Second try
/// -----------------------------------------------------------------------

TestFunctLPAUTROT_S::TestFunctLPAUTROT_S( NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z ) :
	ZZ(Z),
	AHAT( NDIM*(NDEG*NINT+1), 0, 3, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	A_p( NDIM*(NDEG*NINT+1) ),
	A_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	rhs( NDIM*(NDEG*NINT+1) ),
	uu1( NDIM*(NDEG*NINT+1) ),
	vv1( NDIM*(NDEG*NINT+1) ),
	uu2( NDIM*(NDEG*NINT+1) ),
	vv2( NDIM*(NDEG*NINT+1) ),
	gg2(2),
	hh2(2),
	one2(2),
	uu3( NDIM*(NDEG*NINT+1) ),
	vv3( NDIM*(NDEG*NINT+1) ),
	gg3(3),
	hh3(3),
	one3(3),
	Re(CRe), Im(CIm),
	vv3Data( NDEG*NINT, NDIM, 2*NTAU+1 )
{
	first = true;
}

TestFunctLPAUTROT_S::~TestFunctLPAUTROT_S() 
{
}

void inline rotbord( Vector& V, NColloc& col, const JagMatrix3D& solData, Array1D<int>& Re, Array1D<int>& Im )
{
	V.Clear();
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			for( int k = 0; k < Re.Size(); k++ )
			{
				V( NDIM + Re(k) + NDIM*(j+NDEG*i) ) = solData( idx )( Im(k), 0 );
				V( NDIM + Im(k) + NDIM*(j+NDEG*i) ) = - solData( idx )( Re(k), 0 );
			}
		}
	}
}

template<bool trans> void inline rotbord( Vector& V, NColloc& col, Vector& IN, Array1D<int>& Re, Array1D<int>& Im )
{
	V.Clear();
	for( int i = 0; i < NINT; i++ )  // i: interval; j: which collocation point
	{
		for( int j = 0; j < NDEG; j++ )
		{
			const int idx = j+i*NDEG;
			for( int k = 0; k < Re.Size(); k++ )
			{
				if(trans)
				{
					V( NDIM + Re(k) + NDIM*(j+NDEG*i) ) = - IN( NDIM + Im(k) + NDIM*(j+NDEG*i) );
					V( NDIM + Im(k) + NDIM*(j+NDEG*i) ) = IN( NDIM + Re(k) + NDIM*(j+NDEG*i) );
				}else
				{
					V( NDIM + Re(k) + NDIM*(j+NDEG*i) ) = IN( NDIM + Im(k) + NDIM*(j+NDEG*i) );
					V( NDIM + Im(k) + NDIM*(j+NDEG*i) ) = - IN( NDIM + Re(k) + NDIM*(j+NDEG*i) );
				}
			}
		}
	}
}

void TestFunctLPAUTROT_S::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	// creating the matrix
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
	col.CharJac_phi( vv1, par, solData );
	col.Star( AHAT.getA31(0), vv1 );
	AHAT.getA13(0) = mB * vv1;

	rotbord( vv2, col, solData, Re, Im );
	col.Star( AHAT.getA31(1), vv2 );
	AHAT.getA13(1) = mB * vv2;
	
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
	
// 	AHAT.getA13(0) = mB * vv1;
// 	AHAT.getA13(1) = mB * vv2;
	
	for( int i = 0; i < NKERNITER; i++ )
	{
		AHAT.Solve  ( 3, vv3, gg3, rhs, one3 );
		AHAT.SolveTR( 3, uu3, hh3, rhs, one3 );
		const double nrm_v = (1.0/sqrt(vv3*vv3));
		const double nrm_u = (1.0/sqrt(uu3*uu3));
		AHAT.getA31(2) = nrm_v * vv3;
		AHAT.getA13(2) = nrm_u * uu3;
	}
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";
}

double TestFunctLPAUTROT_S::Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	if( first ){ Init( col, par, solData ); first = false; }
	
	col.CharJac_x( AHAT.getA11(), par, solData, ZZ, true );
	col.CharJac_mB( mB, par, solData, ZZ, true );
	
	col.CharJac_phi( vv1, par, solData );
	col.Star( AHAT.getA31(0), vv1 );
	AHAT.getA13(0) = mB * vv1;

	rotbord( vv2, col, solData, Re, Im );
	col.Star( AHAT.getA31(1), vv2 );
	AHAT.getA13(1) = mB * vv2;
	
	AHAT.Solve  ( 3, vv3, gg3, rhs, one3 );
	AHAT.SolveTR( 3, uu3, hh3, rhs, one3 );
	const double nrm3_v = (1.0/sqrt(vv3*vv3));
	const double nrm3_u = (1.0/sqrt(uu3*uu3));
	AHAT.getA31(2) = nrm3_v * vv3;
	AHAT.getA13(2) = nrm3_u * uu3;
	std::cout<<"TF: "<<gg3(2)<<", "<<hh3(2)<<"\n";

	if( gg3(2) > 0.0 ) std::cout<<"\t+++\n";
	else               std::cout<<"\t---\n";
	return gg3(2);
}

double TestFunctLPAUTROT_S::Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha )
{
	col.Interpolate( vv3Data, vv3 );
	col.CharJac_x_p( A_p, par, solData, vv3Data, ZZ, alpha );
	double gg_p = (uu3 * A_p);
	std::cout<<"\nGP0 "<<alpha<<":"<<gg_p;
	
	Vector phiG0pLAMxG1( NDIM*(NDEG*NINT+1) );
	phiG0pLAMxG1 =  gg3(0) * vv1;/*phi*/
	phiG0pLAMxG1 += gg3(1) * vv2;/*LAMBDA * v*/
	JagMatrix3D phiG0pLAMxG1Data( NDEG*NINT, NDIM, 2*NTAU+1 );
	col.Interpolate( phiG0pLAMxG1Data, phiG0pLAMxG1 );
	Vector DxmBphiG0pLAMxG1( NDIM*(NDEG*NINT+1) );
	col.CharJac_mB_p( DxmBphiG0pLAMxG1, par, solData, phiG0pLAMxG1Data, ZZ, alpha );
	gg_p += uu3 * DxmBphiG0pLAMxG1;
	std::cout<<"\nGP1 "<<alpha<<":"<<gg_p;
	
	Vector wSmB( NDIM*(NDEG*NINT+1) );
	wSmB = !mB * uu3;
	Vector wSmBDxphi( NDIM*(NDEG*NINT+1) );
	col.CharJac_phi_p( wSmBDxphi, par, solData, alpha );
	gg_p += gg3(0) * (wSmBDxphi * wSmB);
	std::cout<<"\nGP2 "<<alpha<<":"<<gg_p;
	
	Vector Dxphiv( NDIM*(NDEG*NINT+1) );
	col.CharJac_phi_p( Dxphiv, par, solData, alpha );
	/*if (alpha!=0)*/ gg_p += hh3(0) * col.Integrate(Dxphiv, vv3);

	std::cout<<"\nGP "<<alpha<<":"<<gg_p;

	if(alpha==0) return 22*gg_p;
	else return gg_p;
}

void   TestFunctLPAUTROT_S::Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	col.Interpolate( vv3Data, vv3 );
	col.CharJac_x_x( A_x, par, solData, vv3Data, ZZ );
	func = !A_x * uu3;
	
	Vector phiG0pLAMxG1( NDIM*(NDEG*NINT+1) );
	phiG0pLAMxG1 =  gg3(0) * vv1;/*phi*/
	phiG0pLAMxG1 += gg3(1) * vv2;/*LAMBDA * v*/
	JagMatrix3D phiG0pLAMxG1Data( NDEG*NINT, NDIM, 2*NTAU+1 );
	col.Interpolate( phiG0pLAMxG1Data, phiG0pLAMxG1 );
	SpMatrix DxmBphiG0pLAMxG1( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) );
	col.CharJac_mB_x( DxmBphiG0pLAMxG1, par, solData, phiG0pLAMxG1Data, ZZ );
	phiG0pLAMxG1 = !DxmBphiG0pLAMxG1 * uu3;
	func += phiG0pLAMxG1;
	
	Vector wSmB( NDIM*(NDEG*NINT+1) );
	wSmB = !mB * uu3;
	Vector wSmBDxphi( NDIM*(NDEG*NINT+1) );
	col.CharJac_phi_x<true>( wSmBDxphi, par, solData, wSmB );
	func += gg3(0)*wSmBDxphi;
	
	Vector wSmBLAM( NDIM*(NDEG*NINT+1) );
	rotbord<true>( wSmBLAM, col, wSmB, Re, Im );
	func += gg3(1)*wSmBLAM;
	
	Vector Dxphiv( NDIM*(NDEG*NINT+1) );
	Vector DxphivS( NDIM*(NDEG*NINT+1) );
	col.CharJac_phi_x<true>( Dxphiv, par, solData, vv3 );
	col.Star( DxphivS, Dxphiv );
	func += hh3(0) * DxphivS;
	
	Vector LAMv( NDIM*(NDEG*NINT+1) );
	Vector LAMvS( NDIM*(NDEG*NINT+1) );
	rotbord<true/*false*/>( LAMv, col, vv3, Re, Im );
	col.Star( LAMvS, LAMv );
	func += hh3(1) * LAMvS;
}
