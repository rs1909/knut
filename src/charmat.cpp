// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "matrix.h"
#include "spmatrix.h"
#include "ncolloc.h"

#include "charmat.h"

#include <cmath>

#define NDIM (col.Ndim())
#define NDEG (col.Ndeg())
#define NINT (col.Nint())
#define NTAU (col.Ntau())
#define NPAR (col.Npar())



CharMat::CharMat( NColloc& col ) : 
	phi( NDIM*(NDEG*NINT+1) ),
	AmzB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	AmzB_x( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	AmzB_p( NDIM*(NDEG*NINT+1) ), 
	phiData( NDEG*NINT, NDIM, 2*NTAU+1 )  // it is not required
{
	ZZ = 1.0;
}

void CharMat::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z )
{
	ZZ = Z;
	col.CharJac_x( AmzB, par, solData, ZZ );
}

void CharMat::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector& qq )
{
	ZZ = Z;
	col.CharJac_x( AmzB, par, solData, ZZ );
	
	Vector shift( phi.Size() );
	
	phi.Clear();
	for( int i = 0; i < NDIM; i++ ) shift( i ) = qq(i);
	AmzB.Solve( phi, shift );
	
	col.Interpolate( phiData, phi );
}

// creates a new delta, by overwriting it

void CharMat::Delta( MatFact& delta, NColloc& col )
{
	Matrix shift( phi.Size(), NDIM );
	Matrix out( phi.Size(), NDIM );
  
	delta.New();
	for( int i = 0; i < NDIM; i++ ) shift( i, i ) = 1.0;
	AmzB.Solve( out, shift );
	// delta need not to be NDIM*NDIM
	for( int i = 0; i < NDIM; i++ )
	{
		for( int j = 0; j < NDIM; j++ )  // L-zM
		{
			delta( i, j ) = out( i, j ) - ZZ*out( i + NDIM*NDEG*NINT, j );
			// std::cout<<out(i,j)<<","<<out( i + NDIM*NDEG*NINT, j )<<" ";
		}
	}
}

void CharMat::Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	Matrix LmzM( phi.Size(), NDIM );
	Matrix temp_x( phi.Size(), NDIM );
	
// 	col.Interpolate( phiData, phi );
	col.CharJac_x_x( AmzB_x, par, solData, phiData, ZZ );
	
	for( int i = 0; i < NDIM; i++ ) // L - zM
	{ 
		LmzM( i, i ) = -1.0;
		LmzM( i + NDIM*NDEG*NINT, i ) = 1.0*ZZ;
	}
	AmzB.Solve( temp_x, LmzM, true );
	//AmzB_u.AvRT( delta_x, temp_x );
	// AmzB_x.Swap();
	// delta needs to be NDIM*NDIM
	AmzB_x.AX( delta_x, temp_x, 1.0, true );
	//delta_x.Print();
}

void CharMat::Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p )
{
	Vector temp_p( phi.Size() );
	
	col.CharJac_x_p( AmzB_p, par, solData, phiData, ZZ, p );
	
	AmzB.Solve( temp_p, AmzB_p );
	// delta need not to be NDIM*NDIM
	for( int i = 0; i < NDIM; i++ )
	{
		delta_p(i) = ZZ*temp_p( i + NDIM*NDEG*NINT ) - temp_p( i );
	}
	// std::cout<<" DELTA_p: "<<delta_p*delta_p<<" AmzB_p: "<<AmzB_p*AmzB_p<<"\n";
}

// colloc is not needed
void CharMat::Switch( Vector& tan, NColloc& /*col*/ )
{
// 	std::cout<<"Charmat::Switch entering\n";
	if( tan.Size() != phi.Size() )
	{
		std::cout<<"Charmat::Switch: Bad sizes\n";
		throw(-12);
	}
	tan = phi;
}

// complex routines

CharMatCPLX::CharMatCPLX( NColloc& col ) : 
	phi( 2*NDIM*(NDEG*NINT+1) ),
	AmzB( 'R', 2*NDIM*(NDEG*NINT+1), 4*NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	AmzB_x( 'R', 2*NDIM*(NDEG*NINT+1), 4*NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	AmzB_p( 2*NDIM*(NDEG*NINT+1) ),
	AmzB_z( 2*NDIM*(NDEG*NINT+1) ),
	phiData( 2*NDEG*NINT, NDIM, NTAU+1 )
{
	zRe = 1.0;
	zIm = 0.0;
}

void CharMatCPLX::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Re, double Im )
{
	zRe = Re;
	zIm = Im;
	col.CharJac_x( AmzB, par, solData, zRe, zIm );
}

void CharMatCPLX::Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Re, double Im, Vector& qq )
{
	zRe = Re;
	zIm = Im;
	col.CharJac_x( AmzB, par, solData, zRe, zIm );
	
	Vector shift( phi.Size() );
	
	phi.Clear();
	for( int i = 0; i < 2*NDIM; i++ ) shift( i ) = qq(i);
	AmzB.Solve( phi, shift );
	
	col.InterpolateCPLX( phiData, phi );
}

void CharMatCPLX::Delta( MatFact& delta, NColloc& col )
{
	Matrix shift( phi.Size(), 2*NDIM );
	Matrix out( phi.Size(), 2*NDIM );

	delta.New();
	for( int i = 0; i < 2*NDIM; i++ ) shift( i, i ) = 1.0;
	AmzB.Solve( out, shift );
	for( int i = 0; i < NDIM; i++ )
	{
		for( int j = 0; j < NDIM; j++ )  // L-zM
		{
			delta( 2*i, 2*j )     = out( 2*i, 2*j )     - zRe*out( 2*(i + NDIM*NDEG*NINT), 2*j )   + zIm*out( 2*(i + NDIM*NDEG*NINT)+1, 2*j );
			delta( 2*i, 2*j+1 )   = out( 2*i, 2*j+1 )   - zRe*out( 2*(i + NDIM*NDEG*NINT), 2*j+1 ) + zIm*out( 2*(i + NDIM*NDEG*NINT)+1, 2*j+1 );
			delta( 2*i+1, 2*j )   = out( 2*i+1, 2*j )   - zIm*out( 2*(i + NDIM*NDEG*NINT), 2*j )   - zRe*out( 2*(i + NDIM*NDEG*NINT)+1, 2*j );
			delta( 2*i+1, 2*j+1 ) = out( 2*i+1, 2*j+1 ) - zIm*out( 2*(i + NDIM*NDEG*NINT), 2*j+1 ) - zRe*out( 2*(i + NDIM*NDEG*NINT)+1, 2*j+1 );
		}
	}
}

void CharMatCPLX::Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	Matrix LmzM( phi.Size(), 2*NDIM );
	Matrix temp_x( phi.Size(), 2*NDIM );
	
	col.CharJac_x_x( AmzB_x, par, solData, phiData, zRe, zIm );
	
	for( int i = 0; i < NDIM; i++ ) // L - zM
	{ 
		LmzM( 2*i, 2*i ) = -1.0;
		LmzM( 2*i+1, 2*i+1 ) = -1.0;
		LmzM( 2*(i + NDIM*NDEG*NINT), 2*i )     = 1.0*zRe;
		LmzM( 2*(i + NDIM*NDEG*NINT), 2*i+1 )   = 1.0*zIm; // because it is transposed
		LmzM( 2*(i + NDIM*NDEG*NINT)+1, 2*i )   = -1.0*zIm;
		LmzM( 2*(i + NDIM*NDEG*NINT)+1, 2*i+1 ) = 1.0*zRe;
	}
	AmzB.Solve( temp_x, LmzM, true );
	//AmzB_u.AvRT( delta_x, temp_x );
	// AmzB_x.Swap();
	AmzB_x.AX( delta_x, temp_x, 1.0, true );
	//delta_x.Print();
}

void CharMatCPLX::Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p )
{
	Vector temp_p( phi.Size() );
	
	col.CharJac_x_p( AmzB_p, par, solData, phiData, zRe, zIm, p );
	
	AmzB.Solve( temp_p, AmzB_p );
	for( int i = 0; i < NDIM; i++ )
	{
		delta_p(2*i)   = zRe*temp_p( 2*(i + NDIM*NDEG*NINT) ) - zIm*temp_p( 2*(i + NDIM*NDEG*NINT)+1 ) - temp_p( 2*i );
		delta_p(2*i+1) = zIm*temp_p( 2*(i + NDIM*NDEG*NINT) ) + zRe*temp_p( 2*(i + NDIM*NDEG*NINT)+1 ) - temp_p( 2*i+1 );
	}
	// std::cout<<" DELTA_p: "<<delta_p*delta_p<<" AmzB_p: "<<AmzB_p*AmzB_p<<"\n";
}

void CharMatCPLX::Delta_z( Vector& delta_z, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	Vector temp_z( phi.Size() );
	
	col.CharJac_x_z( AmzB_z, par, solData, phiData, zRe, zIm );
	
	AmzB.Solve( temp_z, AmzB_z );
	for( int i = 0; i < NDIM; i++ )
	{
		delta_z(2*i)   = -2.0*zIm*zRe*temp_z( 2*(i + NDIM*NDEG*NINT) ) + (zIm*zIm-zRe*zRe)*temp_z( 2*(i + NDIM*NDEG*NINT)+1 ) + 
		                 (  zIm*( temp_z( 2*i ) + phi( 2*(i + NDIM*NDEG*NINT) ) ) + zRe*( temp_z( 2*i+1 ) + phi( 2*(i + NDIM*NDEG*NINT)+1 ) ) );
		delta_z(2*i+1) = (zRe*zRe-zIm*zIm)*temp_z( 2*(i + NDIM*NDEG*NINT) ) - 2.0*zIm*zRe*temp_z( 2*(i + NDIM*NDEG*NINT)+1 ) +
		                 ( -zRe*( temp_z( 2*i ) + phi( 2*(i + NDIM*NDEG*NINT) ) ) + zIm*( temp_z( 2*i+1 ) + phi( 2*(i + NDIM*NDEG*NINT)+1 ) ) );
	}
	// std::cout<<" DELTA_z: "<<delta_z*delta_z<<" AmzB_z: "<<AmzB_z*AmzB_z<<"\n";
}

// colloc is not needed (just for compatibility reasons)
void CharMatCPLX::Switch( Vector& Re, Vector& Im, double& alpha, NColloc& /*col*/ )
{
	std::cout<<"CharmatCPLX::Switch entering\n";
	if( (2*Re.Size() != phi.Size())||(2*Im.Size() != phi.Size()) )
	{
		std::cout<<"CharmatCPLX::Switch: Bad sizes\n";
		throw(-12);
	}
	for( int i=0; i<Re.Size(); i++ )
	{
		Re(i) = phi( 2*i );
		Im(i) = phi( 2*i + 1 );
	}
	if( zRe > 0.0 )
	{
		alpha = atan( fabs(zIm/zRe) );
	}
	else
	{
		alpha = atan( fabs(zRe/zIm) ) + M_PI/2.0;
	}
	std::cout<<"alpha "<<alpha<<"\n";
// 	throw(-1);
}

CharMatLPAUT::CharMatLPAUT( NColloc& col ) :
	Q3( NDIM ),
	Phi1( NDIM*(NDEG*NINT+1) ),
	Psi1( NDIM*(NDEG*NINT+1) ),
	Phi3( NDIM*(NDEG*NINT+1) ),
	Sigma( NDIM*(NDEG*NINT+1) ),
	
	DxAmzBPhi1( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	DpAmzBPhi1( NDIM*(NDEG*NINT+1) ),
	
	MAmzBinv( NDIM*(NDEG*NINT+1), NDIM ),           // (NDIM+1)*(NDIM*(NDEG*NINT+1))
	
	AmzB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	mB( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	
	DxAmzBSigma( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	DpAmzBSigma( NDIM*(NDEG*NINT+1) ),
	
	DxmBPhi1( 'R', NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) ),
	DpmBPhi1( NDIM*(NDEG*NINT+1) ),
	
	Phi1Data( NDEG*NINT, NDIM, NTAU+1 ),
	SigmaData( NDEG*NINT, NDIM, NTAU+1 ),
	tempAx( NDIM*(NDEG*NINT+1), NDIM ), tempBx( NDIM*(NDEG*NINT+1), NDIM ),
	tempCx( NDIM*(NDEG*NINT+1), NDIM ), tempDx( NDIM*(NDEG*NINT+1), NDIM ),
	tempAp( NDIM*(NDEG*NINT+1) ), tempBp( NDIM*(NDEG*NINT+1) ),
	stmpAp( NDIM ), stmpBp( NDIM ),
	stmpCp( NDIM ), stmpDp( NDIM )
{
	ZZ = 1.0;
}

// CharMatLPAUTROT::CharMatLPAUTROT( NColloc& col, Array1D<int>& Re_, Array1D<int>& Im_ ) :
// 	CharMatLPAUT( col ),
// 	Re( Re_ ), Im( Im_ ),
// 	Q2( NDIM )
// {
// 
// }

void CharMatLPAUT::getQs( NColloc& col, const Vector& par, const JagMatrix3D& solData, Vector* qq )
{

}

// void CharMatLPAUTROT::getQs( NColloc& col, const Vector& par, const JagMatrix3D& solData, Vector* qq )
// {
// 	col.CharJac_phi( Q1, par, solData );
// 	// fill up Q2
// 	Q2.Clear();
// 	for( int i=0; i<Re.Size(); i++ )
// 	{
// 		Q2(Re(i)) = -solData(0)(Im(i),0);
// 		Q2(Im(i)) = solData(0)(Re(i),0);
// 	}
// 	Q1 /= sqrt(Q1*Q1);
// 	Q2 /= sqrt(Q2*Q2);
// 	if( qq )
// 	{
// 		Vector& Q = *qq;
// 		
// 		alpha1 = Q(0);
// 		for( int i = 0; i < NDIM; i++ )
// 		{
// 			Q3(i) = Q(1+i);
// 			QN(i) = alpha1 * Q1(i);// + alpha2 * Q2(i);
// 		}
// 	}
// }

void CharMatLPAUT::init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector* qq )
{
	Matrix shift( Phi1.Size(), NDIM );
	Vector temp( Phi1.Size() );
	
	ZZ = Z;
	
	col.CharJac_phi( Psi1, par, solData );
	double norm = 0.0;
	for( int i = 0; i < NDIM; i++ ) { norm += Psi1(NDIM*NDEG*NINT + i)*Psi1(NDIM*NDEG*NINT + i); }
	norm = sqrt(norm);
	
	col.CharJac_x( AmzB, par, solData, ZZ );
	col.CharJac_mB( mB, par, solData, ZZ );
	
	// MAmzBinv
	for( int i = 0; i < NDIM; i++ ) shift( NDIM*NDEG*NINT + i, i ) = 1.0;
	AmzB.Solve( MAmzBinv, shift, true );
	
	// Phi1
	for( int i = 0; i < NDIM; i++ ) { temp( i ) = Psi1( i )/norm; }
	AmzB.Solve( Phi1, temp, false );
// 	Psi1 -= Phi1;
// 	std::cout<<"phidiff "<<sqrt(col.Integrate(Psi1,Psi1))<<" ";
	// Psi1
	mB.AX( temp, Phi1, 1.0, false );
	AmzB.Solve( Psi1, temp, false );
	
	if( qq )
	{
		Vector& Q = *qq;
		alpha = Q(0);
		for( int i = 0; i < NDIM; i++ ) Q3(i) = Q(1+i);

		// Phi3
		temp.Clear();
		for( int i = 0; i < NDIM; i++ ) temp( i ) = Q3(i);
		AmzB.Solve( Phi3, temp );
		// Sigma
		for( int i=0; i<Phi1.Size(); i++ )
		{
			Sigma(i) = ZZ*( Phi3(i) - alpha*Psi1(i) );
		}
		col.InterpolateREAL( Phi1Data, Phi1 );
		col.InterpolateREAL( SigmaData, Sigma );
	}
}

void CharMatLPAUT::Delta( MatFact& delta, NColloc& col )
{
	delta.New(); // in order to factorize it again...
	delta.Clear();
	
	// delta(NDIM+1,NDIM+1)
	
	for( int i = 0; i < NDIM; i++ )
	{
		// 1,1
		delta( i, 0 ) = - Phi1( NDIM*NDEG*NINT + i ) + ZZ*Psi1( NDIM*NDEG*NINT + i );
		// 2,2
		delta( NDIM, 1+i ) = Phi1( NDIM*NDEG*NINT + i ); // M * Phi1
		//
		for( int j = 0; j < NDIM; j++ )  // L-zM
		{
			// 1,2
			if( i == j ) delta( i, 1+j ) = 1.0 - ZZ*MAmzBinv( j, i );
			else delta( i, 1+j ) = - ZZ*MAmzBinv( j, i );
		}
	}
}

// void CharMatLPAUTROT::Delta( MatFact& delta, NColloc& col )
// {
// 	delta.New(); // in order to factorize it again...
// 	delta.Clear();
// 	
// 	// delta(NDIM+1,NDIM+1)
// 	Vector DDQ1(NDIM);
// 	for( int i = 0; i < NDIM; i++ )
// 	{
// 		DDQ1(i) = 0.0;
// 		for( int j = 0; j < NDIM; j++ )
// 		{
// 			DDQ1(i) += DzDelta(i,j) * Q1(j);
// 		}
// 	}
// 	delta(0,0) = DDQ1 * DDQ1;
// 	
// 	for( int i = 0; i < NDIM; i++ )
// 	{
// 		for( int j = 0; j < NDIM; j++ )
// 		{
// 			delta( i+1, 0 ) += DDQ1(j) * NDelta(j,i);
// 			delta( 0, i+1 ) += DDQ1(j) * NDelta(j,i);
// 		}
// 		for( int j = 0; j < NDIM; j++ )
// 		{
// 			delta(1+i,1+j) += Q1(i) * Q1(j) + Q2(i) * Q2(j);
// 			for( int k = 0; k < NDIM; k++ )
// 			{
// 				delta(1+i,1+j) += NDelta(k,i) * NDelta(k,j);
// 			}
// 		}
// 	}
// }

inline void CharMatLPAUT::mktmp_x( NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	col.CharJac_x_x( DxAmzBPhi1, par, solData, Phi1Data, ZZ );
	
	// Note:  DxPhi1 = - AmzBinv * DxAmzB * { Phi1, . }
	
	// 4. - M Dx Phi1 = - M * ( - AmzBinv * DxAmzB * { Phi1, . } )             ==>> tempDx
	DxAmzBPhi1.AX( tempDx, MAmzBinv, 1.0, true );
	
	// o.k.!
	// 1. \mu * alpha * MAmzBinv * mB * Dx Phi1 = 
	// =  \mu * alpha * MAmzBinv * mB * ( - AmzBinv * DxAmzB * { Phi1, . } )   ==>>   tempBx
	mB.AX( tempBx, MAmzBinv, -alpha*ZZ, true );
	AmzB.Solve( tempAx, tempBx, true );
	DxAmzBPhi1.AX( tempBx, tempAx, 1.0, true );
	
	// o.k.!
	// 2. MAmzBinv * DxAmzB * Sigma              ==>>   tempAx
	col.CharJac_x_x( DxAmzBSigma, par, solData, SigmaData, ZZ );
	DxAmzBSigma.AX( tempAx, MAmzBinv, 1.0, true );
	
	// o.k.!
	// 3. \mu * alpha * MAmzBinv * DxmB * Phi1   ==>>   tempCx
	col.CharJac_mB_x( DxmBPhi1, par, solData, Phi1Data, ZZ );
	DxmBPhi1.AX( tempCx, MAmzBinv, ZZ*alpha, true );
}

void CharMatLPAUT::Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData )
{
	mktmp_x( col, par, solData );
	//. summing up
	for( int i=0; i<Phi1.Size(); i++ )
	{
		delta_x(i,NDIM) = 0.0;
		for( int j=0; j<NDIM; j++ )
		{
			delta_x(i,j) = tempAx(i,j) + tempBx(i,j) + tempCx(i,j) + tempDx(i,j);
// 			delta_x(i,NDIM) -= tempDx(i,j) * Q3(j);
		}
	}
}

// void CharMatLPAUTROT::Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData )
// {
// 	mktmp_x( col, par, solData );
// 	
// 	Vector DDQ1(NDIM);
// 	for( int i = 0; i < NDIM; i++ )
// 	{
// 		DDQ1(i) = 0.0;
// 		for( int j = 0; j < NDIM; j++ )
// 		{
// 			DDQ1(i) += DzDelta(i,j) * Q1(j);
// 		}
// 	}
// 	
// 	delta_x.Clear();
// 	for( int i=0; i<PhiN.Size(); i++ )
// 	{
// 		for( int j=0; j<NDIM; j++ )
// 		{
// 			delta_x(i,0) += DDQ1(j) * ( tempAx(i,j) + tempBx(i,j) + tempCx(i,j) );
// 			for( int k=0; k<NDIM; k++ )
// 			{
// 				delta_x(i,1+j) += NDelta(k,j) * ( tempAx(i,k) + tempBx(i,k) + tempCx(i,k) );
// 			}
// 		}
// 	}
// }


inline void CharMatLPAUT::mktmp_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int p )
{
	// making the derivative of the trivial eigenvector
	col.CharJac_x_p( DpAmzBPhi1, par, solData, Phi1Data, ZZ, p );
	
	// Note:  DpPhi1 = - AmzBinv * DpAmzB * Phi1
	
	// o.k.!
	// 1. \mu * alpha * MAmzBinv * mB * DpPhi1   ==>>   stmpAp
	// =  \mu * alpha * MAmzBinv * mB * ( - AmzBinv * DpAmzB * Phi1 )
	AmzB.Solve( tempAp, DpAmzBPhi1 );
	mB.AX( tempBp, tempAp, -alpha*ZZ, false );
	MAmzBinv.AX( stmpAp, tempBp, 1.0, true );
	
	// o.k.!
	// 2. MAmzBinv * DpAmzB * Sigma              ==>>   stmpBp
	col.CharJac_x_p( DpAmzBSigma, par, solData, SigmaData, ZZ, p );
	MAmzBinv.AX( stmpBp, DpAmzBSigma, 1.0, true );
	
	// o.k.!
	// 3. \mu * alpha * MAmzBinv * DpmB * Phi1   ==>>   stmpCp
	col.CharJac_mB_p( DpmBPhi1, par, solData, Phi1Data, ZZ, p );
	MAmzBinv.AX( stmpCp, DpmBPhi1, ZZ*alpha, true );
	
	// 4. -M * DpPhi1                            ==>>   stpmDp
	// =  -M * ( - AmzBinv * DpAmzB * Phi1 )
	MAmzBinv.AX( stmpBp, DpAmzBPhi1, 1.0, true );
}

void CharMatLPAUT::Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p )
{
	mktmp_p( col, par, solData, p );
	//. summing up
	delta_p(NDIM) = 0.0;
	for( int i=0; i<NDIM; i++ )
	{
		delta_p(i) = stmpAp(i) + stmpBp(i) + stmpCp(i) + stmpDp(i);
// 		delta_p(NDIM) -= stmpDp(i) * Q3(i);
	}
}

// void CharMatLPAUTROT::Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p )
// {
// 	mktmp_p( col, par, solData, p );
// 	
// 	Vector DDQ1(NDIM);
// 	for( int i = 0; i < NDIM; i++ )
// 	{
// 		DDQ1(i) = 0.0;
// 		for( int j = 0; j < NDIM; j++ )
// 		{
// 			DDQ1(i) += DzDelta(i,j) * Q1(j);
// 		}
// 	}
// 	
// 	delta_p.Clear();
// 	for( int j=0; j<NDIM; j++ )
// 	{
// 		delta_p(0) += DDQ1(j) * ( stmpAp(j) + stmpBp(j) + stmpCp(j) );
// 		for( int k=0; k<NDIM; k++ )
// 		{
// 			delta_p(1+j) += NDelta(k,j) * ( stmpAp(k) + stmpBp(k) + stmpCp(k) );
// 		}
// 	}
// }
