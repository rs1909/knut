// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "pderror.h"
#include "point.h"
#include "ncolloc.h"
#include "system.h"
#include "matrix.h"
#include "spmatrix.h"
#include "plot.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

// specified in the constsnts file
#define NITER 30
#define ERRTOL 1e-4
#define ERRTOLCONT 1e-4
#define BIFTOL 1e-6  // only for starting bifurcations
// not specified
#define MIN_NSIMAG 1e-4

#define NDIM colloc.Ndim()
#define NTAU colloc.Ntau()
#define NPAR colloc.Npar()
#define NINT colloc.Nint()
#define NDEG colloc.Ndeg()
// #define CARS 9

// private
void Point::FillSol( System& sys_ )
{
	Vector fx( colloc.Ndim() );
	
	sys_.stpar( par );
	
	for( int i = 0; i < colloc.Nint(); i++ )
	{
		for( int d = 0; d <  colloc.Ndeg(); d++ )
		{
			sys_.stsol( fx, colloc.Profile( i, d ) );
			for( int j=0; j < colloc.Ndim(); j++ )
			{
				sol( NDIM*(i*NDEG+d) + j ) = fx( j );
			}
		}
	}
	for( int j=0; j < colloc.Ndim(); j++ )
	{
		sol( NDIM*NDEG*NINT + j ) = sol( j );
	}
}


Point::Point( System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg, int nmul, int nmat ) :
	PointData( sys_.ndim()*(nint*ndeg+1), sys_.npar()+ParEnd, nmul ),  // sol, par, mRe, mIm
	var( var_ ), eqn( eqn_ ), varMap( var_.Size() ), varMapCont( var_.Size() + 1 ),
	solNu( sys_.ndim()*(nint*ndeg+1) ), parNu( sys_.npar()+ParEnd ),
	rotRe( 2 ), rotIm( 2 ),              // !!!only for Hartmut's equation!!!!
	colloc( sys_, nint, ndeg, nmat ),
	solData( ndeg*nint, sys_.ndim(), 2*sys_.ntau()+1 ),
	jacStab( nmat, NDIM*(NDEG*NINT+1), NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) )
{
	RefEps    = ERRTOL;
	ContEps   = ERRTOLCONT;
	StartEps  = BIFTOL;
	RefIter   = NITER;
	ContIter  = NITER;
	StartIter = NITER;
	
	rotRe(0) = 0;
	rotIm(0) = 1;
	rotRe(1) = 3;
	rotIm(1) = 4;
	
	Construct();
	FillSol( sys_ );
	par(NPAR+ParNorm) = 0.0;
	par(NPAR+ParAngle) = 0.0;
}

Point::~Point()
{
	Dispose();
}

// public
// remember that this ereases the state variables except sol, qq, par xxDot
void Point::Reset( Array1D<Eqn>& eqn_, Array1D<Var>& var_ )
{
	Vector* qq_temp=0;
	HyperVector* xxDot_temp=0;
	if( qq ) qq_temp = new Vector( *qq );
	if( xxDot ) xxDot_temp = new HyperVector( *xxDot );
	Dispose();
	eqn.Init( eqn_.Size() );
	eqn = eqn_;
	var.Init( var_.Size() );
	var = var_;
	varMap.Init( var_.Size() );
	varMapCont.Init( var_.Size() + 1 );
	Construct( );
	if( qq_temp && qq ) *qq = *qq_temp;
	if( xxDot_temp && xxDot ) xxDot->getV1() = xxDot_temp->getV1();
	delete xxDot_temp;
	delete qq_temp;
}

// not even member function public
void PtToEqnVar( Array1D<Eqn>& eqnr, Array1D<Var>& varr, PtType Pt, int nparx, const int* parx, int npar_ )
{
	switch( Pt ){
		case Sol:
		case SolPDSW:
		case SolLPSW:
			if( nparx != 0 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 2 );
			varr.Init( 2 );
			eqnr(0) = EqnSol; eqnr(1) = EqnNone;
			varr(0) = VarSol; varr(1) = VarNone;
			break;
		case BifLP:
			if( nparx != 1 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 3 );
			varr.Init( 3 );
			eqnr(0) = EqnSol; eqnr(1) = EqnLPNullSpace; eqnr(2) = EqnNorm;
			varr(0) = VarSol; varr(1) = VarNullSpace;   varr(2) = (Var)(VarPAR0 + parx[0]);
			break;
		case BifPD:
			if( nparx != 1 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 3 );
			varr.Init( 3 );
			eqnr(0) = EqnSol; eqnr(1) = EqnPDNullSpace; eqnr(2) = EqnNorm;
			varr(0) = VarSol; varr(1) = VarNullSpace;   varr(2) = (Var)(VarPAR0 + parx[0]);
			break;
		case BifNS:
			if( nparx != 1 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 4 );
			varr.Init( 4 );
			eqnr(0) = EqnSol; eqnr(1) = EqnCPLXNullSpace; eqnr(2) = EqnCPLXNormRe;                     eqnr(3) = EqnCPLXNormIm;
			varr(0) = VarSol; varr(1) = VarNullSpace;     varr(2) = (Var)(VarPAR0 + npar_ + ParAngle); varr(3) = (Var)(VarPAR0 + parx[0]);
			break;
		case SolAUT:
		case SolAUTPDSW:
		case SolAUTHOPFSW:
			if( nparx != 1 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 3 );
			varr.Init( 3 );
			eqnr(0) = EqnSol; eqnr(1) = EqnNone; eqnr(2) = EqnPhase;
			varr(0) = VarSol; varr(1) = VarNone; varr(2) = (Var)(VarPAR0 + parx[0]);
			break;
		case BifAUTLP:
			if( nparx != 2 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 4 );
			varr.Init( 4 );
			eqnr(0) = EqnSol; eqnr(1) = EqnLPAUTNullSpace; eqnr(2) = EqnPhase;                 eqnr(3) = EqnNorm;
			varr(0) = VarSol; varr(1) = VarNullSpace;      varr(2) = (Var)(VarPAR0 + parx[0]); varr(3) = (Var)(VarPAR0 + parx[1]);
			break;
		case BifAUTPD:
			if( nparx != 2 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 4 );
			varr.Init( 4 );
// 			eqnr(0) = EqnSol; eqnr(1) = EqnPDNullSpace; eqnr(2) = EqnPhase;                 eqnr(3) = EqnNorm;
// 			varr(0) = VarSol; varr(1) = VarNullSpace;   varr(2) = (Var)(VarPAR0 + parx[0]); varr(3) = (Var)(VarPAR0 + parx[1]);
			eqnr(0) = EqnSol; eqnr(1) = EqnPDNullSpace; eqnr(2) = EqnNorm;                  eqnr(3) = EqnPhase;  
			varr(0) = VarSol; varr(1) = VarNullSpace;   varr(2) = (Var)(VarPAR0 + parx[0]); varr(3) = (Var)(VarPAR0 + parx[1]);
			break;
		case BifAUTNS:
			if( nparx != 2 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 5 );
			varr.Init( 5 );
			eqnr(0) = EqnSol; eqnr(1) = EqnCPLXNullSpace; eqnr(2) = EqnPhase;                          eqnr(3) = EqnCPLXNormRe;            eqnr(4) = EqnCPLXNormIm;
			varr(0) = VarSol; varr(1) = VarNullSpace;     varr(2) = (Var)(VarPAR0 + npar_ + ParAngle); varr(3) = (Var)(VarPAR0 + parx[0]); varr(4) = (Var)(VarPAR0 + parx[1]);
			break;
		case SolTor:
			if( nparx != 1 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 2 );
			varr.Init( 2 );
			eqnr(0) = EqnTORSol; eqnr(1) = EqnTORPhase1;
			varr(0) = VarTORSol; varr(1) = (Var)(VarPAR0 + parx[0]);
			break;
		case SolAUTTor:
			if( nparx != 2 ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
			eqnr.Init( 3 );
			varr.Init( 3 );
			eqnr(0) = EqnTORSol; eqnr(1) = EqnTORPhase0;             eqnr(2) = EqnTORPhase1;
			varr(0) = VarTORSol; varr(1) = (Var)(VarPAR0 + parx[0]); varr(2) = (Var)(VarPAR0 + parx[1]);
			break;
		default:
			std::cout<<"No such pointtype\n"; PDError(-1);
			break;
	}
}

// private
// What does it construct?
// qq, qq0, qqR, qqNu   why are all these necessary?
// charMat
// xxDot, xx, rhs, jac
void Point::Construct( )
{	
	if( (eqn.Size() < 2)||(var.Size() < 2)||(eqn.Size() != var.Size()) )
	{
		std::cout<<" bad sizes ";
		PDError(-1);
	}
	else
	{
		dim3 = eqn.Size() - 2;
	}
	
	if( (eqn(0) != EqnSol)||(var(0) != VarSol) )
	{
		std::cout<<" bad EqnSol ";
		PDError(-1);
	}
	else
	{
		dim1 = NDIM*(NINT*NDEG+1);
	}
	
	if( ((eqn(1) != EqnNone)||(var(1) != VarNone))&&
	    ((eqn(1) != EqnLPNullSpace)||(var(1) != VarNullSpace))&&
	    ((eqn(1) != EqnLPAUTNullSpace)||(var(1) != VarNullSpace))&&
	    ((eqn(1) != EqnPDNullSpace)||(var(1) != VarNullSpace))&&
	    ((eqn(1) != EqnCPLXNullSpace)||(var(1) != VarNullSpace))&&
	    ((eqn(1) != EqnLPAUTROTNullSpace)||(var(1) != VarNullSpace)) )
	{
		std::cout<<" bad EqnNullSpace ";
		PDError(-1);
	}
	else
	{
		switch( eqn(1) )
		{
			case EqnLPNullSpace:
			 goto next;
			case EqnPDNullSpace:
			 next:
				qq   = new Vector( NDIM );
				qq0  = new Vector( NDIM );
				qqR  = new Vector( NDIM );
				qqNu = new Vector( NDIM );
				charMat = new CharMat( colloc );
				dim2 = NDIM;
				break;
				
			case EqnCPLXNullSpace:
				qq   = new Vector( 2*NDIM );
				qq0  = new Vector( 2*NDIM );
				qqR  = new Vector( 2*NDIM );
				qqNu = new Vector( 2*NDIM ); // new Vector( 2*NDIM );
				charMat = new CharMatCPLX( colloc );
				dim2 = 2*NDIM;
				break;
			case EqnLPAUTNullSpace:
				dim2 = NDIM+1;
				qq   = new Vector( dim2 );
				qq0  = new Vector( dim2 );
				qqR  = new Vector( dim2 );
				qqNu = new Vector( dim2 ); 
				charMat = new CharMatLPAUT( colloc );
				break;
			case EqnLPAUTROTNullSpace:
				dim2 = NDIM+1;
				qq   = new Vector( dim2 );
				qq0  = new Vector( dim2 );
				qqR  = new Vector( dim2 );
				qqNu = new Vector( dim2 );
				charMat = 0;//new CharMatLPAUTROT( colloc, rotRe, rotIm );
				break;
			case EqnNone:
				qq = 0;
				qq0 = 0;
				qqR = 0;
				qqNu = 0;
				charMat = 0;
				dim2 = 0;
				break;
			default:
				std::cout<<" bad EqnNullSpace ";
				PDError(-1);
				break;
		}
	}
	
	for( int i = 2; i < var.Size(); i++ )
	{
		if( (var(i) - VarPAR0 >= 0)&&(var(i) - VarPAR0 < NPAR + ParEnd) )
		{
			varMap( i ) = var(i) - VarPAR0;
		}
		else
		{
			std::cout<<" bad PARAMETERS1 "<<i<<", "<<var(i)<<"\n";
			PDError(-1);
		}
	}
	for( int i = 0; i < var.Size(); i++ ) varMapCont(i) = varMap(i);
	
	xxDot   = new HyperVector( dim1, dim2, dim3+1 );
	
	xx      = new HyperVector( dim1, dim2, dim3+1 );
	
	rhs     = new HyperVector( dim1, dim2, dim3+1 );
	
	jac     = new HyperMatrix( dim1, dim2, dim3+1, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) );
}

// private
void Point::Dispose()
{
	delete jac;
	
	delete rhs;
	delete xx;
	
	delete xxDot;
	
	delete charMat;
	delete qqNu;
	delete qq0;
	delete qqR;
	delete qq;
}

// private
inline double Point::SolNorm( Vector& sol, Vector* qq, Vector& par )
{
	double Xnorm = colloc.Integrate( sol, sol );
	
	if( qq ) Xnorm += (*qq)*(*qq);
	for( int i = 2; i < varMap.Size(); i++ ) Xnorm += par( varMap(i) )*par( varMap(i) );
	
	return sqrt( Xnorm );
}

// **************************************************************************************************************** //
//
//                           THE JACOBIAN
//
// **************************************************************************************************************** //

void Point::Jacobian(
		HyperMatrix& AA, HyperVector& RHS,                      // output
		Vector& /*parPrev*/, Vector& par,                           // parameters
		Vector& solPrev, Vector& sol, JagMatrix3D& solData,     // solution
		Vector* qPrev,   Vector* q,                             // eigenvector
		Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
		bool cont )                                             // cont: true if continuation
{
// uses also: eqn, var, CharMat, qqR
// ---------------------------------------------------------------------------------------------------------------- //
//
//                           The periodic solution part
//
// ---------------------------------------------------------------------------------------------------------------- //

	if( eqn(0) == EqnSol )
	{
		colloc.RHS_x( AA.getA11(), par, sol, solData );
		colloc.RHS( RHS.getV1(), par, sol, solData );
		
		for( int i=2; i<varMap.Size(); i++ )
		{
			if( varMap(i) < NPAR )
			{
				colloc.RHS_p( AA.getA13(i-2), par, sol, solData, varMap(i) );
			}
			else
			{
				switch( varMap(i)-NPAR )
				{
					case ParNorm:
					case ParAngle:
						AA.getA13(i-2).Clear();
						break;
					default:
						std::cout<<" bad PARAMETERS2 "<<varMap(i)<<"\n";
						PDError(-1);
						break;
				}
			}
		}
	}

// ---------------------------------------------------------------------------------------------------------------- //
//
//                           The Characteristic Matrix part
//
// ---------------------------------------------------------------------------------------------------------------- //
	
	switch( eqn(1) )
	{
		case EqnNone:
			// nothing to do
			break;
		case EqnLPAUTROTNullSpace:
		case EqnLPAUTNullSpace:
		case EqnLPNullSpace:
			// a) initialize
			charMat->Init( colloc, par, solData, 1.0, *q );
			goto skip1;
			break;
		case EqnPDNullSpace:
			// a) initialize
			charMat->Init( colloc, par, solData, -1.0, *q );
		  skip1:
			// b) charMat
			charMat->Delta( AA.getA22(), colloc );
			// c) deriv w.r.t. the solution
			charMat->Delta_x( AA.getA21(), colloc, par, solData );
			// d) w.r.t the parameters
			for( int i=2; i<varMap.Size(); i++ )
			{
				if( varMap(i) < NPAR )
				{
					charMat->Delta_p( AA.getA23(i-2), colloc, par, solData, varMap(i) );
				}else
				{
					switch( varMap(i)-NPAR )
					{
						case ParNorm:
							AA.getA23(i-2) = *qqR;
							break;
						case ParAngle: // not supported
						default:
							std::cout<<" bad PARAMETERS3 "<<varMap(i)<<"\n";
							PDError(-1);
							break;
					}
				}
			}
			// e) writing the RHS
			RHS.getV2().Clear();
			for( int i=0; i < dim2; i++ )
			{
				for( int j=0; j < dim2; j++ )
				{
					RHS.getV2()(i) -= AA.getA22()(i,j) * (*q)(j);
				}
				RHS.getV2()(i) -= par(NPAR+ParNorm)*(*qqR)(i);
			}
			break;
		case EqnCPLXNullSpace:
			// a) initialize
			charMat->Init( colloc, par, solData, cos(par(NPAR+ParAngle)), sin(par(NPAR+ParAngle)), *q );
			// b) charMat
			charMat->Delta( AA.getA22(), colloc );
			// c) deriv w.r.t. the solution
			charMat->Delta_x( AA.getA21(), colloc, par, solData );
			// d) w.r.t the parameters
			for( int i=2; i<varMap.Size(); i++ )
			{
				if( varMap(i) < NPAR )
				{
					charMat->Delta_p( AA.getA23(i-2), colloc, par, solData, varMap(i) );
				}else
				{
					switch( varMap(i)-NPAR )
					{
						case ParNorm:
							AA.getA23(i-2) = *qqR; // xxRM->getV2();
							break;
						case ParAngle:
							charMat->Delta_z( AA.getA23(i-2), colloc, par, solData );
							break;
						default:
							std::cout<<" bad PARAMETERS4 "<<varMap(i)<<"\n";
							PDError(-1);
							break;
					}
				}
			}
			// e) writing the RHS
			RHS.getV2().Clear();
			for( int i=0; i < dim2; i++ )
			{
				for( int j=0; j < dim2; j++ )
				{
					RHS.getV2()(i) -= AA.getA22()(i,j) * (*q)(j);
				}
				RHS.getV2()(i) -= par(NPAR+ParNorm)*(*qqR)(i); // par(NPAR+ParNorm)*(xxRM->getV2()(i));
			}
// 			AA.getA22().Print();
			break;
		default:
			std::cout<<" EqnNullSpace: not supported ";
			PDError(-1);
			break;
	}

// ---------------------------------------------------------------------------------------------------------------- //
//
//                           Other equations
//
// ---------------------------------------------------------------------------------------------------------------- //

	for( int i=2; i<eqn.Size(); i++ )
	{
		switch( eqn(i) )
		{
			case EqnNorm:
				// the norm of the eigenvector remains the same
				AA.getA31(i-2).Clear();
				AA.getA32(i-2) = *qPrev;
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(i) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}else
					{
						switch( varMap(i)-NPAR )
						{
							case ParNorm:
								AA.getA33(i-2,j-2) = 0.0; // xxRM->getV3()(j-2);
								break;
							case ParAngle: // not supported in Real cases
							default:
								std::cout<<" bad PARAMETERS5 "<<varMap(i)<<"\n";
								PDError(-1);
								break;
						}
					}
				}
				// RHS
				RHS.getV3()(i-2) = 1.0 - (*qPrev)*(*q); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			case EqnCPLXNormRe:
				// the norm of the eigenvector remains the same
				AA.getA31(i-2).Clear();
				for( int j = 0; j < NDIM; j++ )
				{
					AA.getA32(i-2)(2*j)   = (*qPrev)(2*j);
					AA.getA32(i-2)(2*j+1) = (*qPrev)(2*j+1);
				}
				// AA.getA32(i-2) = *qPrev; this is the same...
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(i) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}else
					{
						switch( varMap(i)-NPAR )
						{
							case ParNorm:
								AA.getA33(i-2,j-2) = 0.0; // xxRM->getV3()(j-2);
								break;
							case ParAngle: 
								AA.getA33(i-2,j-2) = 0.0;
								break;
							default:
								std::cout<<" bad PARAMETERS6 "<<var(i)<<"\n";
								PDError(-1);
								break;
						}
					}
				}
				// RHS
				RHS.getV3()(i-2) = 1.0 - AA.getA32(i-2)*(*q); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			case EqnCPLXNormIm:
				AA.getA31(i-2).Clear();
				for( int j=0; j<NDIM; j++ )
				{
					AA.getA32(i-2)(2*j)   = -(*qPrev)(2*j+1);
					AA.getA32(i-2)(2*j+1) = (*qPrev)(2*j);
				}
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(i) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}else
					{
						switch( varMap(i)-NPAR )
						{
							case ParNorm:
								AA.getA33(i-2,j-2) = 0.0; // xxRM->getV3()(j-2);
								break;
							case ParAngle:
								AA.getA33(i-2,j-2) = 0.0;
								break;
							default:
								std::cout<<" bad PARAMETERS7 "<<varMap(i)<<"\n";
								PDError(-1);
								break;
						}
					}
				}
				// RHS
				RHS.getV3()(i-2) = 0.0 - AA.getA32(i-2)*(*q); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			// Phase conditions. It is already linear!
			case EqnPhase:
				colloc.PhaseStar( AA.getA31(i-2), solPrev ); // this should be the previous solution!!!
				if( dim2 != 0 ) AA.getA32(i-2).Clear();
				// other variables
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(i) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}else
					{
						switch( varMap(i)-NPAR )
						{
							case ParNorm:
								AA.getA33(i-2,j-2) = 0.0; // xxRM->getV3()(j-2);
								break;
							case ParAngle:
								AA.getA33(i-2,j-2) = 0.0;
								break;
							default:
								std::cout<<" bad PARAMETERS8 "<<varMap(i)<<"\n";
								PDError(-1);
								break;
						}
					}
				}
				// RHS
				RHS.getV3()(i-2) = -( AA.getA31(i-2)*sol ); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			case EqnPhaseRot:
				colloc.PhaseRotStar( AA.getA31(i-2), solPrev, rotRe, rotIm ); // this should be the previous solution!!!
				if( dim2 != 0 ) AA.getA32(i-2).Clear();
				// other variables
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(i) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}else
					{
						switch( varMap(i)-NPAR )
						{
							case ParNorm:
								AA.getA33(i-2,j-2) = 0.0; // xxRM->getV3()(j-2);
								break;
							case ParAngle:
								AA.getA33(i-2,j-2) = 0.0;
								break;
							default:
								std::cout<<" bad PARAMETERS8 "<<varMap(i)<<"\n";
								PDError(-1);
								break;
						}
					}
				}
				// RHS
				RHS.getV3()(i-2) = -( AA.getA31(i-2)*sol ); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			default:
				std::cout<<" EqnOther: not supported ";
				PDError(-1);
				break;
		}
	}
	if( cont )
	{
		// copying the tangent
		if( dim1 != 0 ) colloc.Star( AA.getA31(dim3), xxDot->getV1() );
		if( dim2 != 0 ) AA.getA32(dim3) = xxDot->getV2();
		for( int i = 0; i < dim3; i++ ) AA.getA33()(dim3,i) = xxDot->getV3()(i);
		AA.getA33()(dim3,dim3) = p1Dot;
		
		// the arclength difference
// 		double R = ds - p1Dot * (par(p1) - parPrev(p1)) - colloc.IntegrateCont( xxDot->getV1(), sol, solPrev );
// 		for( int i = 2; i < varMap.Size()-1; i++ ) R -= xxDot->getV3()(i-2) * ( par(varMap(i)) - parPrev(varMap(i)) );
// 		if( q ) for( int i = 0; i < qPrev->Size(); i++ ) R -= xxDot->getV2()(i) * ( (*q)(i) - (*qPrev)(i) );
		RHS.getV3()(dim3) = 0.0;
	}
}

// **************************************************************************************************************** //
//
//                          !END! THE JACOBIAN !END!
//
// **************************************************************************************************************** //

inline void Point::Update( HyperVector& X )
{
	sol += X.getV1();
	if( qq ) *qq += X.getV2();
	for( int i = 2; i < varMap.Size(); i++ ) par( varMap(i) ) += X.getV3()(i-2);
}

inline void Point::ContUpdate( HyperVector& X )
{
	solNu += X.getV1();
	if( qqNu ) *qqNu += X.getV2();
	for( int i = 2; i < varMap.Size(); i++ )
	{
		parNu( varMap(i) ) += X.getV3()(i-2);
	}
	parNu( p1 ) += X.getV3()(dim3);
}

// rewrite with the secant method! The derivative is not reliable.
int Point::Start(  )
{
	colloc.Init( par, sol );
	colloc.Interpolate( solData, sol );
	
	// if there is a characteristic matrix to compute its kernel
	if( eqn(1) != EqnNone )
	{
		Vector wr( dim2 ), wi( dim2 );
		Matrix vr( dim2, dim2 ), vl( dim2, dim2 );
//!!!!!!!!!!!!!!!!
// making charmat...
//!!!!!!!!!!!!!!!!
		switch( eqn(1) )
		{
			case EqnLPAUTROTNullSpace:
			case EqnLPAUTNullSpace:
			case EqnLPNullSpace:
				charMat->Init( colloc, par, solData, 1.0 );
				charMat->Delta( jac->getA22(), colloc );
				break;
			case EqnPDNullSpace:
				charMat->Init( colloc, par, solData, -1.0 );
				charMat->Delta( jac->getA22(), colloc );
				break;
			case EqnCPLXNullSpace:
				{
				// searching for the nearest multiplier to the unit circle
					Stability();
					double dmin = 10.0;
					int imin=0;
					for( int i=0; i< mRe.Size(); i++ )
					{
						if( dmin > fabs(sqrt(mRe(i)*mRe(i) + mIm(i)*mIm(i)) - 1.0) )
						{
							if( fabs(mIm(i)) > MIN_NSIMAG )
							{
								dmin = fabs(sqrt(mRe(i)*mRe(i) + mIm(i)*mIm(i)) - 1.0);
								imin = i;
							}
						}
					}
					double zRe = mRe(imin), zIm = fabs(mIm(imin));
					double nrm = sqrt( zRe*zRe + zIm*zIm );
					
					std::cout<<"mRe(imin) "<<zRe<<" mIm(imin) "<<zIm<<" nrm: "<<nrm<<"\n";
					
					charMat->Init( colloc, par, solData, zRe/nrm, zIm/nrm );
					charMat->Delta( jac->getA22(), colloc );
					if( zRe > 0 )
					{
						par(NPAR+ParAngle) = atan( zIm/zRe );
					}
					else
					{
						par(NPAR+ParAngle) = atan( fabs(zRe/zIm) ) + M_PI/2.0;
					}
// 					std::cout<<"COS "<<cos(par(NPAR+ParAngle))<<" SIN "<<sin(par(NPAR+ParAngle))<<"\n";
				}
				break;
			default:
				std::cout<<" EqnNullSpace: not supported ";
				PDError(-1);
				break;
		}
//!!!!!!!!!!!!!!!!
// END making charmat...
//!!!!!!!!!!!!!!!!
// 		std::cout<<"DELTA:\n";
// 		jac->getA22().Print();
		// getting the eigenvalues and eigenvectors
		jac->getA22().Eigval( wr, wi, vl, vr );
// 		std::cout<<"Eigenvalues of DELTA:\n";
// 		wr.Print();
// 		wi.Print();
// 		jac->getA22().Print();
		
		// finding the minimal eigenvalueone = 1.;
		double absval, min = DBL_MAX;
		int imin=-1;
		for( int i = 0; i < dim2; i++ )
		{
			absval = sqrt(wr(i)*wr(i) + wi(i)*wi(i));
			if( min > absval ){
				min = absval;
				imin = i;
			}
		}
		if( imin != -1 )
		{
			for( int i = 0; i < dim2; i++ )
			{
				(*qq)(i) = vr(i,imin);
			}
		}else
		{
			std::cout<<" Point::Start: haven't found a nullvector ";
			PDError(-1);
		}
		*qq /= sqrt( (*qq)*(*qq) );
		*qq0 = *qq;

//!!!!!!!!!!!!!!!!
// setting up things for continuation
//!!!!!!!!!!!!!!!!
		
		double a = 0.0, b=0.0;
		charMat->Delta( jac->getA22(), colloc );
		switch( eqn(1) )
		{
			case EqnLPAUTROTNullSpace:
			case EqnLPAUTNullSpace:
			case EqnLPNullSpace:
			case EqnPDNullSpace:
				*qqR = jac->getA22() * (*qq0);
// 				jac->getA22().AX( *qqR, *qq0, 1.0, false );
				par(NPAR + ParNorm) = -((*qq0)*(*qqR));
// 				std::cout<<" par(NPAR + ParNorm) "<<par(NPAR + ParNorm)<<"\n";
				*qqR = *qq0;
				break;
			case EqnCPLXNullSpace:
				*qqR = jac->getA22() * (*qq0);
// 				jac->getA22().AX( *qqR, *qq0, 1.0, false );
// 				std::cout<<"\n qqR_e: "; qqR->Print();
// 				std::cout<<" qq: "; qq->Print();
				a = (*qq0)*(*qqR);
				for( int i=0; i<NDIM; i++ )
				{ 
					b += -(*qq0)(2*i+1)*(*qqR)(2*i) + (*qq0)(2*i)*(*qqR)(2*i+1);
				}
// 				std::cout<<" a: "<<a<<" b: "<<b<<"\n";
				
				par(NPAR + ParNorm) = -sqrt(a*a+b*b);
// 				std::cout<<" par(NPAR + ParNorm)CPLX "<<par(NPAR + ParNorm)<<"\n";
				for( int i=0; i<NDIM; i++ )
				{ 
					(*qqR)(2*i)   = (a*(*qq0)(2*i) - b*(*qq0)(2*i+1))/sqrt(a*a+b*b);
					(*qqR)(2*i+1) = (a*(*qq0)(2*i+1) + b*(*qq0)(2*i))/sqrt(a*a+b*b);
				}
// 				std::cout<<" qqR_v: "; qqR->Print();
				break;
			default: 
				std::cout<<" EqnNullSpace: not supported ";
				PDError(-1);
				break;
		}
//!!!!!!!!!!!!!!!!
// END setting up things for continuation
//!!!!!!!!!!!!!!!!
// 		wr.Print();
// 		wi.Print();
// 		qq->Print();
// 		qqR->Print();
		
		p1 = NPAR+ParNorm;
		Tangent();
		
		std::cout<<"\nLABEL\t"<<"   NORM\t\t"<<parType( NPAR, p1 )<<parNum( NPAR, p1 )<<"\t";
		for( int j = 2; j < varMap.Size(); j++ ) std::cout<<"\t"<<parType( NPAR, varMap(j) )<<parNum( NPAR, varMap(j) )<<"\t";
		std::cout<<"\tIT\n";
		
		int cit=0,it=0;
		do
		{
			cit = Continue( -par(p1)/p1Dot/1.5 );
			std::cout<<" s"<<it+1<<"\t"<<Norm()<<"\t"<<par(p1);
			for( int j = 2; j < varMap.Size(); j++ ) std::cout<<"\t"<<par(varMap(j));
			std::cout<<"\t"<<cit<<"\n";
			std::cout.flush();
		}while( (fabs( par(p1) ) > StartEps)&&(++it < StartIter) );
		par( p1 ) = 0.0;
		qqR->Clear();
	}else
	{
		Refine();
	}
	
	return 0;
}

int Point::Refine( )
{
	int it=0;
	double Xnorm, Dnorm;
	
	if( qq ) *qq0 = *qq;
	solNu = sol; // here solNu is the previous solution
	parNu = par; // here parNu is the previous parameter

	std::cout<<"\nIT\tERR\t\tSOLnorm\t\tDIFFnorm\n";
	
	do
	{
		colloc.Init( par, sol );
		colloc.Interpolate( solData, sol );
	
		Jacobian( *jac, *rhs, parNu, par, solNu, sol, solData, qq0, qq, varMap, false );
		jac->Solve( *xx, *rhs, dim3 );
		
		Update( *xx );
		
		// computing norms to determine convergence
		Xnorm = sqrt( colloc.Integrate( sol, sol ) );
		Dnorm = sqrt( colloc.Integrate( xx->getV1(), xx->getV1() ) + (xx->getV2())*(xx->getV2()) + (xx->getV3())*(xx->getV3()) );
		std::cout<<" "<<it<<"\t"<<Dnorm/(1.0+Xnorm)<<"\t"<<Xnorm<<"\t"<<Dnorm<<'\n';
		std::cout.flush();
	}
	while( (Dnorm/(1.0+Xnorm) >= RefEps)&&(it++ < RefIter) );
	
	return it;
}


// this is just for periodic systems in case of pitchfork bifurcation
void Point::SwitchLP( double ds )
{
	Vector qq_(NDIM);
	MatFact DD( NDIM, NDIM );
	Vector wr( NDIM ), wi( NDIM );
	Matrix vr( NDIM, NDIM ), vl( NDIM, NDIM );

	std::cout<<"Point::SwitchLP\n";
	xxDot->getV1().Clear();
	xxDot->getV2().Clear();
	xxDot->getV3().Clear();
	p1Dot = 0.0;
	
	// a) finding the eigenvector
	CharMat *chm = new CharMat( colloc );
	colloc.Init( par, sol );
	colloc.Interpolate( solData, sol );
	chm->Init( colloc, par, solData, 1.0 );
	chm->Delta( DD, colloc );
	DD.Eigval( wr, wi, vl, vr );

	// b) finding q
	// finding the minimal eigenvalueone = 1.;
	double absval, min=10.0;
	int imin=-1;
	for( int i = 0; i < NDIM; i++ )
	{
		absval = sqrt(wr(i)*wr(i) + wi(i)*wi(i));
		if( min > absval ){
			min = absval;
			imin = i;
		}
	}
	if( imin != -1 )
	{
		for( int i = 0; i < NDIM; i++ )
		{
			qq_(i) = vr(i,imin);
		}
	}else
	{
		std::cout<<" Point::SwitchLP: haven't found a nullvector ";
		PDError(-1);
	}
	qq_ /= sqrt( qq_*qq_ );
	
	// c) computing the eigenvector
	chm->Init( colloc, par, solData, -1.0, qq_ );
	chm->Switch( xxDot->getV1(), colloc );
	delete chm;
	
	double norm = sqrt(colloc.Integrate( xxDot->getV1(), xxDot->getV1() ));
	xxDot->getV1() /= norm;
	xxDot->getV2().Clear();
	xxDot->getV3().Clear();
	sol += ds * xxDot->getV1();
}

void Point::SwitchPD( double ds )
{
	Vector qq_(NDIM);
	MatFact DD( NDIM, NDIM );
	Vector wr( NDIM ), wi( NDIM );
	Matrix vr( NDIM, NDIM ), vl( NDIM, NDIM );

	xxDot->getV1().Clear();
	xxDot->getV2().Clear();
	xxDot->getV3().Clear();
	p1Dot = 0.0;
	
	// a) finding the eigenvector
	CharMat *chm = new CharMat( colloc );
	colloc.Init( par, sol );
	colloc.Interpolate( solData, sol );
	chm->Init( colloc, par, solData, -1.0 );
	chm->Delta( DD, colloc );
	DD.Eigval( wr, wi, vl, vr );

	// b) finding q
	// finding the minimal eigenvalueone = 1.;
	double absval, min=1.0;
	int imin=-1;
	for( int i = 0; i < NDIM; i++ )
	{
		absval = sqrt(wr(i)*wr(i) + wi(i)*wi(i));
		if( min > absval ){
			min = absval;
			imin = i;
		}
	}
	if( imin != -1 )
	{
// 		std::cout<<"The smallest eigenvalue: "<<sqrt(wr(imin)*wr(imin) + wi(imin)*wi(imin))<<"\n";
		for( int i = 0; i < NDIM; i++ )
		{
			qq_(i) = vr(i,imin);
		}
	}else
	{
		std::cout<<" Point::SwitchPD: haven't found a nullvector ";
		PDError(-1);
	}
	qq_ /= sqrt( qq_*qq_ );
	
	// c) computing the eigenvector
	Vector tan(xxDot->getV1().Size());
	chm->Init( colloc, par, solData, -1.0, qq_ );
	chm->Switch( tan, colloc );
	delete chm;
	
	// setting the period two double
	par(0) *= 2.0;
	//	only in case of period doubling and equidistant mesh
	{
		Vector tmp(sol);
		for( int i=0; i<(NDEG*NINT+1)/2; i++ )
		{
			for( int j=0; j<NDIM; j++ )
			{
				xxDot->getV1()(NDIM*i+j)                   =  tan(NDIM*2*i+j);
				xxDot->getV1()(NDIM*((NDEG*NINT+1)/2+i)+j) = -tan(NDIM*2*i+j);
				sol(NDIM*i+j)                   = tmp(NDIM*2*i+j);
				sol(NDIM*((NDEG*NINT+1)/2+i)+j) = tmp(NDIM*2*i+j);
			}
		}
		if( (NDEG*NINT+1)%2 != 0 )
		{
			for( int j=0; j<NDIM; j++ )
			{
				xxDot->getV1()(NDIM*NDEG*NINT+j) = xxDot->getV1()(NDIM*(NDEG*NINT-1)+j);
				sol(NDIM*NDEG*NINT+j) = sol(j);
			}
		}
	}
	double norm = sqrt(colloc.Integrate( xxDot->getV1(), xxDot->getV1() ));
	xxDot->getV1() /= norm;
	xxDot->getV2().Clear();
	xxDot->getV3().Clear();
	sol += ds * xxDot->getV1();
}

// it should be moved to NColloc, because we need the mesh to construct tangent
// Here, it is NOT assumed that the mesh is equidistant
// Moreover, ds should be controlled from the constants file. Currently it is fixed at 0.2,
// which seems to be working for the traffic model
void Point::SwitchHOPF( double ds )
{
	if( par(NPAR+ParAngle) < 0 ) par(NPAR+ParAngle) *= -1.0;
	if( par(NPAR+ParAngle) < M_PI ) par(0) *= 2*M_PI/par(NPAR+ParAngle);
	else par(0) *= 2*M_PI/(par(NPAR+ParAngle) - M_PI);
	std::cout<<"SwitchHOPF: T = "<<par(0)<<", angle/2PI = "<<par(NPAR+ParAngle)/(2*M_PI)<<"\n";
	
	for( int i = 0; i < NINT; i++ )
	{
		for( int j = 0; j < NDEG+1; j++ )
		{
			const double t = colloc.Profile( i, j );
			for( int p = 0; p < NDIM; p++ )
			{
				xxDot->getV1()( p + (j+i*NDEG)*NDIM ) = sin(2.0*M_PI*t) * (*qq)(2*p+1) + cos(2.0*M_PI*t) * (*qq)(2*p);
			}
		}
	}
	double norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) );
	xxDot->getV1() /= norm;
	xxDot->getV2().Clear();
	xxDot->getV3().Clear();
// 	std::cout<<"HOPFtanNorm "<<norm<<"\n";
	
	for( int i = 0; i < NDIM*(NINT*NDEG+1); i++ )
	{
		sol( i ) += ds*xxDot->getV1()(i);
	}
}

void Point::Tangent( )
{
	double norm;
	
	if( qq )
	{
		*qq0 = *qq;
	}
	
	colloc.Init( par, sol );
	colloc.Interpolate( solData, sol );
	
	varMapCont( varMap.Size() ) = p1; // not necessary
	// az RHS-t feleslegesen szamolja ki && the first qq should be qq0
	Jacobian( *jac, *rhs, par, par, sol, sol, solData, qq0, qq, varMapCont, false );
	
	rhs->getV1() = jac->getA13( dim3 );
	if( dim2 != 0 ) rhs->getV2() = jac->getA23( dim3 );
	for( int i = 0; i < dim3; i++ ) rhs->getV3()(i) = jac->getA33( i, dim3 );
	
	jac->Solve( *xxDot, *rhs, dim3 );
	xxDot->getV3()(dim3) = -1.0;
	
	norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV2())*(xxDot->getV2()) + (xxDot->getV3())*(xxDot->getV3()) );
	
	xxDot->getV1() /= -norm;
	xxDot->getV2() /= -norm;
	xxDot->getV3() /= -norm;
	p1Dot = 1.0/norm;
	xxDot->getV3()(dim3) = p1Dot;
	
// 	std::cout<<"Tangent Cnorm: "<<norm<<'\n';
// 	double Pnorm = sqrt(p1Dot*p1Dot), Qnorm = sqrt( (xxDot->getV2())*(xxDot->getV2()) ); 
// 	double Xnorm = sqrt(colloc.Integrate( xxDot->getV1(), xxDot->getV1() )), Onorm = sqrt( (xxDot->getV3())*(xxDot->getV3()) );
// 	std::cout<<"Pnorm: "<<Pnorm<<" Qnorm: "<<Qnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm<<'\n';
}

int Point::Continue( double ds )
{
	double Xnorm, Dnorm, Rnorm, Tnorm;
// 	double Dnorm1, Dnorm2, Dnorm3;
	
	parNu = par;
	for( int i=0; i< solNu.Size(); i++ )     solNu(i)           = sol(i)           + ds * xxDot->getV1()(i);
	if( qqNu ) for( int i=0; i< qqNu->Size(); i++ ) (*qqNu)(i)  = (*qq)(i)         + ds * xxDot->getV2()(i);
	for( int i = 2; i < varMap.Size(); i++ ) parNu( varMap(i) ) = par( varMap(i) ) + ds * xxDot->getV3()(i-2);
	parNu(p1) = par(p1) + ds*p1Dot;
	
	if( qq )
	{
		*qq0 = *qq;
	}
	
	varMapCont( varMap.Size() ) = p1;

	int  it=0;
	bool conv;
	do
	{
		colloc.Init( parNu, solNu );
		colloc.Interpolate( solData, solNu );
		
		Jacobian( *jac, *rhs, par, parNu, sol, solNu, solData, qq, qqNu, varMapCont, true );
		
		jac->Solve( *xx, *rhs );
		
		ContUpdate( *xx );
		
		Rnorm = sqrt( colloc.Integrate( rhs->getV1(), rhs->getV1() ) + 
		              (rhs->getV2())*(rhs->getV2()) + (rhs->getV3())*(rhs->getV3()) );
// 		if( qq) std::cout<<"RN: "<<Rnorm<<" r "<<sqrt((rhs->getV2())*(rhs->getV2()))<<" q "<<sqrt((*qqNu)*(*qqNu))<<' ';
		
		Xnorm = sqrt( colloc.Integrate( solNu, solNu ) );//SolNorm( solNu, qqNu, parNu );
// 		std::cout<<" XN: "<<Xnorm;
		Dnorm = sqrt( colloc.Integrate( xx->getV1(), xx->getV1() ) + 
		              (xx->getV2())*(xx->getV2()) + (xx->getV3())*(xx->getV3()) );
// 		Dnorm1 = sqrt( colloc.Integrate( xx->getV1(), xx->getV1() ) );
// 		Dnorm2 = sqrt( (xx->getV2())*(xx->getV2()) );
// 		Dnorm3 = sqrt( (xx->getV3())*(xx->getV3()) );
// 		std::cout<<"Cont:";
// 		for( int k=0;k<parNu.Size();k++) std::cout<<" par("<<k<<"): "<<par(k);
// 		std::cout<<" Xnorm:"<<Xnorm<<" Dnorm:"<<Dnorm<<"\nDnorm1:"<<Dnorm1<<" Dnorm2:"<<Dnorm2<<" Dnorm3:"<<Dnorm3<<" "<<it<<'\n';
// 		std::cout<<"Par   "; par.Print(); //std::cout<<"\n";
// 		std::cout<<"ParNu "; parNu.Print(); std::cout<<"\n";
// 		if( qqNu ) { std::cout<<"Dq "; rhs->getV2().Print(); std::cout<<"\n"; }
		conv = (Dnorm/(1.0+Xnorm) >= ContEps) || (Rnorm >= 10.0*ContEps);
		
		// updating the tangent
		jac->AX( *rhs, *xxDot );
		rhs->getV3()(dim3) = 0.0;
		jac->Solve( *xx, *rhs );
		xxDot->getV1() -= xx->getV1();
		xxDot->getV2() -= xx->getV2();
		xxDot->getV3() -= xx->getV3();
		Tnorm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV2())*(xxDot->getV2()) + (xxDot->getV3())*(xxDot->getV3()) );
		xxDot->getV1() /= Tnorm;
		xxDot->getV2() /= Tnorm;
		xxDot->getV3() /= Tnorm;
		p1Dot = xxDot->getV3()(dim3);
		// end updating tangent

	}
	while( conv /*&& (Dnorm/(1.0+Xnorm) < 1.0)*/&&(++it < ContIter) );
	if( !conv )
	{
		
		
// 		std::cout<<"Tan1 "<<colloc.Integrate( jacCont->getA13(dim3), jacCont->getA13(dim3) );
// 		if( qq ) std::cout<<" Tan2 "<<(jacCont->getA23(dim3))*(jacCont->getA23(dim3));
// 		for( int i=0; i<dim3; i++ ) std::cout<<" t3 "<<jacCont->getA33(i,dim3);
// 		std::cout<<"\n";
		
		// checking the tangent and the secant
// 		double Pnorm = sqrt(p1Dot*p1Dot), Qnorm = sqrt( (xxDot->getV2())*(xxDot->getV2()) ); 
// 		double Xnorm = sqrt(colloc.Integrate( xxDot->getV1(), xxDot->getV1() )), Onorm = sqrt( (xxDot->getV3())*(xxDot->getV3()) );
// 		std::cout<<"Cnorm: "<<Tnorm<<"\nDot Pnorm: "<<Pnorm<<" Qnorm: "<<Qnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm<<'\n';
// 		
// 		xx->getV1() = solNu;
// 		xx->getV1() -= sol;
// 		if( qq ){ xx->getV2() = *qqNu; xx->getV2() -= *qq; }
// 		for( int i = 2; i < varMap.Size(); i++ ) xx->getV3()(i-2) = parNu( varMap(i) ) - par( varMap(i) );
// 		xx->getV3()(dim3) = parNu(p1) - par(p1);
		
		// replacing xxDot with the secant... it may improve...
// 		xxDot->getV1() = xxCont->getV1(); 
// 		xxDot->getV2() = xxCont->getV2(); 
// 		for( int i=0; i<dim3; i++ ) xxDot->getV3()(i) = xxCont->getV3()(i);
// 		p1Dot = xxCont->getV3()(dim3);
// 		
// 		norm = sqrt( p1Dot*p1Dot + colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV2())*(xxDot->getV2()) + (xxDot->getV3())*(xxDot->getV3()) );
// 		
// 		xxDot->getV1() /= norm;
// 		xxDot->getV2() /= norm;
// 		xxDot->getV3() /= norm;
// 		p1Dot /= norm;
		// end of replacing!!!		
		
// 		Pnorm = sqrt( xx->getV3()(dim3)*xx->getV3()(dim3) )/ds;
// 		Qnorm = sqrt( (xx->getV2())*(xx->getV2()) )/ds; 
// 		Xnorm = sqrt(colloc.Integrate( xx->getV1(), xx->getV1() ))/ds; 
// 		Onorm = 0;
// 		for( int i=0; i<dim3+1; i++ ) Onorm += (xx->getV3()(i))*(xx->getV3()(i));
// 		Onorm = sqrt(Onorm)/ds;
// 		std::cout<<"Dif Pnorm: "<<Pnorm<<" Qnorm: "<<Qnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm<<'\n';
		
		// copying back the solution
		sol = solNu;
		par = parNu;
		if( qq ) *qq = *qqNu;
// 		if( qq ){ std::cout<<" Cont end qq: "; qq->Print(); }
		// renorming qq // 
		if( qq )
		{
			double norm = sqrt((*qq)*(*qq));
			(*qq) /= norm;
		}
		
		
	}else
	{
// 		sol = solNu;
// 		par = parNu;
// 		if( qq ) *qq = *qqNu;
		
		std::cout<<"\n\n\n ------------------- NO CONVERGENCE -------------------\n\n\n\n";
 		// PDError(12);
	}
	
	return it;
}

void Point::Stability( )
{
	mRe.Clear(); 
	mIm.Clear();
	
	colloc.Init( par, sol );
	colloc.Interpolate( solData, sol );
	
	colloc.StabJac( jacStab, par, solData );
	
	jacStab.Eigval( mRe, mIm );
// 	mRe.Print();
// 	mIm.Print();
}

static inline int instab( Vector& mRe, Vector& mIm, int aut )
{
	double dmin1 = DBL_MAX, dmin2 = DBL_MAX;
	int imin1 = -1, imin2 = -1;
	if( aut > 0 )
	{
		for( int i=0; i< mRe.Size(); i++ )
		{
			const double mabs = fabs(sqrt((mRe(i)-1.0)*(mRe(i)-1.0) + mIm(i)*mIm(i)));
			if( dmin1 > mabs ) { dmin1 = mabs; imin1 = i; }
		}
	}
	if( aut > 1 )
	{
		for( int i=0; i< mRe.Size(); i++ )
		{
			const double mabs = fabs(sqrt((mRe(i)-1.0)*(mRe(i)-1.0) + mIm(i)*mIm(i)));
			if( (dmin2 > mabs)&&(i != imin1) ) { dmin2 = mabs; imin2 = i; }
		}
	}
	
	int ustab = 0;
	for( int i=0; i<mRe.Size(); i++)
	{
		if( (mRe(i)*mRe(i) + mIm(i)*mIm(i) >= 1.0)&&(i != imin1)&&(i != imin2) ) ++ustab;
	}
	return ustab;
}

static inline PtType testbif( Vector& mRe, Vector& mIm, int aut )
{
	double dmin1 = 0.1, dmin2 = 0.1;
	int imin1 = -1, imin2 = -1;
	if( aut > 0 )
	{
		for( int i = 0; i < mRe.Size(); i++ )
		{
			const double mabs = fabs(sqrt((mRe(i)-1.0)*(mRe(i)-1.0) + mIm(i)*mIm(i)));
			if( dmin1 > mabs ) { dmin1 = mabs; imin1 = i; }
		}
		if( imin1 == -1 ) std::cout<<"testbif-1: Accuracy problem?\n";
	}
	if( aut > 1 )
	{
		for( int i = 0; i < mRe.Size(); i++ )
		{
			const double mabs = fabs(sqrt((mRe(i)-1.0)*(mRe(i)-1.0) + mIm(i)*mIm(i)));
			if( (dmin2 > mabs)&&(i != imin1) ) { dmin2 = mabs; imin2 = i; }
		}
		if( imin2 == -1 ) std::cout<<"testbif-2: Accuracy problem?\n";
	}

	double dminLP = DBL_MAX, dminPD = DBL_MAX, dminNS = DBL_MAX;
	int iminLP = -1, iminPD = -1, iminNS = -1;
	for( int i = 0; i < mRe.Size(); i++ )
	{
		if( (i != imin1)&&(i != imin2) )
		{
			const double LPabs = fabs( sqrt((mRe(i)-1.0)*(mRe(i)-1.0)) );
			const double PDabs = fabs( sqrt((mRe(i)+1.0)*(mRe(i)+1.0)) );
			const double NSabs = fabs( sqrt(mRe(i)*mRe(i) + mIm(i)*mIm(i)) - 1.0 );
			if( (dminLP > LPabs) && (mIm(i) == 0.0) ) { dminLP = LPabs; iminLP = i; }
			if( (dminPD > PDabs) && (mIm(i) == 0.0) ) { dminPD = PDabs; iminPD = i; }
			if( (dminNS > NSabs) && (mIm(i) != 0.0) ) { dminNS = NSabs; iminNS = i; }
		}
	}
	if( (dminLP < dminPD) && (dminLP < dminNS) ) return BifLP;
	else if( (dminPD < dminLP) && (dminPD < dminNS) ) return BifPD;
	else if( (dminNS < dminPD) && (dminNS < dminLP) ) return BifNS;
	else return Sol;
}

int  Point::UStab()
{
	return instab( mRe, mIm, 0 );
}

int  Point::UStabAUT()
{
	return instab( mRe, mIm, 1 );
}

int  Point::UStabAUTRot()
{
	return instab( mRe, mIm, 2 );
}

PtType Point::testBif()
{
	return testbif( mRe, mIm, 0 );
}

PtType Point::testBifAUT()
{
	return testbif( mRe, mIm, 1 );
}

PtType Point::testBifAUTRot()
{
	return testbif( mRe, mIm, 2 );
}

void Point::Plot( GnuPlot& pl )
{
	pl.SetPointSize( 0.8 );
	
	for( int i = 0; i < NDIM; i++ )
	{
		pl.Plot(i, "with lines");
		for( int j = 0; j < NINT*NDEG+1; j++ )
		{
			double t = (double)j/((double)( NINT*NDEG ));
			pl.AddData( i, t, sol( i + NDIM*j ) ); 
		}
	}
	pl.Show();
}

void Point::Write( std::ofstream& file )
{
	Vector msh( NDEG*NINT+1 );
	colloc.getMesh( msh );
	
	file<<NPAR<<"\t";
	for( int i=0; i<NPAR; i++ ) file<<par(i)<<"\t";
	
	file<<mRe.Size()<<"\t";
	for( int i=0; i<mRe.Size(); i++ ) file<<mRe(i)<<"\t"<<mIm(i)<<"\t";
	
	file<<NDIM<<"\t";
	file<<NINT<<"\t";
	file<<NDEG<<"\t";
	for( int i=0; i<NDEG*NINT+1; i++ ) file<<msh(i)<<"\t";
	for( int i=0; i<NDIM*(NINT*NDEG+1); i++ ) file<<sol(i)<<"\t";
	file<<"\n";
	file.flush();
}

void Point::Read( std::ifstream& file, bool tan )
{
	int npar_, nmul_, ndim_, nint_, ndeg_;
	file>>npar_;
	if( NPAR != npar_ ) { std::cout<<"Not compatible file (NPAR) "<<nmul_<<"\n"; PDError(-1); }
	for( int i=0; i<NPAR; i++ ) file>>par(i);
	
	file>>nmul_;
	if( mRe.Size() < nmul_ ) { std::cout<<"Not compatible file (NMUL) "<<nmul_<<"\n"; PDError(-1); }
	for( int i=0; i<nmul_; i++ ) { file>>mRe(i); file>>mIm(i); }
	
	file>>ndim_;
	file>>nint_;
	file>>ndeg_;
	
	if( NDIM != ndim_ ) { std::cout<<"Not compatible file (NDIM) "<<ndim_<<"\n"; PDError(-1); }
// 	if( NINT != nint_ ) { std::cout<<"Not compatible file (NINT) "<<nint_<<"\n"; PDError(-1); }
// 	if( NDEG != ndeg_ ) { std::cout<<"Not compatible file (NDEG) "<<ndeg_<<"\n"; PDError(-1); }
	
	Vector msh( ndeg_*nint_ + 1 );
	for( int i=0; i<ndeg_*nint_+1; i++ ) file>>msh(i);
	
	if( (NINT == nint_)&&(NDEG == ndeg_) )
	{
		colloc.setMesh( msh );
		for( int i=0; i<NDIM*(nint_*ndeg_+1); i++ ) file>>sol(i);
	}
	else
	{
		Vector in( NDIM*(nint_*ndeg_+1) );
	
		for( int i=0; i<NDIM*(nint_*ndeg_+1); i++ ) file>>in(i);
		colloc.Import( sol, in, msh, ndeg_ );
	}
	if( tan )
	{
		solNu = sol;
		parNu = par;
		Read( file, false );
		xxDot->getV1() = sol;
		xxDot->getV1() -= solNu;
		if( qq ) xxDot->getV2().Clear();
		for( int i = 2; i < varMapCont.Size(); i++ ) xxDot->getV3()(i-2) = par( varMapCont(i) ) - parNu( varMapCont(i) );
		
		// renorming everything
		p1Dot = xxDot->getV3()(dim3);
		
		double norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV2())*(xxDot->getV2()) + (xxDot->getV3())*(xxDot->getV3()) );
		
		xxDot->getV1() /= norm;
		xxDot->getV2() /= norm;
		xxDot->getV3() /= norm;
		p1Dot /= norm;
		std::cout<<"norm "<<norm<<" P!DOT :: "<<p1Dot<<"\n";
		
		sol = solNu;
		par = parNu;
	}
}

void Point::ReadNull( std::ifstream& file )
{
	double tmp;
	int npar_, nmul_, ndim_, nint_, ndeg_;
	file>>npar_;
	if( NPAR != npar_ ) { std::cout<<"RN:Not compatible file (NPAR) "<<nmul_<<"\n"; PDError(-1); }
	for( int i=0; i<NPAR; i++ ) file>>tmp;
	
	file>>nmul_;
	if( mRe.Size() < nmul_ ) { std::cout<<"RN:Not compatible file (NMUL) "<<nmul_<<"\n"; PDError(-1); }
	for( int i=0; i<nmul_; i++ ) { file>>tmp; file>>tmp; }
	
	file>>ndim_;
	file>>nint_;
	file>>ndeg_;
	
	if( NDIM != ndim_ ) { std::cout<<"RN:Not compatible file (NDIM) "<<ndim_<<"\n"; PDError(-1); }
	
	for( int i=0; i<ndeg_*nint_+1; i++ ) file>>tmp;
	for( int i=0; i<NDIM*(nint_*ndeg_+1); i++ ) file>>tmp;
}

void Point::SwitchTRSol( Vector& Sol, const Vector& mshint, const Vector& mshdeg ) // starting data for tori: solution
{
	colloc.Export( Sol, mshint, mshdeg, sol );
}

void Point::SwitchTRTan( Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg ) // starting data for tori: tangent
{
	Vector TRe(sol.Size()), TIm(sol.Size());
	CharMatCPLX* chm;
	chm = static_cast< CharMatCPLX* >(charMat);
	if( chm )
	{
		chm->Switch( TRe, TIm, alpha, colloc );
		colloc.Export( Re, mshint, mshdeg, TRe );
		colloc.Export( Im, mshint, mshdeg, TIm );
	}else
	{
		std::cout<<"Not the complex characteristic matrix\n";
		PDError(-1);
	}
}
