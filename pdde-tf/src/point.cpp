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
	par(NPAR+ParAngle) = 0.0;
	par(NPAR+ParRot) = 0.0;
}

Point::~Point()
{
	Dispose();
}

// public
// remember that this ereases the state variables except sol, qq, par xxDot
void Point::Reset( Array1D<Eqn>& eqn_, Array1D<Var>& var_ )
{
	HyperVector* xxDot_temp=0;
	if( xxDot ) xxDot_temp = new HyperVector( *xxDot );
	Dispose();
	eqn.Init( eqn_.Size() );
	eqn = eqn_;
	var.Init( var_.Size() );
	var = var_;
	varMap.Init( var_.Size() );
	varMapCont.Init( var_.Size() + 1 );
	Construct( );
	if( xxDot_temp && xxDot ) xxDot->getV1() = xxDot_temp->getV1();
	delete xxDot_temp;
}

struct PtTab
{
	BranchSW sw;
	int      neqn;
	int      nparx;
	Eqn      eqns[5];
	Var      vars[5];
};

// not even member function public
BranchSW PtToEqnVar( Array1D<Eqn>& eqnr, Array1D<Var>& varr, PtType Pt, int nparx, const int* parx, int npar_ )
{
	PtTab tab;
	const Var PANGLE = (Var)(VarPAR0+npar_+ParAngle);
	const Var PX0 = (Var)(VarPAR0 + parx[0]);
	const Var PX1 = (Var)(VarPAR0 + parx[1]);
	switch( Pt )
	{
	/// TIME-PERIODIC TEST-FUNCTIONAL
		case SolTF: // same as "case Sol:"
			tab = (PtTab){ NOSwitch,   2, 0,
			 {EqnSol, EqnNone },
			 {VarSol, VarNone } }; break;
		case SolTFBRSW:
			tab = (PtTab){ TFBRSwitch, 2, 0,
			 {EqnSol, EqnNone },
			 {VarSol, VarNone } }; break;
		case SolTFPDSW:
			tab = (PtTab){ TFPDSwitch, 2, 0,
			 {EqnSol, EqnNone },
			 {VarSol, VarNone } }; break;
		case BifTFLP:
			tab = (PtTab){ NOSwitch,   3, 1,
			 {EqnSol, EqnNone,           EqnTFLP },
			 {VarSol, VarNone,           PX0 } }; break;
		case BifTFPD:
			tab = (PtTab){ NOSwitch,   3, 1,
			 {EqnSol, EqnNone,           EqnTFPD },
			 {VarSol, VarNone,           PX0 } }; break;
		case BifTFNS:
			tab = (PtTab){ NOSwitch,   4, 1,
			 {EqnSol, EqnNone,           EqnTFCPLX_RE,  EqnTFCPLX_IM },
			 {VarSol, VarNone,           PANGLE,        PX0 } };        break;
	/// AUTONOMOUS TEST-FUNCTIONAL
		case SolTFAUT: // same as "case SolAUT:"
			tab = (PtTab){ NOSwitch,   3, 1,
			 {EqnSol, EqnNone,           EqnPhase },
			 {VarSol, VarNone,           PX0 } };    break;
		case SolTFAUTBRSW:
			tab = (PtTab){ TFBRSwitch, 3, 1,
			 {EqnSol, EqnNone,           EqnPhase },
			 {VarSol, VarNone,           PX0 } };    break;
		case SolTFAUTPDSW:
			tab = (PtTab){ TFPDSwitch, 3, 1,
			 {EqnSol, EqnNone,           EqnPhase },
			 {VarSol, VarNone,           PX0 } };    break;
		case SolTFAUTHBSW:
			tab = (PtTab){ TFHBSwitch, 3, 1,
			 {EqnSol, EqnNone,           EqnPhase },
			 {VarSol, VarNone,           PX0 } };    break;
		case BifTFAUTLP:
			tab = (PtTab){ NOSwitch,   4, 2,
			 {EqnSol, EqnNone,           EqnPhase,      EqnTFLP },
			 {VarSol, VarNone,           PX0,           PX1 } };   break;
		case BifTFAUTPD:
			tab = (PtTab){ NOSwitch,   4, 2,
			 {EqnSol, EqnNone,           EqnPhase,      EqnTFPD },
			 {VarSol, VarNone,           PX0,           PX1 } };   break;
		case BifTFAUTNS:
			tab = (PtTab){ NOSwitch,   5, 2,
			 {EqnSol, EqnNone,           EqnPhase,      EqnTFCPLX_RE,  EqnTFCPLX_IM },
			 {VarSol, VarNone,           PANGLE,        PX0,           PX1 } };        break;
	/// TORUS
		case SolTor:
			tab = (PtTab){ TFTRSwitch, 2, 1,
			 {EqnTORSol, EqnTORPhase1 },
			 {VarTORSol, PX0 } };        break;
		case SolAUTTor:
			tab = (PtTab){ TFTRSwitch, 2, 2,
			 {EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
			 {VarTORSol, PX0,          PX1 } };        break;
		default:
			tab = (PtTab){ NOSwitch,  0, 0, { EqnNone }, { VarNone } };
			std::cout<<"No such pointtype\n"; PDError(-1);
			break;
	}
	if( tab.nparx != nparx ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
	eqnr.Init( tab.neqn );
	varr.Init( tab.neqn );
	for( int i = 0; i < tab.neqn; i++ )
	{
		eqnr(i) = tab.eqns[i];
		varr(i) = tab.vars[i];
	}
	return tab.sw;
}

// private
// What does it construct?
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
	
	if( ((eqn(1) != EqnNone)||(var(1) != VarNone)) )
	{
		std::cout<<" bad EqnNullSpace ";
		PDError(-1);
	}
	
	testFunct = 0;
	for( int i = 2; i < eqn.Size(); i++ )
	{
		switch( eqn(i) ) 
		{
			case EqnTFLP:
				if( testFunct == 0 ) testFunct = new TestFunct( colloc, 1.0 );
				else PDError(-1);
				break;
			case EqnTFPD:
				if( testFunct == 0 ) testFunct = new TestFunct( colloc, -1.0 );
				else PDError(-1);
				break;
			case EqnTFLPAUT:
				if( testFunct == 0 ) testFunct = new TestFunctLPAUT( colloc, 1.0 );
				else PDError(-1);
				break;
			case EqnTFLPAUTROT:
				if( testFunct == 0 ) testFunct = new TestFunctLPAUTROT( colloc, rotRe, rotIm, 1.0 );
				else PDError(-1);
				break;
			case EqnTFCPLX_RE:
				if( eqn(i+1) != EqnTFCPLX_IM ) { std::cout<<"EqnTFCPLX_RE is not paired\n"; PDError(-1); }
				if( testFunct == 0 ) testFunct = new TestFunctCPLX( colloc );
				else PDError(-1);
				break;
			case EqnTFCPLX_IM:
				if( eqn(i-1) != EqnTFCPLX_RE ) { std::cout<<"EqnTFCPLX_IM is not paired\n"; PDError(-1); }
			default:
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
	
	xxDot   = new HyperVector( dim1, 0, dim3+1 );
	
	xx      = new HyperVector( dim1, 0, dim3+1 );
	
	rhs     = new HyperVector( dim1, 0, dim3+1 );
	
	jac     = new HyperMatrix( dim1, 0, dim3+1, NDIM*(NDEG*NINT+1)*NTAU*NDIM*(NDEG+1) );
}

// private
void Point::Dispose()
{
	delete testFunct;
	delete jac;
	
	delete rhs;
	delete xx;
	
	delete xxDot;
}

// private
// inline double Point::SolNorm( Vector& sol, Vector& par )
// {
// 	double Xnorm = colloc.Integrate( sol, sol );
// 	
// 	for( int i = 2; i < varMap.Size(); i++ ) Xnorm += par( varMap(i) )*par( varMap(i) );
// 	
// 	return sqrt( Xnorm );
// }

// **************************************************************************************************************** //
//
//                           THE JACOBIAN
//
// **************************************************************************************************************** //

void Point::Jacobian(
		HyperMatrix& AA, HyperVector& RHS,                      // output
		Vector& /*parPrev*/, Vector& par,                           // parameters
		Vector& solPrev, Vector& sol, JagMatrix3D& solData,     // solution
		Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
		bool cont )                                             // cont: true if continuation
{
// uses also: eqn, var
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
			else if( varMap(i)-NPAR == ParAngle )
			{
				AA.getA13(i-2).Clear();
			} else
			{
				std::cout<<" bad PARAMETERS2 "<<varMap(i)<<"\n";
				PDError(-1);
			}
		}
	}

// ---------------------------------------------------------------------------------------------------------------- //
//
//                           The Characteristic Matrix part
//
// ---------------------------------------------------------------------------------------------------------------- //
	
	switch( eqn(1) != EqnNone )
	{
		std::cout<<" EqnNullSpace: not supported ";
		PDError(-1);
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
			// Phase conditions. It is already linear!
			case EqnPhase:
				colloc.PhaseStar( AA.getA31(i-2), solPrev ); // this should be the previous solution!!!
				// other variables
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						AA.getA33(i-2,j-2) = 0.0;
					} else
					{
						std::cout<<" bad PARAMETERS8 "<<varMap(j)<<"\n";
						PDError(-1);
					}
				}
				// RHS
				RHS.getV3()(i-2) = -( AA.getA31(i-2)*sol ); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			case EqnPhaseRot:
				colloc.PhaseRotStar( AA.getA31(i-2), solPrev, rotRe, rotIm ); // this should be the previous solution!!!
				// other variables
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						AA.getA33()(i-2,j-2) = 0.0;
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						AA.getA33(i-2,j-2) = 0.0;
					} else
					{
						std::cout<<" bad PARAMETERS9 "<<varMap(j)<<"\n";
						PDError(-1);
					}
				}
				// RHS
				RHS.getV3()(i-2) = -( AA.getA31(i-2)*sol ); // -par(NPAR+ParNorm)*(xxRM->getV3()(i-2))
				break;
			case EqnTFLP:
			case EqnTFPD:
			case EqnTFLPAUT:
			case EqnTFLPAUTROT:
				RHS.getV3()(i-2) = testFunct->Funct( colloc, par, sol, solData );
				testFunct->Funct_x( AA.getA31(i-2), colloc, par, sol, solData );
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						AA.getA33()(i-2,j-2) = testFunct->Funct_p( colloc, par, sol, solData, varMap(j) );
// 						if( varMap(j) == 0 ) std::cout<<"TF-T "<<AA.getA33()(i-2,j-2)<<" ";
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						AA.getA33(i-2,j-2) = 0.0;
					} else
					{
						std::cout<<" bad PARAMETERS10 "<<varMap(j)<<"\n";
						PDError(-1);
					}
				}
				break;
			case EqnTFCPLX_RE:
				if( eqn(i+1) != EqnTFCPLX_IM ) { std::cout<<"EqnTFCPLX_RE is not paired "<<varMap(i)<<"\n"; PDError(-1); }
				testFunct->Funct( RHS.getV3()(i-2), RHS.getV3()(i-1), colloc, par, sol, solData, 
				                  cos(par(NPAR+ParAngle)), sin(par(NPAR+ParAngle)) );
				testFunct->Funct_x( AA.getA31(i-2), AA.getA31(i-1), colloc, par, sol, solData );
				for( int j=2; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						testFunct->Funct_p( AA.getA33()(i-2,j-2), AA.getA33()(i-1,j-2), colloc, par, sol, solData, varMap(j) );
// 						if( varMap(j) == 0 ) std::cout<<"TF-T "<<AA.getA33()(i-2,j-2)<<" ";
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						testFunct->Funct_z( AA.getA33()(i-2,j-2), AA.getA33()(i-1,j-2), colloc, par, sol, solData );
					} else
					{
						std::cout<<" bad PARAMETERS11 "<<varMap(j)<<"\n";
						PDError(-1);
					}
				}
				break;
			case EqnTFCPLX_IM:
				if( eqn(i-1) != EqnTFCPLX_RE )
				{
					std::cout<<"EqnTFCPLX_IM is not paired "<<varMap(i)<<"\n"; PDError(-1);
				}
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
		for( int i = 0; i < dim3; i++ ) AA.getA33()(dim3,i) = xxDot->getV3()(i);
		AA.getA33()(dim3,dim3) = p1Dot;
		
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
	for( int i = 2; i < varMap.Size(); i++ ) par( varMap(i) ) += X.getV3()(i-2);
}

inline void Point::ContUpdate( HyperVector& X )
{
	solNu += X.getV1();
	for( int i = 2; i < varMap.Size(); i++ )
	{
		parNu( varMap(i) ) += X.getV3()(i-2);
	}
	parNu( p1 ) += X.getV3()(dim3);
}

/// --------------------------------------------------------------
/// Starting bifurcation continuation using TEST FUNCTIONS
/// --------------------------------------------------------------

/// It only computes the critical characteristic multiplier and refines the solution

int Point::StartTF( Eqn FN )
{
	if( (FN == EqnTFCPLX_RE) || (FN == EqnTFCPLX_IM) )
	{
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
		double zRe = mRe(imin);
		double zIm = fabs(mIm(imin));
		double nrm = sqrt( zRe*zRe + zIm*zIm );
		std::cout<<"mRe(imin) "<<zRe<<" mIm(imin) "<<zIm<<" nrm: "<<nrm<<"\n";

		if( zRe > 0 )
		{
			par(NPAR+ParAngle) = atan( zIm/zRe );
		}
		else
		{
			par(NPAR+ParAngle) = atan( fabs(zRe/zIm) ) + M_PI/2.0;
		}
	}
	return Refine();
}

int Point::Refine( )
{
	int it=0;
	double Xnorm, Dnorm;
	
	solNu = sol; // here solNu is the previous solution
	parNu = par; // here parNu is the previous parameter

	std::cout<<"\nIT\tERR\t\tSOLnorm\t\tDIFFnorm\n";
	
	do
	{
		colloc.Init( par, sol );
		colloc.Interpolate( solData, sol );
	
		Jacobian( *jac, *rhs, parNu, par, solNu, sol, solData, varMap, false );
		jac->Solve( *xx, *rhs, dim3 );
		
		Update( *xx );
		
		// computing norms to determine convergence
		Xnorm = sqrt( colloc.Integrate( sol, sol ) );
		Dnorm = sqrt( colloc.Integrate( xx->getV1(), xx->getV1() ) + (xx->getV3())*(xx->getV3()) );
		std::cout<<" "<<it<<"\t"<<Dnorm/(1.0+Xnorm)<<"\t"<<Xnorm<<"\t"<<Dnorm<<'\n';
		std::cout.flush();
	}
	while( (Dnorm/(1.0+Xnorm) >= RefEps)&&(it++ < RefIter) );
	
	return it;
}

void Point::SwitchTFHB( double ds )
{
	if( par(NPAR+ParAngle) < 0 ) par(NPAR+ParAngle) *= -1.0;
	if( par(NPAR+ParAngle) < M_PI ) par(0) *= 2*M_PI/par(NPAR+ParAngle);
	else par(0) *= 2*M_PI/(par(NPAR+ParAngle) - M_PI);
	std::cout<<"SwitchHOPF: T = "<<par(0)<<", angle/2PI = "<<par(NPAR+ParAngle)/(2*M_PI)<<"\n";
	
	TestFunctCPLX *tf = static_cast<TestFunctCPLX*>(testFunct);
	Vector QRE(NDIM), QIM(NDIM);
	tf->SwitchHB( QRE, QIM, colloc );
	
	for( int i = 0; i < NINT; i++ )
	{
		for( int j = 0; j < NDEG+1; j++ )
		{
			const double t = colloc.Profile( i, j );
			for( int p = 0; p < NDIM; p++ )
			{
				xxDot->getV1()( p + (j+i*NDEG)*NDIM ) = sin(2.0*M_PI*t) * QIM(p) + cos(2.0*M_PI*t) * QRE(p);
			}
		}
	}

	double norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) );
	xxDot->getV1() /= norm;
	xxDot->getV3().Clear();
// 	std::cout<<"HOPFtanNorm "<<norm<<"\n";
	
	for( int i = 0; i < NDIM*(NINT*NDEG+1); i++ )
	{
		sol( i ) += ds*xxDot->getV1()(i);
	}
}

/// Switching with the test functionals!!!

void Point::SwitchTFLP( double ds )
{
	Vector qq_(NDIM);
	MatFact DD( NDIM, NDIM );
	Vector wr( NDIM ), wi( NDIM );
	Matrix vr( NDIM, NDIM ), vl( NDIM, NDIM );

	std::cout<<"Point::SwitchTF_LP\n";
	xxDot->getV1().Clear();
	xxDot->getV3().Clear();
	p1Dot = 0.0;
	
	TestFunct* tf = new TestFunct( colloc, 1.0 );
	tf->Funct( colloc, par, sol, solData );
	tf->Switch( xxDot->getV1() );
	delete tf;
	
	sol += ds * xxDot->getV1();
}

void Point::SwitchTFPD( double ds )
{
	Vector qq_(NDIM);
	MatFact DD( NDIM, NDIM );
	Vector wr( NDIM ), wi( NDIM );
	Matrix vr( NDIM, NDIM ), vl( NDIM, NDIM );

	std::cout<<"Point::SwitchTF_PD\n";
	xxDot->getV1().Clear();
	xxDot->getV3().Clear();
	p1Dot = 0.0;
	
	Vector tan(xxDot->getV1().Size());
	TestFunct* tf = new TestFunct( colloc, -1.0 );
	tf->Funct( colloc, par, sol, solData );
	tf->Switch( tan );
	delete tf;
	
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
	xxDot->getV3().Clear();
	sol += ds * xxDot->getV1();
}

void Point::Tangent( )
{
	double norm;
	
	colloc.Init( par, sol );
	colloc.Interpolate( solData, sol );
	
	varMapCont( varMap.Size() ) = p1; // not necessary
	// az RHS-t feleslegesen szamolja ki && the first qq should be qq0
	Jacobian( *jac, *rhs, par, par, sol, sol, solData, varMapCont, false );
	
	rhs->getV1() = jac->getA13( dim3 );
	for( int i = 0; i < dim3; i++ ) rhs->getV3()(i) = jac->getA33( i, dim3 );
	
	jac->Solve( *xxDot, *rhs, dim3 );
	xxDot->getV3()(dim3) = -1.0;
	
	norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV3())*(xxDot->getV3()) );
	
	xxDot->getV1() /= -norm;
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
	
	parNu = par;
	for( int i=0; i< solNu.Size(); i++ )     solNu(i)           = sol(i)           + ds * xxDot->getV1()(i);
	for( int i = 2; i < varMap.Size(); i++ ) parNu( varMap(i) ) = par( varMap(i) ) + ds * xxDot->getV3()(i-2);
	parNu(p1) = par(p1) + ds*p1Dot;
	
	varMapCont( varMap.Size() ) = p1;

	int  it=0;
	bool conv;
	do
	{
		colloc.Init( parNu, solNu );
		colloc.Interpolate( solData, solNu );
		
		Jacobian( *jac, *rhs, par, parNu, sol, solNu, solData, varMapCont, true );
		
		jac->Solve( *xx, *rhs );
		
		ContUpdate( *xx );
		
		Rnorm = sqrt( colloc.Integrate( rhs->getV1(), rhs->getV1() ) + (rhs->getV3())*(rhs->getV3()) );
		Xnorm = sqrt( colloc.Integrate( solNu, solNu ) );
		Dnorm = sqrt( colloc.Integrate( xx->getV1(), xx->getV1() ) + (xx->getV3())*(xx->getV3()) );
		conv = (Dnorm/(1.0+Xnorm) >= ContEps) || (Rnorm/(1.0+Xnorm) >= ContEps);
		
		// updating the tangent
		jac->AX( *rhs, *xxDot );
		rhs->getV3()(dim3) = 0.0;
		jac->Solve( *xx, *rhs );
		xxDot->getV1() -= xx->getV1();
		xxDot->getV3() -= xx->getV3();
		Tnorm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV3())*(xxDot->getV3()) );
		xxDot->getV1() /= Tnorm;
		xxDot->getV3() /= Tnorm;
		p1Dot = xxDot->getV3()(dim3);
		// end updating tangent

	}
	while( conv /*&& (Dnorm/(1.0+Xnorm) < 1.0)*/&&(++it < ContIter) );
	if( !conv )
	{
		/// checking the tangent and the secant
	#ifdef DEBUG
		double Pnorm = sqrt(p1Dot*p1Dot)
		double Xnorm = sqrt(colloc.Integrate( xxDot->getV1(), xxDot->getV1() )), Onorm = sqrt( (xxDot->getV3())*(xxDot->getV3()) );
		std::cout<<"Cnorm: "<<Tnorm<<"\nDot Pnorm: "<<Pnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm;
		for( int i = 2; i < varMap.Size(); i++ ) std::cout<<" O"<<varMap(i)<<": "<<xxDot->getV3()(i-2);
		std::cout<<'\n';
		
		xx->getV1() = solNu;
		xx->getV1() -= sol;
		for( int i = 2; i < varMap.Size(); i++ ) xx->getV3()(i-2) = parNu( varMap(i) ) - par( varMap(i) );
		xx->getV3()(dim3) = parNu(p1) - par(p1);
		
		Pnorm = sqrt( xx->getV3()(dim3)*xx->getV3()(dim3) )/ds;
		Xnorm = sqrt(colloc.Integrate( xx->getV1(), xx->getV1() ))/ds; 
		Onorm = 0;
		for( int i=0; i<dim3+1; i++ ) Onorm += (xx->getV3()(i))*(xx->getV3()(i));
		Onorm = sqrt(Onorm)/ds;
		std::cout<<"Dif Pnorm: "<<Pnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm;
		for( int i = 2; i < varMap.Size(); i++ ) std::cout<<" O"<<varMap(i)<<": "<<xx->getV3()(i-2)/ds;
		std::cout<<'\n';
	#endif
		/// END OF CHECKING
		
		// copying back the solution
		sol = solNu;
		par = parNu;
		
	}else
	{
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
	if( (dminLP < dminPD) && (dminLP < dminNS) ) return BifTFLP;
	else if( (dminPD < dminLP) && (dminPD < dminNS) ) return BifTFPD;
	else if( (dminNS < dminPD) && (dminNS < dminLP) ) return BifTFNS;
	else return SolTF;
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
	if( NPAR != npar_ ) { std::cout<<"Not compatible file (NPAR) "<<npar_<<"\n"; PDError(-1); }
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
		for( int i = 2; i < varMapCont.Size(); i++ ) xxDot->getV3()(i-2) = par( varMapCont(i) ) - parNu( varMapCont(i) );
		
		// renorming everything
		p1Dot = xxDot->getV3()(dim3);
		
		double norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + (xxDot->getV3())*(xxDot->getV3()) );
		
		xxDot->getV1() /= norm;
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
	if( NPAR != npar_ ) { std::cout<<"RN:Not compatible file (NPAR) "<<npar_<<"\n"; PDError(-1); }
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

void Point::SwitchTFTRTan( Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg ) // starting data for tori: tangent
{
	Vector TRe(sol.Size()), TIm(sol.Size());
	TestFunctCPLX* tf = static_cast< TestFunctCPLX* >(testFunct);
	if( tf )
	{
		tf->Switch( TRe, TIm, alpha );
		colloc.Export( Re, mshint, mshdeg, TRe );
		colloc.Export( Im, mshint, mshdeg, TIm );
	}else
	{
		std::cout<<"Not the complex test functional\n";
		PDError(-1);
	}
}