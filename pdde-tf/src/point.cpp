// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "config.h"

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
	var( var_ ), eqn( eqn_ ), varMap( var_.Size() ), varMapCont( var_.Size() + 1 ),
	sol( sys_.ndim()*(nint*ndeg+1) ), par( sys_.npar()+ParEnd ),
	solNu( sys_.ndim()*(nint*ndeg+1) ), parNu( sys_.npar()+ParEnd ),
	mRe(nmul), mIm(nmul),
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
	PtType   type;
	BranchSW sw;
	int      neqn;
	int      nparx;
	Eqn      eqns[5];
	Var      vars[5];
};

// not even member function public
BranchSW PtToEqnVar( Array1D<Eqn>& eqnr, Array1D<Var>& varr, PtType Pt, Array1D<Var> parx, int npar_ )
{
	PtTab tab;
	const Var PANGLE = (Var)(VarPAR0+npar_+ParAngle);
	const Var PX0 = parx(0);
	const Var PX1 = parx(1);
	switch( Pt )
	{
	/// TIME-PERIODIC TEST-FUNCTIONAL
		case SolTF:
			{ PtTab tmp = { SolTF, NOSwitch,   1, 0,
			 { EqnSol },
			 { VarSol } }; tab = tmp; } break;
		case SolTFBRSW:
			{ PtTab tmp = { SolTFBRSW, TFBRSwitch, 1, 0,
			 { EqnSol },
			 { VarSol } }; tab = tmp; } break;
		case SolTFPDSW:
			{ PtTab tmp = { SolTFPDSW, TFPDSwitch, 1, 0,
			 { EqnSol },
			 { VarSol } }; tab = tmp; } break;
		case BifTFLP:
			{ PtTab tmp = { BifTFLP, NOSwitch,   2, 1,
			 { EqnSol, EqnTFLP },
			 { VarSol, PX0 } }; tab = tmp; } break;
		case BifTFPD:
			{ PtTab tmp = { BifTFPD, NOSwitch,   2, 1,
			 { EqnSol, EqnTFPD },
			 { VarSol, PX0 } }; tab = tmp; } break;
		case BifTFNS:
			{ PtTab tmp = { BifTFNS, NOSwitch,   3, 1,
			 { EqnSol, EqnTFCPLX_RE,  EqnTFCPLX_IM },
			 { VarSol, PANGLE,        PX0 } }; tab = tmp; } break;
	/// AUTONOMOUS TEST-FUNCTIONAL
		case SolTFAUT:
			{ PtTab tmp = { SolTFAUT, NOSwitch,   2, 1,
			 { EqnSol, EqnPhase },
			 { VarSol, PX0 } }; tab = tmp; } break;
		case SolTFAUTBRSW:
			{ PtTab tmp = { SolTFAUTBRSW, TFBRSwitch, 2, 1,
			 { EqnSol, EqnPhase },
			 { VarSol, PX0 } }; tab = tmp; } break;
		case SolTFAUTPDSW:
			{ PtTab tmp = { SolTFAUTPDSW, TFPDSwitch, 2, 1,
			 { EqnSol, EqnPhase },
			 { VarSol, PX0 } }; tab = tmp; } break;
		case SolTFAUTHBSW:
			{ PtTab tmp = { SolTFAUTHBSW, TFHBSwitch, 2, 1,
			 { EqnSol, EqnPhase },
			 { VarSol, PX0 } }; tab = tmp; } break;
		case BifTFAUTLP:
			{ PtTab tmp = { BifTFAUTLP, NOSwitch,   3, 2,
			 { EqnSol, EqnPhase,      EqnTFLPAUT },
			 { VarSol, PX0,           PX1 } }; tab = tmp; } break;
		case BifTFAUTPD:
			{ PtTab tmp = { BifTFAUTPD, NOSwitch,   3, 2,
			 { EqnSol, EqnPhase,      EqnTFPD },
			 { VarSol, PX0,           PX1 } }; tab = tmp; } break;
		case BifTFAUTNS:
			{ PtTab tmp = { BifTFAUTNS, NOSwitch,   4, 2,
			 { EqnSol, EqnPhase,      EqnTFCPLX_RE,  EqnTFCPLX_IM },
			 { VarSol, PANGLE,        PX0,           PX1 } }; tab = tmp; } break;
	/// TORUS
		case SolTor:
			{ PtTab tmp = { SolTor, TFTRSwitch, 2, 1,
			 { EqnTORSol, EqnTORPhase1 },
			 { VarTORSol, PX0 } }; tab = tmp; } break;
		case SolAUTTor:
			{ PtTab tmp = { SolAUTTor, TFTRSwitch, 2, 2,
			 { EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
			 { VarTORSol, PX0,          PX1 } }; tab = tmp; } break;
		default:
			{ PtTab tmp = { SolUser, NOSwitch,  0, 0, { EqnNone }, { VarNone } }; tab = tmp; }
			std::cout<<"No such pointtype\n"; PDError(-1);
			break;
	}
	if( tab.nparx != parx.Size() ) { std::cout<<"Error: wrong number of parameters\n"; PDError(1); }
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
	if( (eqn.Size() == 0)||(var.Size() == 0)||(eqn.Size() != var.Size()) )
	{
		std::cout<<"Bad equation and variable sizes!";
		PDError(-1);
	}
	else
	{
		dim3 = eqn.Size() - 1;
	}
	
	if( (eqn(0) != EqnSol)||(var(0) != VarSol) )
	{
		std::cout<<"Wrong first equation!";
		PDError(-1);
	}
	else
	{
		dim1 = NDIM*(NINT*NDEG+1);
	}
	
	testFunct = 0;
	for( int i = 1; i < eqn.Size(); i++ )
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
	
	for( int i = 1; i < var.Size(); i++ )
	{
		if( (var(i)-VarPAR0 >= 0)&&(var(i)-VarPAR0 < NPAR+ParEnd) )
		{
			varMap( i ) = var(i)-VarPAR0;
		}
		else
		{
			std::cout<<"{1} Non-existing parameter P"<<var(i)-VarPAR0<<" was specified at position "<<i<<".\n";
			PDError(-1);
		}
	}
	for( int i = 0; i < var.Size(); i++ ) varMapCont(i) = varMap(i);
	varMapCont( varMap.Size() ) = p1;
	
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

// **************************************************************************************************************** //
//
//                           THE JACOBIAN
//
// **************************************************************************************************************** //

void Point::Jacobian(
		HyperMatrix& AA, HyperVector& RHS,                      // output
		Vector& parPrev, Vector& par,                           // parameters
		Vector& solPrev, Vector& sol, JagMatrix3D& solData,     // solution
		Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
		double ds, bool cont )                                             // cont: true if continuation
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
		
		for( int i = 1; i < varMap.Size(); i++ )
		{
			if( varMap(i) < NPAR )
			{
				colloc.RHS_p( AA.getA13(i-1), par, sol, solData, varMap(i) );
			}
			else if( varMap(i)-NPAR == ParAngle )
			{
				AA.getA13(i-1).Clear();
			} else
			{
				std::cout<<"{2} Non-existing parameter P"<<varMap(i)<<" was specified at position "<<i<<".\n";
				PDError(-1);
			}
		}
	}

// ---------------------------------------------------------------------------------------------------------------- //
//
//                           Other equations
//
// ---------------------------------------------------------------------------------------------------------------- //

	for( int i=1; i<eqn.Size(); i++ )
	{
		switch( eqn(i) )
		{
			// Phase conditions.
			case EqnPhase:
				colloc.PhaseStar( AA.getA31(i-1), solPrev ); // this should be the previous solution!!!
				// other variables
				for( int j = 1; j < varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						AA.getA33()(i-1,j-1) = 0.0;
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						AA.getA33(i-1,j-1) = 0.0;
					} else
					{
						std::cout<<"{2} Non-existing parameter P"<<varMap(j)<<" was specified at position "<<j<<".\n";
						PDError(-1);
					}
				}
				RHS.getV3()(i-1) = -( AA.getA31(i-1)*sol );
				break;
			case EqnPhaseRot:
				colloc.PhaseRotStar( AA.getA31(i-1), solPrev, rotRe, rotIm ); // this should be the previous solution!!!
				// other variables
				for( int j = 1; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						AA.getA33()(i-1,j-1) = 0.0;
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						AA.getA33(i-1,j-1) = 0.0;
					} else
					{
						std::cout<<"{3} Non-existing parameter P"<<varMap(j)<<" was specified at position "<<j<<".\n";
						PDError(-1);
					}
				}
				RHS.getV3()(i-1) = -( AA.getA31(i-1)*sol );
				break;
			case EqnTFLP:
			case EqnTFPD:
			case EqnTFLPAUT:
			case EqnTFLPAUTROT:
				RHS.getV3()(i-1) = testFunct->Funct( colloc, par, sol, solData );
				testFunct->Funct_x( AA.getA31(i-1), colloc, par, sol, solData );
				for( int j = 1; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						AA.getA33()(i-1,j-1) = testFunct->Funct_p( colloc, par, sol, solData, varMap(j) );
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						AA.getA33(i-1,j-1) = 0.0;
					} else
					{
						std::cout<<"{4} Non-existing parameter P"<<varMap(j)<<" was specified at position "<<j<<".\n";
						PDError(-1);
					}
				}
				break;
			case EqnTFCPLX_RE:
				if( eqn(i+1) != EqnTFCPLX_IM )
				{
					std::cout<<"EqnTFCPLX_RE is not paired."<<varMap(i)<<"\n";
					PDError(-1);
				}
				testFunct->Funct( RHS.getV3()(i-1), RHS.getV3()(i), colloc, par, sol, solData, 
				                  cos(par(NPAR+ParAngle)), sin(par(NPAR+ParAngle)) );
				testFunct->Funct_x( AA.getA31(i-1), AA.getA31(i), colloc, par, sol, solData );
				for( int j = 1; j<varMap.Size(); j++ )
				{
					if( varMap(j) < NPAR )
					{
						testFunct->Funct_p( AA.getA33()(i-1,j-1), AA.getA33()(i,j-1), colloc, par, sol, solData, varMap(j) );
					}
					else if( varMap(j)-NPAR == ParAngle )
					{
						testFunct->Funct_z( AA.getA33()(i-1,j-1), AA.getA33()(i,j-1), colloc, par, sol, solData );
					} else
					{
						std::cout<<"{5} Non-existing parameter P"<<varMap(j)<<" was specified at position "<<j<<".\n";
						PDError(-1);
					}
				}
				break;
			case EqnTFCPLX_IM:
				if( eqn(i-1) != EqnTFCPLX_RE )
				{
					std::cout<<"EqnTFCPLX_IM is not paired."<<varMap(i)<<"\n"; PDError(-1);
				}
				break;
			default:
				std::cout<<"Unknown equation type encountered.";
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
		
		if( ds != 0.0 )
		{
			RHS.getV3()(dim3) = ds - colloc.IntegrateCont( xxDot->getV1(), sol, solPrev );
			for( int j = 1; j < varMap.Size(); ++j ) RHS.getV3()(dim3) -= xxDot->getV3()(j)*(par(varMap(j))-parPrev(varMap(j)));
		}else
		{
			RHS.getV3()(dim3) = 0.0;
		}
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
	for( int i = 1; i < varMap.Size(); i++ ) par( varMap(i) ) += X.getV3()(i-1);
}

inline void Point::ContUpdate( HyperVector& X )
{
	solNu += X.getV1();
	for( int i = 1; i < varMapCont.Size(); i++ ) parNu( varMapCont(i) ) += X.getV3()(i-1);
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
	
		Jacobian( *jac, *rhs, parNu, par, solNu, sol, solData, varMap, 0.0, false );
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
	TestFunctCPLX *tf = static_cast<TestFunctCPLX*>(testFunct);
	Vector QRE(NDIM), QIM(NDIM);
	tf->SwitchHB( QRE, QIM, colloc, par );

	if( par(NPAR+ParAngle) < 0 ) par(NPAR+ParAngle) *= -1.0;
	if( par(NPAR+ParAngle) < M_PI ) par(0) *= 2*M_PI/par(NPAR+ParAngle);
	else par(0) *= 2*M_PI/(par(NPAR+ParAngle) - M_PI);
	std::cout<<"SwitchHOPF: T = "<<par(0)<<", angle/2PI = "<<par(NPAR+ParAngle)/(2*M_PI)<<"\n";

#ifdef DEBUG
	std::ofstream file( "neweigenvec" );
	file<<std::scientific;
	file.precision(12);
#endif

	for( int i = 0; i < NINT; i++ )
	{
		for( int j = 0; j < NDEG+1; j++ )
		{
			const double t = colloc.Profile( i, j );
			for( int p = 0; p < NDIM; p++ )
			{
				xxDot->getV1()( p + (j+i*NDEG)*NDIM ) = cos(2.0*M_PI*t) * QRE(p) + sin(2.0*M_PI*t) * QIM(p);
			#ifdef DEBUG
				file<<xxDot->getV1()( p + (j+i*NDEG)*NDIM )<<"\t";
			#endif
			}
		#ifdef DEBUG
			file<<par(0)*t<<"\n";
		#endif
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
	Jacobian( *jac, *rhs, par, par, sol, sol, solData, varMapCont, 0.0, false );
	
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

int Point::Continue( double ds, bool jacstep )
{
	double Xnorm, Dnorm, Rnorm, Tnorm;
	
	parNu = par;
	for( int i = 0; i < solNu.Size(); i++ )  solNu(i)           = sol(i)           + ds * xxDot->getV1()(i);
	for( int i = 1; i < varMap.Size(); i++ ) parNu( varMap(i) ) = par( varMap(i) ) + ds * xxDot->getV3()(i-1);
	parNu(p1) = par(p1) + ds*p1Dot;
	
	varMapCont( varMap.Size() ) = p1;

	int  it=0;
	bool conv;
	do
	{
		colloc.Init( parNu, solNu );
		colloc.Interpolate( solData, solNu );
		
		if( jacstep ) Jacobian( *jac, *rhs, par, parNu, sol, solNu, solData, varMapCont, ds, true );
		else          Jacobian( *jac, *rhs, par, parNu, sol, solNu, solData, varMapCont, 0.0, true );
		
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
	#ifdef DEBUG
		/// checking the tangent and the secant
		double Pnorm = sqrt(p1Dot*p1Dot);
		double Xnorm = sqrt(colloc.Integrate( xxDot->getV1(), xxDot->getV1() )), Onorm = sqrt( (xxDot->getV3())*(xxDot->getV3()) );
		std::cout<<"Cnorm: "<<Tnorm<<"\nDot Pnorm: "<<Pnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm;
		for( int i = 1; i < varMap.Size(); i++ ) std::cout<<" O"<<varMap(i)<<": "<<xxDot->getV3()(i-1);
		std::cout<<'\n';
		
		xx->getV1() = solNu;
		xx->getV1() -= sol;
		for( int i = 1; i < varMap.Size(); i++ ) xx->getV3()(i-1) = parNu( varMap(i) ) - par( varMap(i) );
		xx->getV3()(dim3) = parNu(p1) - par(p1);
		
		Pnorm = sqrt( xx->getV3()(dim3)*xx->getV3()(dim3) )/ds;
		Xnorm = sqrt(colloc.Integrate( xx->getV1(), xx->getV1() ))/ds; 
		Onorm = 0;
		for( int i = 0; i < dim3+1; i++ ) Onorm += (xx->getV3()(i))*(xx->getV3()(i));
		Onorm = sqrt(Onorm)/ds;
		std::cout<<"Dif Pnorm: "<<Pnorm<<" Xnorm: "<<Xnorm<<" Onorm: "<<Onorm;
		for( int i = 1; i < varMap.Size(); i++ ) std::cout<<" O"<<varMap(i)<<": "<<xx->getV3()(i-2)/ds;
		std::cout<<'\n';
		/// END OF CHECKING
	#endif
		
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

///--------------------------------
/// INPUT & OUTPUT
///--------------------------------

void Point::Write( std::ofstream& file )
{
	Vector msh( NDEG*NINT+1 );
	colloc.getMesh( msh );
	
	file<<NPAR+ParEnd<<"\t";
	for( int i = 0; i < NPAR+ParEnd; i++ ) file<<par(i)<<"\t";
	
	file<<mRe.Size()<<"\t";
	for( int i = 0; i < mRe.Size(); i++ ) file<<mRe(i)<<"\t"<<mIm(i)<<"\t";
	
	file<<NDIM<<"\t";
	file<<NINT<<"\t";
	file<<NDEG<<"\t";
	for( int i = 0; i < NDEG*NINT+1; i++ ) file<<msh(i)<<"\t";
	for( int i = 0; i < NDIM*(NINT*NDEG+1); i++ ) file<<sol(i)<<"\t";
	file<<"\n";
	file.flush();
}

void Point::Read( std::ifstream& file, bool tan )
{
	int npar_, nmul_, ndim_, nint_, ndeg_;
	file>>npar_;
	if( NPAR+ParEnd != npar_ ) { std::cout<<"Not compatible file (NPAR) "<<npar_<<"\n"; PDError(-1); }
	for( int i = 0; i < NPAR+ParEnd; i++ ) file>>par(i);
	
	file>>nmul_;
	if( mRe.Size() < nmul_ ) { std::cout<<"Not compatible file (NMUL) "<<nmul_<<"\n"; PDError(-1); }
	for( int i = 0; i < nmul_; i++ ) { file>>mRe(i); file>>mIm(i); }
	
	file>>ndim_;
	file>>nint_;
	file>>ndeg_;
	
	if( NDIM != ndim_ ) { std::cout<<"Not compatible file (NDIM) "<<ndim_<<"\n"; PDError(-1); }
// 	if( NINT != nint_ ) { std::cout<<"Not compatible file (NINT) "<<nint_<<"\n"; PDError(-1); }
// 	if( NDEG != ndeg_ ) { std::cout<<"Not compatible file (NDEG) "<<ndeg_<<"\n"; PDError(-1); }
	
	Vector msh( ndeg_*nint_ + 1 );
	for( int i = 0; i < ndeg_*nint_+1; i++ ) file>>msh(i);
	
	if( (NINT == nint_)&&(NDEG == ndeg_) )
	{
		colloc.setMesh( msh );
		for( int i = 0; i < NDIM*(nint_*ndeg_+1); i++ ) file>>sol(i);
	}
	else
	{
		Vector in( NDIM*(nint_*ndeg_+1) );
	
		for( int i = 0; i < NDIM*(nint_*ndeg_+1); i++ ) file>>in(i);
		colloc.Import( sol, in, msh, ndeg_ );
	}
	if( tan )
	{
		solNu = sol;
		parNu = par;
		Read( file, false );
		xxDot->getV1() = sol;
		xxDot->getV1() -= solNu;
		for( int i = 1; i < varMapCont.Size(); i++ ) xxDot->getV3()(i-1) = par( varMapCont(i) ) - parNu( varMapCont(i) );
		
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
	if( NPAR+ParEnd != npar_ ) { std::cout<<"RN:Not compatible file (NPAR) "<<npar_<<"\n"; PDError(-1); }
	for( int i = 0; i < NPAR+ParEnd; i++ ) file>>tmp;
	
	file>>nmul_;
	if( mRe.Size() < nmul_ ) { std::cout<<"RN:Not compatible file (NMUL) "<<nmul_<<"\n"; PDError(-1); }
	for( int i = 0; i < nmul_; i++ ) { file>>tmp; file>>tmp; }
	
	file>>ndim_;
	file>>nint_;
	file>>ndeg_;
	
	if( NDIM != ndim_ ) { std::cout<<"RN:Not compatible file (NDIM) "<<ndim_<<"\n"; PDError(-1); }
	
	for( int i = 0; i < ndeg_*nint_+1; i++ ) file>>tmp;
	for( int i = 0; i < NDIM*(nint_*ndeg_+1); i++ ) file>>tmp;
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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <sys/mman.h>

void Point::BinaryWrite( mat4Data& data, int n )
{
	Vector msh( NDEG*NINT+1 );
	colloc.getMesh( msh );
	
	data.setPar( n, par );
	data.setMul( n, mRe, mIm );
	data.setMesh( n, msh );
	data.setProfile( n, sol );
}

void Point::BinaryRead( mat4Data& data, int n )
{
	Vector msh( data.getNInt()*data.getNDeg()+1 );
	if( data.getNPar() == (NPAR+ParEnd) )
	{
		data.getPar( n, par );
	}else
	{
		std::cout<<"Wrong number of parameters\n";
		PDError(-1);
	}
	data.getMul( n, mRe, mIm );
	data.getMesh( n, msh );
	if( data.getNDim() == NDIM )
	{
		if( data.getNInt() == NINT && data.getNDeg() == NDEG )
		{
			colloc.setMesh( msh );
			data.getProfile( n, sol );
			std::cout<<"not converted\n";
		}else
		{
			Vector tmp( data.getNDim()*(data.getNDeg()*data.getNInt()+1) );
			data.getProfile( n, tmp );
			colloc.Import( sol, tmp, msh, data.getNDeg() );
		}
	}else
	{
		std::cout<<"binaryread failed";
		PDError(-1);
	}
	par.Print();
}

int mat4Data::findMatrix( const char* name, mat4Data::header* found, int *sz )
{
	struct header hd;
	int cur_off = 0;
	int cur_size;
	do{
		memcpy( &hd, (char*)address + cur_off, sizeof(struct header) );
		if( hd.type != 0 )
		  { std::cout<<"not a double matrix\n"; PDError(-1); }
		if( hd.imagf == 0 )
		  cur_size = sizeof(struct header) + hd.namelen*sizeof(char) + hd.mrows*hd.ncols*sizeof(double);
		else
		  cur_size = sizeof(struct header) + hd.namelen*sizeof(char) + 2*hd.mrows*hd.ncols*sizeof(double);
		if( strncmp( name, (char*)address + cur_off + sizeof(struct header), 20 ) == 0 )
		{
			memcpy( found, &hd, sizeof(struct header) );
			*sz = cur_size;
			return cur_off;
		}
		cur_off += cur_size;
	}while( cur_off < size );
	return -1;
}

mat4Data::mat4Data( const std::string& fileName, int steps_, int ndim_, int npar_, int nint_, int ndeg_, int nmul_ )
{
	wperm = true;
	
	ncols = steps_;
	ndim = ndim_;
	npar = npar_;
	nint = nint_;
	ndeg = ndeg_;
	nmul = nmul_;
	
	par_offset = 0;
	strncpy( par_name, "pdde_par", 20 );
	par_header.type = 0;
	par_header.mrows = npar;
	par_header.ncols = ncols;
	par_header.imagf = 0;
	par_header.namelen = ((strlen(par_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	par_size = sizeof(struct header) + par_header.namelen + ncols * npar * sizeof(double);
	
	mul_offset = par_offset + par_size;
	strncpy( mul_name, "pdde_mul", 20 );
	mul_header.type = 0;
	mul_header.mrows = nmul;
	mul_header.ncols = ncols;
	mul_header.imagf = 1;
	mul_header.namelen = ((strlen(mul_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	mul_size = sizeof(struct header) + mul_header.namelen + ncols * 2 * nmul * sizeof(double);
	
	ndim_offset = mul_offset + mul_size;
	strncpy( ndim_name, "pdde_ndim", 20 );
	ndim_header.type = 0;
	ndim_header.mrows = 1;
	ndim_header.ncols = ncols;
	ndim_header.imagf = 0;
	ndim_header.namelen = ((strlen(ndim_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	ndim_size = sizeof(struct header) + ndim_header.namelen + ncols * 1 * sizeof(double);

	nint_offset = ndim_offset + ndim_size;
	strncpy( nint_name, "pdde_nint", 20 );
	nint_header.type = 0;
	nint_header.mrows = 1;
	nint_header.ncols = ncols;
	nint_header.imagf = 0;
	nint_header.namelen = ((strlen(nint_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	nint_size = sizeof(struct header) + nint_header.namelen + ncols * 1 * sizeof(double);

	ndeg_offset = nint_offset + nint_size;
	strncpy( ndeg_name, "pdde_ndeg", 20 );
	ndeg_header.type = 0;
	ndeg_header.mrows = 1;
	ndeg_header.ncols = ncols;
	ndeg_header.imagf = 0;
	ndeg_header.namelen = ((strlen(ndeg_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	ndeg_size = sizeof(struct header) + ndeg_header.namelen + ncols * 1 * sizeof(double);

	mesh_offset = ndeg_offset + ndeg_size;
	strncpy( mesh_name, "pdde_mesh", 20 );
	mesh_header.type = 0;
	mesh_header.mrows = ndeg*nint+1;
	mesh_header.ncols = ncols;
	mesh_header.imagf = 0;
	mesh_header.namelen = ((strlen(mesh_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	mesh_size = sizeof(struct header) + mesh_header.namelen + ncols * (ndeg*nint+1) * sizeof(double);

	prof_offset = mesh_offset + mesh_size;
	strncpy( prof_name, "pdde_prof", 20 );
	prof_header.type = 0;
	prof_header.mrows = ndim*(ndeg*nint+1);
	prof_header.ncols = ncols;
	prof_header.imagf = 0;
	prof_header.namelen = ((strlen(prof_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	prof_size = sizeof(struct header) + prof_header.namelen + ncols * ndim*(ndeg*nint+1) * sizeof(double);
	size = prof_offset + prof_size;
	
	// creating the mmapped file
	if( ( file = open( fileName.c_str(), O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR ) ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to open file\n"); throw(-1); }
	
	if( lseek( file, size-1, SEEK_SET ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to seek file\n"); throw(-1); }
	
	if( write( file, "\0", 1 ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to write file\n"); throw(-1); }
	
	if( ( address = mmap( 0, size, PROT_WRITE, MAP_SHARED, file, 0 ) ) == MAP_FAILED )
	{ perror("mmappedPointData::mmappedPointData: unable to mmap file\n"); throw(-1); }

	*((struct header*)( (char*)address + par_offset) ) = par_header;
	*((struct header*)( (char*)address + mul_offset) ) = mul_header;
	*((struct header*)( (char*)address + ndim_offset) ) = ndim_header;
	*((struct header*)( (char*)address + nint_offset) ) = nint_header;
	*((struct header*)( (char*)address + ndeg_offset) ) = ndeg_header;
	*((struct header*)( (char*)address + mesh_offset) ) = mesh_header;
	*((struct header*)( (char*)address + prof_offset) ) = prof_header;
	strncpy( (char*)address + par_offset + sizeof(struct header), par_name, 20 );
	strncpy( (char*)address + mul_offset + sizeof(struct header), mul_name, 20 );
	strncpy( (char*)address + ndim_offset + sizeof(struct header), ndim_name, 20 );
	strncpy( (char*)address + nint_offset + sizeof(struct header), nint_name, 20 );
	strncpy( (char*)address + ndeg_offset + sizeof(struct header), ndeg_name, 20 );
	strncpy( (char*)address + mesh_offset + sizeof(struct header), mesh_name, 20 );
	strncpy( (char*)address + prof_offset + sizeof(struct header), prof_name, 20 );

// 	std::cout<<"NDIM "<<ndim<<" NINT "<<nint<<" NDEG "<<ndeg<<" NPAR "<<npar<<" NMUL "<<nmul<<"\n";
}

mat4Data::mat4Data( const std::string& fileName )
{
	wperm = false;
	this->openReadOnly( fileName );
}

void mat4Data::openReadOnly( const std::string& fileName )
{
	if( ( file = open( fileName.c_str(), O_RDONLY ) ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to open file\n"); throw(-1); }
	
	struct stat filestat;
	if( fstat( file, &filestat ) != 0 )
	{ perror("mmappedPointData::mmappedPointData: unable to stat file\n"); throw(-1); }
	filesize = filestat.st_size;
	size = filesize;
	
	if( ( address = mmap( 0, filesize, PROT_READ, MAP_PRIVATE, file, 0 ) ) == MAP_FAILED )
	{ perror("mmappedPointData::mmappedPointData: unable to mmap file\n"); throw(-1); }
	
	if( (par_offset = findMatrix( "pdde_par", &par_header, &par_size )) == -1 ) { std::cout<<"err1"; PDError(-1); }
	npar = par_header.mrows; std::cout<<"NPAR "<<npar<<"\n";
	ncols = par_header.ncols;
	if( par_header.imagf != 0 ) { std::cout<<"err2"; PDError(-1); }

	if( (mul_offset = findMatrix( "pdde_mul", &mul_header, &mul_size )) == -1 ) { std::cout<<"err3"; PDError(-1); }
	nmul = mul_header.mrows; std::cout<<"NMUL "<<nmul<<"\n";
	if( mul_header.ncols != ncols ) { std::cout<<"err4"; PDError(-1); }
	if( mul_header.imagf == 0 ) { std::cout<<"err5"; PDError(-1); }
	
	if( (ndim_offset = findMatrix( "pdde_ndim", &ndim_header, &ndim_size )) == -1 ) { std::cout<<"err6"; PDError(-1); }
	if( ndim_header.mrows != 1 ) { std::cout<<"err7 "<<ndim_header.mrows<<"\n"; PDError(-1); }
	if( ndim_header.ncols != ncols ) { std::cout<<"err9"; PDError(-1); }
	if( ndim_header.imagf != 0 ) { std::cout<<"err9"; PDError(-1); }
	ndim = *((double*)((char*)address + ndim_offset + ndim_header.col_off(0)));
// 	std::cout<<"NDIM "<<ndim<<"\n";
	
	if( (nint_offset = findMatrix( "pdde_nint", &nint_header, &nint_size )) == -1 ) { std::cout<<"err10"; PDError(-1); }
	if( nint_header.mrows != 1 ) { std::cout<<"err11"; PDError(-1); }
	if( nint_header.ncols != ncols ) { std::cout<<"err12"; PDError(-1); }
	if( nint_header.imagf != 0 ) { std::cout<<"err13"; PDError(-1); }
	nint = *((double*)((char*)address + nint_offset + nint_header.col_off(0)));
// 	std::cout<<"NINT "<<nint<<"\n";
	
	if( (ndeg_offset = findMatrix( "pdde_ndeg", &ndeg_header, &ndeg_size )) == -1 ) { std::cout<<"err14"; PDError(-1); }
	if( ndeg_header.mrows != 1 ) { std::cout<<"err15"; PDError(-1); }
	if( ndeg_header.ncols != ncols ) { std::cout<<"err16"; PDError(-1); }
	if( ndeg_header.imagf != 0 ) { std::cout<<"err17"; PDError(-1); }
	ndeg = *((double*)((char*)address + ndeg_offset + ndeg_header.col_off(0)));
// 	std::cout<<"NDEG "<<ndeg<<"\n";

	if( (mesh_offset = findMatrix( "pdde_mesh", &mesh_header, &mesh_size )) == -1 ) { std::cout<<"err18"; PDError(-1); }
	if( mesh_header.mrows != ndeg*nint+1 ) { std::cout<<"err19 "<<mesh_header.mrows<<"\n"; PDError(-1); }
	if( mesh_header.ncols != ncols ) { std::cout<<"err20"; PDError(-1); }
	if( mesh_header.imagf != 0 ) { std::cout<<"err21"; PDError(-1); }
	
	if( (prof_offset = findMatrix( "pdde_prof", &prof_header, &prof_size )) == -1 ) { std::cout<<"err22"; PDError(-1); }
	if( prof_header.mrows != ndim*(ndeg*nint+1) ) { std::cout<<"err23"; PDError(-1); }
	if( prof_header.ncols != ncols ) { std::cout<<"err24"; PDError(-1); }
	if( prof_header.imagf != 0 ) { std::cout<<"err25"; PDError(-1); }
	
// 	std::cout<<"NDIM "<<ndim<<" NINT "<<nint<<" NDEG "<<ndeg<<" NPAR "<<npar<<" NMUL "<<nmul<<"\n";
}

mat4Data::~mat4Data()
{
	if( munmap( address, size ) != 0 )
	{ perror("mmappedPointData::~mmappedPointData: unable to munmap file\n"); throw(-1); }
	if( close( file ) != 0 )
	{ perror("mmappedPointData::~mmappedPointData: unable to close file\n"); throw(-1); }
}

void mat4Data::setPar( int n, const Vector& par )
{
	if( wperm && n < ncols )
	{
		if( par.Size() <= npar )
		{
			for( int i = 0; i < par.Size(); ++i)
				((double*)( (char*)address + par_offset + par_header.col_off(n) ))[i] = par(i);
			for( int i = par.Size(); i < npar; ++i)
				((double*)( (char*)address + par_offset + par_header.col_off(n) ))[i] = 0.0;
		}else
		{
			std::cout<<"setPar 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setPar 2";
		PDError(-1);
	}
}

void mat4Data::setMul( int n, const Vector& re, const Vector& im )
{
	if( wperm && n < ncols && re.Size() == im.Size() )
	{
		if( re.Size() <= nmul )
		{
			for( int i = 0; i < re.Size(); ++i )
			{
				((double*)( (char*)address + mul_offset + mul_header.col_off(n) ))[i] = re(i);
				((double*)( (char*)address + mul_offset + mul_header.col_off_im(n) ))[i] = im(i);
			}
			for( int i = re.Size(); i < nmul; ++i )
			{
				((double*)( (char*)address + mul_offset + mul_header.col_off(n) ))[i] = 0.0;
				((double*)( (char*)address + mul_offset + mul_header.col_off_im(n) ))[i] = 0.0;
			}
		}else
		{
			std::cout<<"setMul 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setMul 2";
		PDError(-1);
	}
}

void mat4Data::setMesh( int n, const Vector& mesh )
{
	if( wperm && n < ncols )
	{
		if( mesh.Size() == ndeg*nint+1 )
		{
			for( int i = 0; i < ndeg*nint+1; ++i)
				((double*)( (char*)address + mesh_offset + mesh_header.col_off(n) ))[i] = mesh(i);
		}else
		{
			std::cout<<"setMesh 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setMesh 2";
		PDError(-1);
	}
}

void mat4Data::setProfile( int n, const Vector& prof )
{
	if( wperm && n < ncols )
	{
		if( prof.Size() == ndim*(ndeg*nint+1) )
		{
			((double*)((char*)address + ndim_offset + ndim_header.col_off(n)))[0] = ndim;
			((double*)((char*)address + nint_offset + nint_header.col_off(n)))[0] = nint;
			((double*)((char*)address + ndeg_offset + ndeg_header.col_off(n)))[0] = ndeg;
			for( int i = 0; i < ndim*(ndeg*nint+1); ++i)
				((double*)( (char*)address + prof_offset + prof_header.col_off(n) ))[i] = prof(i);
		}else
		{
			std::cout<<"setProf 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setProf 2";
		PDError(-1);
	}
}

void mat4Data::getPar( int n, Vector& par )
{
	if( n < ncols )
	{
		if( par.Size() <= npar )
		{
			for( int i = 0; i < par.Size(); ++i)
				par(i) = ((double*)( (char*)address + par_offset + par_header.col_off(n) ))[i];
		}else
		{
			std::cout<<"getPar 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getPar 2";
		PDError(-1);
	}
}

void mat4Data::getMul( int n, Vector& re, Vector& im )
{
	if( n < ncols && re.Size() == im.Size() )
	{
		if( re.Size() <= nmul )
		{
			for( int i = 0; i < re.Size(); ++i )
			{
				re(i) = ((double*)( (char*)address + mul_offset + mul_header.col_off(n) ))[i];
				im(i) = ((double*)( (char*)address + mul_offset + mul_header.col_off_im(n) ))[i];
			}
		}else
		{
			std::cout<<"getMul 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getMul 2";
		PDError(-1);
	}
}

void mat4Data::getMesh( int n, Vector& mesh )
{
	if( n < ncols )
	{
		if( mesh.Size() == ndeg*nint+1 )
		{
			for( int i = 0; i < ndeg*nint+1; ++i)
				mesh(i) = ((double*)( (char*)address + mesh_offset + mesh_header.col_off(n) ))[i];
		}else
		{
			std::cout<<"getMesh 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getMesh 2";
		PDError(-1);
	}
}

void mat4Data::getProfile( int n, Vector& prof )
{
	if( n < ncols )
	{
		if( prof.Size() == ndim*(ndeg*nint+1) )
		{
			for( int i = 0; i < ndim*(ndeg*nint+1); ++i)
				prof(i) = ((double*)( (char*)address + prof_offset + prof_header.col_off(n) ))[i];
		}else
		{
			std::cout<<"getProf 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getProf 2";
		PDError(-1);
	}
}
