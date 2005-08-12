// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cmath>

#include "system.h"
#include "torpoint.h"

PointTR::PointTR( System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int ndeg1_, int ndeg2_, int nint1_, int nint2_ ) :
	var(var_), eqn(eqn_),
	colloc( sys_, ndeg1_, ndeg2_, nint1_, nint2_ ),
	par( sys_.npar() + ParEnd ),
	parNu( sys_.npar() + ParEnd ) //,
//	omega(2), omegaNu(2)
{
	Construct();
}

PointTR::~PointTR()
{
	Destruct();
}

#define NITER 10
#define ERRTOL 1e-4
#define ERRTOLCONT 1e-4
#define BIFTOL 1e-6  // only for starting bifurcations

#define NDEG1 colloc.Ndeg1()
#define NDEG2 colloc.Ndeg2()
#define NINT1 colloc.Nint1()
#define NINT2 colloc.Nint2()
#define SYS   colloc.Sys()

#define RHO (SYS.npar()+ParRot)

void PointTR::Construct()
{
	if( eqn(0) != EqnTORSol ){ std::cout<<"Bad variables\n"; throw(-1); }

	RefEps = ERRTOL;
	ContEps = ERRTOLCONT;
	StartEps = BIFTOL;
	Iter = NITER;
	
	jac   = new HyperMatrix( SYS.ndim()*((NDEG1*NINT1))*((NDEG2*NINT2)), 0, eqn.Size(),
	               SYS.ndim()*(NDEG1*NINT1)*(NDEG2*NINT2) * SYS.ndim()*(SYS.ntau()+1)*(NDEG1+1)*(NDEG2+1) );
	xx    = new HyperVector( SYS.ndim()*((NDEG1*NINT1))*((NDEG2*NINT2)), 0, eqn.Size() );
	xxNu  = new HyperVector( SYS.ndim()*((NDEG1*NINT1))*((NDEG2*NINT2)), 0, eqn.Size() );
	Dxx   = new HyperVector( SYS.ndim()*((NDEG1*NINT1))*((NDEG2*NINT2)), 0, eqn.Size() );
	rhs   = new HyperVector( SYS.ndim()*((NDEG1*NINT1))*((NDEG2*NINT2)), 0, eqn.Size() );
	xxDot = new HyperVector( SYS.ndim()*((NDEG1*NINT1))*((NDEG2*NINT2)), 0, eqn.Size() );
}

void PointTR::Destruct()
{
	delete xxDot;
	delete rhs;
	delete Dxx;
	delete xxNu;
	delete xx;
	delete jac;
}

//
// here starts the real computation
//

void PointTR::ToJacVar( Array1D<Vector*>& A13, Array1D<int>& JacVar, bool cont )
{
	// it is just for the first equation i.e. EqnTORSol
	const int nvar = var.Size()-1; // number of border lines
	for( int i = 1; i < var.Size(); i++ )
	{
		JacVar(i-1) = var(i) - VarPAR0;
		A13(i-1) = & jac->getA13(i-1);
	}
	if( cont )
	{
		JacVar( nvar ) = p1;
		A13( nvar ) = & jac->getA13( nvar );
	}
}

void PointTR::JacobianFixed( Vector& /*sol*/, Vector& presol, Vector& /*par*/, Vector& /*prepar*/, bool cont )
{
	// first the phase conditions
	int ph0 = -1, ph1 = -1;
	for( int i = 1; i < eqn.Size(); i++ )
	{
		if( eqn(i) == EqnTORPhase0 ){ if( ph0==(-1) ) { ph0 = i-1; } else { throw(-1); } }
		if( eqn(i) == EqnTORPhase1 ){ if( ph1==(-1) ) { ph1 = i-1; } else { throw(-1); } }
	}
	if( ph0!=(-1) && ph1!=(-1) )
	{
		colloc.PhaseBOTH( jac->getA31(ph0), jac->getA31(ph1), presol ); 
		rhs->getV3()(ph0) = 0.0;
		rhs->getV3()(ph1) = 0.0;
	}
	else
	{
		if( ph1 != (-1) ){ colloc.PhaseONE( jac->getA31(ph1), presol ); rhs->getV3()(ph1) = 0.0; }
		if( ph0 != (-1) ) throw(-1);
	}
	
	// second the tangent
	if( cont )
	{
		const int nvar = var.Size()-1;
		colloc.Star( jac->getA31(nvar), xxDot->getV1() );
		jac->getA33().Clear();
		for( int i = 0; i < nvar; i++ ) jac->getA33(nvar, i) = xxDot->getV3()(i);
		jac->getA33(nvar,nvar) = p1Dot;
	}
}

void PointTR::JacobianVarying( Array1D<Vector*> A13, Array1D<int>& JacVar, 
                          Vector& sol, Vector& presol, Vector& par, Vector& prepar, bool cont, double ds )
{

	const int nvar = var.Size() - 1;

	// jacobian, derivatives of the right-hand side
	colloc.Jacobian( jac->getA11(), A13, rhs->getV1(), par, sol, JacVar );

	// the other equations
	// Currently: none
	for( int i = 1; i<eqn.Size(); i++ )
	{
		switch( eqn(i) )
		{
			case EqnTORPhase0:
			case EqnTORPhase1:
				// nothing happens
				break;
			default:
				std::cout<<"Torpoint: No such equation."<<eqn(i)<<"\n";
				break;
		}
	}
	
	// rhs from the tangent
	if( cont )
	{
		rhs->getV3()(nvar) = ds - (jac->getA31(nvar)*sol - jac->getA31(nvar)*presol);
		rhs->getV3()(nvar) -= (par(p1) - prepar(p1)) * p1Dot;
	}
}

// it does not work from here

// int PointTR::Refine( )
// {
// 	Array1D<int>       var( nvar );
// 	Array1D< Vector* > A13( nvar );
// 
// 	for( int i=0; i<nvar; i++ ) 
// 	{
// 		var(i) = i;
// 		A13(i) = & jac.getA13(i);
// 	}
// 	if( nvar == 1 ) var(0) = 1;
// 	int it = 0;
// 	xxDot.getV1() = xx.getV1();
// 	do
// 	{
// 		colloc.Jacobian( jac.getA11(), A13, rhs.getV1(), par, omega, xx.getV1(), var );
// 		if( nvar == 1 ) colloc.PhaseONE( jac.getA31(0), xxDot.getV1() );
// 		else colloc.PhaseBOTH( jac.getA31(0), jac.getA31(1), xxDot.getV1() );
// 		for( int r=0; r<nvar; r++ ) rhs.getV3()(r) = 0.0;//-(jac.getA31(r)*xx.getV1());
// 		jac.Solve( Dxx, rhs, nvar );
// 		xx.getV1() += Dxx.getV1();
// 		for( int r=0; r<nvar; r++ ) omega(var(r)) += Dxx.getV3()(r);
// 		std::cout<<"DNORM "<<sqrt(Dxx.getV1()*Dxx.getV1()+Dxx.getV3()*Dxx.getV3())<<" NORM "<<sqrt(xx.getV1()*xx.getV1())<<" omega "<<omega(0)<<","<<omega(1)<<"\n";
// 	}
// 	while( (sqrt(Dxx.getV1()*Dxx.getV1()+Dxx.getV3()*Dxx.getV3()) > 1e-5)&&(it++<5) );
// 	return it;
// }
// 
// void PointTR::Tangent( )
// {
// 	Array1D<int>       var( nvar + 1 );
// 	Array1D< Vector* > A13( nvar + 1 );
// 
// 	for( int i=0; i<nvar; i++ ) 
// 	{
// 		var(i) = i;
// 		A13(i) = & jac.getA13(i);
// 	}
// 	if( nvar == 1 ) var(0) = 1;
// 
// 	var( nvar ) = p1;
// 	A13( nvar ) = & xxDot.getV1();
// 
// 	rhs.getV3().Clear();
// 	colloc.Jacobian( jac.getA11(), A13, rhs.getV1(), par, omega, xx.getV1(), var );
// 	if( nvar == 1 ) colloc.PhaseONE( jac.getA31(0), xx.getV1() );
// 	else colloc.PhaseBOTH( jac.getA31(0), jac.getA31(1), xx.getV1() );
// 	std::cout<<"SOLVing\n";
// 	jac.Solve( rhs, xxDot, nvar );
// 	std::cout<<"NORMing\n";
// 	double norm = sqrt(1.0 + colloc.Integrate( rhs.getV1(), rhs.getV1() ));
// 	std::cout<<"NORMing2\n";
// 	p1Dot = 1.0/norm;
// 	xxDot.getV1() = rhs.getV1();
// 	xxDot.getV1() /= -norm;
// 	for( int i=0; i<nvar; i++ ) xxDot.getV3()(i) = -rhs.getV3()(i)/norm;
// }
// 
// int PointTR::Continue( double ds, bool first )
// {
// 	Array1D<int>       var( nvar + 1 );
// 	Array1D< Vector* > A13( nvar + 1 );
// 
// 	for( int i=0; i<nvar+1; i++ ) 
// 	{
// 		var(i) = i;
// 		A13(i) = & jac.getA13(i);
// 	}
// 	if( nvar == 1 ) var(0) = 1;
// 
// 	var( nvar ) = p1;
// 
// 	std::cout<<"initializing\n";
// 	// initializing xxNu, parNu
// 	xxNu.getV1() = xxDot.getV1();
// 	xxNu.getV3() = xxDot.getV3();
// 	xxNu.getV1() *= ds;
// 	xxNu.getV3() *= ds;
// 	xxNu.getV1() += xx.getV1();
// 	xxNu.getV3() += xx.getV3();
// 	parNu = par;
// 	omegaNu = omega;
// 	for( int i=0; i<nvar; i++ ) omegaNu(var(i)) += ds*xx.getV3()(i);
// 	parNu(p1) += ds*p1Dot;
// 
// 	// copying the phase conditions
// 		if( nvar == 1 ) colloc.PhaseONE( jac.getA31(0), xx.getV1() );
// 		else colloc.PhaseBOTH( jac.getA31(0), jac.getA31(1), xx.getV1() );
// 
// 	// copying the tangent
// 	colloc.Star( jac.getA31(nvar), xxDot.getV1() );
// 	for( int i=0; i<nvar; i++ ) jac.getA33(nvar, i) = xxDot.getV3()(i);
// 	jac.getA33().Clear();
// 	jac.getA33(nvar,nvar) = p1Dot;
// 	
// 	std::cout<<"entering cycle\n";
// 	// approximating the
// 	int it = 0;
// 	do
// 	{
// 		colloc.Jacobian( jac.getA11(), A13, rhs.getV1(), parNu, omegaNu, xxNu.getV1(), var );
// 		// computing the end of rhs
// 		rhs.getV3().Clear();
// 		// right hnad side of the phase conditions
// 		// for( int i=0; i<nvar; i++ ) rhs.getV3()(i) = - ( jac.getA31(i) * xxNu.getV1() );
// 		//rhs.getV3()(nvar) = ds - colloc.IntegrateDIFF( xxNu.getV1(), xx.getV1(), xxDot.getV1() );
// 		rhs.getV3()(nvar) = ds - (jac.getA31(nvar)*xxNu.getV1() - jac.getA31(nvar)*xx.getV1());
// 		for( int i=0; i<nvar; i++ ) rhs.getV3()(nvar) -= (omegaNu(var(i)) - omega(var(i))) * xxDot.getV3()(i);
// 		rhs.getV3()(nvar) -= (parNu(p1) - par(p1)) * p1Dot;
// 		// solving
// 		jac.Solve( Dxx, rhs );
// 		xxNu.getV1() += Dxx.getV1();
// 		for( int r=0; r<nvar; r++ ) omegaNu(var(r)) += Dxx.getV3()(r);
// 		parNu(p1) += Dxx.getV3()(nvar);
// 		std::cout<<"DNORM "<<sqrt(Dxx.getV1()*Dxx.getV1())<<" DNRM3 "<<sqrt(Dxx.getV3()*Dxx.getV3())<<" par(p1) "<<parNu(p1)
// 					<<" NORM "<<sqrt(xxNu.getV1()*xxNu.getV1())<<" omega "<<omegaNu(0)<<","<<omegaNu(1)<<","<<omegaNu(1)/omegaNu(0)<<"\n";
// 	}
// 	while( (sqrt(Dxx.getV1()*Dxx.getV1()/*+Dxx.getV3()*Dxx.getV3()*/) > 1e-4)&&(it++<5) );
// 	// check that the tangent is close to the result
// 	Vector DXX(xxNu.getV1());
// 	DXX -= xx.getV1();
// 	DXX /= ds;
// 	Vector DX3(xxNu.getV3());
// 	DX3 -= xx.getV3();
// 	DX3 /= ds;
// 	double NRM1 = sqrt(colloc.Integrate( DXX, DXX ));
// 	double NRM3 = sqrt(DX3*DX3);
// 	std::cout<<"D-DXXNRM "<<NRM1<<" DX3NRM "<<NRM3<<"\n";
// 	NRM1 = colloc.Integrate( xxDot.getV1(), xxDot.getV1() );
// 	NRM3 = sqrt( xxDot.getV3() * xxDot.getV3() );
// 	std::cout<<"T-DXXNRM "<<NRM1<<" DX3NRM "<<NRM3<<"\n";
// 	// update tangent
// 	rhs.getV1().Clear();
// 	rhs.getV3().Clear();
// 	rhs.getV3()(nvar) = 1.0;
// 	jac.Solve( xxDot, rhs );
// 	double norm = sqrt( colloc.Integrate( xxDot.getV1(), xxDot.getV1() ) + xxDot.getV3() * xxDot.getV3() );
// 	std::cout<<"nrm: "<<norm<<"\n";
// 	xxDot.getV1() /= norm;
// 	xxDot.getV3() /= norm;
// 	p1Dot = xxDot.getV3()(nvar) / norm;
// 	// update esolution
// 	xx.getV1() = xxNu.getV1();
// 	xx.getV3() = xxNu.getV3();
// 	par = parNu;
// 	omega = omegaNu;
// 	return it;
// }

int PointTR::Continue( double ds, bool /*first*/ )
{
// 	std::cout<<"initializing\n";
	// initializing xxNu, parNu

	const int        nvar = var.Size() - 1;
	Array1D<Vector*> A13( nvar+1 );
	Array1D<int>     JacVar( nvar+1 );

	ToJacVar( A13, JacVar, true );

// 	std::cout<<"initializing: approx nextsol\n";
	xxNu->getV1() = xxDot->getV1();
	xxNu->getV3() = xxDot->getV3();
	xxNu->getV1() *= ds;
	xxNu->getV3() *= ds;
	xxNu->getV1() += xx->getV1();
	xxNu->getV3() += xx->getV3();
	parNu = par; //contains the rotation number as well

// 	std::cout<<"initializing: approx omegaNu\n";
	for( int i = 0; i < nvar; i++ )
	{
// 		std::cout<<" JacVar(i) "<<JacVar(i); 
		parNu(JacVar(i)) += ds*xx->getV3()(i);
	}
	parNu(p1) += ds*p1Dot;

// 	std::cout<<"initializing: jacfixed\n";
	JacobianFixed( xxNu->getV1(), xx->getV1(), parNu, par, true );
	
// 	std::cout<<"entering cycle\n";
	// approximating the
	int it = 0;
	do
	{
		JacobianVarying( A13, JacVar, xxNu->getV1(), xx->getV1(), parNu, par, true, ds );
		jac->Solve( *Dxx, *rhs );
		xxNu->getV1() += Dxx->getV1();
		for( int r = 0; r < nvar; r++ ) parNu(JacVar(r)) += Dxx->getV3()(r);
		parNu(p1) += Dxx->getV3()(nvar);
		std::cout<<"DNORM "<<sqrt(Dxx->getV1()*Dxx->getV1())<<" DNRM3 "<<sqrt(Dxx->getV3()*Dxx->getV3())<<" par(p1) "<<parNu(p1)
					<<" NORM "<<sqrt(xxNu->getV1()*xxNu->getV1())<<" rho "<<parNu(RHO)<<"\n";
	}
	while( (sqrt(Dxx->getV1()*Dxx->getV1()+Dxx->getV3()*Dxx->getV3()) > ContEps)&&(it++<Iter) );
	// check that the tangent is close to the result
// 	Vector DXX(xxNu->getV1());
// 	DXX -= xx->getV1();
// 	DXX /= ds;
// 	Vector DX3(xxNu->getV3());
// 	DX3 -= xx->getV3();
// 	DX3 /= ds;
// 	double NRM1 = sqrt(colloc.Integrate( DXX, DXX ));
// 	double NRM3 = sqrt(DX3*DX3);
// 	std::cout<<"D-DXXNRM "<<NRM1<<" DX3NRM "<<NRM3<<"\n";
// 	NRM1 = colloc.Integrate( xxDot->getV1(), xxDot->getV1() );
// 	NRM3 = sqrt( xxDot->getV3() * xxDot->getV3() );
// 	std::cout<<"T-DXXNRM "<<NRM1<<" DX3NRM "<<NRM3<<"\n";
	// update tangent
	rhs->getV1().Clear();
	rhs->getV3().Clear();
	rhs->getV3()(nvar) = 1.0;
	jac->Solve( *xxDot, *rhs );
	double norm = sqrt( colloc.Integrate( xxDot->getV1(), xxDot->getV1() ) + xxDot->getV3() * xxDot->getV3() );
	std::cout<<"nrm: "<<norm<<"\n";
	xxDot->getV1() /= norm;
	xxDot->getV3() /= norm;
	p1Dot = xxDot->getV3()(nvar) / norm;
	// update esolution
	xx->getV1() = xxNu->getV1();
	xx->getV3() = xxNu->getV3();
	par = parNu;
	return it;
}
