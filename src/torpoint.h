// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef TORPOINT_H
#define TORPOINT_H

class System;

#include "matrix.h"
#include "hypermatrix.h"
#include "torcolloc.h"
#include "pointtype.h"

class PointTR
{
	public:
		PointTR( System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int ndeg1_, int ndeg2_, int nint1_, int nint2_ );
		~PointTR();
		void Tangent();
		int  Refine();
		int  Continue( double ds, bool first );

		// supplementary
		void    setCont( int i ) { p1 = i; }
		void    setSol( Vector& s ) { xx->getV1() = s; }
		Vector& getSol() { return xx->getV1(); }
		Vector& getTan() { return xxDot->getV1(); }
		void    setPar( Vector& p ) { for( int i=0; (i<p.Size())&&(i<par.Size()); i++ ) par(i) = p(i); }
		void    setOmega( double om0, double om1 ) 
		{
			par( colloc.Sys().npar() + ParRot ) = om1/om0/2.0;
		}
		Vector& getPar() { return par; }
		
		void    setIter( int i ) { Iter = i; }
		void    setRefEps( double d ) { RefEps = d; }
		void    setContEps( double d ) { ContEps = d; }
		void    setStartEps( double d ) { StartEps = d; }

		int     idxmap( int j1, int j2, int i1, int i2 ) { return colloc.idxmap( j1, j2, i1, i2 ); }

		void    ImportSol( Vector& Sol ) { colloc.ImportSol( xx->getV1(), Sol ); }
		void    ImportTan( Vector& Re, Vector& Im, double alpha ) 
		{
			colloc.ImportTan( xxDot->getV1(), Re, Im, alpha ); 
			xxDot->getV3().Clear(); 
			p1Dot = 0.0;
		}
		void    Start( double ds )
		{
			for( int i = 0; i < xx->getV1().Size(); i++ ) xx->getV1()(i) += ds*( xxDot->getV1()(i) );
		}
		double  Norm( ) { return sqrt( colloc.Integrate( xx->getV1(), xx->getV1() ) ); }
		void    SaveSol( char* dat, char* idx ) { colloc.Save( dat, idx, xx->getV1() ); }
		void    SaveTan( char* dat, char* idx ) { colloc.Save( dat, idx, xxDot->getV1() ); }
		void    SaveAll( char* dat, char* idx ) 
		{
			Vector tmp(xxDot->getV1()); 
			tmp *= 0.1;
			tmp += xx->getV1();
			colloc.Save( dat, idx, tmp ); 
		}
		int     npar() { return colloc.Sys().npar(); }

	private:
		void Construct();
		void Destruct();

		void ToJacVar( Array1D<Vector*>& A13, Array1D<int>& JacVar, bool cont );
		void JacobianFixed( Vector& sol, Vector& presol, Vector& par, Vector& prepar, bool cont );
		void JacobianVarying( Array1D<Vector*> A13, Array1D<int>& JacVar, 
                            Vector& sol, Vector& presol, Vector& par, Vector& prepar, bool cont, double ds );
		
		Array1D<Var> var;
		Array1D<Eqn> eqn;

		// convergence parameters
		double       RefEps;
		double       ContEps;
		double       StartEps;
		int          Iter;

		CollocTR     colloc;
		int          p1;
		Vector       par;
		Vector       parNu;
		HyperMatrix* jac;
		HyperVector* xx;
		HyperVector* xxNu;
		HyperVector* Dxx;
		HyperVector* rhs;
		HyperVector* xxDot;
		double       p1Dot; // actually this is equal to xxDot.getV3()(end);
};

#endif