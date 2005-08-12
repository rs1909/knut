// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef NPOINT_H
#define NPOINT_H

#include "matrix.h"
#include "hypermatrix.h"
#include "charmat.h"

#include "ncolloc.h"

#include "plot.h"

#include <fstream>
#include <list>
#include <cmath>
 
#define MULTIPLIERS 8
#define BMATRICES 1

#include "pointtype.h"

class PointData
{
	public:
		PointData( int nsol, int npar, int nev ) : sol(nsol), par(npar), mRe(nev), mIm(nev) {}

		void Save( ofstream& file )
		{
			for( int i=0; i < sol.Size(); i++ ) file<<sol(i)<<"\n";
		}
		void Load( ifstream& file );

	protected:

		// solution
		Vector sol;
		Vector par;

		// multipliers
		Vector mRe;
		Vector mIm;
};

class Point : public PointData
{
	public:

		// constructor
		Point( System& sys, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg, int nmul = MULTIPLIERS, int nmat = BMATRICES );
		~Point();

		// algorithms
		void    Reset( Array1D<Eqn>& eqn_, Array1D<Var>& var_ );
		void    SwitchLP( ); // switches branch
		void    SwitchPD( ); // switches branch
		void    SwitchHOPF( double ds );
		void    SwitchTRSol( Vector& Sol, const Vector& mshint, const Vector& mshdeg ); // starting data for tori: solution
		void    SwitchTRTan( Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg ); // starting data for tori: tangent
		int     Start( ); // approximating qq
		int     Refine( );
		void    Tangent( );
		int     Continue( double ds );
		void    Stability( );

		// supplementary
		void    setSol( Vector& s ) { sol = s; }
		Vector& getSol() { return sol; }
		
		void    setPar( Vector& p ) { for( int i=0; (i<p.Size())&&(i<par.Size()); i++ ) par(i) = p(i); }
		Vector& getPar() { return par; }
		
		void    setCont( int i );
		
		void    setIter( int i ) { Iter = i; }
		void    setRefEps( double d ) { RefEps = d; }
		void    setContEps( double d ) { ContEps = d; }
		void    setStartEps( double d ) { StartEps = d; }
		
		void    getqq( Vector& qq_ ) { if( qq ) qq_ = *qq; }
		double  Norm() 
		{
			// double e=0;
			// for( int j=0; j<colloc.Ndim(); j++ ) e += sol(j)*sol(j); // does not matter that headpoint or tailpoint, its periodic...
			return sqrt( colloc.Integrate( sol, sol ) );
		}
		double  NormMX() 
		{ 
			double max=0.0, min=1.0e32;
			for( int i=0; i<colloc.Nint()*colloc.Ndeg()+1; i++ )
			{
				double e=0.0;
				for( int j=0; j<colloc.Ndim(); j++ ) e += sqrt(sol(colloc.Ndim()*i + j)*sol(colloc.Ndim()*i + j));
				if( e > max ) max = e;
				if( e < min ) min = e;
			}
			return max-min;
		}
		int     UStab();
		int     UStabAUT();
		int     UStabAUTRot();
		PtType  testBif();
		PtType  testBifAUT();
		PtType  testBifAUTRot();
		
		void    Plot( GnuPlot& pl );
		void    Print( char* file )
		{
			ofstream ff(file);
			for( int i=0; i<colloc.Nint()*colloc.Ndeg()+1; i++ )
			{
				for( int j=0; j<colloc.Ndim(); j++ ) ff<<sol(colloc.Ndim()*i + j)<<"\t";
				ff<<"\n";
			}
		}
		void ReadNull( ifstream& file );
		void Read( ifstream& file, bool tan=false );
		void Write( ofstream& file );

	private:
		
		void    Construct();
		void    Dispose();
		void    FillSol( System& sys_ );
	
		// internal member-functions
		double inline SolNorm( Vector& sol, Vector* qq, Vector& par );
		
		void /*inline*/   Jacobian( HyperMatrix& AA, HyperVector& RHS, Vector& par, 
		                        Vector& solPrev, Vector& sol, JagMatrix3D& solData, 
		                        Vector* q0, Vector* q );
		void Point::Jacobian
		(
			HyperMatrix& AA, HyperVector& RHS,                      // output
			Vector& parPrev, Vector& par,                           // parameters
			Vector& solPrev, Vector& sol, JagMatrix3D& solData,     // solution
			Vector* qPrev,   Vector* q,                             // eigenvector
			Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
			double ds,       bool cont                              // the arclength; cont: true if continuation
		);
		
		void inline   Update( HyperVector& X );                          // sol,   qq,   par            += X
		void inline   ContUpdate( HyperVector& X );                      // solNu, qqNu, parNu, par(p1) += X

		// convergence parameters
		double       RefEps;
		double       ContEps;
		double       StartEps;
		int          Iter;
		
		// variables and equations
		Array1D<Var> var;
		Array1D<Eqn> eqn;
		Array1D<int> varMap;
		Array1D<int> varMapCont;
		
		int          dim1;
		int          dim2;
		int          dim3;
		
		// parameter to continue
		int          p1;
		
		// solutions
//		Vector       sol;
//		Vector       par;
		Vector*      qq0; // the eigenvector for the norm
		Vector*      qqR; // remainder eigenvector
		Vector*      qq;

		// solution in the continuation
		Vector       solNu;
		Vector       parNu;
		Vector*      qqNu;
		
		// tangents
		double       p1Dot;
		HyperVector* xxDot;  // -> solDot, qqDot, parDot
		
		// multipliers
//		Vector       mRe;
//		Vector       mIm;
		
		// for the rotation phase conditions
		Array1D<int> rotRe;
		Array1D<int> rotIm;
		
		// int the Newton iterations
		HyperVector* xx;     // -> sol,   --- this will be the solution in the iterations i.e., Dxx
		HyperVector* rhs;
//		HyperVector* xxRM; // it retires qqR.
		
		// store the Jacobian
		HyperMatrix* jac;
		
		// for collocation
		NColloc      colloc;
		
		// the solutions interpolated values
		JagMatrix3D  solData;
		
		// characteristic matrix
		baseCharMat* charMat;
// 		SpFact       jacStabA;
// 		SpMatrix     jacStabB;
		StabMatrix   jacStab;
};


#endif
