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
#include "testfunct.h"

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
		inline PointData( int nsol, int npar, int nev ) : sol(nsol), par(npar), mRe(nev), mIm(nev) {}

		inline void Save( std::ofstream& file )
		{
			for( int i=0; i < sol.Size(); i++ ) file<<sol(i)<<"\n";
		}
		void Load( std::ifstream& file );

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
		void    SwitchTFLP( double ds ); // switches branch with testFunct
		void    SwitchTFPD( double ds ); // switches branch with testFunct
		void    SwitchTFHB( double ds ); // switches branch with testFunct
		void    SwitchTRSol( Vector& Sol, const Vector& mshint, const Vector& mshdeg )
			{ colloc.Export( Sol, mshint, mshdeg, sol ); } // starting data for tori: solution
		void    SwitchTFTRTan( Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg ); // starting data for tori: tangent WITH testFunct
		int     StartTF( Eqn FN );
		int     Refine( );
		void    Tangent( );
		int     Continue( double ds );
		void    Stability( );

		// supplementary
		inline void    setSol( Vector& s ) { sol = s; }
		inline Vector& getSol() { return sol; }
		
		inline void    setPar( Vector& p ) { for( int i=0; (i<p.Size())&&(i<par.Size()); i++ ) par(i) = p(i); }
		inline Vector& getPar() { return par; }
		
		inline void    setCont( int p ) { p1 = p; varMapCont( varMap.Size() ) = p1; }
		inline void    setSym( int n, int* sRe, int* sIm )
		{
			rotRe.Init(n); rotIm.Init(n);
			for( int i=0; i<n; i++ ) { rotRe(i) = sRe[i]; rotIm(i) = sIm[i]; }
		}
		
		inline void    setRefIter( int i ) { RefIter = i; }
		inline void    setContIter( int i ) { ContIter = i; }
		inline void    setStartIter( int i ) { StartIter = i; }
		inline void    setRefEps( double d ) { RefEps = d; }
		inline void    setContEps( double d ) { ContEps = d; }
		inline void    setStartEps( double d ) { StartEps = d; }
		
		inline double  Norm() 
		{
			// double e=0;
			// for( int j=0; j<colloc.Ndim(); j++ ) e += sol(j)*sol(j); // does not matter that headpoint or tailpoint, its periodic...
			return sqrt( colloc.Integrate( sol, sol ) );
		}
		inline double  NormMX() 
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
		inline void    Print( char* file )
		{
			std::ofstream ff(file);
			for( int i=0; i<colloc.Nint()*colloc.Ndeg()+1; i++ )
			{
				for( int j=0; j<colloc.Ndim(); j++ ) ff<<sol(colloc.Ndim()*i + j)<<"\t";
				ff<<"\n";
			}
		}
		void ReadNull( std::ifstream& file );
		void Read( std::ifstream& file, bool tan=false );
		void Write( std::ofstream& file );

	private:
		
		void    Construct();
		void    Dispose();
		void    FillSol( System& sys_ );
	
		// internal member-functions
		inline double SolNorm( Vector& sol, Vector& par );
		
// 		void    Jacobian( HyperMatrix& AA, HyperVector& RHS, Vector& par, 
// 		                        Vector& solPrev, Vector& sol, JagMatrix3D& solData, 
// 		                        Vector* q0, Vector* q );
		void Jacobian
		(
			HyperMatrix& AA, HyperVector& RHS,                      // output
			Vector& parPrev, Vector& par,                           // parameters
			Vector& solPrev, Vector& sol, JagMatrix3D& solData,     // solution
			Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
			bool cont                                               // cont: true if continuation
		);
		
		inline void   Update( HyperVector& X );                          // sol,   qq,   par            += X
		inline void   ContUpdate( HyperVector& X );                      // solNu, qqNu, parNu, par(p1) += X

		// convergence parameters
		double       RefEps;
		double       ContEps;
		double       StartEps;
		int          RefIter;
		int          ContIter;
		int          StartIter;
		
		// variables and equations
		Array1D<Var> var;
		Array1D<Eqn> eqn;
		Array1D<int> varMap;
		Array1D<int> varMapCont;
		
		int          dim1;
		int          dim3;
		
		// parameter to continue
		int          p1;
		
		// solutions
//		Vector       sol;
//		Vector       par;

		// solution in the continuation
		Vector       solNu;
		Vector       parNu;
		
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
		
		// store the Jacobian
		HyperMatrix* jac;
		
		// for collocation
		NColloc      colloc;
		
		// the solutions interpolated values
		JagMatrix3D  solData;
		
		// stability matrix
		StabMatrix   jacStab;
		
		// test functionals
		baseTestFunct* testFunct;
};

#endif
