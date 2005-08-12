// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef NCOLLOC_H
#define NCOLLOC_H

#include "system.h"
#include "matrix.h"
#include "spmatrix.h"

class NColloc
{
	public:

		NColloc( System& _sys, const int nint, const int ndeg, int nmat );      // computes mesh, metric, metricD
//		NColloc( System& _sys, const Vector& Tmesh, const int ndeg ); // computes metric, metricD

		~NColloc() {}

		void Init( const Vector& par, const Vector& sol ); // computes time, kk, ee, dd, rr ...

		void Interpolate( JagMatrix3D& out, const Vector& sol );
		void InterpolateREAL( JagMatrix3D& out, const Vector& sol );
		void InterpolateCPLX( JagMatrix3D& out, const Vector& sol );
			// out(NDIM,NTAU+1,NDEG*NINT)
			// out( ., NTAU, . ) contains the derivative

		void   Star( Vector& out, Vector& sol );
		double Integrate( Vector& v1, Vector& v2 );
		double IntegrateDerivative( Vector& v1, Vector& v2 );
		double IntegrateCont( Vector& v1, Vector& v2, Vector& v3 );
		
		void   PhaseStar( Vector& V1, Vector& V2 );
		void   PhaseRotStar( Vector& V1, Vector& V2, Array1D<int>& Re, Array1D<int>& Im );
		
		void   Import( Vector& out, const Vector& in, const Vector& mesh, int deg_ );
		void   Export( Vector& out, const Vector& mshint, const Vector& mshdeg, const Vector& in );
		
		// computing the Jacobians, right-hand sides, characteristic matrices		

		// continuation of solution
	
		void RHS( Vector& rhs, const Vector& par, const Vector& sol, const JagMatrix3D& solData );          // 
		void RHS_p( Vector& rhs, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int p ); // sol is currently not needed
		void RHS_x( SpMatrix& A, const Vector& par, const Vector& sol, const JagMatrix3D& solData );        // sol is currently not needed

		// for stability computation
		
// 		void StabJac( SpMatrix& A, SpMatrix& B, const Vector& par, const JagMatrix3D& solData );
		void StabJac( StabMatrix& AB, const Vector& par, const JagMatrix3D& solData );
		// this computes its own matrix structures, because it is nowhere else needed: kkSI, eeSI, rrSI, ddSI, etc.
		// However, the variables will be contained within the NColloc class

		// continuation of bifurcations -> characteristic matrices

		void CharJac_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, double Z );
		void CharJac_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, double ZRe, double ZIm );
		
		void CharJac_x_p( Vector& V, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z, int p );
		void CharJac_x_p( Vector& V, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double ZRe, double ZIm, int p );
		
		void CharJac_x_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z );
		void CharJac_x_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Re, double Im );
		
		void CharJac_x_z( Vector& V, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Re, double Im );
		// we do need Re, Im
		
		// for the autonomous FOLD bifurcation
		void CharJac_mB(   SpMatrix& A, const Vector& par, const JagMatrix3D& solData, double Z );
		void CharJac_mB_p( Vector& V,   const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z, int alpha );
		void CharJac_mB_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData, const JagMatrix3D& phiData, double Z );
		
		// autonoumous trivial eigenvector
		void CharJac_phi( Vector& V, const Vector& par, const JagMatrix3D& solData );
		void CharJac_phi_x( SpMatrix& A, const Vector& par, const JagMatrix3D& solData );
		void CharJac_phi_p( Vector& V, const Vector& par, const JagMatrix3D& solData, int alpha );

		// supplementary
		int Ndim() const { return ndim; }
		int Npar() const { return npar; }
		int Ntau() const { return ntau; }
		int Nint() const { return nint; }
		int Ndeg() const { return ndeg; }
		
		void setMesh( const Vector& msh );
		void getMesh( Vector& msh );
		double Profile( int i, int d ) { return mesh(i) + meshINT(d)*(mesh(i+1)-mesh(i)); }

	private:

		// helper functions

		int& WRIDX( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int cTAU, int cDEG, int cDIM )
		{ 
			return A.WrLi(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM, cDIM + ndim*( cDEG - dd(cTAU,idx)+rr(cTAU,idx)*(ndeg+1)) );
		}

		double& WRDAT( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int cTAU, int cDEG, int cDIM )
		{ 
			return A.WrLx(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM, cDIM + ndim*( cDEG - dd(cTAU,idx)+rr(cTAU,idx)*(ndeg+1)) );
		}

		// for CharJac_x
		int& WRIDXCPLX( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int rIM, int cTAU, int cDEG, int cDIM, int cIM )
		{ 
			return A.WrLi( rIM + 2*(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM), cIM + 2*(cDIM + ndim*( cDEG - dd(cTAU,idx)+rr(cTAU,idx)*(ndeg+1))) );
		}

		double& WRDATCPLX( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int rIM, int cTAU, int cDEG, int cDIM, int cIM )
		{ 
			return A.WrLx( rIM + 2*(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM), cIM + 2*(cDIM + ndim*( cDEG - dd(cTAU,idx)+rr(cTAU,idx)*(ndeg+1))) );
		}
		
		// for CharJac_x_x
		int& WRIDXCPLXM( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int rIM, int cTAU, int cDEG, int cDIM )
		{ 
			return A.WrLi( rIM + 2*(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM), cDIM + ndim*( cDEG - dd(cTAU,idx)+rr(cTAU,idx)*(ndeg+1)) );
		}

		double& WRDATCPLXM( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int rIM, int cTAU, int cDEG, int cDIM )
		{ 
			return A.WrLx( rIM + 2*(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM), cDIM + ndim*( cDEG - dd(cTAU,idx)+rr(cTAU,idx)*(ndeg+1)) );
		}
		
		int& WRIDXS( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int cTAU, int cDEG, int cDIM )
		{
			return A.WrLi(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM, cDIM + ndim*( cDEG - ddS(cTAU,idx)+rrS(cTAU,idx)*(ndeg+1)) );
		}

		double& WRDATS( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int cTAU, int cDEG, int cDIM )
		{
			return A.WrLx(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM, cDIM + ndim*( cDEG - ddS(cTAU,idx)+rrS(cTAU,idx)*(ndeg+1)) );
		}

		int& WRIDXI( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int cTAU, int cDEG, int cDIM )
		{
			return A.WrLi(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM, cDIM + ndim*( cDEG - ddI(cTAU,idx)+rrI(cTAU,idx)*(ndeg+1)) );
		}

		double& WRDATI( SpMatrix& A, int idx, int rINT, int rDEG, int rDIM, int cTAU, int cDEG, int cDIM )
		{
			return A.WrLx(ndim + ndim*( rDEG + ndeg*rINT ) + rDIM, cDIM + ndim*( cDEG - ddI(cTAU,idx)+rrI(cTAU,idx)*(ndeg+1)) );
		}
		
		// the equations 
		System* sys;

		int ndim;
		int npar;
		int ntau;
		
		int nint;
		int ndeg;
		
		int nmat;

		Vector mesh;
		Vector time;

		// these store the structure of the sparse matrix NO SEPARATION
		Array2D<int> kk;   // dim(NTAU+1,NDEG*NINT) : which delay to which interval 
		Array2D<int> ee;   // dim(NTAU+1,NDEG*NINT) : the ordering of kk
		Array2D<int> dd;	 // dim(NTAU+1,NDEG*NINT) : how many neighbouring intervals up to its index
		Array2D<int> rr;	 // dim(NTAU+1,NDEG*NINT) : how many non overlapping intervals up to its index

		// these store the structure of the sparse matrix WITH SEPARATION
		Array2D<int> kkS;  // same as kk, but not folded back
		Array2D<int> eeS;  //  ...
		Array2D<int> rrS;  //  ...
		Array2D<int> ddS;  //  ...
		Array2D<int> mmS;  // which matrix are we in
		Array2D<int> szS;  // szI(mmI( . ,idx),idx) size of the line within the matrix mmS
		
		// For the stability matrices
		Array2D<int> kkI;
		Array2D<int> eeI;
		Array2D<int> rrI;
		Array2D<int> ddI;
		Array2D<int> mmI;
		Array2D<int> szI;
		
		// it stores all the collocation data
		Array3D<double> tt;

		// matrix for integration
		Matrix metric;
		// integration with derivatives (for phase conditions)
		Matrix metricPhase;

		// internal use for the initialization
		Vector col;
		Vector out;
		Vector meshINT;
		
		// for rhs, and derivatives
		// Matrix fx, dfx, t_dfx, dummy
};

#endif
