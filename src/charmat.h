// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef CHARMAT_H
#define CHARMAT_H

#include "matrix.h"
#include "spmatrix.h"
#include "ncolloc.h"
#include "error.h"

class baseCharMat
{
	public:
		virtual ~baseCharMat() {}
		
		virtual void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z ) = 0;
		virtual void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector& q ) = 0;
		
		virtual void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Re, double Im ) = 0;
		virtual void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Re, double Im, Vector& q ) = 0;
		
		virtual void Delta( MatFact& delta, NColloc& col ) = 0;
		virtual void Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData ) = 0;
		virtual void Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p ) = 0;
		
		virtual void Delta_z( Vector& delta_z, NColloc& col, const Vector& par, const JagMatrix3D& solData ) = 0;
		virtual void Switch( Vector& sol, NColloc& col ) = 0;
};

class CharMat : public baseCharMat
{
	public:
		CharMat( NColloc& col );
		virtual ~CharMat() { }
		
		void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z );
		void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector& q );
		
		void Init( NColloc& , const Vector&, const JagMatrix3D&, double, double ) { }
		void Init( NColloc& , const Vector&, const JagMatrix3D&, double, double, Vector& ) { }
		
		void Delta( MatFact& delta, NColloc& col );
		void Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData );
		void Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p );
		
		void Delta_z( Vector&, NColloc&, const Vector&, const JagMatrix3D& ) { }
		
		void Switch( Vector& tan, NColloc& col );
		
	private:
	
		double      ZZ;
		Vector      phi;
		
		SpFact      AmzB;    // A - zB
		SpMatrix    AmzB_x;
		Vector      AmzB_p;
		
		JagMatrix3D phiData;
};

class CharMatCPLX : public baseCharMat
{
	public:
		CharMatCPLX( NColloc& col );
		virtual ~CharMatCPLX() { }
		
		void Init( NColloc& , const Vector& , const JagMatrix3D& , double ) { PDError(-1); }
		void Init( NColloc& , const Vector& , const JagMatrix3D& , double , Vector&  ) { PDError(-1); }
		
		void Init( NColloc& , const Vector&, const JagMatrix3D&, double zzRe, double zzIm );
		void Init( NColloc& , const Vector&, const JagMatrix3D&, double zzRe, double zzIm, Vector& q );
		
		void Delta( MatFact& delta, NColloc& col );
		void Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData );
		void Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p );
		
		void Delta_z( Vector& delta_z, NColloc& col, const Vector& par, const JagMatrix3D& solData );
		
		void Switch( Vector& /*tan*/, NColloc& /*col*/ ) { PDError(-1); }
		
		void Switch( Vector& Re, Vector& Im, double& alpha, NColloc& col );
	private:
	
		double      zRe, zIm;
		Vector      phi;
		
		SpFact      AmzB;    // A - zB
		SpMatrix    AmzB_x;
		Vector      AmzB_p;
		Vector      AmzB_z;
		
		JagMatrix3D phiData;
};

class CharMatLPAUT : public baseCharMat
{
	public:
		CharMatLPAUT( NColloc& col );
		virtual ~CharMatLPAUT() { }
		
		void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z ) 
			{ getQs( col, par, solData, 0 ); init( col, par, solData, Z, 0 ); }
		void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector& q )
			{ getQs( col, par, solData, &q ); init( col, par, solData, Z, &q ); }
		
		void Init( NColloc& , const Vector&, const JagMatrix3D&, double, double ) { }
		void Init( NColloc& , const Vector&, const JagMatrix3D&, double, double, Vector& ) { }
		
		void Delta( MatFact& delta, NColloc& col );
		void Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData );
		void Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p );
		
		void Delta_z( Vector&, NColloc&, const Vector&, const JagMatrix3D& ) { }
		
		void Switch( Vector& /*tan*/, NColloc& /*col*/ ) { }
		
	protected:
		void init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector* q );
		void getQs( NColloc& col, const Vector& par, const JagMatrix3D& solData, Vector* qq );
		
		inline void mktmp_x( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		inline void mktmp_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int p );
		
		double      ZZ;
		double      alpha;
		Vector      Q3;
		Vector      Phi1, Psi1, Phi3, Sigma;
		SpMatrix    DxAmzBPhi1;
		Vector      DpAmzBPhi1;
		
		Matrix      MAmzBinv;           // (NDIM+1)*(NDIM*(NDEG*NINT+1))
		
		SpFact      AmzB;
		SpMatrix    mB;
		
		SpMatrix    DxAmzBSigma;
		Vector      DpAmzBSigma;
		
		SpMatrix    DxmBPhi1;
		Vector      DpmBPhi1;
		
		JagMatrix3D Phi1Data;
		JagMatrix3D SigmaData;
		
		//temporary storage
		Matrix tempAx, tempBx, tempCx, tempDx;
		Vector tempAp, tempBp;
		Vector stmpAp, stmpBp, stmpCp, stmpDp;
};

// class CharMatLPAUTROT : public CharMatLPAUT
// {
// 	public:
// 		CharMatLPAUTROT( NColloc& col, Array1D<int>& Re_, Array1D<int>& Im_ );
// 		virtual ~CharMatLPAUTROT() { }
// 		
// 		void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z ) 
// 			{ getQs( col, par, solData, 0 ); init( col, par, solData, Z, 0 ); }
// 		void Init( NColloc& col, const Vector& par, const JagMatrix3D& solData, double Z, Vector& q )
// 			{ getQs( col, par, solData, &q ); init( col, par, solData, Z, &q ); }
// 		
// 		void Init( NColloc& , const Vector&, const JagMatrix3D&, double, double ) { }
// 		void Init( NColloc& , const Vector&, const JagMatrix3D&, double, double, Vector& ) { }
// 		
// 		void Delta( MatFact& delta, NColloc& col );
// 		void Delta_x( Matrix& delta_x, NColloc& col, const Vector& par, const JagMatrix3D& solData );
// 		void Delta_p( Vector& delta_p, NColloc& col, const Vector& par, const JagMatrix3D& solData, int p );
// 		
// 		void Delta_z( Vector&, NColloc&, const Vector&, const JagMatrix3D& ) { }
// 		void Switch( Vector& /*tan*/, NColloc& /*col*/ ) { }
// 	
// 	private:
// 		void getQs( NColloc& col, const Vector& par, const JagMatrix3D& solData, Vector* qq );
// 		
// 		Array1D<int> Re, Im;
// // 		double       alpha2;
// 		Vector       Q2;
// // 		Vector       DpQ2;
// // 		Matrix       DxQ2;
// };

#endif
