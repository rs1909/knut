// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef TESTFUNCT_H
#define TESTFUNCT_H

#include <hypermatrix.h>

// forward declarations
class NColloc;

class baseTestFunct
{
	public:
		virtual        ~baseTestFunct() {}
		virtual void   Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData ) = 0;
		virtual double Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData ) = 0;
		virtual double Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha ) = 0;
		virtual void   Funct_x( Vector& func, NColloc& col, const Vector& sol, const Vector& par, const JagMatrix3D& solData ) = 0;
};

class TestFunct : public baseTestFunct
{
	public:
		TestFunct( NColloc& col, double Z );
		~TestFunct( );
		void   Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
	
	private:
		bool        first;
		double      ZZ;
		HyperMatrix AHAT;
		Vector      A_p;
		SpMatrix    A_x;
		Vector      rhs;
		Vector      uu;
		Vector      vv;
		JagMatrix3D vvData;
};

class TestFunctLPAUT : public baseTestFunct
{
	public:
		TestFunctLPAUT( NColloc& col, double Z );
		~TestFunctLPAUT( );
		void   Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& sol, const Vector& par, const JagMatrix3D& solData );
	
	private:
		bool        first;
		double      ZZ;
		HyperMatrix AHAT;
		Vector      A_p;
		SpMatrix    A_x;
		SpMatrix    mB;
		Vector      mB_p;
		SpMatrix    mB_x;
		Vector      rhs;   // this is zero all the time
		Vector      phi;
		Vector      temp;
		Vector      DpPhi;
		Vector      DxPhi;
		Vector      uu2;   // for computing the generalized eigenvector
		Vector      vv2;
		Vector      gg2;
		Vector      hh2;
		Vector      one2;  // this is ( 0.0, 1.0 )
		JagMatrix3D phiData;
		JagMatrix3D vv2Data;
		JagMatrix3D solMSHData;
};

class TestFunctLPAUTROT : public baseTestFunct
{
	public:
		TestFunctLPAUTROT( NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z );
		~TestFunctLPAUTROT( );
		void   Init( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& par, const Vector& sol, const JagMatrix3D& solData );
	
	private:
		bool         first;
		double       ZZ;
		Array1D<int> Re, Im;
		HyperMatrix  AHAT;
		Vector       A_p;
		SpMatrix     A_x;
		SpMatrix     mB;
		Vector       mB_p;
		SpMatrix     mB_x;
		Vector       phi;
		Vector       DpPhi;
		Vector       DxPhi;
		Vector       LAM;
		Vector       DxLAM;
		Vector       uu3;   // for computing the generalized eigenvector
		Vector       vv3;
		Vector       gg3;
		Vector       hh3;
		Vector       rhs3;
		Vector       one3;  // this is ( 0.0, 1.0 )
		Vector       temp;
		JagMatrix3D  phiLAMData;
		JagMatrix3D  vv3Data;
		JagMatrix3D  solMSHData;
};

#endif
