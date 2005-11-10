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
		virtual void   Init( NColloc& col, const Vector& par, const JagMatrix3D& solData ) = 0;
		virtual double Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData ) = 0;
		virtual double Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha ) = 0;
		virtual void   Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData ) = 0;
};

class TestFunct : public baseTestFunct
{
	public:
		TestFunct( NColloc& col, double Z );
		~TestFunct( );
		void   Init( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData );
	
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
		void   Init( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData );
	
	private:
		bool        first;
		double      ZZ;
		HyperMatrix AHAT;
		Vector      A_p;
		SpMatrix    A_x;
		SpMatrix    mB;
		Vector      rhs;   // this is zero all the time
		Vector      uu1;   // for computing the trivial kernel
		Vector      vv1;
		Vector      uu2;   // for computing the generalized eigenvector
		Vector      vv2;
		Vector      gg2;
		Vector      hh2;
		Vector      one2;  // this is ( 0.0, 1.0 )
		JagMatrix3D vv2Data;
};


class TestFunctLPAUTROT : public baseTestFunct
{
	public:
		TestFunctLPAUTROT( NColloc& col, double Z );
		~TestFunctLPAUTROT( );
		void   Init( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData );
	
	private:
		bool        first;
		double      ZZ;
		HyperMatrix AHAT;
		Vector      A_p;
		SpMatrix    A_x;
		SpMatrix    mB;
		Vector      rhs;   // this is zero all the time
		Vector      uu1;   // for computing the trivial kernel
		Vector      vv1;
		Vector      uu2;   // for computing the second trivial kernel
		Vector      vv2;
		Vector      gg2;
		Vector      hh2;
		Vector      one2;
		Vector      uu3;   // for computing the generalized eigenvector
		Vector      vv3;
		Vector      gg3;
		Vector      hh3;
		Vector      one3;  // this is ( 0.0, 1.0 )
		JagMatrix3D vv3Data;
};

class TestFunctLPAUTROT_S : public baseTestFunct
{
	public:
		TestFunctLPAUTROT_S( NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z );
		~TestFunctLPAUTROT_S( );
		void   Init( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct( NColloc& col, const Vector& par, const JagMatrix3D& solData );
		double Funct_p( NColloc& col, const Vector& par, const JagMatrix3D& solData, int alpha );
		void   Funct_x( Vector& func, NColloc& col, const Vector& par, const JagMatrix3D& solData );
	
	private:
		bool        first;
		double      ZZ;
		HyperMatrix AHAT;
		Vector      A_p;
		SpMatrix    A_x;
		SpMatrix    mB;
		Vector      rhs;   // this is zero all the time
		Vector      uu1;   // for computing the trivial kernel
		Vector      vv1;
		Vector      uu2;   // for computing the second trivial kernel
		Vector      vv2;
		Vector      gg2;
		Vector      hh2;
		Vector      one2;
		Vector      uu3;   // for computing the generalized eigenvector
		Vector      vv3;
		Vector      gg3;
		Vector      hh3;
		Vector      one3;  // this is ( 0.0, 1.0 )
		Array1D<int> Re, Im;
		JagMatrix3D vv3Data;
};

#endif
