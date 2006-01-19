// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef TORCOLLOC_H
#define TORCOLLOC_H

#include "matrix.h"

class System;
class SpMatrix;

class CollocTR
{
	public:
		CollocTR( System& sys_, int ndeg1_, int ndeg2_, int nint1_, int nint2_ );
		// this provides the jacobian, the right hand side and the derivatives w.r.t. var
		// the difficulty is with the derivative w.r.t. the frequencies
		void Jacobian( SpMatrix& A, Array1D< Vector* > Avar, Vector& rhs, Vector& par, Vector& sol, Array1D<int>& var ); 
		// not yet implemented
		// these are easy
		void PhaseONE( Vector& ph, Vector& presol );                                     // implemented
		void PhaseBOTH( Vector& ph0, Vector& ph1, Vector& presol );                      // implemented
		double Integrate( Vector& ph1, Vector& ph2 );
		double IntegrateDIFF( Vector& ph1, Vector& ph2, Vector& ph3 );
		void Star( Vector& ph1, Vector& ph2 );

		void ImportSol( Vector& out, Vector& in );
		void ImportTan( Vector& out, Vector& Re, Vector& Im, double alpha );
		void Save( const char* dat, const char* idx, const Vector& in );

		inline int Ndeg1( ) { return ndeg1; }
		inline int Ndeg2( ) { return ndeg2; }
		inline int Nint1( ) { return nint1; }
		inline int Nint2( ) { return nint2; }

		//utils
		inline int idxmap( int j1, int j2, int i1, int i2 );
		inline int idxmapSep( int j1, int j2, int i1, int i2 );
		inline int idxkk( int j1, int j2, int k );
		inline System& Sys() { return *sys; }

	private:
		System* sys;
		int ndeg1, ndeg2;
		int nint1, nint2;
		Vector col1,  col2;
		Vector mesh1, mesh2;
		Array1D< Array1D<double> > lgr1, lgr2; // 1. meshpoint 2. polynom
		Array1D< Array1D<double> > dlg1, dlg2; // 1. meshpoint 2. polynom
		
		Matrix I1, ID1, I2, ID2;
		Vector mlg1, mlg2, mlgd1, mlgd2, ilg1, ilg2, ilgd1, ilgd2;
		// functions
		void index( int* kk, int* ee, int* rr, Vector& par, Vector& tau, double* t1, double* t2 );
		void indexSep( int* kk, int* ee, int* rr, Vector& par, Vector& tau, double* t1, double* t2 );
		void interpolate( Matrix& xx, Array1D<double>& sol, int* kk, double* t1, double* t2 );
};

inline int CollocTR::idxmap( int j1, int j2, int i1, int i2 )
{
	if( j1 < ndeg1 ) 
	{
		if( j2 < ndeg2 )
		{
			// j1 < NDEG1
			// j2 < NDEG2
			return j1 + j2*ndeg1 + i1*ndeg1*ndeg2 + i2*ndeg1*ndeg2*nint1;
		}else
		{
			// j1 < NDEG1
			// j2 = NDEG2
			return j1 + 0*ndeg1 + i1*ndeg1*ndeg2 + ((i2+1)%nint2)*ndeg1*ndeg2*nint1;
		}
	}else
	{
		if( j2 < ndeg2 )
		{
			// j1 = NDEG1
			// j2 < NDEG2
			return 0 + j2*ndeg1 + ((i1+1)%nint1)*ndeg1*ndeg2 + i2*ndeg1*ndeg2*nint1;
		}else
		{
			// j1 = NDEG1
			// j2 = NDEG2
			if( (i1 == nint1-1)&&(i2 == nint2-1) ) return 0;
			else return 0 + 0*ndeg1 + ((i1+1)%nint1)*ndeg1*ndeg2 + ((i2+1)%nint2)*ndeg1*ndeg2*nint1;
		}
	}
}

inline int CollocTR::idxmapSep( int j1, int j2, int i1, int i2 )
{
	if( j1 < ndeg1 )
	{
		if( j2 < ndeg2 )
		{
			// j1 < NDEG1
			// j2 < NDEG2
			return j1 + j2*ndeg1 + i1*ndeg1*ndeg2 + i2*ndeg1*ndeg2*nint1;
		}else
		{
			// j1 < NDEG1
			// j2 = NDEG2
			return j1 + 0*ndeg1 + i1*ndeg1*ndeg2 + ((i2+1)%nint2)*ndeg1*ndeg2*nint1;
		}
	}else
	{
		if( j2 < ndeg2 )
		{
			// j1 = NDEG1
			// j2 < NDEG2
			return 0 + j2*ndeg1 + (i1+1)*ndeg1*ndeg2 + i2*ndeg1*ndeg2*nint1;
		}else
		{
			// j1 = NDEG1
			// j2 = NDEG2
			if( (i1 == nint1-1)&&(i2 == nint2-1) ) return 0;
			else return 0 + 0*ndeg1 + (i1+1)*ndeg1*ndeg2 + ((i2+1)%nint2)*ndeg1*ndeg2*nint1;
		}
	}
}

inline int CollocTR::idxkk( int j1, int j2, int k )
{
	return j1 + j2*(ndeg1+1) + k*(ndeg1+1)*(ndeg2+1);
}

#endif
