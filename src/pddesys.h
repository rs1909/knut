// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef PDDESYS_H
#define PDDESYS_H

#include "matrix.h"

extern "C"
{
	int    sys_ndim();
	int    sys_npar();
	int    sys_ntau();
	int    sys_nderi();
	void   sys_tau( Vector& out, double t, const Vector& par );
	void   sys_dtau( Vector& out, double t, const Vector& par, int vp );
	void   sys_rhs( Vector& out, double t, const Matrix& x, const Vector& par );
	void   sys_deri( Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v );
	void   sys_stpar( Vector& par );
	void   sys_stsol( Vector& out, double t );
};

#endif
