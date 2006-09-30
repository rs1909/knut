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

#ifdef _MSC_VER
# define SYSEXPORT __declspec( dllexport )
#else
# define SYSEXPORT
#endif

extern "C"
{
	SYSEXPORT int    sys_ndim();
	SYSEXPORT int    sys_npar();
	SYSEXPORT int    sys_ntau();
	SYSEXPORT int    sys_nderi();
	SYSEXPORT void   sys_tau( Vector& out, double t, const Vector& par );
	SYSEXPORT void   sys_dtau( Vector& out, double t, const Vector& par, int vp );
	SYSEXPORT void   sys_rhs( Vector& out, double t, const Matrix& x, const Vector& par );
	SYSEXPORT void   sys_deri( Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v );
	SYSEXPORT void   sys_stpar( Vector& par );
	SYSEXPORT void   sys_stsol( Vector& out, double t );
};

#endif
