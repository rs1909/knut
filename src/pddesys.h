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

namespace Sys
{
	int    ndim();
	int    npar();
	int    ntau();
	int    nderi();
	void   tau( Vector& out, double t, const Vector& par );
	void   dtau( Vector& out, double t, const Vector& par, int vp );
	void   rhs( Vector& out, double t, const Matrix& x, const Vector& par );
	void   deri( Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v );
	void   stpar( Vector& par );
	void   stsol( Vector& out, double t );
};

#endif
