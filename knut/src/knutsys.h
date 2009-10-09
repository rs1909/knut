// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef KNUTSYS_H
#define KNUTSYS_H

#include "matrix.h"
#include <string>
#include <vector>

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
  SYSEXPORT int    sys_nevent();
  SYSEXPORT void   sys_tau(Vector& out, double t, const Vector& par);
  SYSEXPORT void   sys_dtau(Vector& out, double t, const Vector& par, int vp);
  SYSEXPORT void   sys_rhs(Vector& out, double t, const Matrix& x, const Vector& par);
  SYSEXPORT void   sys_deri(Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v);
  SYSEXPORT void   sys_stpar(Vector& par);
  SYSEXPORT void   sys_stsol(Vector& out, double t);
  // Vectorized version of each function
  SYSEXPORT void   sys_p_tau( Array2D<double>& out, const Array1D<double>& time, const Array1D<double>& par );
  SYSEXPORT void   sys_p_dtau( Array2D<double>& out, const Array1D<double>& time, const Array1D<double>& par, int vp );
  SYSEXPORT void   sys_p_rhs( Array2D<double>& out, const Array1D<double>& time, const Array3D<double>& x, const Array1D<double>& par, int sel );
  SYSEXPORT void   sys_p_deri( Array3D<double>& out, const Array1D<double>& time, const Array3D<double>& x, const Array1D<double>& par, int sel, int nx, const int* vx, int np, const int* vp, const Array3D<double>& vv );
  SYSEXPORT void   sys_p_event( Array2D<double>& out, const Array1D<double>& time, const Array3D<double>& x, const Array1D<double>& par );
  SYSEXPORT void   sys_parnames( std::vector<std::string>& out );
};

#endif
