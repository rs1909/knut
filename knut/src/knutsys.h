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

#include <cstdlib>
#include "matrix.h"

#ifdef _MSC_VER
# define SYSEXPORT __declspec( dllexport )
#else
# define SYSEXPORT
#endif

extern "C"
{
  SYSEXPORT size_t sys_ndim();
  SYSEXPORT size_t sys_npar();
  SYSEXPORT size_t sys_ntau();
  SYSEXPORT size_t sys_nderi();
  SYSEXPORT size_t sys_nevent();
  SYSEXPORT void   sys_tau(KNVector& out, double t, const KNVector& par);
  SYSEXPORT void   sys_dtau(KNVector& out, double t, const KNVector& par, size_t vp);
  SYSEXPORT void   sys_mass(KNArray1D<double>& out);
  SYSEXPORT void   sys_rhs(KNVector& out, double t, const KNMatrix& x, const KNVector& par);
  SYSEXPORT void   sys_deri(KNMatrix& out, double t, const KNMatrix& x, const KNVector& par, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNMatrix& v);
  SYSEXPORT void   sys_stpar(KNVector& par);
  SYSEXPORT void   sys_stsol(KNVector& out, double t);
  // Vectorized version of each function
  SYSEXPORT void   sys_p_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par );
  SYSEXPORT void   sys_p_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par, size_t vp );
  SYSEXPORT void   sys_p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel );
  SYSEXPORT void   sys_p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
  SYSEXPORT void   sys_p_event( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par );
  SYSEXPORT void   sys_parnames( const char *names[] );
};

#endif
