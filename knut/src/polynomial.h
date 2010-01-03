// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <iostream>
#include "matrix.h"

void poly_lgr(const KNVector& t, KNVector &out, double c);
void poly_dlg(const KNVector& t, KNVector& out, double c);
void poly_d2lg(const KNVector& t, KNVector& out, double c);

/// returns (pp(0) + pp(1)*x + ... + pp(n)*x^n) * (aa + bb*x)
///
inline void poly_linmul(KNArray1D<double>& pp, double aa, double bb)
{
  P_ASSERT_X1(pp(pp.size() - 1) == 0.0, "poly_linmul: truncating the highest order term!");
  for (int i = pp.size() - 1; i > 0; --i)
  {
    pp(i) = aa * pp(i) + bb * pp(i - 1);
  }
  pp(0) = aa * pp(0);
}

void poly_coeff_lgr(KNArray1D<double>& out, const KNArray1D<double>& t, int i);

void poly_coeff_mul(KNArray1D<double>& out, KNArray1D<double>& in1, KNArray1D<double>& in2);

void poly_coeff_int(KNArray1D<double>& out, KNArray1D<double>& in);

void poly_coeff_diff(KNArray1D<double>& out, KNArray1D<double>& in);

double poly_eval(const KNArray1D<double>& in, double c);

void lobatto(KNArray1D<double>& mesh);

void gauss(KNArray1D<double>& mesh);

double poly_lgr_eval( const KNArray1D<double>& t, int i, double c);

double poly_dlg_eval( const KNArray1D<double>& t, int i, double c);
