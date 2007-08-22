#include <iostream>
#include "matrix.h"

/// returns (pp(0) + pp(1)*x + ... + pp(n)*x^n) * (aa + bb*x)
///
inline void poly_linmul(Array1D<double>& pp, double aa, double bb)
{
  if (pp(pp.Size() - 1) != 0.0) std::cout << "Warning: \"poly_linmul\" is truncating the highest order term!\n";
  for (int i = pp.Size() - 1; i > 0; --i)
  {
    pp(i) = aa * pp(i) + bb * pp(i - 1);
  }
  pp(0) = aa * pp(0);
}

void poly_coeff_lgr(Array1D<double>& out, const Array1D<double>& t, int i);

void poly_coeff_mul(Array1D<double>& out, Array1D<double>& in1, Array1D<double>& in2);

void poly_coeff_int(Array1D<double>& out, Array1D<double>& in);

void poly_coeff_diff(Array1D<double>& out, Array1D<double>& in);

double poly_eval(const Array1D<double>& in, double c);

void lobatto(Array1D<double>& mesh);

void gauss(Array1D<double>& mesh);

double poly_lgr_eval( const Array1D<double>& t, int i, double c);

double poly_dlg_eval( const Array1D<double>& t, int i, double c);
