#include "polynomial.h"
#include <cmath>

/// returns the coefficient of the lagrangian interpolation of a
/// function at the i-th point on the mesh t
///
void poly_coeff_lgr(Array1D<double>& out, const Array1D<double>& t, int i)
{
  out.Clear();
  out(0) = 1.0;
  for (int j = 0; j < t.Size(); j++)
  {
    if (j != i)
    {
      poly_linmul(out, -t(j) / (t(i) - t(j)), 1.0 / (t(i) - t(j)));
    }
  }
}

/// multiplies the two polynomials in1, in2 in out
/// out = (in1(0) + in1(1)*x + ... + in1(n)*x^n) * (in2(0) + in2(1)*x + ... + in2(n)*x^n)
///
void poly_coeff_mul(Array1D<double>& out, Array1D<double>& in1, Array1D<double>& in2)
{
  if (out.Size() < (in1.Size() - 1) + (in2.Size() - 1) + 1) std::cout << "Error in \"poly_mul\"\n";
  out.Clear();
  for (int i = 0; i < in1.Size(); i++)
  {
    for (int j = 0; j < in2.Size(); j++)
    {
      out(i + j) = in1(i) * in2(j);
    }
  }
}

/// integrates the polynomial in and writes the result in out
///
void poly_coeff_int(Array1D<double>& out, Array1D<double>& in)
{
  out.Clear();
  out(0) = 0.0;
  for (int i = 0; i < in.Size() - 1; i++)
  {
    out(i + 1) += in(i) / (i + 1);
  }
}

/// differentiate the polynomial in and write to out
///
void poly_coeff_diff(Array1D<double>& out, Array1D<double>& in)
{
  out.Clear();
  out(in.Size() - 1) = 0.0;
  for (int i = 1; i < in.Size(); i++)
  {
    out(i - 1) += in(i) * i;
  }
}

/// evaluate the polynomial in at point c
///
double poly_eval(const Array1D<double>& in, double c)
{
  double tmp = in(in.Size() - 1);
  for (int i = in.Size() - 1; i > 0; --i)
  {
    tmp = in(i - 1) + c * tmp;
  }
  return tmp;
}

/// computes the lobatto mesh on [0,1] 
/// this is defined on the boundary
///
void lobatto(Array1D<double>& mesh)
{
  const int N = mesh.Size() - 1;
  if (mesh.Size() > 1) for (int i = 0; i < mesh.Size(); i++) mesh(i) = 1.0 / 2.0 * (1 - cos((i * M_PI) / N));
}

/// computes the gauss mesh on [0,1] 
/// this is defined on the boundary
///
void gauss(Array1D<double>& mesh)
{
  const int N = mesh.Size();
  if (mesh.Size() > 1) for (int i = 0; i < mesh.Size(); i++) mesh(i) = 1.0 / 2.0 * (1 - cos((2.0 * i + 1.0) / (2.0 * N + 2) * M_PI));
}
