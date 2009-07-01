// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "polynomial.h"
#include <cmath>

void poly_lgr(const Vector& t, Vector &out, double c)
{

  P_ASSERT_X1(t.Size() == out.Size(), "poly_lgr: wrong dimensions");
  for (int i = 0; i < t.Size(); i++)
  {
    out(i) = 1.0;
    for (int j = 0; j < t.Size(); j++)
    {
      if (i != j)
      {
        out(i) *= (c - t(j)) / (t(i) - t(j));
      }
    }
  }
}

void poly_dlg(const Vector& t, Vector& out, double c)
{
  int j, k, l;
  double f;

  P_ASSERT_X1(t.Size() == out.Size(), "poly_dlg: wrong dimensions");

  for (j = 0; j < t.Size(); j++)
  {
    out(j) = 0.0;
    for (k = 0; k < t.Size(); k++)
    {
      if (k != j)
      {
        f = 1.0;
        for (l = 0; l < t.Size(); l++)
        {
          if ((l != k) && (l != j)) f *= (c - t(l)) / (t(j) - t(l));
        }
        out(j) += f / (t(j) - t(k));
      }
    }
  }
}

void poly_d2lg(const Vector& t, Vector& out, double c)
{
  int i, j, k, l;
  double f;

  P_ASSERT_X1(t.Size() == out.Size(), "poly_d2lg: wrong dimensions.");

  for (i = 0; i < t.Size(); i++)
  {
    out(i) = 0.0;
    for (l = 0; l < t.Size(); l++)
    {
      if (l==i) continue;
      for (k = 0; k < t.Size(); k++)
      {
        if (k==l) continue;
        if (k==i) continue;
        f = 1.0;
        for (j = 0; j < t.Size(); j++)
        {
          if (j==k) continue;
          if (j==l) continue;
          if (j==i) continue;
          f *= (c - t(j)) / (t(i) - t(j));
        }
        out(i) += f / (t(i) - t(k)) / (t(i) - t(l));
      }
    }
  }
}

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
  P_ASSERT_X1(out.Size() > (in1.Size() - 1) + (in2.Size() - 1), "poly_coeff_mul: wrong dimensions.");
  out.Clear();
  for (int i = 0; i < in1.Size(); i++)
  {
    for (int j = 0; j < in2.Size(); j++)
    {
      out(i + j) += in1(i) * in2(j);
    }
  }
}

/// integrates the polynomial in and writes the result in out
///
void poly_coeff_int(Array1D<double>& out, Array1D<double>& in)
{
  out.Clear();
  for (int i = 0; i < in.Size() - 1; i++)
  {
    out(i + 1) = in(i) / (i + 1);
  }
  P_ASSERT_X1(in(in.Size() - 1) == 0.0, "poly_coeff_int: overflow.");
}

/// differentiate the polynomial in and write to out
///
void poly_coeff_diff(Array1D<double>& out, Array1D<double>& in)
{
  out.Clear();
  for (int i = 1; i < in.Size(); i++)
  {
    out(i - 1) = in(i) * i;
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

double poly_lgr_eval( const Array1D<double>& t, int i, double c)
{
  double res = 1.0;
  for (int j = 0; j < t.Size(); j++)
  {
    if (j != i)
    {
      res *= (c-t(j))/(t(i)-t(j));
    }
  }
  return res;
}

double poly_dlg_eval( const Array1D<double>& t, int i, double c)
{
  double res = 0.0;
  for (int k = 0; k < t.Size(); ++k)
  {
    if (k != i)
    {
      double prod = 1.0/(t(i)-t(k));
      for (int j = 0; j < t.Size(); ++j)
      {
        if ((j != i)&&(j != k))
        {
          prod *= (c-t(j))/(t(i)-t(j));
        }
      }
      res += prod;
    }
  }
  return res;
}
