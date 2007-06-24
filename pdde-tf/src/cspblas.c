/* ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages packages root directory
//
// ------------------------------------------------------------------------- */

#include "cspblas.h"

static
inline void SpMulR(const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, const double* in, const double alpha, const int incIN)
{
  register int j, k;
  register double r0, r1, r2, r3, r4, r5, r6, r7, t0, t1;

  for (j = 0; j < n; j += 1)
  {
    out[j] = 0.0;
    for (k = Ap[j]; k + 7 < Ap[j+1]; k += 8)
    {
      r0 = Ax[k] * in[incIN*Ai[k]];
      r1 = Ax[k+1] * in[incIN*Ai[k+1]];
      r2 = Ax[k+2] * in[incIN*Ai[k+2]];
      r3 = Ax[k+3] * in[incIN*Ai[k+3]];
      t0 = r0 + r1 + r2 + r3;
      r4 = Ax[k+4] * in[incIN*Ai[k+4]];
      r5 = Ax[k+5] * in[incIN*Ai[k+5]];
      r6 = Ax[k+6] * in[incIN*Ai[k+6]];
      r7 = Ax[k+7] * in[incIN*Ai[k+7]];
      t1 = r4 + r5 + r6 + r7;
      out[j] += t0 + t1;
    }
    for (; k < Ap[j+1]; k++)
    {
      out[j] += Ax[k] * in[incIN*Ai[k]];
    }
    out[j] *= alpha;
  }
}

static
inline void SpMulRpY(const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                     double* out, const double* in, const double alpha,
                     const double* Y, double beta, const int incIN, const int incY)
{
  register int j, k;
  register double r0, r1, r2, r3, r4, r5, r6, r7, t0, t1;

  for (j = 0; j < n; j += 1)
  {
    out[j] = 0.0;
    for (k = Ap[j]; k + 7 < Ap[j+1]; k += 8)
    {
      r0 = Ax[k] * in[incIN*Ai[k]];
      r1 = Ax[k+1] * in[incIN*Ai[k+1]];
      r2 = Ax[k+2] * in[incIN*Ai[k+2]];
      r3 = Ax[k+3] * in[incIN*Ai[k+3]];
      t0 = r0 + r1 + r2 + r3;
      r4 = Ax[k+4] * in[incIN*Ai[k+4]];
      r5 = Ax[k+5] * in[incIN*Ai[k+5]];
      r6 = Ax[k+6] * in[incIN*Ai[k+6]];
      r7 = Ax[k+7] * in[incIN*Ai[k+7]];
      t1 = r4 + r5 + r6 + r7;
      out[j] += t0 + t1;
    }
    for (; k < Ap[j+1]; k++)
    {
      out[j] += Ax[k] * in[incIN*Ai[k]];
    }
    out[j] = out[j] * alpha + Y[incY*j] * beta;
  }
}

static
inline void SpMulC(const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, const double* in, const double alpha, const int incIN)
{
  register int i, j, k;
  register double t0;

  for (i = 0; i + 7 < m; i += 8)
  {
    out[i] = 0.0;
    out[i+1] = 0.0;
    out[i+2] = 0.0;
    out[i+3] = 0.0;
    out[i+4] = 0.0;
    out[i+5] = 0.0;
    out[i+6] = 0.0;
    out[i+7] = 0.0;
  }
  for (; i < m; i += 1)
  {
    out[i] = 0.0;
  }
  for (j = 0; j < n; j += 1)
  {
    t0 = alpha * in[incIN*j];
    k = Ap[j];
    for (k = Ap[j]; k + 7 < Ap[j+1]; k += 8)
    {
      out[Ai[k]] += Ax[k] * t0;
      out[Ai[k+1]] += Ax[k+1] * t0;
      out[Ai[k+2]] += Ax[k+2] * t0;
      out[Ai[k+3]] += Ax[k+3] * t0;
      out[Ai[k+4]] += Ax[k+4] * t0;
      out[Ai[k+5]] += Ax[k+5] * t0;
      out[Ai[k+6]] += Ax[k+6] * t0;
      out[Ai[k+7]] += Ax[k+7] * t0;
    }
    for (; k < Ap[j+1]; k += 1)
    {
      out[Ai[k]] += Ax[k] * t0;
    }
  }
}

static
inline void SpMulCpY(const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                     double* out, const double* in, const double alpha,
                     const double* Y, const double beta, const int incIN, const int incY)
{
  register int i, j, k;
  register double t0;

  for (i = 0; i + 7 < m; i += 8)
  {
    out[i] = Y[incY*i] * beta;
    out[i+1] = Y[incY*(i+1)] * beta;
    out[i+2] = Y[incY*(i+2)] * beta;
    out[i+3] = Y[incY*(i+3)] * beta;
    out[i+4] = Y[incY*(i+4)] * beta;
    out[i+5] = Y[incY*(i+5)] * beta;
    out[i+6] = Y[incY*(i+6)] * beta;
    out[i+7] = Y[incY*(i+7)] * beta;
  }
  for (i = 0; i < m; i += 1)
  {
    out[i] = Y[incY*i] * beta;
  }
  for (j = 0; j < n; j += 1)
  {
    t0 = alpha * in[incIN*j];
    for (k = Ap[j]; k + 7 < Ap[j+1]; k += 8)
    {
      out[Ai[k]] += Ax[k] * t0;
      out[Ai[k+1]] += Ax[k+1] * t0;
      out[Ai[k+2]] += Ax[k+2] * t0;
      out[Ai[k+3]] += Ax[k+3] * t0;
      out[Ai[k+4]] += Ax[k+4] * t0;
      out[Ai[k+5]] += Ax[k+5] * t0;
      out[Ai[k+6]] += Ax[k+6] * t0;
      out[Ai[k+7]] += Ax[k+7] * t0;
    }
    for (; k < Ap[j+1]; k += 1)
    {
      out[Ai[k]] += Ax[k] * t0;
    }
  }
}

void cspblas_mmx(char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                 double* out, const double* in, double alpha)
{
  if (alpha == 1.0)
  {
    if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
    {
      SpMulC(n, m, Ap, Ai, Ax, out, in, 1.0, 1);
    }
    else
    {
      SpMulR(n, m, Ap, Ai, Ax, out, in, 1.0, 1);
    }
  }
  else
    if (alpha == -1.0)
    {
      if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
      {
        SpMulC(n, m, Ap, Ai, Ax, out, in, -1.0, 1);
      }
      else
      {
        SpMulR(n, m, Ap, Ai, Ax, out, in, -1.0, 1);
      }
    }
    else
    {
      if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
      {
        SpMulC(n, m, Ap, Ai, Ax, out, in, alpha, 1);
      }
      else
      {
        SpMulR(n, m, Ap, Ai, Ax, out, in, alpha, 1);
      }
    }
}

void cspblas_mmxpy(char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, const double* in, double alpha, const double* C, double beta)
{
  if (alpha == 1.0)
  {
    if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
    {
      SpMulCpY(n, m, Ap, Ai, Ax, out, in, 1.0, C, beta, 1, 1);
    }
    else
    {
      SpMulRpY(n, m, Ap, Ai, Ax, out, in, 1.0, C, beta, 1, 1);
    }
  }
  else
    if (alpha == -1.0)
    {
      if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
      {
        SpMulCpY(n, m, Ap, Ai, Ax, out, in, -1.0, C, beta, 1, 1);
      }
      else
      {
        SpMulRpY(n, m, Ap, Ai, Ax, out, in, -1.0, C, beta, 1, 1);
      }
    }
    else
    {
      if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
      {
        SpMulCpY(n, m, Ap, Ai, Ax, out, in, alpha, C, beta, 1, 1);
      }
      else
      {
        SpMulRpY(n, m, Ap, Ai, Ax, out, in, alpha, C, beta, 1, 1);
      }
    }
}

void cspblas_mmxm(char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                  double* out, int ldout, const double* in, int ldin, double alpha, int nrhs)
{
  int c;
  if (alpha == 1.0)
  {
    if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulC(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, 1);
      }
    }
    else
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulR(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, 1);
      }
    }
  }
  else
  {
    if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulC(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, 1);
      }
    }
    else
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulR(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, 1);
      }
    }
  }
}

void cspblas_mmxmpym(char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                     double* out, int ldout,
                     const double* in, int ldin, double alpha,
                     const double* Y, int ldY, double beta, int nrhs)
{
  int c;
  if (alpha == 1.0)
  {
    if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulCpY(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, &Y[c*ldY], beta, 1, 1);
      }
    }
    else
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulRpY(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, &Y[c*ldY], beta, 1, 1);
      }
    }
  }
  else
  {
    if (((format == 'C') && (trans == NoTrans)) || ((format == 'R') && (trans == Trans)))
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulCpY(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, &Y[c*ldY], beta, 1, 1);
      }
    }
    else
    {
      for (c = 0; c < nrhs; c++)
      {
        SpMulRpY(n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, &Y[c*ldY], beta, 1, 1);
      }
    }
  }
}
