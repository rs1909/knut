// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cmath>

#include "knerror.h"
#include "system.h"
#include "torpoint.h"
#include "mat4data.h"

PointTR::PointTR(System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int ndeg1_, int ndeg2_, int nint1_, int nint2_)
  : BasePoint(sys_, eqn_, var_, sys_.ndim()*((ndeg1_*nint1_))*((ndeg2_*nint2_)),
    sys_.ndim()*(ndeg1_*nint1_)*(ndeg2_*nint2_)*sys_.ndim()*(sys_.ntau()+1)*(ndeg1_+1)*(ndeg2_+1))
{
  colloc = new CollocTR(sys_, ndeg1_, ndeg2_, nint1_, nint2_);
  basecolloc = colloc;

  P_ERROR_X1((eqn(0) == EqnTORSol) && (var(0) == VarTORSol), "The first equation must the boundary value problem of the periodic solution.");

  BasePoint::Construct();
}

PointTR::~PointTR()
{
  BasePoint::Destruct();
  delete colloc;
}

void PointTR::Jacobian
(
  HyperMatrix& jac, HyperVector& rhs, // output
  Vector& parPrev,  Vector& par,      // parameters
  Vector& solPrev,  Vector& sol,      // solution
  Array1D<int>&     varMap,           // contains the variables. If cont => contains the P0 too.
  double ds, bool cont                // ds stepsize, cont: true if continuation
)
{
  // uses also: eqn, var
  const int        nvar = dim3;
  Array1D<Vector*> A13(nvar + 1);
  Array1D<int>     JacVar(nvar + 1);

  // it is just for the first equation i.e. EqnTORSol
  for (int i = 1; i < varMap.Size(); i++)
  {
    JacVar(i-1) = varMap(i);
    A13(i-1) = & jac.getA13(i-1);
  }

  // first the phase conditions
  int ph0 = -1, ph1 = -1;
  for (int i = 1; i < eqn.Size(); i++)
  {
    if (eqn(i) == EqnTORPhase0)
    {
      P_ERROR_X1(ph0 == -1, "Too many phase conditions.");
      ph0 = i - 1;
    }
    if (eqn(i) == EqnTORPhase1)
    {
      P_ERROR_X1(ph1 == -1, "Too many phase conditions.");
      ph1 = i - 1;
    }
  }
  if (ph0 != (-1) && ph1 != (-1))
  {
    colloc->PhaseBOTH(jac.getA31(ph0), jac.getA31(ph1), solPrev);
    rhs.getV3()(ph0) = 0.0;
    rhs.getV3()(ph1) = 0.0;
  }
  else
  {
    P_ERROR_X1(ph0 == -1, "There is a first phase condition, but no second.");
    if (ph1 != (-1))
    {
      colloc->PhaseONE(jac.getA31(ph1), sol/*Prev*/);
      rhs.getV3()(ph1) = 0.0;
    }
  }

  // jacobian, derivatives of the right-hand side
  colloc->Jacobian(jac.getA11(), A13, rhs.getV1(), par, sol, JacVar);
  // the other equations
  // Currently: none
  for (int i = 1; i < eqn.Size(); i++)
  {
    switch (eqn(i))
    {
      case EqnTORPhase0:
      case EqnTORPhase1:
        // nothing happens
        break;
      default:
        P_MESSAGE3("An unknown type of equation is encountered: ", eqn(i), ".");
        break;
    }
  }

  // rhs from the tangent
  if (cont)
  {
    // copying the tangent
    if (dim1 != 0) colloc->Star(jac.getA31(dim3), xxDot->getV1());
    for (int i = 0; i < xxDot->getV3().Size(); i++) jac.getA33()(dim3, i) = xxDot->getV3()(i);
    if (ds != 0.0)
    {
      rhs.getV3()(dim3) = ds - (jac.getA31(dim3) * sol - jac.getA31(dim3) * solPrev);
      for (int j = 1; j < varMap.Size(); ++j) rhs.getV3()(dim3) -= xxDot->getV3()(j - 1) * (par(varMap(j)) - parPrev(varMap(j)));
    }
    else
    {
      rhs.getV3()(dim3) = 0.0;
    }
  }
}

#define NDEG1 colloc->Ndeg1()
#define NDEG2 colloc->Ndeg2()
#define NINT1 colloc->Nint1()
#define NINT2 colloc->Nint2()

void PointTR::ReadBinary(mat4Data& data, int n)
{
  data.getBlanket(n, sol);
  data.getPar(n, par);
}

void PointTR::WriteBinary(mat4Data& data, int n)
{
  data.setPar(n, par);
  for (int i = 0; i < NINT1; ++i)
  {
    for (int j = 0; j < NDEG1; ++j)
    {
      data.setMesh1(n, i*NDEG1 + j, (getMesh1()(j) + i) / ((double)NINT1));
    }
  }
  for (int i = 0; i < NINT2; ++i)
  {
    for (int j = 0; j < NDEG2; ++j)
    {
      data.setMesh2(n, i*NDEG2 + j, (getMesh2()(j) + i) / ((double)NINT2));
    }
  }
  data.setBlanket(n, sol);
}

