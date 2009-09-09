// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "odepoint.h"
#include "odecolloc.h"
#include "matrix.h"
#include "hypermatrix.h"

#define NDIM colloc->Ndim()
#define NPAR colloc->Npar()
#define NINT colloc->Nint()
#define NDEG colloc->Ndeg()

ODEPoint::ODEPoint(System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg)
 : PerSolPoint(sys_, eqn_, var_, sys_.ndim()*(ndeg*nint + 1), sys_.ndim()*(ndeg*nint + 1)*sys_.ndim()*(ndeg+1), sys_.ndim())
{
  colloc = new ODEColloc(sys_, nint, ndeg);
  basecolloc = colloc;
  persolcolloc = colloc;
  Construct();
  FillSol(sys_);
  par(NPAR + ParAngle) = 0.0;
  par(NPAR + ParRot) = 0.0;
}

ODEPoint::~ODEPoint()
{
  delete colloc;
  Destruct();
}

// its unmodified from point.cpp. Only some stuff is removed.
void ODEPoint::Jacobian
(
  HyperMatrix& AA, HyperVector& RHS, // output
  Vector& parPrev, Vector& par,      // parameters
  Vector& solPrev, Vector& sol,      // solution
  Array1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
  double ds, bool cont               // ds stepsize, cont: true if continuation
)
{
  if (eqn(0) == EqnODESol)
  {
    colloc->RHS_x(AA.getA11(), par, sol);
    colloc->RHS(RHS.getV1(), par, sol);
    for (int i = 1; i < varMap.Size(); i++)
    {
      if (varMap(i) < NPAR)
      {
        colloc->RHS_p(AA.getA13(i - 1), par, sol, varMap(i));
      }
      else if (varMap(i) - NPAR == ParAngle)
      {
        AA.getA13(i - 1).Clear();
      }
      else
      {
        P_MESSAGE5("Non-existing parameter P", varMap(i), " was specified at position ", i, ".");
      }
    }
  } else P_MESSAGE3("Invalid first equation, that is not EqnODESol: ", eqn(0), ".");
  // ADDITIONAL EQUATIONS
  // ONLY THE PHASE CONDITION AT THE MOMENT
  for (int i = 1; i < eqn.Size(); i++)
  {
    switch (eqn(i))
    {
        // Phase conditions.
      case EqnPhase:
        colloc->PhaseStar(AA.getA31(i - 1), solPrev); // this should be the previous solution!!!
        // other variables
        for (int j = 1; j < varMap.Size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            AA.getA33()(i - 1, j - 1) = 0.0;
          }
          else if (varMap(j) - NPAR == ParAngle)
          {
            AA.getA33(i - 1, j - 1) = 0.0;
          }
          else
          {
            P_MESSAGE5("Non-existing parameter P", varMap(j), " was specified at position ", j, ".");
          }
        }
        RHS.getV3()(i - 1) = -(AA.getA31(i - 1) * sol);
        break;
      default:
        P_MESSAGE5("Unknown equation type ", eqn(i), " is encountered at position ", varMap(i), ".");
        break;
    }
  }
  if (cont)
  {
    // copying the tangent
    if (dim1 != 0) colloc->Star(AA.getA31(dim3), xxDot->getV1());
    for (int i = 0; i < xxDot->getV3().Size(); i++) AA.getA33()(dim3, i) = xxDot->getV3()(i);
    if (ds != 0.0)
    {
      RHS.getV3()(dim3) = ds - colloc->IntegrateCont(xxDot->getV1(), sol, solPrev);
      for (int j = 1; j < varMapCont.Size(); ++j) RHS.getV3()(dim3) -= xxDot->getV3()(j - 1) * (par(varMap(j)) - parPrev(varMap(j)));
    }
    else
    {
      RHS.getV3()(dim3) = 0.0;
    }
  }
}