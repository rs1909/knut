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

KNOdePeriodicSolution::KNOdePeriodicSolution(KNAbstractContinuation* cnt, KNSystem& sys_, 
  KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_, int nint, int ndeg)
  : KNAbstractPeriodicSolution(cnt, sys_, eqn_, var_, sys_.ndim()*(ndeg*nint + 1), sys_.ndim()*(ndeg*nint + 1)*sys_.ndim()*(ndeg+1), sys_.ndim()),
   jacStab('R', sys_.ndim()*(ndeg*nint + 1), sys_.ndim()*(ndeg*nint + 1)*sys_.ndim()*(ndeg+1)),
   matrixInitialCondition(sys_.ndim()*(ndeg*nint+1), sys_.ndim()),
   matrixSolution(sys_.ndim()*(ndeg*nint+1), sys_.ndim()),
   monodromyMatrix(sys_.ndim(), sys_.ndim())
{
  colloc = new KNOdeBvpCollocation(sys_, nint, ndeg);
  basecolloc = colloc;
  persolcolloc = colloc;
  construct();
  FillSol(sys_);
  par(colloc->nPar() + ParAngle) = 0.0;
  par(colloc->nPar() + ParRot) = 0.0;
}

KNOdePeriodicSolution::~KNOdePeriodicSolution()
{
  delete colloc;
  destruct();
}

#define NDIM (colloc->nDim())
#define NPAR (colloc->nPar())
#define NINT (colloc->nInt())
#define NDEG (colloc->nDeg())

// its unmodified from point.cpp. Only some stuff is removed.
void KNOdePeriodicSolution::jacobian(
  KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
  KNVector& parPrev, KNVector& par,      // parameters
  KNVector& solPrev, KNVector& sol,      // solution
  KNArray1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
  double ds, bool cont               // ds stepsize, cont: true if continuation
)
{
  if (eqn(0) == EqnODESol)
  {
    colloc->rightHandSide_x(AA.getA11(), par, sol);
    colloc->rightHandSide(RHS.getV1(), par, sol);
    for (int i = 1; i < varMap.size(); i++)
    {
      if (varMap(i) < NPAR)
      {
        colloc->rightHandSide_p(AA.getA13(i - 1), par, sol, varMap(i));
      }
      else if (varMap(i) - NPAR == ParAngle)
      {
        AA.getA13(i - 1).clear();
      }
      else
      {
        P_MESSAGE5("Non-existing parameter P", varMap(i), " was specified at position ", i, ".");
      }
    }
  } else P_MESSAGE3("Invalid first equation, that is not EqnODESol: ", eqn(0), ".");
  // ADDITIONAL EQUATIONS
  // ONLY THE PHASE CONDITION AT THE MOMENT
  for (int i = 1; i < eqn.size(); i++)
  {
    switch (eqn(i))
    {
        // Phase conditions.
      case EqnPhase:
        colloc->phaseStar(AA.getA31(i - 1), solPrev); // this should be the previous solution!!!
        // other variables
        for (int j = 1; j < varMap.size(); j++)
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
    if (dim1 != 0) colloc->star(AA.getA31(dim3), xxDot->getV1());
    for (int i = 0; i < xxDot->getV3().size(); i++) AA.getA33()(dim3, i) = xxDot->getV3()(i);
    if (ds != 0.0)
    {
      RHS.getV3()(dim3) = ds - colloc->integrateWithCp(xxDot->getV1(), sol, solPrev);
      for (int j = 1; j < varMapCont.size(); ++j) RHS.getV3()(dim3) -= xxDot->getV3()(j - 1) * (par(varMap(j)) - parPrev(varMap(j)));
    }
    else
    {
      RHS.getV3()(dim3) = 0.0;
    }
  }
}

void KNOdePeriodicSolution::Stability(bool init)
{
  for (int i = 0; i < NDIM; ++i)
  {
    matrixInitialCondition(i,i) = 1.0;
  }
  colloc->jacobianOfStability(jacStab, par);
  jacStab.solve(matrixSolution, matrixInitialCondition);
  for (int i = 0; i < NDIM; ++i)
  {
    for (int j = 0; j< NDIM; ++j)
    {  
      monodromyMatrix(i,j) = matrixSolution(NDIM*NDEG*NINT+i,j);
    }
  }
  monodromyMatrix.eigenvalues(mRe, mIm);
}
