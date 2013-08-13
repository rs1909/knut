// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "stpoint.h"
#include "stcolloc.h"
#include "exprsystem.h"
#include "mat4data.h"

#define NDIM (colloc->nDim())
#define NPAR (colloc->nPar())

class KNAbstractContinuation;
          
KNSteadyStateSolution::KNSteadyStateSolution(KNAbstractContinuation* cnt, KNExprSystem& sys_, 
  KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_)
  : KNAbstractPoint(cnt, sys_, eqn_, var_, sys_.ndim(), sys_.ndim()*sys_.ndim()),
   jacStab('R', sys_.ndim(), sys_.ndim()*sys_.ndim())
{
  colloc = new KNSteadyStateJacobian(sys_);
  par(VarToIndex(VarAngle,colloc->nPar())) = 0.0;
  par(VarToIndex(VarRot,colloc->nPar())) = 0.0;
  par(VarToIndex(VarPeriod,NPAR)) = 1.0;
  basecolloc = colloc;
  construct();

  // filling solution
  KNArray1D<double> tim((size_t)1);
  tim(0) = 0.0;
  KNArray2D<double> tmp(NDIM, 1);
  sys_.stsol (tmp, tim);
  for (size_t q=0; q < NDIM; q++) sol(q) = tmp(q,0);
  sys_.stpar(par);
}

KNSteadyStateSolution::~KNSteadyStateSolution()
{
  destruct();
  delete colloc;
}

// its unmodified from point.cpp. Only some stuff is removed.
void KNSteadyStateSolution::jacobian(
  KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
  KNVector& parPrev, KNVector& par,      // parameters
  KNVector& solPrev, KNVector& sol,      // solution
  KNArray1D<size_t>& varMap,           // contains the variables. If cont => contains the P0 too.
  double ds, bool cont               // ds stepsize, cont: true if continuation
)
{
  if (eqn(0) == EqnSteady)
  {
    colloc->rightHandSide_x(AA.getA11(), par, sol);
    colloc->rightHandSide(RHS.getV1(), par, sol);
    for (size_t i = 1; i < varMap.size(); i++)
    {
      if (varMap(i) < NPAR)
      {
        colloc->rightHandSide_p(AA.getA13(i - 1), par, sol, varMap(i));
      }
      else if (IndexToVar(varMap(i),NPAR) == VarAngle)
      {
        AA.getA13(i - 1).clear();
      }
      else
      {
        P_MESSAGE5("Non-existing parameter P", varMap(i), " was specified at position ", i, ".");
      }
    }
  } else P_MESSAGE3("Invalid first equation, that is not EqnSteady: ", eqn(0), ".");
  // ADDITIONAL EQUATIONS
  // ONLY THE PHASE CONDITION AT THE MOMENT
  for (size_t i = 1; i < eqn.size(); i++)
  {
    switch (eqn(i))
    {
      // NO ADDITIONAL EQUATIONS
      default:
        P_MESSAGE5("Unknown equation type ", eqn(i), " is encountered at position ", varMap(i), ".");
        break;
    }
  }
  if (cont)
  {
    // copying the tangent
    if (dim1 != 0) colloc->star(AA.getA31(dim3), xxDot->getV1());
    for (size_t i = 0; i < xxDot->getV3().size(); i++) AA.getA33()(dim3, i) = xxDot->getV3()(i);
    if (ds != 0.0)
    {
      RHS.getV3()(dim3) = ds - colloc->integrateWithCp(xxDot->getV1(), sol, solPrev);
      for (size_t j = 1; j < varMapCont.size(); ++j) RHS.getV3()(dim3) -= xxDot->getV3()(j - 1) * (par(varMap(j)) - parPrev(varMap(j)));
    }
    else
    {
      RHS.getV3()(dim3) = 0.0;
    }
//    AA.getA31(dim3).print();
//    AA.getA13(dim3).print();
//    AA.getA33().print();
  }
//  sol.print();
//  par.print();
}

void KNSteadyStateSolution::BinaryRead(KNDataFile& data, size_t n)
{
  data.lockRead();
  P_ERROR_X1(data.getNPar() == VarToIndex(VarEnd,NPAR), "Wrong number of parameters in the input MAT file.");
  data.getPar(n, par);
//  data.getMul(n, mRe, mIm);
  P_ERROR_X1(data.getNDim() == NDIM, "Wrong number of dimensions in the input MAT file.");
  P_ERROR_X1(data.getNInt() == 0, "Nonzero intervals."); 
  P_ERROR_X1(data.getNDeg() == 0, "Nonzero degree.");
  data.getProfile(n, sol);
  data.unlock();
}

void KNSteadyStateSolution::BinaryWrite(KNDataFile& data, BifType bif, size_t n)
{
        //  std::cout << "DAT " << &data << "\n";
  data.lockWrite();
//  data.setNTrivMul(0, nTrivMulLP);
//  data.setNTrivMul(1, nTrivMulPD);
//  data.setNTrivMul(2, nTrivMulNS);
  //
  data.setMagic(n, bif);
  data.setPar(n, par);
//  data.setMul(n, mRe, mIm);
//  data.setElem(n, persolcolloc->getElem());
//  data.setMesh(n, persolcolloc->getMesh());
  data.setProfile(n, sol);
  data.unlock();
}
