// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "config.h"

#include "knerror.h"
#include "point.h"
#include "ncolloc.h"
#include "exprsystem.h"
#include "matrix.h"
#include "spmatrix.h"
#include "plot.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

#define NDIM colloc->nDim()
#define NTAU colloc->nTau()
#define NPAR colloc->nPar()
#define NINT colloc->nInt()
#define NDEG colloc->nDeg()


KNDdePeriodicSolution::KNDdePeriodicSolution(KNAbstractContinuation* cnt, KNExprSystem& sys_, 
  KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_, size_t nint, size_t ndeg, size_t nmul) 
  : KNAbstractPeriodicSolution(cnt, sys_, eqn_, var_, sys_.ndim()*(ndeg*nint + 1), sys_.ndim()*(ndeg*nint + 1)*sys_.ntau()*sys_.ndim()*(ndeg+1), nmul),
    jacStab(0)
{
  colloc = new KNDdeBvpCollocation(sys_, nint, ndeg);
  basecolloc = colloc;
  persolcolloc = colloc;

  construct();
  FillSol(sys_);
  par(VarToIndex(VarAngle,NPAR)) = 0.0;
  par(VarToIndex(VarRot,NPAR)) = 0.0;
}

KNDdePeriodicSolution::~KNDdePeriodicSolution()
{
  destruct();
  delete jacStab;
  delete colloc;
}

// private
// What does it construct?
// charMat
// xxDot, xx, rhs, jac
void KNDdePeriodicSolution::construct()
{
  P_ERROR_X1((eqn(0) == EqnSol) && (var(0) == VarSol), "The first equation must be the boundary value problem of the periodic solution.");

  KNAbstractPoint::construct();
  
  // a) setting the test functionals b) determining the number of trivial multipliers
  testFunct.init(eqn.size());
  nTrivMulLP = 0;
  nTrivMulPD = 0;
  nTrivMulNS = 0;
  for (size_t i = 1; i < eqn.size(); i++)
  {
    switch (eqn(i))
    {
      case EqnTFLP:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new KNTestFunctional(*colloc, 1.0);
        ++nTrivMulLP;
        break;
      case EqnTFPD:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new KNTestFunctional(*colloc, -1.0);
        ++nTrivMulPD;
        break;
      case EqnTFLPAUT:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new KNLpAutTestFunctional(*colloc, 1.0);
        ++nTrivMulLP;
        break;
      case EqnTFLPAUTROT:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new KNLpAutRotTestFunctional(*colloc, rotRe, rotIm, 1.0);
        ++nTrivMulLP;
        break;
      case EqnTFCPLX_RE:
        P_ERROR_X1(eqn(i + 1) == EqnTFCPLX_IM,
                   "The real and imaginary parts of the complex test functional are not paired.");
        P_ERROR(testFunct(i) == 0);
        testFunct(i) = new KNComplexTestFunctional(*colloc);
        nTrivMulNS += 2;
        break;
      case EqnTFCPLX_IM:
        P_ERROR_X1(eqn(i - 1) == EqnTFCPLX_RE,
                   "The real and imaginary parts of the complex test functional are not paired.");
        break;
      case EqnPhase:
        ++nTrivMulLP;
        break;
      case EqnPhaseRot:
        ++nTrivMulLP;
        break;
      default:
        break;
    }
    if (testFunct(i)) testFunct(i)->setKernelTolerance(KernEps, KernIter);
  }
}

// private
void KNDdePeriodicSolution::destruct()
{
  for (size_t i=0; i<testFunct.size(); ++i) delete testFunct(i);
  KNAbstractPoint::destruct();
}

// **************************************************************************************************************** //
//
//                           THE JACOBIAN
//
// **************************************************************************************************************** //

void KNDdePeriodicSolution::jacobian(
  KNSparseBlockMatrix& AA, KNBlockVector& RHS,                      // output
  KNVector& parPrev, KNVector& par,                           // parameters
  KNVector& solPrev, KNVector& sol,                           // solution
  KNArray1D<size_t>& varMap,                                // contains the variables. If cont => contains the P0 too.
  double ds, bool cont)                                              // cont: true if continuation
{
// uses also: eqn, var, varMapCont
// ---------------------------------------------------------------------------------------------------------------- //
//
//                           The periodic solution part
//
// ---------------------------------------------------------------------------------------------------------------- //

  if (eqn(0) == EqnSol)
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
  }

// ---------------------------------------------------------------------------------------------------------------- //
//
//                           Other equations
//
// ---------------------------------------------------------------------------------------------------------------- //

  for (size_t i = 1; i < eqn.size(); i++)
  {
    switch (eqn(i))
    {
        // Phase conditions.
      case EqnPhase:
        colloc->phaseStar(AA.getA31(i - 1), solPrev); // this should be the previous solution!!!
        // other variables
        for (size_t j = 1; j < varMap.size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            AA.getA33()(i - 1, j - 1) = 0.0;
          }
          else if (IndexToVar(varMap(j),NPAR) == VarAngle)
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
      case EqnPhaseRot:
        colloc->phaseRotationStar(AA.getA31(i - 1), solPrev, rotRe, rotIm); // this should be the previous solution!!!
        // other variables
        for (size_t j = 1; j < varMap.size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            AA.getA33()(i - 1, j - 1) = 0.0;
          }
          else if (IndexToVar(varMap(j),NPAR) == VarAngle)
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
      case EqnTFLP:
      case EqnTFPD:
      case EqnTFLPAUT:
      case EqnTFLPAUTROT:
        RHS.getV3()(i - 1) = testFunct(i)->funct(*colloc, par, sol);
        testFunct(i)->funct_x(AA.getA31(i - 1), *colloc, par, sol);
        for (size_t j = 1; j < varMap.size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            AA.getA33()(i - 1, j - 1) = testFunct(i)->funct_p(*colloc, par, sol, varMap(j));
          }
          else if (IndexToVar(varMap(j),NPAR) == VarAngle)
          {
            AA.getA33(i - 1, j - 1) = 0.0;
          }
          else
          {
            P_MESSAGE5("Non-existing parameter P", varMap(j), " was specified at position ", j, ".");
          }
        }
        break;
      case EqnTFCPLX_RE:
        P_ERROR_X3(eqn(i + 1) == EqnTFCPLX_IM,
                   "The real and imaginary parts of the complex test functional are not paired at position ", varMap(i), ".");
        testFunct(i)->funct(RHS.getV3()(i - 1), RHS.getV3()(i), *colloc, par, sol,
                         cos(par(VarToIndex(VarAngle,NPAR))), sin(par(VarToIndex(VarAngle,NPAR))));
        testFunct(i)->funct_x(AA.getA31(i - 1), AA.getA31(i), *colloc, par, sol);
        for (size_t j = 1; j < varMap.size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            testFunct(i)->funct_p(AA.getA33()(i - 1, j - 1), AA.getA33()(i, j - 1), *colloc, par, sol, varMap(j));
          }
          else if (IndexToVar(varMap(j),NPAR) == VarAngle)
          {
            testFunct(i)->funct_z(AA.getA33()(i - 1, j - 1), AA.getA33()(i, j - 1), *colloc, par, sol);
          }
          else
          {
            P_MESSAGE5("Non-existing parameter P", varMap(j), " was specified at position ", j, ".");
          }
        }
        break;
      case EqnTFCPLX_IM:
        P_ERROR_X3(eqn(i - 1) == EqnTFCPLX_RE,
                   "The real and imaginary parts of the complex test functional are not paired at position ", varMap(i), ".");
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
  }
}

// **************************************************************************************************************** //
//
//                          !END! THE JACOBIAN !END!
//
// **************************************************************************************************************** //

void KNDdePeriodicSolution::postProcess()
{
  // update test functional
  for (size_t idx=0; idx < testFunct.size(); ++idx)
  {
    if (testFunct(idx)) testFunct(idx)->initStep(*colloc);
  }
}

/// --------------------------------------------------------------
/// Starting bifurcation continuation using TEST FUNCTIONS
/// --------------------------------------------------------------

void KNDdePeriodicSolution::SwitchTFHB(double ds)
{
  std::ostream& out = outStream();
  size_t idx = eqn.size();
  for (size_t i=0; i<eqn.size(); ++i) if (eqn(i) == EqnTFCPLX_RE) { idx = i; break; }
  P_ERROR_X1(idx != eqn.size(), "No test functional was selected for Hopf bifurcation switch!" );
  KNComplexTestFunctional *tf = static_cast<KNComplexTestFunctional*>(testFunct(idx));
  KNVector QRE(NDIM), QIM(NDIM);

  double newperiod = par(0);
  double dist = tf->kernelComplex(newperiod, QRE, QIM, *colloc, par);
  par(0) = newperiod;
  
  out << "    T = " << par(0) << ", arg(Z) = " << par(VarToIndex(VarAngle,NPAR)) / (2*M_PI) << " * 2Pi," << " DST=" << dist << "\n";
//   std::cerr << "Distance between the predicted and real critial eigenfunctions: " << dist << "\n";
  P_ERROR_X3(dist < 0.01, "Cannot switch branches."
    " The predicted and the computed eigenfunctions are not close enough."
    " Check your selection of P0, perhaps it is not small enough. The distance (dist/norm1) is ", dist, ".");

#ifdef DEBUG
  std::cout << "Printing: neweigenvec\n";
  std::ofstream file("neweigenvec");
  file << std::scientific;
  file.precision(12);
#endif

  for (size_t i = 0; i < NINT; i++)
  {
    for (size_t j = 0; j < NDEG + 1; j++)
    {
      const double t = colloc->Profile(i, j);
      for (size_t p = 0; p < NDIM; p++)
      {
        xxDot->getV1()(p + (j + i*NDEG)*NDIM) = cos(2.0 * M_PI * t) * QRE(p) + sin(2.0 * M_PI * t) * QIM(p);
#ifdef DEBUG
        file << xxDot->getV1()(p + (j + i*NDEG)*NDIM) << "\t";
#endif
      }
#ifdef DEBUG
      file << par(0)*t << "\n";
#endif
    }
  }
  const double norm = sqrt(colloc->integrate(xxDot->getV1(), xxDot->getV1()));
  xxDot->getV1() /= norm;
  xxDot->getV3().clear();
  KNVector eql(NDIM);
  for (size_t p = 0; p < NDIM; ++p) eql(p) = sol(p);
  for (size_t i = 0; i < NDEG*NINT + 1; ++i)
  {
    for (size_t p = 0; p < NDIM; ++p)
    {
      sol(p + i*NDIM) = eql(p) + ds * xxDot->getV1()(p + i * NDIM);
    }
  }
}

/// Switching with the test functionals!!!

void KNDdePeriodicSolution::SwitchTFLP(BranchSW type, double ds)
{
  xxDot->getV1().clear();
  xxDot->getV3().clear();

  KNAbstractTestFunctional* tf = 0;
  switch (type)
  {
    case TFBRSwitch:
      tf = static_cast<KNAbstractTestFunctional*>(new KNTestFunctional(*colloc, 1.0));
      break;
    case TFBRAUTSwitch:
      tf = static_cast<KNAbstractTestFunctional*>(new KNLpAutTestFunctional(*colloc, 1.0));
      break;
    case TFBRAUTROTSwitch:
      tf = static_cast<KNAbstractTestFunctional*>(new KNLpAutRotTestFunctional(*colloc, rotRe, rotIm, 1.0));
      break;
    default:
      return;
      break;
  }
  tf->setKernelTolerance(KernEps, KernIter);
  tf->funct(*colloc, par, sol);
  tf->kernel(xxDot->getV1());
  delete tf;
  double norm = sqrt(colloc->integrate(xxDot->getV1(), xxDot->getV1()));
  xxDot->getV1() /= norm;
  xxDot->getV3().clear();

  sol += ds * xxDot->getV1();
}

void KNDdePeriodicSolution::SwitchTFPD(double ds)
{
  KNVector tan(xxDot->getV1().size());
  KNTestFunctional* tf = new KNTestFunctional(*colloc, -1.0);
  tf->setKernelTolerance(KernEps, KernIter);
  tf->funct(*colloc, par, sol);
  tf->kernel(tan);
  delete tf;

  // setting the period two double
  par(0) *= 2.0;
  par(VarToIndex(VarPeriod,NPAR)) *= 2.0;

  solNu = sol;
  colloc->pdMeshConvert(sol, xxDot->getV1(), solNu, tan);

  double norm = sqrt(colloc->integrate(xxDot->getV1(), xxDot->getV1()));
  xxDot->getV1() /= norm;
  xxDot->getV3().clear();
  sol += ds * xxDot->getV1();
}


void KNDdePeriodicSolution::Stability(bool init)
{
  std::ostream& out = outStream();
  mRe.clear();
  mIm.clear();

  // it is not always necessary to initialize 
  // especially not after a continuation step
  if (init) colloc->init(sol, par);
  size_t nmat = colloc->nMat();

  if (jacStab)
  {
  	if (jacStab->nmat() != nmat)
  	{
      out << "(NMAT:" << jacStab->nmat() << "->" << nmat << ")";
      printStream();
  	  delete jacStab;
  	  jacStab = new KNSparseMatrixPolynomial(nmat, NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1));
  	}
  } else
  {
  	out << "(NMAT=" << nmat << ")";
  	printStream();
  	jacStab = new KNSparseMatrixPolynomial(nmat, NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1));
  }

  colloc->jacobianOfStability(*jacStab, par);

  jacStab->eigenvalues(mRe, mIm);
//  mRe.print();
//  mIm.print();
}

void KNDdePeriodicSolution::Plot(GnuPlot& pl)
{
  pl.SetPointSize(0.8);

  for (size_t i = 0; i < NDIM; i++)
  {
    pl.Plot(static_cast<unsigned int>(i), "with lines");
    for (size_t j = 0; j < NINT*NDEG + 1; j++)
    {
      double t = (double)j / ((double)(NINT * NDEG));
      pl.AddData(static_cast<unsigned int>(i), t, sol(i + NDIM*j));
    }
  }
  pl.Show();
}

///--------------------------------
/// INPUT & OUTPUT
///--------------------------------

void KNDdePeriodicSolution::Write(std::ofstream& file)
{
  KNVector msh(colloc->getMesh());

  file << VarToIndex(VarEnd,NPAR) << "\t";
  for (size_t i = 0; i < VarToIndex(VarEnd,NPAR); i++) file << par(i) << "\t";

  file << mRe.size() << "\t";
  for (size_t i = 0; i < mRe.size(); i++) file << mRe(i) << "\t" << mIm(i) << "\t";

  file << NDIM << "\t";
  file << NINT << "\t";
  file << NDEG << "\t";
  for (size_t i = 0; i < NINT; i++)
  {
    for (size_t j = 0; j < NDEG; j++)
    {
      file << colloc->Profile(i, j) << "\t";
    }
  }
  file << colloc->getMesh()(NINT) << "\t";
  for (size_t i = 0; i < NDIM*(NINT*NDEG + 1); i++) file << sol(i) << "\t";
  file << "\n";
  file.flush();
}

void KNDdePeriodicSolution::Read(std::ifstream& file)
{
  size_t npar_, nmul_, ndim_, nint_, ndeg_;
  file >> npar_;
  P_ERROR_X3(VarToIndex(VarEnd,NPAR) == npar_, "Not compatible input file (NPAR) ", npar_, ".");
  for (size_t i = 0; i < VarToIndex(VarEnd,NPAR); i++) file >> par(i);

  file >> nmul_;
  P_ERROR_X3(mRe.size() >= nmul_, "Not compatible input file (NMUL) ", nmul_, ".");
  for (size_t i = 0; i < nmul_; i++)
  {
    file >> mRe(i);
    file >> mIm(i);
  }

  file >> ndim_;
  file >> nint_;
  file >> ndeg_;

  P_ERROR_X3(NDIM == ndim_, "Not compatible input file (NDIM) ", ndim_, ".");

  KNVector msh(nint_ + 1);
  double t;
  for (size_t i = 0; i < ndeg_*nint_ + 1; i++)
  {
    file >> t;
    if (i % ndeg_ == 0) msh(i / ndeg_) = t;
  }

  if ((NINT == nint_) && (NDEG == ndeg_))
  {
    colloc->setMesh(msh);
    for (size_t i = 0; i < NDIM*(nint_*ndeg_ + 1); i++) file >> sol(i);
  }
  else
  {
    KNVector in(NDIM*(nint_*ndeg_ + 1));

    for (size_t i = 0; i < NDIM*(nint_*ndeg_ + 1); i++) file >> in(i);
    colloc->importProfile(sol, in, msh, ndeg_, true);
  }
}

void KNDdePeriodicSolution::ReadNull(std::ifstream& file)
{
  double tmp;
  size_t npar_, nmul_, ndim_, nint_, ndeg_;
  file >> npar_;
  P_ERROR_X3(VarToIndex(VarEnd,NPAR) == npar_, "Not compatible input file (NPAR) ", npar_, ".");
  for (size_t i = 0; i < VarToIndex(VarEnd,NPAR); i++) file >> tmp;

  file >> nmul_;
  P_ERROR_X3(mRe.size() >= nmul_, "Not compatible input file (NMUL) ", nmul_, ".");
  for (size_t i = 0; i < nmul_; i++)
  {
    file >> tmp;
    file >> tmp;
  }

  file >> ndim_;
  file >> nint_;
  file >> ndeg_;

  P_ERROR_X3(NDIM == ndim_, "Not compatible input file (NDIM) ", ndim_, ".");

  for (size_t i = 0; i < ndeg_*nint_ + 1; i++) file >> tmp;
  for (size_t i = 0; i < NDIM*(nint_*ndeg_ + 1); i++) file >> tmp;
}

void KNDdePeriodicSolution::SwitchTFTRTan(KNVector& Re, KNVector& Im, double& alpha, const KNVector& mshint, const KNVector& mshdeg)   // starting data for tori: tangent
{
  KNVector TRe(sol.size()), TIm(sol.size());
  size_t idx = 0;
  for (size_t i=0; i<eqn.size(); ++i) if (eqn(i) == EqnTFCPLX_RE) { idx = i; break; }
  KNComplexTestFunctional* tf = static_cast< KNComplexTestFunctional* >(testFunct(idx));
  if (tf)
  {
    tf->kernel(TRe, TIm, alpha);
    colloc->exportProfile(Re, mshint, mshdeg, TRe);
    colloc->exportProfile(Im, mshint, mshdeg, TIm);
  }
  else
  {
    P_MESSAGE1("Cannot compute the initial tangent to the torus branch, because no complex test functional was defined.");
  }
}
