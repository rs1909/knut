// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //
#define _USE_MATH_DEFINES
#include <cmath>

#include "basepoint.h"
#include "basecomp.h"
#include "exprsystem.h"
#include "mat4data.h"
#include <sstream>
#include <algorithm>

// specified in the constants file
#define REFEPS    (1E-5)
#define NREFITER  (20)
#define CONTEPS   (1E-5)
#define NCONTITER (20)
#define KERNEPS   (1E-10)
#define NKERNITER (20)
#define NCONTCURVATURE (0.25*M_PI/2.0)

#define MAXEQNS 12

struct PtTab
{
  PtType   type;
  BranchSW sw;
  size_t   neqn;
  size_t   nparx;
  Eqn      eqns[MAXEQNS];
  Var      vars[MAXEQNS];
};

// not even member function public
BranchSW PtToEqnVar(KNArray1D<Eqn>& eqnr, KNArray1D<Var>& varr, PtType Pt, KNArray1D<Var> parx, size_t /*npar_*/)
{
  PtTab tab;
  const Var PANGLE = VarAngle; // IndexToVar(VarAngle,npar_);
  switch (Pt)
  {
    case SolODE:
    {
      PtTab tmp = { SolODE, NOSwitch,   1, 0,
                    { EqnODESol },
                    { VarODESol } };
      tab = tmp;
    }
    break;
    case SolAUTODE:
    {
      PtTab tmp = { SolAUTODE, NOSwitch,   2, 1,
                    { EqnODESol, EqnPhase },
                    { VarODESol, VarNone } };
      tab = tmp;
    }
    break;
    case SolSteady:
    {
      PtTab tmp = { SolSteady, NOSwitch,   1, 0,
                    { EqnSteady },
                    { VarSteady } };
      tab = tmp;
    }
    break;
    /// TIME-PERIODIC TEST-FUNCTIONAL
    case SolTF:
    {
      PtTab tmp = { SolTF, NOSwitch,   1, 0,
                    { EqnSol },
                    { VarSol } };
      tab = tmp;
    }
    break;
    case SolTFBRSW:
    {
      PtTab tmp = { SolTFBRSW, TFBRSwitch, 1, 0,
                    { EqnSol },
                    { VarSol } };
      tab = tmp;
    }
    break;
    case SolTFPDSW:
    {
      PtTab tmp = { SolTFPDSW, TFPDSwitch, 1, 0,
                    { EqnSol },
                    { VarSol } };
      tab = tmp;
    }
    break;
    case BifTFLP:
    {
      PtTab tmp = { BifTFLP, NOSwitch,   2, 1,
                    { EqnSol, EqnTFLP },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    case BifTFPD:
    {
      PtTab tmp = { BifTFPD, NOSwitch,   2, 1,
                    { EqnSol, EqnTFPD },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    case BifTFNS:
    {
      PtTab tmp = { BifTFNS, NOSwitch,   3, 1,
                    { EqnSol, EqnTFCPLX_RE,  EqnTFCPLX_IM },
                    { VarSol, PANGLE,        VarNone } };
      tab = tmp;
    }
    break;
    /// AUTONOMOUS TEST-FUNCTIONAL
    case SolTFAUT:
    {
      PtTab tmp = { SolTFAUT, NOSwitch,   2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    case SolTFAUTBRSW:
    {
      PtTab tmp = { SolTFAUTBRSW, TFBRAUTSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    case SolTFAUTPDSW:
    {
      PtTab tmp = { SolTFAUTPDSW, TFPDSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    // Hopf bifurcation switch
    case SolTFAUTHBSW:
    {
      PtTab tmp = { SolTFAUTHBSW, TFHBSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    // Spectral submanifold switch: Same as Hopf bifurcation
    case SolTFAUTSSMSW:
    {
      PtTab tmp = { SolTFAUTHBSW, TFSSMSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, VarNone } };
      tab = tmp;
    }
    break;
    case BifTFAUTLP:
    {
      PtTab tmp = { BifTFAUTLP, NOSwitch,   3, 2,
                    { EqnSol, EqnTFLPAUT, EqnPhase },
                    { VarSol, VarNone,    VarNone } };
      tab = tmp;
    }
    break;
    case BifTFAUTPD:
    {
      PtTab tmp = { BifTFAUTPD, NOSwitch,   3, 2,
                    { EqnSol, EqnTFPD, EqnPhase },
                    { VarSol, VarNone, VarNone } };
      tab = tmp;
    }
    break;
    case BifTFAUTNS:
    {
      PtTab tmp = { BifTFAUTNS, NOSwitch,   4, 2,
                    { EqnSol, EqnTFCPLX_RE,  EqnTFCPLX_IM, EqnPhase },
                    { VarSol, PANGLE,        VarNone,      VarNone } };
      tab = tmp;
    }
    break;
    /// TORUS
    case SolTor:
    {
      PtTab tmp = { SolTor, NOSwitch, 2, 1,
                    { EqnTORSol, EqnTORPhase1 },
                    { VarTORSol, VarNone } };
      tab = tmp;
    }
    break;
    case SolTorNS:
    {
      PtTab tmp = { SolTor, TFTRSwitch, 2, 1,
                    { EqnTORSol, EqnTORPhase1 },
                    { VarTORSol, VarNone } };
      tab = tmp;
    }
    break;
    case SolAUTTor:
    {
      PtTab tmp = { SolAUTTor, NOSwitch, 3, 2,
                    { EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
                    { VarTORSol, VarNone,      VarNone } };
      tab = tmp;
    }
    break;
    case SolAUTTorNS:
    {
      PtTab tmp = { SolAUTTor, TFTRSwitch, 3, 2,
                    { EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
                    { VarTORSol, VarNone,      VarNone } };
      tab = tmp;
    }
    break;
    default:
    {
      PtTab tmp = { SolUser, NOSwitch,  0, 0, { EqnNone }, { VarNone } };
      tab = tmp;
    }
    P_MESSAGE3("Invalid point type ", (int)Pt, ".");
    break;
  }
  P_ERROR_X4(tab.nparx == parx.size(), "Wrong number of additional continuation parameters (NPARX). ", tab.nparx, "!=", parx.size());
  P_ERROR_X1(tab.neqn < MAXEQNS, "Too many equations and/or parameters.");
  for (size_t i = 0; i < tab.nparx; i++) tab.vars[(tab.neqn-tab.nparx) + i] = parx(i);
  eqnr.init(tab.neqn);
  varr.init(tab.neqn);
  for (size_t i = 0; i < tab.neqn; i++)
  {
    eqnr(i) = tab.eqns[i];
    varr(i) = tab.vars[i];
  }
  return tab.sw;
}

KNAbstractPoint::KNAbstractPoint(KNAbstractContinuation* cnt, KNExprSystem& sys_,
 const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_,
 const size_t solsize, const size_t nz_jac_) :
    var(var_), eqn(eqn_), varMap(var_.size()), varMapCont(var_.size() + 1), npar(sys_.npar()),
    sol(solsize), par(VarToIndex(VarEnd,sys_.npar())),
    solNu(solsize), parNu(VarToIndex(VarEnd,sys_.npar())),
    notifyObject(cnt)
{
  dim1   = solsize;
  nz_jac = nz_jac_;
  RefEps   = REFEPS;
  RefIter  = NREFITER;
  ContEps  = CONTEPS;
  ContIter = NCONTITER;
  KernEps  = KERNEPS;
  KernIter = NKERNITER;
  ContCurvature = NCONTCURVATURE;

  par(VarToIndex(VarAngle,sys_.npar())) = 0.0;
  par(VarToIndex(VarRot,sys_.npar())) = 0.0;
  // DO NOT CALL Construct! It will be called from its children!
}

KNAbstractPoint::~KNAbstractPoint()
{
  // DO NOT CALL Destruct! It will be called from its children!
}

std::ostream& KNAbstractPoint::outStream()
{
  if (notifyObject) return notifyObject->outStream();
  else return std::cout;
}

void KNAbstractPoint::printStream()
{
  if (notifyObject) notifyObject->printStream();
}

// public
// remember that this ereases the state variables except sol, qq, par xxDot
void KNAbstractPoint::reset(const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_)

{
  KNBlockVector* xxDot_temp = nullptr;
  if (xxDot) xxDot_temp = new KNBlockVector(*xxDot);
  destruct();
  eqn.init(eqn_.size());
  eqn = eqn_;
  var.init(var_.size());
  var = var_;
  varMap.init(var_.size());
  varMapCont.init(var_.size() + 1);
  // this calls the virtual construct
  construct();
  if (xxDot_temp != nullptr)
  {
    xxDot->getV1() = xxDot_temp->getV1();
    for (size_t i = 0; i < std::min<size_t>(xxDot_temp->getV3().size(), xxDot->getV3().size()); ++i)
      xxDot->getV3()(i) = xxDot_temp->getV3()(i);
    delete xxDot_temp;
  }
}

void KNAbstractPoint::construct()
{
  P_ERROR_X1((eqn.size() != 0) && (var.size() != 0) && (eqn.size() == var.size()), "Number of equations and variables do not agree.");
  dim3 = eqn.size() - 1;

  for (size_t i = 1; i < var.size(); i++)
  {
//    std::cout << "E " << eqn(i) << ", V " << var(i) << "\n";
    varMap(i) = VarToIndex(var(i),npar);
//    std::cout << varMap(i) << " " << par.size() <<"\n";
    P_ERROR_X5(varMap(i) < par.size(), "Non-existing parameter P", varMap(i), " at position ", i, ".");
  }
  for (size_t i = 0; i < var.size(); i++) varMapCont(i) = varMap(i);

  xxDot   = new KNBlockVector(dim1, 0, dim3 + 1);
  xxDotNu = new KNBlockVector(dim1, 0, dim3 + 1);

  xx      = new KNBlockVector(dim1, 0, dim3 + 1);

  rhs     = new KNBlockVector(dim1, 0, dim3 + 1);

  jac     = new KNSparseBlockMatrix(dim1, 0, dim3 + 1, nz_jac);
}

// private
void KNAbstractPoint::destruct()
{
  delete jac;

  delete rhs;
  delete xx;

  delete xxDot;
  delete xxDotNu;
}

// this will adapt the mesh whenever it is specified
size_t KNAbstractPoint::refine(bool adapt)
{
  std::ostream& out = outStream();
  // here solNu is the previous solution
  if (adapt)
  {
    basecolloc->meshAdapt(solNu, sol, xxDotNu->getV1(), xxDot->getV1());
    sol = solNu;
    xxDot->getV1() = xxDotNu->getV1();
  }
  solNu = sol;
  parNu = par;

  xx->getV3().clear();

  if(!adapt) { out << "IT\tERR\t\tSOLnorm\t\tDIFFnorm\n"; printStream(); }

  size_t it = 0;
  double Xnorm, Dnorm;
  do
  {
    basecolloc->init(sol, par);

    if (!adapt)
    {
      jacobian(*jac, *rhs, parNu, par, solNu, sol, varMap, 0.0, false);
      jac->solve(*xx, *rhs, dim3);
      update(*xx);
    } else
    {
      jacobian(*jac, *rhs, parNu, par, solNu, sol, varMapCont, 0.0, true);
      jac->solve(*xx, *rhs, dim3+1);
      updateWithAdaptation(*xx);
    }
    // computing norms to determine convergence
    Xnorm = sqrt(basecolloc->integrate(sol, sol));
    Dnorm = sqrt(basecolloc->integrate(xx->getV1(), xx->getV1()) + (xx->getV3()) * (xx->getV3()));
    if(!adapt)
    {
      out << " " << it << "\t" << Dnorm / (1.0 + Xnorm) << "\t" << Xnorm << "\t" << Dnorm << '\n';
      printStream();
    }
  }
  while ((Dnorm / (1.0 + Xnorm) >= RefEps) && (it++ < RefIter));
  if (it >= RefIter)
  {
    out << "Warning: refinement did not converge. "
        << "CritNorm: " << Dnorm / (1.0 + Xnorm) << " SolNorm: " << Xnorm << " DiffNorm: " << Dnorm << '\n';
    printStream();
  }

  return it;
}

size_t KNAbstractPoint::tangent(bool adapt)
{
  std::ostream& out = outStream();
  double norm;

  basecolloc->init(sol, par);

  if (!adapt)
  {
    // setting up a random tangent
    xxDot->getV1().random();
    xxDot->getV3().random();
    norm = sqrt(basecolloc->integrate(xxDot->getV1(), xxDot->getV1()) + (xxDot->getV3()) * (xxDot->getV3()));
    xxDot->getV1() /= norm;
    xxDot->getV3() /= norm;
  }
  // az RHS-t feleslegesen szamolja ki && the first qq should be qq0
  jacobian(*jac, *rhs, par, par, sol, sol, varMapCont, 0.0, true);

  double diffnorm = 1.0;
  size_t it = 0;
  do
  {
    jac->multiply<false>(*rhs, *xxDot, dim3 + 1);
    rhs->getV3()(dim3) -= 1.0;
    jac->solve(*xx, *rhs);
    xxDot->getV1() -= xx->getV1();
    xxDot->getV3() -= xx->getV3();
    diffnorm = sqrt(basecolloc->integrate(xx->getV1(), xx->getV1()) + (xx->getV3()) * (xx->getV3()));
    norm = sqrt(basecolloc->integrate(xxDot->getV1(), xxDot->getV1()) + (xxDot->getV3()) * (xxDot->getV3()));
    xxDot->getV1() /= norm;
    xxDot->getV3() /= norm;
    // putting back the tangent...
    if (dim1 != 0) basecolloc->star(jac->getA31(dim3), xxDot->getV1());
    for (size_t i = 0; i < dim3 + 1; i++) jac->getA33()(dim3, i) = xxDot->getV3()(i);
  }
  while ((++it < KernIter) && (diffnorm > KernEps));
  if (diffnorm > KernEps)
  {
    out << "KNDdePeriodicSolution::Tangent: warning: No convergence in finding the singular vector. Residual = " << diffnorm << ", steps " << it << "\n";
    printStream();
  }
  if (!adapt && (xxDot->getV3()(dim3) < 0.0))
  {
    xxDot->getV1() *= -1.0;
    xxDot->getV3() *= -1.0;
  }

  return it;
}

size_t KNAbstractPoint::nextStep(double ds, double& angle, const IterateTangent jacstep)
{
  std::ostream& out = outStream();
  double Xnorm, Dnorm, Rnorm, Tnorm;

  parNu = par;
  for (size_t i = 0; i < solNu.size(); i++)  solNu(i)           = sol(i)           + ds * xxDot->getV1()(i);
  for (size_t i = 1; i < varMapCont.size(); i++) parNu(varMapCont(i)) = par(varMapCont(i)) + ds * xxDot->getV3()(i - 1);
  xxDotNu->getV1() = xxDot->getV1();
  xxDotNu->getV3() = xxDot->getV3();

#ifdef DEBUG
    if (ds != 0.0) out << "Iteration: Tangent (T), updated tangent (U) and actual difference (A)\n";
    else out << "Refine (Cont): Tangent (T), updated tangent (U) and actual difference (A)\n";
    const double T1norm = sqrt(basecolloc->integrate(xxDot->getV1(), xxDot->getV1()));

    out << "TX = " << T1norm;
    for (size_t i = 1; i < varMapCont.size(); i++) out << " T" << varMapCont(i) << " = " << xxDot->getV3()(i - 1);
    out << " TC = " << xxDot->getV3()(dim3) << '\n';
    printStream();
#endif /*DEBUG*/
  size_t  it = 0;
  bool conv;
  do
  {
    basecolloc->init(solNu, parNu);

    jacobian(*jac, *rhs, par, parNu, sol, solNu, varMapCont, 0.0, true);

    jac->solve(*xx, *rhs);

    updateWithCp(*xx);

    Rnorm = sqrt(basecolloc->integrate(rhs->getV1(), rhs->getV1()) + (rhs->getV3()) * (rhs->getV3()));
    Xnorm = sqrt(basecolloc->integrate(solNu, solNu));
    Dnorm = sqrt(basecolloc->integrate(xx->getV1(), xx->getV1()) + (xx->getV3()) * (xx->getV3()));
    conv = (Dnorm / (1.0 + Xnorm) < ContEps) && (Rnorm / (1.0 + Xnorm) < ContEps);

#ifdef DEBUG
    const double A1norm = sqrt(basecolloc->integrate(xx->getV1(), xx->getV1()));
    const double dsdiv = (ds == 0.0) ? 1.0 : ds;
    out << "AX = " << A1norm/dsdiv;
    for (size_t i = 1; i < varMapCont.size(); i++) out << " A" << varMapCont(i) << " = " << xx->getV3()(i - 1)/dsdiv;
    out << '\n';
    printStream();
#endif /*DEBUG*/
    // updating the tangent if converged or jacstep == false
    if ((jacstep == IterateTangent::yes) || conv)
    {
      // tangent
      jac->multiply<false>(*rhs, *xxDotNu, dim3 + 1);
      rhs->getV3()(dim3) -= 1.0;
      jac->solve(*xx, *rhs);
      xxDotNu->getV1() -= xx->getV1();
      xxDotNu->getV3() -= xx->getV3();
      Tnorm = sqrt(basecolloc->integrate(xxDotNu->getV1(), xxDotNu->getV1()) + (xxDotNu->getV3()) * (xxDotNu->getV3()));
      xxDotNu->getV1() /= Tnorm;
      xxDotNu->getV3() /= Tnorm;
    }
    // end updating tangent
    if (ds != 0.0) out << '.';
    else out << 'o';
    printStream();
  }
  while (!conv /*&& (Dnorm/(1.0+Xnorm) < 1.0)*/ && (++it < ContIter));
  if (conv)
  {
#ifdef DEBUG
    // xx strores the actual difference
    xx->getV1() = solNu;
    xx->getV1() -= sol;
    for (size_t i = 1; i < varMapCont.size(); i++) xx->getV3()(i - 1) = parNu(varMapCont(i)) - par(varMapCont(i));

    // checking the tangent and the secant
    if (ds != 0.0) out << "Finish: tangent (T), updated tangent (U) and actual difference (A)\n";
    else out << "Refine finish: tangent (T), updated tangent (U) and actual difference (A)\n";
    const double T1norm = sqrt(basecolloc->integrate(xxDot->getV1(), xxDot->getV1()));
    const double U1norm = sqrt(basecolloc->integrate(xxDotNu->getV1(), xxDotNu->getV1()));
    const double A1norm = sqrt(basecolloc->integrate(xx->getV1(), xx->getV1()));

    out << "TX = " << T1norm;
    for (size_t i = 1; i < varMapCont.size(); i++) out << " T" << varMapCont(i) << " = " << xxDot->getV3()(i - 1);
    out << '\n';
    out << "UX = " << U1norm;
    for (size_t i = 1; i < varMapCont.size(); i++) out << " U" << varMapCont(i) << " = " << xxDotNu->getV3()(i - 1);
    out << '\n';
    const double dsdiv = (ds == 0.0) ? 1.0 : ds;
    out << "AX = " << A1norm/dsdiv;
    for (size_t i = 1; i < varMapCont.size(); i++) out << " A" << varMapCont(i) << " = " << xx->getV3()(i - 1)/dsdiv;
    out << '\n';
    printStream();
    /// END OF CHECKING
#endif

    // find out how far was from the original solution
    xx->getV1() = solNu;
    xx->getV1() -= sol;
    xx->getV1() -= ds*xxDot->getV1();
    double XnormSQ = basecolloc->integrate(xx->getV1(), xx->getV1());
    for (size_t i = 1; i < varMapCont.size(); i++)
      XnormSQ += pow(par(varMapCont(i)) + ds * xxDot->getV3()(i - 1) - parNu(varMapCont(i)),2);
    // Sometimes ds == 0
    if (fabs(ds) > 0.0) angle = atan(sqrt(XnormSQ)/fabs(ds))/(M_PI/2.0);
    else angle = 0;
    if ((angle < ContCurvature) || (jacstep == IterateTangent::no)) // dont check curvature at the begining
    {
      // copying back the solution
      sol = solNu;
      par = parNu;
      xxDot->getV1() = xxDotNu->getV1();
      xxDot->getV3() = xxDotNu->getV3();
      // test functional
      postProcess();
    }
  }
//  else
//  {
//    std::cout << "\n\n\n ------------------- NO CONVERGENCE -------------------\n\n\n\n";
    // P_MESSAGE1("");
//  }

  out << '|';
  printStream();
  return it;
}

#define NDIM persolcolloc->nDim()
#define NTAU persolcolloc->nTau()
#define NPAR persolcolloc->nPar()
#define NINT persolcolloc->nInt()
#define NDEG persolcolloc->nDeg()

// private
void KNAbstractPeriodicSolution::FillSol(KNExprSystem& sys_)
{
  persolcolloc->resetProfile(sys_, sol);
  sys_.stpar(par);
  par(VarToIndex(VarPeriod,NPAR)) = 1.0;
}

// __not__ specified in the constants file
#define MIN_NSIMAG 1e-4

/// It only computes the critical characteristic multiplier and refines the solution

// id : the index of the characteristic multiplier
void KNAbstractPeriodicSolution::findAngle(const size_t id)
{
  std::ostream& out = outStream();

  Stability(true); // re-initialize the collocation just in case when the period changes
  std::vector<std::pair<double,size_t>> aind;
  aind.resize(mRe.size());
  for (size_t i = 0; i < mRe.size(); i++)
  {
    if (fabs(mIm(i)) > MIN_NSIMAG) aind[i] = std::pair<double,size_t>(fabs(sqrt(mRe(i)*mRe(i) + mIm(i)*mIm(i)) - 1.0), i);
    else aind[i] = std::pair<double,size_t>(DBL_MAX, i);
  }
  std::sort(aind.begin(), aind.end(),
     [](std::pair<double,size_t> a, std::pair<double,size_t> b) { return a.first < b.first; });
//   for(auto& it : aind)
//   {
//     std::cout << it.first << ", " << it.second << "\n";
//   }
//  double dmin = 10.0;
//  size_t imin = 0;
//  for (size_t i = 0; i < mRe.size(); i++)
//  {
//    if (dmin > fabs(sqrt(mRe(i)*mRe(i) + mIm(i)*mIm(i)) - 1.0))
//    {
//      if (fabs(mIm(i)) > MIN_NSIMAG)
//      {
//        dmin = fabs(sqrt(mRe(i) * mRe(i) + mIm(i) * mIm(i)) - 1.0);
//        imin = i;
//      }
//    }
//  }
  auto zRe = mRe(aind[id].second);
  auto zIm = fabs(mIm(aind[id].second));
  auto nrm = sqrt(zRe * zRe + zIm * zIm);

  if (zRe > 0)
  {
    par(VarToIndex(VarAngle,NPAR)) = atan(zIm / zRe);
  }
  else
  {
    par(VarToIndex(VarAngle,NPAR)) = M_PI + atan(zIm/ zRe);
  }
  out << id << "-th closest root " << "Z = " << zRe << " + I*" << zIm << "\n"
         "Z = " << nrm << " * " << "EXP( " << par(VarToIndex(VarAngle,NPAR)) / (2*M_PI) << " * I*2Pi )\n";
  printStream();
}

void KNAbstractPeriodicSolution::BinaryWrite(KNAbstractData& data, BifType bif, size_t n)
{
//  std::cout << "DAT " << &data << "\n";
  data.lockWrite();
  data.setNTrivMul(0, nTrivMulLP);
  data.setNTrivMul(1, nTrivMulPD);
  data.setNTrivMul(2, nTrivMulNS);
  //
  data.setMagic(n, bif);
  data.setPar(n, par);
  data.setMul(n, mRe, mIm);
  data.setElem(n, persolcolloc->getElem());
  data.setMesh(n, persolcolloc->getMesh());
  data.setProfile(n, sol);
  data.unlock();
}

void KNAbstractPeriodicSolution::BinaryRead(const KNAbstractData& data, size_t n)
{
  std::ostream& out = outStream();
  data.lockRead();
  KNVector msh(data.getNInt() + 1);
  P_ERROR_X1(data.getNPar() == VarToIndex(VarEnd,NPAR), "Wrong number of parameters in the input MAT file.");
  data.getPar(n, par);
  data.getMul(n, mRe, mIm);
  data.getMesh(n, msh);
  P_ERROR_X1(data.getNDim() == NDIM, "Wrong number of dimensions in the input MAT file.");
  if (data.getNInt() == NINT && data.getNDeg() == NDEG)
  {
    persolcolloc->setMesh(msh);
    data.getProfile(n, sol);
  }
  else
  {
    KNVector tmp(data.getNDim()*(data.getNDeg()*data.getNInt() + 1));
    data.getProfile(n, tmp);
    const double amp = Amplitude(tmp, data.getNDim(), data.getNDeg(), data.getNInt());
    if (amp < 1e-6)
    {
      out << "Warning: importing without mesh adaptation.\n";
      printStream();
      persolcolloc->importProfile(sol, tmp, msh, data.getNDeg(), false);
    } else
    {
      persolcolloc->importProfile(sol, tmp, msh, data.getNDeg(), true);
    }
  }
  data.unlock();
}
