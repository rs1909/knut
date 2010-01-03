// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "basepoint.h"
#include "system.h"
#include "mat4data.h"

// specified in the constants file
#define REFEPS    (1E-5)
#define NREFITER  (20)
#define CONTEPS   (1E-5)
#define NCONTITER (20)
#define KERNEPS   (1E-10)
#define NKERNITER (20)

struct PtTab
{
  PtType   type;
  BranchSW sw;
  int      neqn;
  int      nparx;
  Eqn      eqns[5];
  Var      vars[5];
};

// not even member function public
BranchSW PtToEqnVar(KNArray1D<Eqn>& eqnr, KNArray1D<Var>& varr, PtType Pt, KNArray1D<Var> parx, int npar_)
{
  PtTab tab;
  const Var PANGLE = (Var)(VarPAR0 + npar_ + ParAngle);
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
                    { VarODESol, parx(0) } };
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
                    { VarSol, parx(0) } };
      tab = tmp;
    }
    break;
    case BifTFPD:
    {
      PtTab tmp = { BifTFPD, NOSwitch,   2, 1,
                    { EqnSol, EqnTFPD },
                    { VarSol, parx(0) } };
      tab = tmp;
    }
    break;
    case BifTFNS:
    {
      PtTab tmp = { BifTFNS, NOSwitch,   3, 1,
                    { EqnSol, EqnTFCPLX_RE,  EqnTFCPLX_IM },
                    { VarSol, PANGLE,        parx(0) } };
      tab = tmp;
    }
    break;
    /// AUTONOMOUS TEST-FUNCTIONAL
    case SolTFAUT:
    {
      PtTab tmp = { SolTFAUT, NOSwitch,   2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, parx(0) } };
      tab = tmp;
    }
    break;
    case SolTFAUTBRSW:
    {
      PtTab tmp = { SolTFAUTBRSW, TFBRSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, parx(0) } };
      tab = tmp;
    }
    break;
    case SolTFAUTPDSW:
    {
      PtTab tmp = { SolTFAUTPDSW, TFPDSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, parx(0) } };
      tab = tmp;
    }
    break;
    case SolTFAUTHBSW:
    {
      PtTab tmp = { SolTFAUTHBSW, TFHBSwitch, 2, 1,
                    { EqnSol, EqnPhase },
                    { VarSol, parx(0) } };
      tab = tmp;
    }
    break;
    case BifTFAUTLP:
    {
      PtTab tmp = { BifTFAUTLP, NOSwitch,   3, 2,
                    { EqnSol, EqnTFLPAUT, EqnPhase },
                    { VarSol, parx(0),    parx(1) } };
      tab = tmp;
    }
    break;
    case BifTFAUTPD:
    {
      PtTab tmp = { BifTFAUTPD, NOSwitch,   3, 2,
                    { EqnSol, EqnTFPD, EqnPhase },
                    { VarSol, parx(0), parx(1) } };
      tab = tmp;
    }
    break;
    case BifTFAUTNS:
    {
      PtTab tmp = { BifTFAUTNS, NOSwitch,   4, 2,
                    { EqnSol, EqnTFCPLX_RE,  EqnTFCPLX_IM, EqnPhase },
                    { VarSol, PANGLE,        parx(0),      parx(1) } };
      tab = tmp;
    }
    break;
    /// TORUS
    case SolTor:
    {
      PtTab tmp = { SolTor, NOSwitch, 2, 1,
                    { EqnTORSol, EqnTORPhase1 },
                    { VarTORSol, parx(0) } };
      tab = tmp;
    }
    break;
    case SolTorNS:
    {
      PtTab tmp = { SolTor, TFTRSwitch, 2, 1,
                    { EqnTORSol, EqnTORPhase1 },
                    { VarTORSol, parx(0) } };
      tab = tmp;
    }
    break;
    case SolAUTTor:
    {
      PtTab tmp = { SolAUTTor, NOSwitch, 3, 2,
                    { EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
                    { VarTORSol, parx(0),      parx(1) } };
      tab = tmp;
    }
    break;
    case SolAUTTorNS:
    {
      PtTab tmp = { SolAUTTor, TFTRSwitch, 3, 2,
                    { EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
                    { VarTORSol, parx(0),      parx(1) } };
      tab = tmp;
    }
    break;
    default:
    {
      PtTab tmp = { SolUser, NOSwitch,  0, 0, { EqnNone }, { VarNone } };
      tab = tmp;
    }
    P_MESSAGE3("Invalid point type ", Pt, ".");
    break;
  }
  P_ERROR_X4(tab.nparx == parx.size(), "Wrong number of additional continuation parameters (NPARX). ", tab.nparx, "!=", parx.size());
  eqnr.init(tab.neqn);
  varr.init(tab.neqn);
  for (int i = 0; i < tab.neqn; i++)
  {
    eqnr(i) = tab.eqns[i];
    varr(i) = tab.vars[i];
  }
  return tab.sw;
}

KNAbstractPoint::KNAbstractPoint(KNSystem& sys_, const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_, 
          const int solsize, const int nz_jac_) :
    var(var_), eqn(eqn_), varMap(var_.size()), varMapCont(var_.size() + 1),
    sol(solsize), par(sys_.npar() + ParEnd),
    solNu(solsize), parNu(sys_.npar() + ParEnd)
{
  dim1   = solsize;
  nz_jac = nz_jac_;
  RefEps   = REFEPS;
  RefIter  = NREFITER;
  ContEps  = CONTEPS;
  ContIter = NCONTITER;
  KernEps  = KERNEPS;
  KernIter = NKERNITER;

  par(sys_.npar() + ParAngle) = 0.0;
  par(sys_.npar() + ParRot) = 0.0;
  // DO NOT CALL Construct! It will be called from its children!
}

KNAbstractPoint::~KNAbstractPoint()
{
  // DO NOT CALL Destruct! It will be called from its children!
}

// public
// remember that this ereases the state variables except sol, qq, par xxDot
void KNAbstractPoint::reset(const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_)

{
  KNBlockVector* xxDot_temp = 0;
  if (xxDot) xxDot_temp = new KNBlockVector(*xxDot);
  destruct();
  eqn.init(eqn_.size());
  eqn = eqn_;
  var.init(var_.size());
  var = var_;
  varMap.init(var_.size());
  varMapCont.init(var_.size() + 1);
  construct();
  xxDot->getV1() = xxDot_temp->getV1();
  for (int i = 0; i < std::min<int>(xxDot_temp->getV3().size(), xxDot->getV3().size()); ++i)
    xxDot->getV3()(i) = xxDot_temp->getV3()(i);
  delete xxDot_temp;
}

void KNAbstractPoint::construct()
{
  P_ERROR_X1((eqn.size() != 0) && (var.size() != 0) && (eqn.size() == var.size()), "Number of equations and variables do not agree.");
  dim3 = eqn.size() - 1;

  for (int i = 1; i < var.size(); i++)
  {
    P_ERROR_X5((var(i) - VarPAR0 >= 0) && (var(i) - VarPAR0 < par.size()), "Non-existing parameter P", var(i) - VarPAR0, " at position ", i, ".");
    varMap(i) = var(i) - VarPAR0;
  }
  for (int i = 0; i < var.size(); i++) varMapCont(i) = varMap(i);

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
int KNAbstractPoint::refine(std::ostream& out, bool adapt)
{
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

  if(!adapt) out << "IT\tERR\t\tSOLnorm\t\tDIFFnorm\n";

  int it = 0;
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
      out.flush();
    }
  }
  while ((Dnorm / (1.0 + Xnorm) >= RefEps) && (it++ < RefIter));
  if (it >= RefIter) std::cerr << "Warning: refinement did not converge. "
    << "CritNorm: " << Dnorm / (1.0 + Xnorm) << " SolNorm: " << Xnorm << " DiffNorm: " << Dnorm << '\n';

  return it;
}

int KNAbstractPoint::tangent(bool adapt)
{
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
  int it = 0;
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
    for (int i = 0; i < dim3 + 1; i++) jac->getA33()(dim3, i) = xxDot->getV3()(i);
  }
  while ((++it < KernIter) && (diffnorm > KernEps));
  if (diffnorm > KernEps) std::cerr << "KNDdePeriodicSolution::Tangent: warning: No convergence in finding the singular vector. Residual = " << diffnorm << ", steps " << it << "\n";
  if (!adapt && (xxDot->getV3()(dim3) < 0.0))
  {
    xxDot->getV1() *= -1.0;
    xxDot->getV3() *= -1.0;
  }

  return it;
}

int KNAbstractPoint::nextStep(double ds, bool jacstep)
{
  double Xnorm, Dnorm, Rnorm, Tnorm;

  parNu = par;
  for (int i = 0; i < solNu.size(); i++)  solNu(i)           = sol(i)           + ds * xxDot->getV1()(i);
  for (int i = 1; i < varMapCont.size(); i++) parNu(varMapCont(i)) = par(varMapCont(i)) + ds * xxDot->getV3()(i - 1);
  xxDotNu->getV1() = xxDot->getV1();
  xxDotNu->getV3() = xxDot->getV3();

  int  it = 0;
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
    conv = (Dnorm / (1.0 + Xnorm) >= ContEps) || (Rnorm / (1.0 + Xnorm) >= ContEps);

#ifdef DEBUG
    std::cout << "Dnorm: " << Dnorm << " Rnorm: " << Rnorm << " Xnorm: " << Xnorm << "\n";
    std::cout.flush();
#endif /*DEBUG*/
    // updating the tangent
    if (!jacstep)
    {
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
  }
  while (conv /*&& (Dnorm/(1.0+Xnorm) < 1.0)*/ && (++it < ContIter));
  if (!conv)
  {
#ifdef DEBUG
    /// checking the tangent and the secant
    double Pnorm = sqrt(xxDotNu->getV3()(dim3) * xxDotNu->getV3()(dim3));
    double Xnorm = sqrt(basecolloc->integrate(xxDotNu->getV1(), xxDotNu->getV1())), Onorm = sqrt((xxDotNu->getV3()) * (xxDotNu->getV3()));
    std::cout << "Cnorm: " << Tnorm << "\nDot Pnorm: " << Pnorm << " Xnorm: " << Xnorm << " Onorm: " << Onorm;
    for (int i = 1; i < varMap.size(); i++) std::cout << " O" << varMap(i) << ": " << xxDotNu->getV3()(i - 1);
    std::cout << '\n';

    xx->getV1() = solNu;
    xx->getV1() -= sol;
    for (int i = 1; i < varMapCont.size(); i++) xx->getV3()(i - 1) = parNu(varMapCont(i)) - par(varMapCont(i));

    Pnorm = sqrt(xx->getV3()(dim3) * xx->getV3()(dim3)) / ds;
    Xnorm = sqrt(basecolloc->integrate(xx->getV1(), xx->getV1())) / ds;
    Onorm = 0;
    for (int i = 0; i < dim3 + 1; i++) Onorm += (xx->getV3()(i)) * (xx->getV3()(i));
    Onorm = sqrt(Onorm) / ds;
    std::cout << "Dif Pnorm: " << Pnorm << " Xnorm: " << Xnorm << " Onorm: " << Onorm;
    for (int i = 1; i < varMap.size(); i++) std::cout << " O" << varMap(i) << ": " << xx->getV3()(i - 1) / ds;
    std::cout << '\n';
    /// END OF CHECKING
#endif

    // copying back the solution
    sol = solNu;
    par = parNu;
    xxDot->getV1() = xxDotNu->getV1();
    xxDot->getV3() = xxDotNu->getV3();
  }
//  else
//  {
//    std::cout << "\n\n\n ------------------- NO CONVERGENCE -------------------\n\n\n\n";
    // P_MESSAGE1("");
//  }

  return it;
}

#define NDIM persolcolloc->nDim()
#define NTAU persolcolloc->nTau()
#define NPAR persolcolloc->nPar()
#define NINT persolcolloc->nInt()
#define NDEG persolcolloc->nDeg()

// private
void KNAbstractPeriodicSolution::FillSol(KNSystem& sys_)
{
  KNVector fx(persolcolloc->nDim());

  sys_.stpar(par);
  par(NPAR+ParPeriod) = 1.0;

  for (int i = 0; i < persolcolloc->nInt(); i++)
  {
    for (int d = 0; d <  persolcolloc->nDeg(); d++)
    {
      sys_.stsol(fx, persolcolloc->Profile(i, d));
      for (int j = 0; j < persolcolloc->nDim(); j++)
      {
        sol(NDIM*(i*NDEG + d) + j) = fx(j);
      }
    }
  }
  for (int j = 0; j < persolcolloc->nDim(); j++)
  {
    sol(NDIM*NDEG*NINT + j) = sol(j);
  }
}

// __not__ specified in the constants file
#define MIN_NSIMAG 1e-4

/// It only computes the critical characteristic multiplier and refines the solution

int KNAbstractPeriodicSolution::StartTF(bool findangle, std::ostream& out)
{
  if (findangle)
  {
    Stability();
    double dmin = 10.0;
    int imin = 0;
    for (int i = 0; i < mRe.size(); i++)
    {
      if (dmin > fabs(sqrt(mRe(i)*mRe(i) + mIm(i)*mIm(i)) - 1.0))
      {
        if (fabs(mIm(i)) > MIN_NSIMAG)
        {
          dmin = fabs(sqrt(mRe(i) * mRe(i) + mIm(i) * mIm(i)) - 1.0);
          imin = i;
        }
      }
    }
    double zRe = mRe(imin);
    double zIm = fabs(mIm(imin));
    double nrm = sqrt(zRe * zRe + zIm * zIm);

    if (zRe > 0)
    {
      par(NPAR + ParAngle) = atan(zIm / zRe);
    }
    else
    {
      par(NPAR + ParAngle) = atan(fabs(zRe / zIm)) + M_PI / 2.0;
    }
    out << "    Z = " << zRe << " + I*" << zIm << "\n" 
           "    Z = " << nrm << " * " << "EXP( " << par(NPAR + ParAngle) / (2*M_PI) << " * I*2Pi )\n";
  }
  return refine(out);
}

void KNAbstractPeriodicSolution::BinaryWrite(KNDataFile& data, int n)
{
  data.lock();
  data.setNTrivMul(0, nTrivMulLP);
  data.setNTrivMul(1, nTrivMulPD);
  data.setNTrivMul(2, nTrivMulNS);
  //
  data.setPar(n, par);
  data.setMul(n, mRe, mIm);
  data.setElem(n, persolcolloc->getElem());
  data.setMesh(n, persolcolloc->getMesh());
  data.setProfile(n, sol);
  data.unlock();
}

void KNAbstractPeriodicSolution::BinaryRead(KNDataFile& data, int n)
{
  data.lock();
  KNVector msh(data.getNInt() + 1);
  P_ERROR_X1(data.getNPar() == (NPAR + ParEnd), "Wrong number of parameters in the input MAT file.");
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
      std::cerr << "Warning: importing without mesh adaptation.\n";
      persolcolloc->importProfile(sol, tmp, msh, data.getNDeg(), false);
    } else
    {
      persolcolloc->importProfile(sol, tmp, msh, data.getNDeg(), true);
    }
  }
  data.unlock();
}
