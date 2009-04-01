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
#include "system.h"
#include "matrix.h"
#include "spmatrix.h"
#include "plot.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

// __not__ specified in the constants file
#define MIN_NSIMAG 1e-4

#define NDIM colloc->Ndim()
#define NTAU colloc->Ntau()
#define NPAR colloc->Npar()
#define NINT colloc->Nint()
#define NDEG colloc->Ndeg()

// private
void Point::FillSol(System& sys_)
{
  Vector fx(colloc->Ndim());

  sys_.stpar(par);
  par(NPAR+ParPeriod) = 1.0;

  for (int i = 0; i < colloc->Nint(); i++)
  {
    for (int d = 0; d <  colloc->Ndeg(); d++)
    {
      sys_.stsol(fx, colloc->Profile(i, d));
      for (int j = 0; j < colloc->Ndim(); j++)
      {
        sol(NDIM*(i*NDEG + d) + j) = fx(j);
      }
    }
  }
  for (int j = 0; j < colloc->Ndim(); j++)
  {
    sol(NDIM*NDEG*NINT + j) = sol(j);
  }
}


Point::Point(System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg, int nmul, int nmat) : BasePoint(sys_, eqn_, var_, sys_.ndim()*(ndeg*nint + 1), sys_.ndim()*(ndeg*nint + 1)*sys_.ntau()*sys_.ndim()*(ndeg+1)),
    mRe(nmul), mIm(nmul),
    rotRe(2), rotIm(2),                  // !!!only for Hartmut's equation!!!!
//     colloc(sys_, nint, ndeg, nmat),
    jacStab(nmat, sys_.ndim()*(ndeg*nint + 1), sys_.ndim()*(ndeg*nint + 1)*sys_.ntau()*sys_.ndim()*(ndeg + 1))
{
  colloc = new NColloc(sys_, nint, ndeg, nmat);
  basecolloc = colloc;
  rotRe(0) = 0;
  rotIm(0) = 1;
  rotRe(1) = 3;
  rotIm(1) = 4;

  Construct();
  FillSol(sys_);
  par(NPAR + ParAngle) = 0.0;
  par(NPAR + ParRot) = 0.0;
}

Point::~Point()
{
  Destruct();
  delete colloc;
}

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
BranchSW PtToEqnVar(Array1D<Eqn>& eqnr, Array1D<Var>& varr, PtType Pt, Array1D<Var> parx, int npar_)
{
  PtTab tab;
  const Var PANGLE = (Var)(VarPAR0 + npar_ + ParAngle);
  switch (Pt)
  {
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
                    { EqnSol, EqnPhase,      EqnTFLPAUT },
                    { VarSol, parx(0),       parx(1) } };
      tab = tmp;
    }
    break;
    case BifTFAUTPD:
    {
      PtTab tmp = { BifTFAUTPD, NOSwitch,   3, 2,
                    { EqnSol, EqnPhase,      EqnTFPD },
                    { VarSol, parx(0),       parx(1) } };
      tab = tmp;
    }
    break;
    case BifTFAUTNS:
    {
      PtTab tmp = { BifTFAUTNS, NOSwitch,   4, 2,
                    { EqnSol, EqnPhase,      EqnTFCPLX_RE,  EqnTFCPLX_IM },
                    { VarSol, PANGLE,        parx(0),       parx(1) } };
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
      PtTab tmp = { SolAUTTor, NOSwitch, 2, 2,
                    { EqnTORSol, EqnTORPhase0, EqnTORPhase1 },
                    { VarTORSol, parx(0),      parx(1) } };
      tab = tmp;
    }
    break;
    case SolAUTTorNS:
    {
      PtTab tmp = { SolAUTTor, TFTRSwitch, 2, 2,
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
  P_ERROR_X1(tab.nparx == parx.Size(), "Wrong number of additional continuation parameters (NPARX).");
  eqnr.Init(tab.neqn);
  varr.Init(tab.neqn);
  for (int i = 0; i < tab.neqn; i++)
  {
    eqnr(i) = tab.eqns[i];
    varr(i) = tab.vars[i];
  }
  return tab.sw;
}

// private
// What does it construct?
// charMat
// xxDot, xx, rhs, jac
void Point::Construct()
{
  P_ERROR_X1((eqn(0) == EqnSol) && (var(0) == VarSol), "The first equation must be the boundary value problem of the periodic solution.");

  BasePoint::Construct();

  // a) setting the test functionals b) determining the number of trivial multipliers
  testFunct.Init(eqn.Size());
  nTrivMulLP = 0;
  nTrivMulPD = 0;
  nTrivMulNS = 0;
  for (int i = 1; i < eqn.Size(); i++)
  {
    switch (eqn(i))
    {
      case EqnTFLP:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new TestFunct(*colloc, 1.0);
        ++nTrivMulLP;
        break;
      case EqnTFPD:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new TestFunct(*colloc, -1.0);
        ++nTrivMulPD;
        break;
      case EqnTFLPAUT:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new TestFunctLPAUT(*colloc, 1.0);
        ++nTrivMulLP;
        break;
      case EqnTFLPAUTROT:
        P_ERROR_X1(testFunct(i) == 0, "The test functional already exist.");
        testFunct(i) = new TestFunctLPAUTROT(*colloc, rotRe, rotIm, 1.0);
        ++nTrivMulLP;
        break;
      case EqnTFCPLX_RE:
        P_ERROR_X1(eqn(i + 1) == EqnTFCPLX_IM,
                   "The real and imaginary parts of the complex test functional are not paired.");
        P_ERROR(testFunct(i) == 0);
        testFunct(i) = new TestFunctCPLX(*colloc);
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
void Point::Destruct()
{
  for (int i=0; i<testFunct.Size(); ++i) delete testFunct(i);
  BasePoint::Destruct();
}

// **************************************************************************************************************** //
//
//                           THE JACOBIAN
//
// **************************************************************************************************************** //

void Point::Jacobian(
  HyperMatrix& AA, HyperVector& RHS,                      // output
  Vector& parPrev, Vector& par,                           // parameters
  Vector& solPrev, Vector& sol,                           // solution
  Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
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
  }

// ---------------------------------------------------------------------------------------------------------------- //
//
//                           Other equations
//
// ---------------------------------------------------------------------------------------------------------------- //

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
      case EqnPhaseRot:
        colloc->PhaseRotStar(AA.getA31(i - 1), solPrev, rotRe, rotIm); // this should be the previous solution!!!
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
      case EqnTFLP:
      case EqnTFPD:
      case EqnTFLPAUT:
      case EqnTFLPAUTROT:
        RHS.getV3()(i - 1) = testFunct(i)->Funct(*colloc, par, sol);
        testFunct(i)->Funct_x(AA.getA31(i - 1), *colloc, par, sol);
        for (int j = 1; j < varMap.Size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            AA.getA33()(i - 1, j - 1) = testFunct(i)->Funct_p(*colloc, par, sol, varMap(j));
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
        break;
      case EqnTFCPLX_RE:
        P_ERROR_X3(eqn(i + 1) == EqnTFCPLX_IM,
                   "The real and imaginary parts of the complex test functional are not paired at position ", varMap(i), ".");
        testFunct(i)->Funct(RHS.getV3()(i - 1), RHS.getV3()(i), *colloc, par, sol,
                         cos(par(NPAR + ParAngle)), sin(par(NPAR + ParAngle)));
        testFunct(i)->Funct_x(AA.getA31(i - 1), AA.getA31(i), *colloc, par, sol);
        for (int j = 1; j < varMap.Size(); j++)
        {
          if (varMap(j) < NPAR)
          {
            testFunct(i)->Funct_p(AA.getA33()(i - 1, j - 1), AA.getA33()(i, j - 1), *colloc, par, sol, varMap(j));
          }
          else if (varMap(j) - NPAR == ParAngle)
          {
            testFunct(i)->Funct_z(AA.getA33()(i - 1, j - 1), AA.getA33()(i, j - 1), *colloc, par, sol);
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

// **************************************************************************************************************** //
//
//                          !END! THE JACOBIAN !END!
//
// **************************************************************************************************************** //

/// --------------------------------------------------------------
/// Starting bifurcation continuation using TEST FUNCTIONS
/// --------------------------------------------------------------

/// It only computes the critical characteristic multiplier and refines the solution

int Point::StartTF(Eqn FN, std::ostream& out)
{
  if ((FN == EqnTFCPLX_RE) || (FN == EqnTFCPLX_IM))
  {
    Stability();
    double dmin = 10.0;
    int imin = 0;
    for (int i = 0; i < mRe.Size(); i++)
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
  return Refine(out);
}

void Point::SwitchTFHB(double ds, std::ostream& out)
{
  int idx = -1;
  for (int i=0; i<eqn.Size(); ++i) if (eqn(i) == EqnTFCPLX_RE) { idx = i; break; }
  P_ERROR_X1 (idx != -1, "No test functional was selected for Hopf bifurcation switch!" );
  TestFunctCPLX *tf = static_cast<TestFunctCPLX*>(testFunct(idx));
  Vector QRE(NDIM), QIM(NDIM);

  par(0) = tf->SwitchHB(QRE, QIM, *colloc, par);
  out << "    T = " << par(0) << ", arg(Z) = " << par(NPAR + ParAngle) / (2*M_PI) << " * 2Pi\n";

#ifdef DEBUG
  std::cout << "Printing: neweigenvec\n";
  std::ofstream file("neweigenvec");
  file << std::scientific;
  file.precision(12);
#endif

  for (int i = 0; i < NINT; i++)
  {
    for (int j = 0; j < NDEG + 1; j++)
    {
      const double t = colloc->Profile(i, j);
      for (int p = 0; p < NDIM; p++)
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
  const double norm = sqrt(colloc->Integrate(xxDot->getV1(), xxDot->getV1()));
  xxDot->getV1() /= norm;
  xxDot->getV3().Clear();
  Vector eql(NDIM);
  for (int p = 0; p < NDIM; ++p) eql(p) = sol(p);
  for (int i = 0; i < NDEG*NINT + 1; ++i)
  {
    for (int p = 0; p < NDIM; ++p)
    {
      sol(p + i*NDIM) = eql(p) + ds * xxDot->getV1()(p + i * NDIM);
    }
  }
}

/// Switching with the test functionals!!!

void Point::SwitchTFLP(BranchSW type, double ds)
{
  xxDot->getV1().Clear();
  xxDot->getV3().Clear();

  baseTestFunct* tf = 0;
  switch (type)
  {
    case TFBRSwitch:
      tf = static_cast<baseTestFunct*>(new TestFunct(*colloc, 1.0));
      break;
    case TFBRAUTSwitch:
      tf = static_cast<baseTestFunct*>(new TestFunctLPAUT(*colloc, 1.0));
      break;
    case TFBRAUTROTSwitch:
      tf = static_cast<baseTestFunct*>(new TestFunctLPAUTROT(*colloc, rotRe, rotIm, 1.0));
      break;
    default:
      return;
      break;
  }
  tf->setKernelTolerance(KernEps, KernIter);
  tf->Funct(*colloc, par, sol);
  tf->Switch(xxDot->getV1());
  delete tf;
  double norm = sqrt(colloc->Integrate(xxDot->getV1(), xxDot->getV1()));
  xxDot->getV1() /= norm;
  xxDot->getV3().Clear();

  sol += ds * xxDot->getV1();
}

void Point::SwitchTFPD(double ds)
{
  Vector tan(xxDot->getV1().Size());
  TestFunct* tf = new TestFunct(*colloc, -1.0);
  tf->setKernelTolerance(KernEps, KernIter);
  tf->Funct(*colloc, par, sol);
  tf->Switch(tan);
  delete tf;

  // setting the period two double
  par(0) *= 2.0;
  par(NPAR+ParPeriod) *= 2.0;

  solNu = sol;
  colloc->pdMeshConvert(sol, xxDot->getV1(), solNu, tan);

  double norm = sqrt(colloc->Integrate(xxDot->getV1(), xxDot->getV1()));
  xxDot->getV1() /= norm;
  xxDot->getV3().Clear();
  sol += ds * xxDot->getV1();
}


void Point::Stability()
{
  mRe.Clear();
  mIm.Clear();

  colloc->Init(sol, par);

  colloc->StabJac(jacStab, par);

  jacStab.Eigval(mRe, mIm);
//  mRe.Print();
//  mIm.Print();
}

void Point::Plot(GnuPlot& pl)
{
  pl.SetPointSize(0.8);

  for (int i = 0; i < NDIM; i++)
  {
    pl.Plot(i, "with lines");
    for (int j = 0; j < NINT*NDEG + 1; j++)
    {
      double t = (double)j / ((double)(NINT * NDEG));
      pl.AddData(i, t, sol(i + NDIM*j));
    }
  }
  pl.Show();
}

///--------------------------------
/// INPUT & OUTPUT
///--------------------------------

void Point::Write(std::ofstream& file)
{
  Vector msh(colloc->getMesh());

  file << NPAR + ParEnd << "\t";
  for (int i = 0; i < NPAR + ParEnd; i++) file << par(i) << "\t";

  file << mRe.Size() << "\t";
  for (int i = 0; i < mRe.Size(); i++) file << mRe(i) << "\t" << mIm(i) << "\t";

  file << NDIM << "\t";
  file << NINT << "\t";
  file << NDEG << "\t";
  for (int i = 0; i < NINT; i++)
  {
    for (int j = 0; j < NDEG; j++)
    {
      file << colloc->Profile(i, j) << "\t";
    }
  }
  file << colloc->getMesh()(NINT) << "\t";
  for (int i = 0; i < NDIM*(NINT*NDEG + 1); i++) file << sol(i) << "\t";
  file << "\n";
  file.flush();
}

void Point::Read(std::ifstream& file)
{
  int npar_, nmul_, ndim_, nint_, ndeg_;
  file >> npar_;
  P_ERROR_X3(NPAR + ParEnd == npar_, "Not compatible input file (NPAR) ", npar_, ".");
  for (int i = 0; i < NPAR + ParEnd; i++) file >> par(i);

  file >> nmul_;
  P_ERROR_X3(mRe.Size() >= nmul_, "Not compatible input file (NMUL) ", nmul_, ".");
  for (int i = 0; i < nmul_; i++)
  {
    file >> mRe(i);
    file >> mIm(i);
  }

  file >> ndim_;
  file >> nint_;
  file >> ndeg_;

  P_ERROR_X3(NDIM == ndim_, "Not compatible input file (NDIM) ", ndim_, ".");

  Vector msh(nint_ + 1);
  double t;
  for (int i = 0; i < ndeg_*nint_ + 1; i++)
  {
    file >> t;
    if (i % ndeg_ == 0) msh(i / ndeg_) = t;
  }

  if ((NINT == nint_) && (NDEG == ndeg_))
  {
    colloc->setMesh(msh);
    for (int i = 0; i < NDIM*(nint_*ndeg_ + 1); i++) file >> sol(i);
  }
  else
  {
    Vector in(NDIM*(nint_*ndeg_ + 1));

    for (int i = 0; i < NDIM*(nint_*ndeg_ + 1); i++) file >> in(i);
    colloc->Import(sol, in, msh, ndeg_);
  }
}

void Point::ReadNull(std::ifstream& file)
{
  double tmp;
  int npar_, nmul_, ndim_, nint_, ndeg_;
  file >> npar_;
  P_ERROR_X3(NPAR + ParEnd == npar_, "Not compatible input file (NPAR) ", npar_, ".");
  for (int i = 0; i < NPAR + ParEnd; i++) file >> tmp;

  file >> nmul_;
  P_ERROR_X3(mRe.Size() >= nmul_, "Not compatible input file (NMUL) ", nmul_, ".");
  for (int i = 0; i < nmul_; i++)
  {
    file >> tmp;
    file >> tmp;
  }

  file >> ndim_;
  file >> nint_;
  file >> ndeg_;

  P_ERROR_X3(NDIM == ndim_, "Not compatible input file (NDIM) ", ndim_, ".");

  for (int i = 0; i < ndeg_*nint_ + 1; i++) file >> tmp;
  for (int i = 0; i < NDIM*(nint_*ndeg_ + 1); i++) file >> tmp;
}

void Point::SwitchTFTRTan(Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg)   // starting data for tori: tangent
{
  Vector TRe(sol.Size()), TIm(sol.Size());
  int idx = 0;
  for (int i=0; i<eqn.Size(); ++i) if (eqn(i) == EqnTFCPLX_RE) { idx = i; break; }
  TestFunctCPLX* tf = static_cast< TestFunctCPLX* >(testFunct(idx));
  if (tf)
  {
    tf->Switch(TRe, TIm, alpha);
    colloc->Export(Re, mshint, mshdeg, TRe);
    colloc->Export(Im, mshint, mshdeg, TIm);
  }
  else
  {
    P_MESSAGE1("Cannot compute the initial tangent to the torus branch, because no complex test functional was defined.");
  }
}

void Point::BinaryWrite(mat4Data& data, int n)
{
  data.setNTrivMul(0, nTrivMulLP);
  data.setNTrivMul(1, nTrivMulPD);
  data.setNTrivMul(2, nTrivMulNS);
  //
  data.setPar(n, par);
  data.setMul(n, mRe, mIm);
  data.setElem(n, colloc->getElem());
  data.setMesh(n, colloc->getMesh());
  data.setProfile(n, sol);
}

void Point::BinaryRead(mat4Data& data, int n)
{
  Vector msh(data.getNInt() + 1);
  P_ERROR_X1(data.getNPar() == (NPAR + ParEnd), "Wrong number of parameters in the input MAT file.");
  data.getPar(n, par);
  data.getMul(n, mRe, mIm);
  data.getMesh(n, msh);
  P_ERROR_X1(data.getNDim() == NDIM, "Wrong number of dimensions in the input MAT file.");
  if (data.getNInt() == NINT && data.getNDeg() == NDEG)
  {
    colloc->setMesh(msh);
    data.getProfile(n, sol);
  }
  else
  {
    Vector tmp(data.getNDim()*(data.getNDeg()*data.getNInt() + 1));
    data.getProfile(n, tmp);
    colloc->Import(sol, tmp, msh, data.getNDeg());
  }
}
