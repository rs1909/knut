// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "constants.h"
#include <iostream>
#include <fstream>

inline bool NConstants::inputAssert(std::istream& is)
{
  if (is.rdstate() & (std::istream::eofbit | std::istream::failbit | std::istream::badbit))
  {
    switch (is.rdstate())
    {
      case std::istream::eofbit:
        P_MESSAGE("Unexpected end of file");
        return true;
        break;
      case std::istream::failbit:
        P_MESSAGE("Input failed\n");
        return true;
        break;
      case std::istream::badbit:
        P_MESSAGE("Bad input\n");
        return true;
        break;
      default:
        P_MESSAGE("Unexpected error");
        return true;
        break;
    }
  }
  return false;
}

void NConstants::loadFile(const std::string &fileName)
{
  std::ifstream file;

  file.open(fileName.c_str());
  P_ERROR_X2(file, "Cannot open ", fileName);

  std::string __sysname__;
  file >> __sysname__;
  inputAssert(file);
  while (file.get() != '\n');
  inputAssert(file);

//  if( __sysname__.find('/') == std::string::npos )
//  {
//   __sysname__ = QFileInfo( fileName.c_str() ).absolutePath( ).toStdString() + "/" + __sysname__;
//  }
  setSysNameText(std::string(__sysname__.c_str()));

  int __ptlabel__;
  file >> __ptlabel__;
  inputAssert(file);
  setLabel(__ptlabel__);
  while (file.get() != '\n');
  inputAssert(file);

  int __pttype__;
  char __cptype__;
  int __cp__;
  file >> __pttype__ >> __cptype__ >> __cp__;
  inputAssert(file);
  P_ERROR_X1((__cptype__ == 'P') || (__cptype__ == 'I'), "Error: CP: Bad parameter type.");
  setPointType((PtType)__pttype__);
  setCp(__cptype__, __cp__);

  bool loadsym = false;
  if (__pttype__ == SolUser)
  {
    int __brswitch__;
    int __neqn__;
    file >> __brswitch__;
    inputAssert(file);
    file >> __neqn__;
    inputAssert(file);
    setNEqns(__neqn__);
    for (int i = 0; i < __neqn__; i++)
    {
      char __eqntype__;
      int  __eqnnum__;
      file >> __eqntype__ >> __eqnnum__;
      inputAssert(file);
      P_ERROR_X1(__eqntype__ == 'E', "Error: EQN: Bad equation type.");
      setEqns(i, 'E', __eqnnum__);
      if (__eqnnum__ == EqnPhaseRot) loadsym = true;
    }
    int __nvar__;
    file >> __nvar__;
    inputAssert(file);
    P_ERROR_X1(__nvar__ == __neqn__, "Error: NVAR: Number of equations and variables are not the same.");
    for (int i = 0; i < __neqn__; i++)
    {
      char __vartype__;
      int  __varnum__;
      file >> __vartype__ >> __varnum__;
      inputAssert(file);
      P_ERROR_X1((__vartype__ == 'S') || (__vartype__ == 'P') || (__vartype__ == 'I'), "Error: VAR: Bad variable/parameter type.");
      setVars(i, __vartype__, __varnum__);
    }
    while (file.get() != '\n');
  }
  else
  {
    int __nvarx__;
    file >> __nvarx__;
    inputAssert(file);
    setNEqns(__nvarx__);
    for (int i = 0; i < __nvarx__; i++)
    {
      char __vartype__;
      int  __varnum__;
      file >> __vartype__ >> __varnum__;
      inputAssert(file);
      P_ERROR_X1((__vartype__ == 'P') || (__vartype__ == 'I'), "Error: PARX: Bad parameter type.");
      setParX(i, __vartype__, __varnum__);
    }
    while (file.get() != '\n');
  }

  int __nint__, __ndeg__, __nmul__, __stab__, __nmat__;
  file >> __nint__ >> __ndeg__ >> __nmul__ >> __stab__ >> __nmat__;
  inputAssert(file);
  while (file.get() != '\n');
  setNInt(__nint__);
  setNDeg(__ndeg__);
  setNMul(__nmul__);

  setStab(__stab__ != 0);
  setNMat(__nmat__);

  int __nint1__, __nint2__, __ndeg1__, __ndeg2__;
  inputAssert(file);
  file >> __nint1__ >> __nint2__ >> __ndeg1__ >> __ndeg2__;
  while (file.get() != '\n');
  setNInt1(__nint1__);
  setNInt2(__nint2__);
  setNDeg1(__ndeg1__);
  setNDeg2(__ndeg2__);

  int __steps__;
  double __cpmin__, __cpmax__;
  file >> __steps__ >> __cpmin__ >> __cpmax__;
  inputAssert(file);
  while (file.get() != '\n');
  setSteps(__steps__);
  setCpMin(__cpmin__);
  setCpMax(__cpmax__);

  double __ds__, __dsmin__, __dsmax__, __dsstart__;
  file >> __ds__ >> __dsmin__ >> __dsmax__ >> __dsstart__;
  inputAssert(file);
  while (file.get() != '\n');
  setDs(__ds__);
  setDsMin(__dsmin__);
  setDsMax(__dsmax__);
  setDsStart(__dsstart__);

  double __epsc__, __epsr__, __epss__;
  file >> __epsc__ >> __epsr__ >> __epss__;
  inputAssert(file);
  while (file.get() != '\n');
  setEpsC(__epsc__);
  setEpsR(__epsr__);
  setEpsS(__epss__);

  int __nitc__, __nitr__, __nits__;
  file >> __nitc__ >> __nitr__ >> __nits__;
  inputAssert(file);
  setNItC(__nitc__);
  setNItR(__nitr__);
  setNItS(__nits__);

  if (loadsym)
  {
    while (file.get() != '\n');
    int __nsym__;
    file >> __nsym__;
    inputAssert(file);
    setNSym(__nsym__);
    int __re__, __im__;
    for (int i = 0; i < __nsym__; ++i)
    {
      file >> __re__ >> __im__;
      inputAssert(file);
      setSymRe(i, __re__);
      setSymIm(i, __im__);
    }
  }
  else
  {
    setNSym(0);
  }
}

void NConstants::saveFile(const std::string &fileName)
{
  std::ofstream file(fileName.c_str());
  printFile(file);
}

void NConstants::printFile(std::ostream& file)
{

  file << getSysName() << " \t\tSYSNAME\n";

  file << getLabel() << " \t\t\tLABEL\n";

  const PtType __pttype__ = getPointType();
  if (__pttype__ != SolUser)
  {
    file << (int)__pttype__ << " " << getCpType() << getCpNum() << " " << getNEqns() << " ";
    for (int i = 0; i < getNEqns(); ++i) file << getParXType(i) << getParXNum(i) << " ";
    file << "\t\tTYPE, CP, NPARX, PARX[NPARX]\n";
  }
  else
  {
    file << (int)__pttype__ << " " << getCpType() << getCpNum() << " "
    << 0 /*switch */ << " " << getNEqns() << " ";
    for (int i = 0; i < getNEqns(); i++) file << getEqnsType(i) << getEqnsNum(i) << " ";
    file << getNEqns() << " ";
    for (int i = 0; i < getNEqns(); i++) file << getVarsType(i) << getVarsNum(i) << " ";
    file << "\tTYPE, CP, SWITCH, NEQN, EQN[NEQN], NVAR, VAR[NVAR]\n";
  }

  file << getNInt() << " "
  << getNDeg() << " "
  << getNMul() << " "
  << getStab() << " "
  << getNMat() << " \t\tNINT, NDEG, NMUL, STAB, NMAT\n";

  file << getNInt1() << " "
  << getNInt2() << " "
  << getNDeg2() << " "
  << getNDeg2() << " \t\tNINT1, NINT2, NDEG1, NDEG2\n";

  file << getSteps() << " "
  << getCpMin() << " "
  << getCpMax() << " \t\tSTEPS, CPMIN, CPMAX\n";

  file << getDs() << " "
  << getDsMin() << " "
  << getDsMax() << " "
  << getDsStart() << " \tDS, DSMIN, DSMAX, DSSTART \n";

  file << getEpsC() << " "
  << getEpsR() << " "
  << getEpsS() << " \tEPSC, EPSR, EPSS\n";

  file << getNItC() << " " << getNItR() << " " << getNItS() << " \t\tNITC, NITR, NITS\n";

  file << getNSym() << " ";
  for (int i = 0; i < getNSym(); ++i) file << getSymRe(i) << " " << getSymIm(i) << " ";
  file << "\t\tNSYM";
  for (int i = 0; i < getNSym(); ++i) file << ", R" << i << ", I" << i;
  file << "\n";

}

int NConstants::toEqnVar(System& sys,
                         Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                         Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                         Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN)
{
  // initializing the equations and variables
  if (getPointType() == SolUser)
  {
    eqn.Init(getNEqns());
    var.Init(getNEqns());
    for (int i = 0; i < getNEqns(); i++)
    {
      eqn(i) = getEqns(i);
      var(i) = getVars(i);
    }
  }
  else
  {
    Array1D<Var> L_PARX(getNEqns());
    for (int i = 0; i < getNEqns(); i++)
    {
      L_PARX(i) = getParX(i);
    }
    setBranchSW(PtToEqnVar(eqn, var, getPointType(), L_PARX, sys.npar()));
  }

  // checking whether it is an autonomous problem or not
  bool aut = false;
  bool phaseRot = false;
  testFN = EqnNone;
  int  testFN_idx = -1;
  for (int i = 0; i < eqn.Size(); i++)
  {
    if (eqn(i) == EqnPhase) aut = true;
    if (eqn(i) == EqnPhaseRot) phaseRot = true;
    if ((eqn(i) == EqnTFPD) ||
        (eqn(i) == EqnTFLP) ||
        (eqn(i) == EqnTFLPAUT) ||
        (eqn(i) == EqnTFLPAUTROT) ||
        (eqn(i) == EqnTFCPLX_RE))
    {
      if (testFN_idx == -1)
      {
        testFN = eqn(i);
        testFN_idx = i;
      }
      else P_MESSAGE("Too many test functionals.");
    }
  }

  // setting up for refinement
  if (aut)
  {
    if (phaseRot)
    {
      // std::cout<<"Phase and PhaseRot\n";
      eqn_refine.Init(3);
      var_refine.Init(3);
      eqn_refine(0) = EqnSol;
      eqn_refine(1) = EqnPhase;
      eqn_refine(2) = EqnPhaseRot;
      var_refine(0) = VarSol;
      var_refine(1) = var(var.Size() - 2);
      var_refine(2) = var(var.Size() - 1);
    }
    else
    {
      // std::cout<<"Phase\n";
      eqn_refine.Init(2);
      var_refine.Init(2);
      eqn_refine(0) = EqnSol;
      eqn_refine(1) = EqnPhase;
      var_refine(0) = VarSol;
      var_refine(1) = var(var.Size() - 1);
    }
  }
  else
  {
    eqn_refine.Init(1);
    var_refine.Init(1);
    eqn_refine(0) = EqnSol;
    var_refine(0) = VarSol;
  }

  if (getBranchSW() == TFHBSwitch)
  {
    Array1D<Var> d(0);
    PtToEqnVar(eqn_refine, var_refine, SolTF, d, sys.npar());
  }// NPARX == 0

  // Here, we set up the branch switching.
  // We suppose that if there is a switch we use one parameter continuation afterwards
  // without using characteristic matrices. This means that we can switch on the characteristic matrix,
  // include the equation for the eigenvector norm before the other equations and
  // add CP to the variables as a normal parameter.
  Eqn eqn_temp;
  switch (getBranchSW())
  {
/// with TEST FUNCTIONALS
    case TFBRSwitch:
      if (aut) eqn_temp = EqnTFLPAUT;
      else eqn_temp = EqnTFLP;
      goto tfskip;
    case TFPDSwitch:
      eqn_temp = EqnTFPD;
      goto tfskip;
tfskip:
      eqn_start.Init(eqn_refine.Size() + 1);
      var_start.Init(var_refine.Size() + 1);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = eqn_temp;
      var_start(var_refine.Size()) = getCp();
      for (int i = 1; i < eqn_refine.Size(); i++)
      {
        eqn_start(i + 1) = eqn_refine(i);
        var_start(i) = var_refine(i);
      }
      testFN = eqn_temp;
      break;
    case TFTRSwitch:
    case TFHBSwitch:
      eqn_start.Init(eqn_refine.Size() + 2);
      var_start.Init(var_refine.Size() + 2);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = EqnTFCPLX_RE;
      eqn_start(2) = EqnTFCPLX_IM;
      var_start(1) = (Var)(VarPAR0 + sys.npar() + ParAngle); // CH
      var_start(var_refine.Size() + 1) = getCp();
      for (int i = 1; i < eqn_refine.Size(); i++)
      {
        eqn_start(i + 2) = eqn_refine(i);
        var_start(i + 1) = var_refine(i);
      }
      testFN = EqnTFCPLX_RE;
      break;
    default:
      eqn_start.Init(eqn.Size());
      var_start.Init(var.Size());
      eqn_start = eqn;
      var_start = var;
      break;
  }
  if (!aut) return 0;
  else if (!phaseRot) return 1;
  else return 2;
}
