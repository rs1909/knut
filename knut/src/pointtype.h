// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef POINTTYPE_H
#define POINTTYPE_H

enum PtType {
  SolUser      = -1,
  /// ODE
  SolODE       = 0,
  SolAUTODE    = 10,
  /// TORUS for DDE
  SolTor       = 20,
  SolTorNS     = 21,
  SolAUTTor    = 30,
  SolAUTTorNS  = 31,
  /// TIME-PERIODIC TEST-FUNCTIONAL for DDE
  SolTF        = 40,
  BifTFLP      = 41,
  BifTFPD      = 42,
  BifTFNS      = 43,
  SolTFBRSW    = 44,
  SolTFPDSW    = 45,
  /// AUTONOMOUS TEST-FUNCTIONAL for DDE
  SolTFAUT     = 50,
  BifTFAUTLP   = 51,
  BifTFAUTPD   = 52,
  BifTFAUTNS   = 53,
  SolTFAUTBRSW = 54,
  SolTFAUTPDSW = 55,
  SolTFAUTHBSW = 56
};

enum Eqn {
  EqnNone           = 0,
  EqnSol            = 1,
  EqnODESol         = 2,
  EqnPhase          = 24,
  EqnPhaseRot       = 25,
  EqnTORSol         = 32,
  EqnTORPhase0      = 40,
  EqnTORPhase1      = 41,
  EqnTFLP           = 64,
  EqnTFPD           = 65,
  EqnTFLPAUT        = 66,
  EqnTFLPAUTROT     = 67,
  EqnTFCPLX_RE      = 68,
  EqnTFCPLX_IM      = 69
};

enum Var {
  VarNone      = 0,
  VarSol       = 1,
  VarODESol    = 2,
  VarTORSol    = 32,
  VarPAR0      = 100,
  VarInternal  = 4000,
  VarAngle     = 4000,
  VarPeriod    = 4001,
  VarRot       = 4002,
  VarEnd
};

enum BranchSW
{
  NOSwitch = 0,
  TFBRSwitch = 1,
  TFPDSwitch = 2,
  TFHBSwitch = 3,
  TFTRSwitch = 4,
  TFBRAUTSwitch = 5,
  TFBRAUTROTSwitch = 6
};

enum BifType
{
  BifNone,
  BifLP,
  BifPD,
  BifNS,
  BifMax,
  BifEndPoint,
  BifNoConvergence
};

template<typename TP> struct TypeTuple
{
  int          index;
  TP           type;
  const char * name;
};

inline int VarToIndex(Var v, int npar)
{
  if (v < VarInternal) return v - VarPAR0;
  else return npar + (v - VarInternal);
}

inline Var IndexToVar(int idx, int npar)
{
  if (idx < npar) return static_cast<Var>(VarPAR0 + idx);
  else return static_cast<Var>(VarInternal + (idx - npar));
}

inline char parType(int npar, int p)
{
  if (p < npar) return 'P';
  else return 'I';
}

inline int parNum(int npar, int p)
{
  if (p < npar) return p;
  else return p - npar;
}

inline const char* EqnToStr(Eqn eqn)
{
  switch (eqn)
  {
    case EqnNone:
      return "EqnNone";
      break;
    case EqnSol:
      return "EqnSol";
      break;
    case EqnODESol:
      return "EqnODESol";
      break;
    case EqnPhase:
      return "EqnPhase";
      break;
    case EqnPhaseRot:
      return "EqnPhaseRot";
      break;
    case EqnTORSol:
      return "EqnTORSol";
      break;
    case EqnTORPhase0:
      return "EqnTORPhase0";
      break;
    case EqnTORPhase1:
      return "EqnTORPhase1";
      break;
    case EqnTFLP:
      return "EqnTFLP";
      break;
    case EqnTFPD:
      return "EqnTFPD";
      break;
    case EqnTFLPAUT:
      return "EqnTFLPAUT";
      break;
    case EqnTFLPAUTROT:
      return "EqnTFLPAUTROT";
      break;
    case EqnTFCPLX_RE:
      return "EqnTFCPLX_RE";
      break;
    case EqnTFCPLX_IM:
      return "EqnTFCPLX_IM";
      break;
    default:
      return "Eqn ?";
      break;
  }
  return "EqnInvalid";
}

inline const char* VarToStr(Var var)
{
  if (var >= VarPAR0) return "VarPAR>0";
  switch (var)
  {
    case VarNone:
      return "VarNone";
      break;
    case VarSol:
      return "VarSol";
      break;
    case VarODESol:
      return "VarODESol";
      break;
    case VarTORSol:
      return "VarTORSol";
      break;
    default:
      return "Var ?";
      break;
  }
}

template< class T > class KNArray1D;

// helper function, implemented in basepoint.cpp
BranchSW PtToEqnVar(KNArray1D<Eqn>& eqnr, KNArray1D<Var>& varr, PtType Pt, KNArray1D<Var> parx, int npar_);

#endif
