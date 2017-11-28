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

#include <cstring>
#include <string>
#include <iostream>

enum class PtType {
  SolUser      = 64,
  /// ODE
  SolODE       = 0,
  SolAUTODE    = 10,
  /// Steady state for DDE
  SolSteady    = 15,
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
  SolTFAUTHBSW = 56,
  SolTFAUTSSMSW = 57
};

enum Eqn {
  EqnNone           = 0,
  EqnSol            = 1,
  EqnODESol         = 2,
  EqnSteady         = 3,
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
  VarSteady    = 3,
  VarTORSol    = 32,
  VarPAR0      = 100,
  VarInternal  = 4000,
  VarAngle     = 4000,
  VarPeriod    = 4001,
  VarRot       = 4002,
  VarEnd
};

enum class BranchSW
{
  NOSwitch = 0,
  TFBRSwitch = 1,
  TFPDSwitch = 2,
  TFHBSwitch = 3,
  TFSSMSwitch = 4,
  TFTRSwitch = 5,
  TFBRAUTSwitch = 6,
  TFBRAUTROTSwitch = 7
};

inline std::ostream& operator<<(std::ostream& os, BranchSW sw)
{  
  os << static_cast<int>(sw);  
  return os;  
}

enum class BifType
{
  BifNone,
  BifLP,
  BifPD,
  BifNS,
  BifUN,
  BifMax,
  BifEndPoint,
  BifNoConvergence
};

inline size_t VarToIndex(Var v, size_t npar)
{
  if (v < VarInternal) return v - VarPAR0;
  else return npar + (v - VarInternal);
}

inline Var IndexToVar(size_t idx, size_t npar)
{
  if (idx < npar) return static_cast<Var>(VarPAR0 + idx);
  else return static_cast<Var>(VarInternal + (idx - npar));
}

template< class T > class KNArray1D;

// helper function, implemented in basepoint.cpp
BranchSW PtToEqnVar(KNArray1D<Eqn>& eqnr, KNArray1D<Var>& varr, PtType Pt, KNArray1D<Var> parx, size_t npar_);

#define SYS_TYPE_SO "so"
#define SYS_TYPE_CXX "cpp"
#define SYS_TYPE_VFC "vf"
#define SYS_TYPE_VF0 "vf0"

template<typename TP> struct TypeTuple
{
  size_t       index;
  TP           type;
  const char * code;
  const char * name;
};

template<typename TP> struct TypeTupleStr
{
  size_t       index;
  TP           type;
  std::string  code;
  std::string  name;
};

template<typename TP> class TypeTupleTabBase
{
  public:
    static const TypeTuple<TP> tabStatic[];
};

template<> const TypeTuple<PtType> TypeTupleTabBase<PtType>::tabStatic[];
template<> const TypeTuple<Eqn> TypeTupleTabBase<Eqn>::tabStatic[];
template<> const TypeTuple<BranchSW> TypeTupleTabBase<BranchSW>::tabStatic[];
template<> const TypeTuple<Var> TypeTupleTabBase<Var>::tabStatic[];
template<> const TypeTuple<BifType> TypeTupleTabBase<BifType>::tabStatic[];

#endif
