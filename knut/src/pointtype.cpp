// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005, 2012 by Robert Szalai
//
// For license, see the file COPYING in the packages package's root directory
//
// ------------------------------------------------------------------------- //

#include "pointtype.h"

template<> const TypeTuple<PtType> TypeTupleTabBase<PtType>::tabStatic[] = {
  {0,  PtType::SolUser,      "USER",             "User defined"}, // User
  {1,  PtType::SolODE,       "ODE",              "ODE limit cycle"}, // ODE limit cycle
  {2,  PtType::SolAUTODE,    "ODE_AUT",          "ODE limit cycle (aut)" }, // ODE limit cycle (aut)
  {3,  PtType::SolSteady,    "DDE_STEADY",       "Steady State"},
  {4,  PtType::SolTF,        "DDE_LC",           "Limit cycle"}, // "Limit cycle
  {5,  PtType::BifTFLP,      "DDE_LP",           "Limit point"}, // Limit point
  {6,  PtType::BifTFPD,      "DDE_PD",           "Period doubling"}, // Period doubling
  {7,  PtType::BifTFNS,      "DDE_NS",           "Neimark-Sacker"}, // Neimark-Sacker
  {8,  PtType::SolTFBRSW,    "DDE_BP_SW",        "Branch switch"}, // Branch switch
  {9,  PtType::SolTFPDSW,    "DDE_PD_SW",        "Period doubling switch"}, // Period doubling switch
  {10, PtType::SolTFAUT,     "DDE_AUT",          "Limit cycle (aut)"}, // Limit cycle (aut)
  {11, PtType::BifTFAUTLP,   "DDE_AUT_LP",       "Limit point (aut)"}, //
  {12, PtType::BifTFAUTPD,   "DDE_AUT_PD",       "Period doubling (aut)"},
  {13, PtType::BifTFAUTNS,   "DDE_AUT_NS",       "Neimark-Sacker (aut)"},
  {14, PtType::SolTFAUTBRSW, "DDE_AUT_BP_SW",    "Branch switch (aut)"},
  {15, PtType::SolTFAUTPDSW, "DDE_AUT_PD_SW",    "Period doubling switch (aut)"},
  {16, PtType::SolTFAUTHBSW, "DDE_AUT_HOPF_SW",  "Hopf switch (aut)"},
  {17, PtType::SolTFAUTSSMSW,"DDE_AUT_SSM_SW",  "Subcenter manifold switch (aut)"},
  {18, PtType::SolTor,       "DDE_TORUS",        "Torus"},
  {19, PtType::SolTorNS,     "DDE_TORUS_SW",     "Torus from NS"},
  {20, PtType::SolAUTTor,    "DDE_AUT_TORUS",    "Torus (aut)"},
  {21, PtType::SolAUTTorNS,  "DDE_AUT_TORUS_SW", "Torus from NS (aut)"},
  {~(size_t)0, PtType::SolUser, "", ""}};

template<> const TypeTuple<Eqn> TypeTupleTabBase<Eqn>::tabStatic[] = {
  {0,  EqnNone,       "NX",                "None" },
  {1,  EqnSol,        "DDE_LC",            "Limit cycle"},
  {2,  EqnODESol,     "ODE_LC",            "ODE Limit cycle"},
  {3,  EqnSteady,     "DDE_STEADY",        "ODE Limit cycle"},
  {4,  EqnPhase,      "PHASE",             "PHCND (transl)"},
  {5,  EqnPhaseRot,   "PHASE_ROT",         "PHCND (rot)"},
  {6,  EqnTORSol,     "DDE_TORUS",         "Torus"},
  {7,  EqnTORPhase0,  "TORUS_PHS1",        "PHCND 1 (torus)"},
  {8,  EqnTORPhase1,  "TORUS_PHS2",        "PHCND 2 (torus)"},
  {9,  EqnTFLP,       "DDE_TF_LP",         "LP TF"},
  {10, EqnTFPD,       "DDE_TF_PD",         "PD TF"},
  {11, EqnTFLPAUT,    "DDE_AUT_TF_LP",     "LP TF (aut)"},
  {12, EqnTFLPAUTROT, "DDE_AUT_ROT_TF_LP", "LP TF (sym)"},
  {13, EqnTFCPLX_RE,  "DDE_TF_NS_RE",      "NS TF (real)"},
  {14, EqnTFCPLX_IM,  "DDE_TF_NS_IM",      "NS TF (imag)"},
  {~(size_t)0, EqnNone, "", ""}};
  
template<> const TypeTuple<BranchSW> TypeTupleTabBase<BranchSW>::tabStatic[] = {
  {0, BranchSW::NOSwitch,         "NX",         "No switch"},
  {1, BranchSW::TFBRSwitch,       "BP",         "Branch"},
  {2, BranchSW::TFPDSwitch,       "PD",         "Period doubling"},
  {3, BranchSW::TFHBSwitch,       "HOPF",       "Hopf (aut)"},
  {4, BranchSW::TFTRSwitch,       "TORUS",      "Torus"},
  {5, BranchSW::TFBRAUTSwitch,    "AUT_BP",     "Branch (aut)"},
  {6, BranchSW::TFBRAUTROTSwitch, "AUT_ROT_BP", "Branch (sym)"},
  {~(size_t)0, BranchSW::NOSwitch, "", ""}};

template<> const TypeTuple<Var> TypeTupleTabBase<Var>::tabStatic[] = {
  {0, VarNone,   "NX",          "None"},
  {1, VarSol,    "DDE_LC",      "Limit cycle"},
  {2, VarODESol, "ODE_LC",      "ODE limit cycle"},
  {3, VarODESol, "DDE_STEADY",  "ODE limit cycle"},
  {4, VarTORSol, "DDE_TORUS",   "Torus"},
  {5, VarPAR0,   "PAR0",        "Par 0"},
  {6, VarAngle,  "PAR_ANGLE",   "Angle"},
  {7, VarPeriod, "PAR_PERIOD",  "Period"},
  {8, VarRot,    "PAR_ROT_NUM", "Rot. num."},
  {~(size_t)0, VarEnd, "", ""}};

template<> const TypeTuple<BifType> TypeTupleTabBase<BifType>::tabStatic[] = {
  {0, BifType::BifNone,          "NX", "None"},
  {1, BifType::BifLP,            "LP", "LP"},
  {2, BifType::BifPD,            "PD", "PD"},
  {3, BifType::BifNS,            "NS", "NS"},
  {4, BifType::BifNS,            "UN", "UN"},
  {5, BifType::BifMax,           "MX", "MX"},
  {6, BifType::BifEndPoint,      "EP", "EP"},
  {7, BifType::BifNoConvergence, "NC", "NC"},
  {~(size_t)0, BifType::BifNone, "", ""}};
