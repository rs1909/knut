// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef POINTTYPE_H
#define POINTTYPE_H

// this applies to torpoint.h too...

enum Eqn {
	EqnNone           = 0,
	EqnSol            = 1,
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

// this applies to torpoint.h too...

enum Var {
	VarNone      = 0,
	VarSol       = 1,
	VarTORSol    = 32,
	VarPAR0      = 64
};

// this applies to torpoint.h too...

enum PtType {
/// TORUS
	SolTor       = 20,
	SolAUTTor    = 30,
/// TIME-PERIODIC TEST-FUNCTIONAL
	SolTF        = 40,
	BifTFLP      = 41,
	BifTFPD      = 42,
	BifTFNS      = 43,
	SolTFBRSW    = 44,
	SolTFPDSW    = 45,
/// AUTONOMOUS TEST-FUNCTIONAL
	SolTFAUT     = 50,
	BifTFAUTLP   = 51,
	BifTFAUTPD   = 52,
	BifTFAUTNS   = 53,
	SolTFAUTBRSW = 54,
	SolTFAUTPDSW = 55,
	SolTFAUTHBSW = 56
};

// internal parameters, placed after the defined parameters

enum Par { ParAngle=0, ParRot=2, ParEnd };

enum BranchSW { NOSwitch=0, TFBRSwitch=1, TFPDSwitch=2, TFHBSwitch=3, TFTRSwitch=4 };

inline char parType( int npar, int p )
{
	if( p < npar ) return 'P';
	else return 'I';
}

inline int parNum( int npar, int p )
{
	if( p < npar ) return p;
	else return p - npar;
}

inline const char* EqnToStr( Eqn eqn )
{
	switch( eqn )
	{
		case EqnNone:           return "EqnNone          "; break;
		case EqnSol:            return "EqnSol           "; break;
		case EqnPhase:          return "EqnPhase         "; break;
		case EqnPhaseRot:       return "EqnPhaseRot      "; break;
		case EqnTORSol:         return "EqnTORSol        "; break;
		case EqnTORPhase0:      return "EqnTORPhase0     "; break;
		case EqnTORPhase1:      return "EqnTORPhase1     "; break;
		default:                return "Eqn ??????????   "; break;
	}
	return "EqnInvalid";
}

inline const char* VarToStr( Var var )
{
	switch( var )
	{
		case VarNone:           return "VarNone          "; break;
		case VarSol:            return "VarSol           "; break;
		case VarTORSol:         return "VarTORSol        "; break;
		default:                return "Var ??????????   "; break;
	}
	if( var >= VarPAR0 ) return "VarPAR>0         ";
	return "VarInvalid";
}

template< class T > class Array1D;

// helper function
BranchSW PtToEqnVar( Array1D<Eqn>& eqn_, Array1D<Var>& var_, PtType Pt, int nparx, const int* parx, int npar_ );

#endif
