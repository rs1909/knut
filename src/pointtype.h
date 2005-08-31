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
	EqnLPNullSpace    = 8,
	EqnPDNullSpace    = 9,
	EqnCPLXNullSpace  = 10,
	EqnLPAUTNullSpace = 11,
	EqnLPAUTROTNullSpace = 12,
	EqnNorm           = 16,
	EqnCPLXNormRe     = 17,
	EqnCPLXNormIm     = 18,
	EqnPhase          = 24,
	EqnPhaseRot       = 25,
	EqnTORSol         = 32,
	EqnTORPhase0      = 40,
	EqnTORPhase1      = 41
};

// this applies to torpoint.h too...

enum Var {
	VarNone      = 0,
	VarSol       = 1,
	VarNullSpace = 8,
	VarTORSol    = 32,
	VarPAR0      = 64
};

// this applies to torpoint.h too...

enum PtType {
	Sol        = 0,
	BifLP      = 1,
	BifPD      = 2,
	BifNS      = 3,
	SolPDSW    = 4,
	SolLPSW    = 5,
	SolAUT     = 10,
	BifAUTLP   = 11,
	BifAUTPD   = 12,
	BifAUTNS   = 13,
	SolAUTPDSW = 14,
	SolAUTHOPFSW = 16,
	SolTor     = 20,
	SolAUTTor  = 30
};

// internal parameters, placed after the defined parameters

enum Par { ParAngle=0, ParNorm, ParRot, ParEnd };

enum BranchSW { NOSwitch=0, LPSwitch, PDSwitch, HOPFSwitch, TORSwitch };

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
		case EqnLPNullSpace:    return "EqnLPNullSpace   "; break;
		case EqnPDNullSpace:    return "EqnPDNullSpace   "; break;
		case EqnCPLXNullSpace:  return "EqnCPLXNullSpace "; break;
		case EqnLPAUTNullSpace: return "EqnLPAUTNullSpace"; break;
		case EqnNorm:           return "EqnNorm          "; break;
		case EqnCPLXNormRe:     return "EqnCPLXNormRe    "; break;
		case EqnCPLXNormIm:     return "EqnCPLXNormIm    "; break;
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
		case VarNullSpace:      return "VarNullSpace     "; break;
		case VarTORSol:         return "VarTORSol        "; break;
		default:                return "Var ??????????   "; break;
	}
	if( var >= VarPAR0 ) return "VarPAR>0         ";
	return "VarInvalid";
}

template< class T > class Array1D;

// helper function
void PtToEqnVar( Array1D<Eqn>& eqn_, Array1D<Var>& var_, PtType Pt, int nparx, const int* parx, int npar_ );

#endif
