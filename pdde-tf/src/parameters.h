#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "pointtype.h"

// pre-declarations
class System;
template< class T > class Array1D;

#define MAX_PARX 8
#define MAX_EQVAR 8
#define MAX_SYM 4

class ConstFile
{
  public:
	// function members
	void     getParams( const char* filename );
	void     printParams( );
	int      toEqnVar( System& sys,
                Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN );
	// data members
	char     SYSNAME[64];   // name of the shared object file
	int      LABEL;         // the restart label in the input file
	int      TYPE;          // type of the solution
	bool     EXTSYS;        // whether we need extended parametrization
	BranchSW SWITCH;        // what type of switch is made
	int      NEQN;          // number of equations
	int      NVAR;          // numger of variables (must be the same as NVAR)
	int      EQN[MAX_EQVAR];  // listing of equations
	int      VAR[MAX_EQVAR];   // listing of variables
	int      VARType[MAX_EQVAR]; // type of the variables
	int      NINT, NDEG;      // for collocation of periodic orbits;
	int      NMUL;            // number of characteristic multipliers
	int      STAB;            // compute stability or not STAB == 0 -> NO else yes
	int      NMAT;            // number of intervals, when computing stability
	int      NINT1, NINT2, NDEG1, NDEG2;
	int      CP;              // the principal continuation parameter
	int      CPType;          // type of the principal continuation parameter
	int      NPARX;           // number of additional parameters
	int      PARX[MAX_PARX];
	int      PARXType[MAX_PARX];
	int      STEPS;
	double   CPMIN;
	double   CPMAX;
	double   DS;
	double   DSMIN;
	double   DSMAX;
	double   DSSTART;
	int      NADAPT;
	double   EPSC;
	double   EPSR;
	double   EPSS;
	int      NITC;
	int      NITR;
	int      NITS;
	int      NSYM;
	int      RESYM[MAX_SYM];
	int      IMSYM[MAX_SYM];
};

#endif
