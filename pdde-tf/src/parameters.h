#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "pointtype.h"
#include <string>
#include <iostream>

// pre-declarations
class System;
template< class T > class Array1D;

#define MAX_PARX 8
#define MAX_EQVAR 8
#define MAX_SYM 4

/*
sys-glass.so                      SYSNAME
0                                 LABEL
{
50 P2 1 P0                        TYPE, CP, NPARX, PARX[NPARX]
or
-1 P2 0 3 E1 E24 E25 3 S1 P0 P11  TYPE, CP, SWITCH, NEQN, EQN[NEQN], NVAR, VAR[NVAR]
}
200 4 5 1 1                       NINT, NDEG, NMUL, STAB, NMAT
12 12 4 4                         NINT1, NINT2, NDEG1, NDEG2
120 -100 100                      STEPS, CPMIN, CPMAX
-0.002 0.002 0.002 -0.002         DS, DSMIN, DSMAX, DSSTART
1e-05 1e-05 1e-06                 EPSC, EPSR, EPSS
12 12 12                          NITC, NITR, NITS
2 0 1 3 4                         NSYM, R0, I0, R1, I1
*/

class ConstFile
{
  public:
	// function members
	void     clear() {}
	void     getParams( const std::string& filename );
	void     printParams( std::ostream& file = std::cout );
	int      toEqnVar( System& sys,
                Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN );
	// data members
	std::string SYSNAME;    // name of the shared object file
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
