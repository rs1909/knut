// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#include "pderror.h"
#include "system.h"
#include "point.h"
#include "torpoint.h"
#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define MAX_PARX 8
#define MAX_EQVAR 8
#define MAX_SYM 4

/*

sys-milltwofull.so   SYSNAME
0							LABEL
0 1 0						TYPE, CP, NPARX, PARX ....
200 4 5					NINT, NDEG, NMUL, STAB, NMAT
12 12 4 4				NINT1, NINT2, NDEG1, NDEG2 (for torus computations only)
126 -100.0 100.0		STEPS, CPMIN, CPMAX
0.1 0.01 0.2 0.1		DS, DSMIN, DSMAX, DSSTART
1e-4 1e-4 1e-4	 	   EPSC, EPSR, EPSS
30 30 30					NITC, NITR, NITS

*/

struct cfile {
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

inline void parNamePrint( Vector& /*par*/, int npar, Array1D<Var>& var )
{
	for( int j = 1; j < var.Size(); j++ ) std::cout<<"\t"<<parType( npar, var(j) - VarPAR0 )<<parNum( npar, var(j) - VarPAR0 )<<"\t";
}

inline void parValuePrint( Vector& par, int /*npar*/, Array1D<Var>& var )
{
	for( int j = 1; j < var.Size(); j++ ) std::cout<<"\t"<<par( var(j) - VarPAR0 );
}

inline void pdioassert( std::istream& is )
{
	if( is.rdstate() & (std::istream::eofbit | std::istream::failbit | std::istream::badbit) )
	{
		switch( is.rdstate() )
		{
			case std::istream::eofbit:
				std::cout<<"Unexpected end of file\n";
				break;
			case std::istream::failbit:
				std::cout<<"Input failed\n";
				break;
			case std::istream::badbit:
				std::cout<<"Bad input\n";
				break;
			default:
				break;
		}
		PDError(-1);
	}
}

void getParams( cfile* pms, const char* filename )
{
	std::ifstream file;
	char pt;
	int ival;
	
	file.open( filename );
	if(!file){ std::cout<<"Cannot open "<<filename<<"\n"; PDError(-1); }
	
	file >> pms->SYSNAME; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->LABEL; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->TYPE >> pt >> pms->CP; pms->CPType = pt; pdioassert(file);
	if( (pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: CP: Bad parameter type."; PDError(-1); }
	
	pms->EXTSYS = (pms->TYPE == -1);
	if( pms->EXTSYS )
	{
		// we need extended specifications
		file >> ival;
		pms->SWITCH = (BranchSW)ival;
		file >> pms->NEQN; pdioassert(file);
		if( pms->NEQN > MAX_EQVAR ) { std::cout<<"Error: NEQN: System limit has been reached."; PDError(-1); }
		for( int i = 0; i < pms->NEQN; i++ )
		{
			file >> pt >> (pms->EQN)[i]; pdioassert(file);
			if( pt != 'E' ) { std::cout<<"Error: EQN: Bad equation type."; PDError(-1); }
		}
		file >> pms->NVAR; pdioassert(file);
		if( pms->NEQN != pms->NVAR ) { std::cout<<"Error: NVAR: Number of equations and variables are not the same."; PDError(-1); }
		for( int i = 0; i < pms->NVAR; i++ )
		{
			file >> pt >> (pms->VAR)[i]; (pms->VARType)[i] = pt; pdioassert(file);
			if( (pt != 'S')&&(pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: VAR: Bad variable/parameter type."; PDError(-1); }
		}
		while( file.get() != '\n' ); pdioassert(file);
	}else
	{
		file >> pms->NPARX; pdioassert(file);
		if( pms->NPARX > MAX_PARX ) { std::cout<<"Error: NPARX: System limit has been reached."; PDError(-1); }
		for( int i = 0; i < pms->NPARX; i++ )
		{
			file >> pt >> (pms->PARX)[i]; (pms->PARXType)[i] = pt; pdioassert(file);
			if( (pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: PARX: Bad parameter type."; PDError(-1); }
		}
		while( file.get() != '\n' ); pdioassert(file);
	}
	
	file >> pms->NINT >> pms->NDEG >> pms->NMUL >> pms->STAB >> pms->NMAT; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->NINT1 >> pms->NINT2 >> pms->NDEG1 >> pms->NDEG2; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->STEPS >> pms->CPMIN >> pms->CPMAX; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->DS >> pms->DSMIN >> pms->DSMAX >> pms->DSSTART; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->EPSC >> pms->EPSR >> pms->EPSS; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> pms->NITC >> pms->NITR >> pms->NITS; pdioassert(file);
	// Check whether we need to load NSYM
	bool isSYM = false;
	if( pms->EXTSYS)
	{
		for( int i = 2; i < pms->NEQN; i++ )
		{
			if( (Eqn)(pms->EQN)[i] == EqnPhaseRot ) isSYM = true;
		}
	}
	if( isSYM )
	{
		// first go to the next line
		while( file.get() != '\n' ); pdioassert(file);
		file >> pms->NSYM; pdioassert(file);
		if( pms->NSYM > MAX_SYM ) { std::cout<<"Error: NSYM: System limit has been reached."; PDError(-1); }
		for( int i = 0; i < pms->NSYM; i++ )
		{
			file >> (pms->RESYM)[i] >> (pms->IMSYM)[i]; pdioassert(file);
		}
	}else pms->NSYM = 0;
	// Ignore the rest of the file
}

void printParams( cfile* pms )
{
	std::cout<<pms->SYSNAME<<" \t\tSYSNAME\n";
	std::cout<<pms->LABEL<<" \t\t\tLABEL\n";
	if( pms->TYPE != -1 )
	{
		std::cout<<pms->TYPE<<" "<<(char)(pms->CPType)<<pms->CP<<" "<<pms->NPARX<<" ";
		for( int i=0; i<pms->NPARX; i++ ) std::cout<<(char)(pms->PARXType)[i]<<(pms->PARX)[i]<<" ";
		std::cout<<"\t\tTYPE, CP, NPARX, PARX[NPARX]\n";
	}else
	{
		std::cout<<pms->TYPE<<" "<<pms->CP<<" "<<pms->SWITCH<<" "<<pms->NEQN<<" ";
		for( int i=0; i<pms->NEQN; i++ ) std::cout<<"E"<<(pms->EQN)[i]<<" ";
		std::cout<<pms->NVAR<<" ";
		for( int i=0; i<pms->NVAR; i++ ) std::cout<<(char)(pms->VARType)[i]<<(pms->VAR)[i]<<" ";
		std::cout<<"\tTYPE, CP, SWITCH, NEQN, EQN[NEQN], NVAR, VAR[NVAR]\n";
	}
	std::cout<<pms->NINT<<" "<<pms->NDEG<<" "<<pms->NMUL<<" "<<pms->STAB<<" "<<pms->NMAT<<" \t\tNINT, NDEG, NMUL, STAB, NMAT\n";
	std::cout<<pms->NINT1<<" "<<pms->NINT2<<" "<<pms->NDEG1<<" "<<pms->NDEG2<<" \t\tNINT1, NINT2, NDEG1, NDEG2\n";
	std::cout<<pms->STEPS<<" "<<pms->CPMIN<<" "<<pms->CPMAX<<" \t\tSTEPS, CPMIN, CPMAX\n";
	std::cout<<pms->DS<<" "<<pms->DSMIN<<" "<<pms->DSMAX<<" "<<pms->DSSTART<<" \tDS, DSMIN, DSMAX, DSSTART \n";
	std::cout<<pms->EPSC<<" "<<pms->EPSR<<" "<<pms->EPSS<<" \tEPSC, EPSR, EPSS\n";
	std::cout<<pms->NITC<<" "<<pms->NITR<<" "<<pms->NITS<<" \t\tNITC, NITR, NITS\n";
	if( pms->NSYM != 0 )
	{
		std::cout<<pms->NSYM<<" ";
		for( int i=0; i<pms->NSYM; i++ ) std::cout<<pms->RESYM[i]<<" "<<pms->IMSYM[i]<<" ";
		std::cout<<"\t\tNSYM";
		for( int i=0; i<pms->NSYM; i++ ) std::cout<<", R"<<i<<", I"<<i;
		std::cout<<"\n";
	}
}

int initEqnVar( System& sys, cfile* params,
                Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN )
{
	// initializing the equations and variables
	if( params->EXTSYS )
	{
		eqn.Init(params->NEQN);
		var.Init(params->NVAR);
		for( int i = 0; i < params->NEQN; i++ )
		{
			eqn(i) = (Eqn)(params->EQN)[i];
			if( (params->VARType)[i] == 'S' ) var(i) = (Var)(params->VAR)[i];
			else if( (params->VARType)[i] == 'P' ) var(i) = (Var)( VarPAR0 + (params->VAR)[i] );
			else if( (params->VARType)[i] == 'I' ) var(i) = (Var)( VarPAR0 + sys.npar() + (params->VAR)[i] );
		}
	}else
	{
		for( int i = 0; i < params->NPARX; i++ )
		{
			if( (params->PARXType)[i] == 'I' ) (params->PARX)[i] += sys.npar();
		}
		params->SWITCH = PtToEqnVar( eqn, var, (PtType)params->TYPE, params->NPARX, params->PARX, sys.npar() );
	}
	// initializing CP
	if( params->CPType == 'I' ) params->CP += sys.npar();
	
	// checking whether it is an autonomous problem or not
	bool aut = false;
	bool phaseRot = false;
	testFN = EqnNone;
	int  testFN_idx = -1;
	for( int i = 0; i < eqn.Size(); i++ ) 
	{
		if( eqn(i) == EqnPhase ) aut = true;
		if( eqn(i) == EqnPhaseRot ) phaseRot = true;
		if( (eqn(i) == EqnTFPD)||
		    (eqn(i) == EqnTFLP)||
		    (eqn(i) == EqnTFLPAUT)||
		    (eqn(i) == EqnTFLPAUTROT)||
		    (eqn(i) == EqnTFCPLX_RE) )
		{
			if( testFN_idx == -1 ) { testFN = eqn(i); testFN_idx = i; }
			else PDError(-1);
		}
	}
	
	// setting up for refinement
	if( aut )
	{
		if( phaseRot )
		{
			// std::cout<<"Phase and PhaseRot\n";
			eqn_refine.Init(3);
			var_refine.Init(3);
			eqn_refine(0) = EqnSol; eqn_refine(1) = EqnPhase;          eqn_refine(2) = EqnPhaseRot;
			var_refine(0) = VarSol; var_refine(1) = var(var.Size()-2); var_refine(2) = var(var.Size()-1);
		}else
		{
			// std::cout<<"Phase\n";
			eqn_refine.Init(2);
			var_refine.Init(2);
			eqn_refine(0) = EqnSol; eqn_refine(1) = EqnPhase;
			var_refine(0) = VarSol; var_refine(1) = var(var.Size()-1);
		}
	}else
	{
		eqn_refine.Init(1);
		var_refine.Init(1);
		eqn_refine(0) = EqnSol;
		var_refine(0) = VarSol;
	}
	
	if( params->SWITCH == TFHBSwitch )
		PtToEqnVar( eqn_refine, var_refine, SolTF, 0, params->PARX, sys.npar() ); // NPARX == 0
	
	// Here, we set up the branch switching.
	// We suppose that if there is a switch we use one parameter continuation afterwards
	// without using characteristic matrices. This means that we can switch on the characteristic matrix,
	// include the equation for the eigenvector norm before the other equations and
	// add CP to the variables as a normal parameter.
	Eqn eqn_temp;
	switch( params->SWITCH )
	{
/// with TEST FUNCTIONALS
		case TFBRSwitch:
			if( aut ) eqn_temp = EqnTFLPAUT; else eqn_temp = EqnTFLP;
			goto tfskip;
		case TFPDSwitch:
			eqn_temp = EqnTFPD;
			goto tfskip;
		tfskip:
			eqn_start.Init( eqn_refine.Size() + 1 );
			var_start.Init( var_refine.Size() + 1 );
			eqn_start(0) = eqn_refine(0);
			var_start(0) = var_refine(0);
			eqn_start(1) = eqn_temp;
			var_start( var_refine.Size() ) = (Var)(VarPAR0 + params->CP);
			for( int i = 1; i < eqn_refine.Size(); i++ )
			{
				eqn_start(i+1) = eqn_refine(i);
				var_start(i) = var_refine(i);
			}
			testFN = eqn_temp;
			break;
		case TFTRSwitch:
		case TFHBSwitch:
			eqn_start.Init( eqn_refine.Size() + 2 );
			var_start.Init( var_refine.Size() + 2 );
			eqn_start(0) = eqn_refine(0);
			var_start(0) = var_refine(0);
			eqn_start(1) = EqnTFCPLX_RE;
			eqn_start(2) = EqnTFCPLX_IM;
			var_start(1) = (Var)(VarPAR0 + sys.npar() + ParAngle); // CH
			var_start( var_refine.Size() + 1 ) = (Var)(VarPAR0 + params->CP);
			for( int i = 1; i < eqn_refine.Size(); i++ )
			{
				eqn_start(i+2) = eqn_refine(i);
				var_start(i+1) = var_refine(i);
			}
			testFN = EqnTFCPLX_RE;
			break;
		default:
			eqn_start.Init( eqn.Size() );
			var_start.Init( var.Size() );
			eqn_start = eqn;
			var_start = var;
			break;
	}
	if( !aut ) return 0; else if( !phaseRot ) return 1; else return 2;
}

int main( int argc, const char** argv )
{
	// parameters
	cfile* params = 0;
	const char*  constFile = 0;
	const char*  outFile = 0;
	const char*  inFile = 0;
	const char*  branchFile = 0;
	
	// argument parsing
	for( int acnt = 1; acnt < argc;  acnt++ )
	{
		if( argv[acnt][0] == '-' )
		{
			switch( argv[acnt][1] )
			{
				case 'c':
					constFile = argv[++acnt];
					params = new cfile;
					getParams( params, constFile );
					printParams( params );
					break;
				case 'i':
					inFile = argv[++acnt];
					break;
				case 'o':
					outFile = argv[++acnt];
					break;
				case 'b':
					branchFile = argv[++acnt];
					break;
				case 'v':
					std::cout<<"This is "<<PACKAGE_NAME<<" version "<<PACKAGE_VERSION<<" ("<<PKG_REV<<")\n";
					exit(0);
					break;
				default:
					std::cout<<"Unexected CL argument.\n"; 
					PDError(-1);
					break;
			}
		}else
		{
			std::cout<<"Unexected CL argument.\n"; 
			PDError(-1);
		}
	}
	

	if( params == 0 )
	{
		std::cout<<"Error: Missing constants file.\n";
		exit(1);
	}
	if( (inFile == 0) && (params->LABEL != 0) )
	{
		std::cout<<"Error: Missing input file.\n";
		exit(1);
	}
	if( outFile == 0)
	{
		outFile ="out.pdde";
		std::cout<<"Warning: Missing output file, using \""<<outFile<<"\" instead.\n";
	}
	if( branchFile == 0)
	{
		branchFile ="branch";
		std::cout<<"Warning: Missing branch file, using \""<<branchFile<<"\" instead.\n";
	}
	
	// **********************************************************************************************************
	
	System sys( params->SYSNAME );

	Vector par(sys.npar()+ParEnd);
	std::ofstream ff( branchFile );
	std::ofstream out( outFile );
	out<<std::scientific;
	out.precision(12);
	ff<<std::scientific;
	ff.precision(12);
	
	//-----------------------------------------------------------------------------------------------------------
	//
	// Initialize the system definition for the three cases of
	// a) refinement b) branch switching c) continuation
	//
	//-----------------------------------------------------------------------------------------------------------
	
	Array1D<Eqn> eqn; // used for continuation
	Array1D<Var> var; // used for continuation
	Array1D<Eqn> eqn_refine;
	Array1D<Var> var_refine;
	Array1D<Eqn> eqn_start;
	Array1D<Var> var_start;
	Eqn          testFN;
	
	int trivial = initEqnVar( sys, params, eqn, var, eqn_refine, var_refine, eqn_start, var_start, testFN );
	const int npar = sys.npar();
	
// 	for( int i=0; i<eqn.Size(); i++ ) std::cout<<EqnToStr( eqn(i) )<<", ";
// 	std::cout<<'\n';
// 	for( int i=0; i<var.Size(); i++ ) std::cout<<VarToStr( var(i) )<<", ";
// 	std::cout<<'\n';
	
	//-----------------------------------------------------------------------------------------------------------
	//
	// END of initialization
	//
	//-----------------------------------------------------------------------------------------------------------
	
	// just a block to contain pt, which eats too much memory
	{
		Point* pt_ptr = new Point( sys, eqn_refine, var_refine, params->NINT, params->NDEG, params->NMUL, params->NMAT );
		Point& pt = *pt_ptr;
		
		pt.setContIter( params->NITC );
		pt.setRefIter( params->NITR );
		pt.setStartIter( params->NITS );
		pt.setRefEps( params->EPSR );
		pt.setContEps( params->EPSC );
		pt.setStartEps( params->EPSS );
		pt.setCont( params->CP );
		if( params->NSYM != 0 ) pt.setSym( params->NSYM, params->RESYM, params->IMSYM );
		
		std::cout<<std::scientific;
		std::cout.precision(6);
		
		// load the initial guess
		if( params->LABEL != 0 )
		{
			std::ifstream istr( inFile );
			for( int i=0; i<params->LABEL-1; i++ )
			{
				pt.ReadNull( istr );
			}
			pt.Read( istr, false );
		}
		
		pt.Refine();
		if( testFN != EqnNone )
		{
			std::cout<<"\n--- Finding the bifurcation point (TF) ---\n";
			pt.Reset( eqn_start, var_start );
			pt.setCont( params->CP );
			pt.StartTF( testFN ); // it only computes the characteristic multiplier refines the solution
		}
		
		// start the continuation!
		if( params->SWITCH != TFTRSwitch )
		{
			std::cout<<"\n--- Starting the continuation ---\n";
			
			for( int j = 0; j < par.Size(); j++ ) par(j) = pt.getPar()(j);
			//
			std::cout<<"\nLABEL\t"<<"   NORM\t\t"<<parType( npar, params->CP )<<parNum( npar, params->CP )<<"\t";
			parNamePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<parType( npar, (params->PARX)[j] )<<parNum( npar, (params->PARX)[j] )<<"\t";
			std::cout<<"\n";
			//
			std::cout<<"  "<<0<<"\t"<<pt.Norm()<<"\t"<<par(params->CP);
			parValuePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<par((params->PARX)[j]);
			std::cout<<"\n";
			
			// making tangents
			if( params->SWITCH == TFPDSwitch )
			{
				std::cout<<"\nSwitching to the period two branch (TF).\n";
				pt.SwitchTFPD( params->DSSTART );
				pt.setCont( params->CP );
			}
			else if( params->SWITCH == TFBRSwitch )
			{
				std::cout<<"\nSwitching to the other branch (TF).\n";
				pt.SwitchTFLP( params->DSSTART );
				pt.setCont( params->CP );
			}
			else if( params->SWITCH == TFHBSwitch )
			{
				std::cout<<"\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
				pt.SwitchTFHB( params->DSSTART );
				pt.setCont( params->CP );
			}
			else
			{
				std::cout<<"\nFinding the tangent.\n";
				pt.setCont( params->CP );
				pt.Tangent();
			}
			pt.Reset( eqn, var );
			
			std::cout<<'\n';
			int ustab = 0, ustabprev = 0;
			double norm = 0.0;
			const int ithist = 5;
			Array1D<int> it( ithist );
			for( int i=0; i < it.Size(); i++ ) it(i) = 3;
			int itpos = 0;
			double ds = params->DS;
			for( int i = 0; i < params->STEPS; i++ ) // 35
			{
				if( i % 24 == 0 )
				{
					std::cout<<"LABEL\t"<<"   NORM\t\t"<<(char)params->CPType<<params->CP<<"\t";
					parNamePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<parType( npar, (params->PARX)[j] )<<parNum( npar, (params->PARX)[j] )<<"\t";
					std::cout<<"\tUSTAB\tIT\n";
				}
				itpos = (itpos+1) % ithist;
				//
				it( itpos ) = pt.Continue( ds );
				//
				if( params->STAB != 0) pt.Stability();
				ustabprev = ustab;
				if( trivial == 0 ) ustab = pt.UStab(); else if( trivial == 1 ) ustab = pt.UStabAUT(); else ustab = pt.UStabAUTRot();
				for( int j = 0; j < par.Size(); j++ ) par(j) = pt.getPar()(j);
				norm = pt.Norm();

				// console output
				std::cout<<"  "<<i+1<<"\t"<<norm<<"\t"<<par(params->CP);
				parValuePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<par((params->PARX)[j]);
				std::cout<<"\t  "<<ustab<<"\t"<<it(itpos)+1;
				if( i != 0  && ustab != ustabprev )
				{
					PtType bif = SolTF;
					if( trivial == 0 ) bif = pt.testBif(); else if( trivial == 1 ) bif = pt.testBifAUT(); else pt.testBifAUTRot();
					switch( bif )
					{
						case BifTFLP:
						case BifTFAUTLP:
							std::cout<<"  LP";
							break;
						case BifTFPD:
						case BifTFAUTPD:
							std::cout<<"  PD";
							break;
						case BifTFNS:
						case BifTFAUTNS:
							std::cout<<"  NS";
							break;
						default:
							std::cout<<"  ??";
							break;
					}
				}
				std::cout<<"\n";
				
				// file output
				pt.Write( out );
				
				// branch output
				for( int j = 0; j < npar; j++ ) ff<<par(j)<<"\t";
				ff<<"\t"<<norm<<"\t"<<pt.NormMX()<<"\t"<<ustab<<"\n";
				ff.flush();
				int itc = it( itpos );
				if( (itc > 3)&&(fabs(ds)/1.414 > params->DSMIN)&&(fabs(ds)/1.414 < params->DSMAX) ) ds /= 1.414;
				if( (itc > 5)&&(fabs(ds)/2.0 > params->DSMIN)&&(fabs(ds)/2.0 < params->DSMAX) ) ds /= 2.0;
				bool decr = true;
				for( int l=0; l < it.Size(); l++ ) if( it(l) > 3 ) decr = false;
				if( decr&&(fabs(ds)*1.414 > params->DSMIN)&&(fabs(ds)*1.414 < params->DSMAX) ) ds *= 1.414;
				if( (itc >= params->NITC)&&(fabs(ds)/2.0 < params->DSMIN) )
				{
					PDError(1);
				}
				std::cout.flush();
			}
		}else
		{
			std::cout<<"ENTERING THE TORUS CODE!\n";
			// construct the initial torus
			double alpha;
			Vector Sol(pt.getSol().Size());
			Vector TRe(pt.getSol().Size()), TIm(pt.getSol().Size());

			// making the mesh for the conversion
			Vector meshint(params->NINT1 + 1), meshdeg( params->NDEG1 + 1 );
			for( int i=0; i<meshint.Size(); i++ ) meshint(i) = (double)i/(params->NINT1);
			for( int i=0; i<meshdeg.Size(); i++ ) meshdeg(i) = (double)i/(params->NDEG1);

			// getting the sol and tangents
			pt.SwitchTRSol( Sol, meshint, meshdeg );
			pt.SwitchTFTRTan( TRe, TIm, alpha, meshint, meshdeg );

			// getting the parameters
			for( int j = 0; j < par.Size(); j++ ) par(j) = pt.getPar()(j);

			// destroy point, construct PointTR
			delete pt_ptr;
			PointTR pttr( sys, eqn, var, params->NDEG1,params->NDEG2, params->NINT1, params->NINT2 );

			// construct the solution tangent from the eigenvectors
			// these next three functions could be only one
			pttr.ImportSol( Sol );
			pttr.ImportTan( TRe, TIm, alpha );
			pttr.Start( params->DSSTART );
			pttr.setPar( par );
			
			pttr.setRho( alpha/(2.0*M_PI) );
			pttr.setCont( params->CP );

			double ds = params->DS;
			std::ostringstream fdata, fidx;
			for( int i = 0; i < params->STEPS; i++ )
			{
				pttr.Continue( ds, false );
				
				// write out the results
				for( int j = 0; j < npar; j++ ) std::cout<<par(j)<<"\t";
				std::cout<<std::endl;
				for( int j = 0; j < npar; j++ ) par(j) = pttr.getPar()(j);
				for( int j = 0; j < npar; j++ ) ff<<par(j)<<"\t";
				ff<<pttr.Norm()<<"\n";
				ff.flush();
				std::ostringstream fdata, fidx;
				fdata << "sol-" << i << ".dat";
				fidx << "sol-" << i << ".idx";
				pttr.SaveSol( fdata.str().c_str(), fidx.str().c_str() );
			}
		}
		// **********************************************************************************************************
	}

	delete params;
}
