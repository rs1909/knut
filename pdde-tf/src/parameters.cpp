#include "parameters.h" // it includes <string>

#include "pderror.h"
#include "matrix.h"
#include "system.h"
#include <iostream>
#include <fstream>

static inline void pdioassert( std::istream& is )
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
		P_MESSAGE("");
	}
}

void ConstFile::getParams( const std::string& filename )
{
	std::ifstream file;
	char pt;
	int ival;
	
	file.open( filename.c_str() );
	if(!file){ std::cout<<"Cannot open "<<filename<<"\n"; P_MESSAGE(""); }
	
	file >> SYSNAME; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> LABEL; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> TYPE >> pt >> CP; CPType = pt; pdioassert(file);
	if( (pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: CP: Bad parameter type."; P_MESSAGE(""); }
	
	EXTSYS = (TYPE == -1);
	if( EXTSYS )
	{
		// we need extended specifications
		file >> ival;
		SWITCH = (BranchSW)ival;
		file >> NEQN; pdioassert(file);
		if( NEQN > MAX_EQVAR ) { std::cout<<"Error: NEQN: System limit has been reached."; P_MESSAGE(""); }
		for( int i = 0; i < NEQN; i++ )
		{
			file >> pt >> (EQN)[i]; pdioassert(file);
			if( pt != 'E' ) { std::cout<<"Error: EQN: Bad equation type."; P_MESSAGE(""); }
		}
		file >> NVAR; pdioassert(file);
		if( NEQN != NVAR ) { std::cout<<"Error: NVAR: Number of equations and variables are not the same."; P_MESSAGE(""); }
		for( int i = 0; i < NVAR; i++ )
		{
			file >> pt >> (VAR)[i]; (VARType)[i] = pt; pdioassert(file);
			if( (pt != 'S')&&(pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: VAR: Bad variable/parameter type."; P_MESSAGE(""); }
		}
		while( file.get() != '\n' ); pdioassert(file);
	}else
	{
		file >> NPARX; pdioassert(file);
		if( NPARX > MAX_PARX ) { std::cout<<"Error: NPARX: System limit has been reached."; P_MESSAGE(""); }
		for( int i = 0; i < NPARX; i++ )
		{
			file >> pt >> (PARX)[i]; (PARXType)[i] = pt; pdioassert(file);
			if( (pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: PARX: Bad parameter type."; P_MESSAGE(""); }
		}
		while( file.get() != '\n' ); pdioassert(file);
	}
	
	file >> NINT >> NDEG >> NMUL >> STAB >> NMAT; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> NINT1 >> NINT2 >> NDEG1 >> NDEG2; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> STEPS >> CPMIN >> CPMAX; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> DS >> DSMIN >> DSMAX >> DSSTART; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> EPSC >> EPSR >> EPSS; pdioassert(file);
	while( file.get() != '\n' ); pdioassert(file);
	
	file >> NITC >> NITR >> NITS; pdioassert(file);
	// Check whether we need to load NSYM
	bool isSYM = false;
	if( EXTSYS)
	{
		for( int i = 2; i < NEQN; i++ )
		{
			if( (Eqn)(EQN)[i] == EqnPhaseRot ) isSYM = true;
		}
	}
	if( isSYM )
	{
		// first go to the next line
		while( file.get() != '\n' ); pdioassert(file);
		file >> NSYM; pdioassert(file);
		if( NSYM > MAX_SYM ) { std::cout<<"Error: NSYM: System limit has been reached."; P_MESSAGE(""); }
		for( int i = 0; i < NSYM; i++ )
		{
			file >> (RESYM)[i] >> (IMSYM)[i]; pdioassert(file);
		}
	}else NSYM = 0;
	// Ignore the rest of the file
}

void ConstFile::printParams( std::ostream& file )
{
	file<<SYSNAME<<" \t\tSYSNAME\n";
	file<<LABEL<<" \t\t\tLABEL\n";
	if( TYPE != -1 )
	{
		file<<TYPE<<" "<<(char)(CPType)<<CP<<" "<<NPARX<<" ";
		for( int i=0; i<NPARX; i++ ) file<<(char)(PARXType)[i]<<(PARX)[i]<<" ";
		file<<"\t\tTYPE, CP, NPARX, PARX[NPARX]\n";
	}else
	{
		file<<TYPE<<" "<<(char)(CPType)<<CP<<" "<<SWITCH<<" "<<NEQN<<" ";
		for( int i=0; i<NEQN; i++ ) file<<"E"<<(EQN)[i]<<" ";
		file<<NVAR<<" ";
		for( int i=0; i<NVAR; i++ ) file<<(char)(VARType)[i]<<(VAR)[i]<<" ";
		file<<"\tTYPE, CP, SWITCH, NEQN, EQN[NEQN], NVAR, VAR[NVAR]\n";
	}
	file<<NINT<<" "<<NDEG<<" "<<NMUL<<" "<<STAB<<" "<<NMAT<<" \t\tNINT, NDEG, NMUL, STAB, NMAT\n";
	file<<NINT1<<" "<<NINT2<<" "<<NDEG1<<" "<<NDEG2<<" \t\tNINT1, NINT2, NDEG1, NDEG2\n";
	file<<STEPS<<" "<<CPMIN<<" "<<CPMAX<<" \t\tSTEPS, CPMIN, CPMAX\n";
	file<<DS<<" "<<DSMIN<<" "<<DSMAX<<" "<<DSSTART<<" \tDS, DSMIN, DSMAX, DSSTART \n";
	file<<EPSC<<" "<<EPSR<<" "<<EPSS<<" \tEPSC, EPSR, EPSS\n";
	file<<NITC<<" "<<NITR<<" "<<NITS<<" \t\tNITC, NITR, NITS\n";
	if( NSYM != 0 )
	{
		file<<NSYM<<" ";
		for( int i=0; i<NSYM; i++ ) file<<RESYM[i]<<" "<<IMSYM[i]<<" ";
		file<<"\t\tNSYM";
		for( int i=0; i<NSYM; i++ ) file<<", R"<<i<<", I"<<i;
		file<<"\n";
	}
}

int ConstFile::toEqnVar( System& sys,
                Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN )
{
	// initializing the equations and variables
	if( EXTSYS )
	{
		eqn.Init(NEQN);
		var.Init(NVAR);
		for( int i = 0; i < NEQN; i++ )
		{
			eqn(i) = (Eqn)(EQN)[i];
			if( (VARType)[i] == 'S' ) var(i) = (Var)(VAR)[i];
			else if( (VARType)[i] == 'P' ) var(i) = (Var)( VarPAR0 + (VAR)[i] );
			else if( (VARType)[i] == 'I' ) var(i) = (Var)( VarPAR0 + sys.npar() + (VAR)[i] );
		}
	}else
	{
		Array1D<Var> L_PARX(NPARX);
		for( int i = 0; i < NPARX; ++i )
		{
			if( (PARXType)[i] == 'P' ) L_PARX(i) = (Var)(VarPAR0 + PARX[i]);
			else if( (PARXType)[i] == 'I' ) L_PARX(i) = (Var)(VarPAR0 + PARX[i]+sys.npar());
		}
		SWITCH = PtToEqnVar( eqn, var, (PtType)TYPE, L_PARX, sys.npar() );
	}
	// initializing CP
	if( CPType == 'I' ) CP += sys.npar();
	
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
			else P_MESSAGE("");
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
	
	if( SWITCH == TFHBSwitch )
	  { Array1D<Var> d(0); PtToEqnVar( eqn_refine, var_refine, SolTF, d, sys.npar() ); }// NPARX == 0
	
	// Here, we set up the branch switching.
	// We suppose that if there is a switch we use one parameter continuation afterwards
	// without using characteristic matrices. This means that we can switch on the characteristic matrix,
	// include the equation for the eigenvector norm before the other equations and
	// add CP to the variables as a normal parameter.
	Eqn eqn_temp;
	switch( SWITCH )
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
			var_start( var_refine.Size() ) = (Var)(VarPAR0 + CP);
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
			var_start( var_refine.Size() + 1 ) = (Var)(VarPAR0 + CP);
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
