#include <iostream>
#include <exception>
#include <fstream>
#include "pointtype.h"
#include "pderror.h"

#define MAX_PARX 8
#define MAX_EQVAR 8

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
};

void getParams( cfile* pms, const char* filename )
{
	std::ifstream file;
	file.exceptions ( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
	try
	{
		char pt;
		int ival;
		file.open( filename );
		file >> pms->SYSNAME;
		while( file.get() != '\n' );
		
		file >> pms->LABEL;
		while( file.get() != '\n' );
		
		file >> pms->TYPE >> pt >> pms->CP; pms->CPType = pt;
		if( (pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: CP: Bad parameter type."; PDError(-1); }
		pms->EXTSYS = (pms->TYPE == -1);
		if( pms->EXTSYS )
		{
			// we need extended specifications
			file >> ival;
			pms->SWITCH = (BranchSW)ival;
			file >> pms->NEQN;
			for( int i = 0; i < pms->NEQN; i++ )
			{
				file >> pt >> (pms->EQN)[i]; 
				if( pt != 'E' ) { std::cout<<"Error: EQN: Bad equation type."; PDError(-1); }
			}
			file >> pms->NVAR;
			if( pms->NEQN != pms->NVAR ) { std::cout<<"Error: NVAR: Number of equations and variables are not the same."; PDError(-1); }
			for( int i = 0; i < pms->NVAR; i++ )
			{
				file >> pt >> (pms->VAR)[i]; (pms->VARType)[i] = pt;
				if( (pt != 'S')&&(pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: VAR: Bad variable/parameter type."; PDError(-1); }
			}
			while( file.get() != '\n' );
		}else
		{
			file >> pms->NPARX;
			for( int i = 0; i < pms->NPARX; i++ )
			{
				file >> pt >> (pms->PARX)[i]; (pms->PARXType)[i] = pt;
				if( (pt != 'P')&&(pt != 'I') ) { std::cout<<"Error: PARX: Bad parameter type."; PDError(-1); }
			}
			while( file.get() != '\n' );
		}
		
		file >> pms->NINT >> pms->NDEG >> pms->NMUL;
		pms->STAB = 1; pms->NMAT = 1;
		while( file.get() != '\n' );
		
		file >> pms->NINT1 >> pms->NINT2 >> pms->NDEG1 >> pms->NDEG2;
		while( file.get() != '\n' );
		
		file >> pms->STEPS >> pms->CPMIN >> pms->CPMAX;
		while( file.get() != '\n' );
		
		file >> pms->DS >> pms->DSMIN >> pms->DSMAX;
		pms->DSSTART = pms->DS;
		while( file.get() != '\n' );
		
		file >> pms->NITC >> pms->EPSC >> pms->EPSR >> pms->EPSS;
		pms->NITR = pms->NITC;
		pms->NITS = pms->NITC;
		// Ignore the rest of the file
	}
	catch (std::ifstream::failure e)
	{
		std::cout << "Exception: opening/reading the constants file.";
		exit(-1);
	}
}

void printParams( std::ofstream& file, cfile* pms )
{
	file<<pms->SYSNAME<<" \t\tSYSNAME\n";
	file<<pms->LABEL<<" \t\t\tLABEL\n";
	if( pms->TYPE != -1 )
	{
		file<<pms->TYPE<<" "<<(char)(pms->CPType)<<pms->CP<<" "<<pms->NPARX<<" ";
		for( int i=0; i<pms->NPARX; i++ ) file<<(char)(pms->PARXType)[i]<<(pms->PARX)[i]<<" ";
		file<<"\t\tTYPE, CP, NPARX, PARX[NPARX]\n";
	}else
	{
		file<<pms->TYPE<<" "<<(char)(pms->CPType)<<pms->CP<<" "<<pms->SWITCH<<" "<<pms->NEQN<<" ";
		for( int i=0; i<pms->NEQN; i++ ) file<<"E"<<(pms->EQN)[i]<<" ";
		file<<pms->NVAR<<" ";
		for( int i=0; i<pms->NVAR; i++ ) file<<(char)(pms->VARType)[i]<<(pms->VAR)[i]<<" ";
		file<<"\t\tTYPE, CP, SWITCH, NEQN, EQN[NEQN], NVAR, VAR[NVAR]\n";
	}
	file<<pms->NINT<<" "<<pms->NDEG<<" "<<pms->NMUL<<" "<<pms->STAB<<" "<<pms->NMAT<<" \t\tNINT, NDEG, NMUL, STAB, NMAT\n";
	file<<pms->NINT1<<" "<<pms->NINT2<<" "<<pms->NDEG1<<" "<<pms->NDEG2<<" \t\tNINT1, NINT2, NDEG1, NDEG2\n";
	file<<pms->STEPS<<" "<<pms->CPMIN<<" "<<pms->CPMAX<<" \t\tSTEPS, CPMIN, CPMAX\n";
	file<<pms->DS<<" "<<pms->DSMIN<<" "<<pms->DSMAX<<" "<<pms->DSSTART<<" \tDS, DSMIN, DSMAX, DSSTART \n";
	file<<pms->EPSC<<" "<<pms->EPSR<<" "<<pms->EPSS<<" \tEPSC, EPSR, EPSS\n";
	file<<pms->NITC<<" "<<pms->NITR<<" "<<pms->NITS<<" \tNITC, NITR, NITS\n";
}

int main( int argc, const char** argv )
{
	struct cfile params;
	if( argc != 2 ) { PDError(-1); }
	else
	{
		getParams( &params, argv[1] );
		std::ofstream fp;
		fp.open( argv[1] );
		printParams( fp, &params );
	}
	return 0;
}
