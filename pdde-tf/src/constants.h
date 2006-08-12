#ifndef CONSTANTS_H
#define CONSTANTS_H

// Includes standard library
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>

// Includes local
#include "pointtype.h"
#include "system.h"

class toNumber
{
 public:
	toNumber( int i ) { int n=snprintf(num,31,"%d", i); }
	toNumber( double d ) { int n=snprintf(num,31,"%lf", d); }
	const char* c_str() const { return num; }
 private:
	char num[32];
};

template<class KEY> class genMap
{
 public:
	genMap() : invalid("Invalid") { }
	int indexof( const KEY& v ) const
	  { for( int i = 0; i < key.size(); ++i ) if( key[i] == v ) return i; return -1; }
	const std::string& string( int i ) const
	  { if( i >= 0 && i < key.size() ) return dsc[i]; else return invalid; }
	const KEY& getkey( int i ) const
	  { if( i >= 0 && i < key.size() ) return key[i]; else return key.at(0); }
	const std::string& map( const KEY& v ) const { return string( indexof(v) ); }
	int size() const { return key.size(); }
 protected:
	std::vector<KEY>         key;
	std::vector<std::string> dsc;
	const std::string        invalid;
};

class brswType : public genMap<BranchSW>
{
 public:
	
	brswType( )
	{
		key.push_back(NOSwitch);   dsc.push_back("No switch");
		key.push_back(TFBRSwitch); dsc.push_back("Branch");
		key.push_back(TFPDSwitch); dsc.push_back("Period doubling");
		key.push_back(TFHBSwitch); dsc.push_back("Hopf (aut)");
		key.push_back(TFTRSwitch); dsc.push_back("Torus");
	}
};

class PointType : public genMap<PtType>
{
 public:
	
	PointType( )
	{
		dsc.push_back( "User" ); key.push_back(SolUser);
		dsc.push_back( "Limit cycle" ); key.push_back(SolTF);
		dsc.push_back( "Limit point" ); key.push_back(BifTFLP);
		dsc.push_back( "Period doubling" ); key.push_back(BifTFPD);
		dsc.push_back( "Neimark-Sacker" ); key.push_back(BifTFNS);
		dsc.push_back( "Branch switch" ); key.push_back(SolTFBRSW);
		dsc.push_back( "Period doubling switch" ); key.push_back(SolTFPDSW);
		dsc.push_back( "Limit cycle (aut)" ); key.push_back(SolTFAUT);
		dsc.push_back( "Limit point (aut)" ); key.push_back(BifTFAUTLP);
		dsc.push_back( "Period doubling (aut)" ); key.push_back(BifTFAUTPD);
		dsc.push_back( "Neimark-Sacker (aut)" ); key.push_back(BifTFAUTNS);
		dsc.push_back( "Branch switch (aut)" ); key.push_back(SolTFAUTBRSW);
		dsc.push_back( "Period doubling switch (aut)" ); key.push_back(SolTFAUTPDSW);
		dsc.push_back( "Hopf switch (aut)" ); key.push_back(SolTFAUTHBSW);
		dsc.push_back( "Torus from NS" ); key.push_back(SolTor);
		dsc.push_back( "Torus from NS (aut)" ); key.push_back(SolAUTTor);
	}
};

class EqnType : public genMap<Eqn>
{
 public:
	EqnType( )
	{
		dsc.push_back( "None" ); key.push_back( EqnNone );
		dsc.push_back( "Limit cycle" ); key.push_back( EqnSol );
		dsc.push_back( "PHCND (transl)" ); key.push_back( EqnPhase );
		dsc.push_back( "PHCND (rot)" ); key.push_back( EqnPhaseRot );
		dsc.push_back( "Torus" ); key.push_back( EqnTORSol );
		dsc.push_back( "PHCND 1 (torus)" ); key.push_back( EqnTORPhase0 );
		dsc.push_back( "PHCND 2 (torus)" ); key.push_back( EqnTORPhase1 );
		dsc.push_back( "LP TF" ); key.push_back( EqnTFLP );
		dsc.push_back( "PD TF" ); key.push_back( EqnTFPD );
		dsc.push_back( "LP TF (aut)" ); key.push_back( EqnTFLPAUT );
		dsc.push_back( "LP TF (sym)" ); key.push_back( EqnTFLPAUTROT );
		dsc.push_back( "NS TF (real)" ); key.push_back( EqnTFCPLX_RE );
		dsc.push_back( "NS TF (imag)" ); key.push_back( EqnTFCPLX_IM );
	}
	char getType( int i ) const
	{
		return 'E';
	}
	
	int  getNum( int i ) const
	{
		return key[i];
	}
};

class VarType : public genMap<Var>
{
 public:
	
	VarType( int npar_, bool full_ = true ) : npar(npar_), sysvar(full_)
	{
		setPar( npar_ );
	}

	void setPar( int npar_ )
	{
		npar = npar_;
		key.clear(); dsc.clear();
		dsc.push_back( "None" ); key.push_back( VarNone );
		if( sysvar )
		{
			dsc.push_back( "Limit cycle" ); key.push_back( VarSol );
			dsc.push_back( "Torus" ); key.push_back( VarTORSol );
		}
		for( int i = 0; i < npar; ++i )
		{
			char _buf[7+1];
			int n = snprintf( _buf, 7, "%d", i );
			dsc.push_back( "P "+std::string(_buf) );
			key.push_back( (Var)(VarPAR0 + i) );
		}
		dsc.push_back( "Angle (NS)" ); key.push_back( (Var)(VarPAR0 + npar + ParAngle) );
		dsc.push_back( "Rotation num." ); key.push_back( (Var)(VarPAR0 + npar + ParRot) );
	}

	void setSysVar( bool b )
	{
		sysvar = b;
		setPar( npar );
	}

	char getType( int i ) const
	{
		if( key[i] < VarPAR0 ) return 'S'; 
		else if( key[i] < VarPAR0+npar ) return 'P';
		else return 'I';
	}
	
	int  getNum( int i ) const
	{
		if( key[i] < VarPAR0 ) return key[i];
		else if( key[i] < VarPAR0+npar ) return key[i]-VarPAR0;
		else return key[i]-VarPAR0-npar;
	}

	Var fromTypeNum( char type, int num ) const
	{
		if( type == 'S' ) return (Var)num;
		else if( type == 'P' ) return (Var)(num+VarPAR0);
		else if( type == 'I' ) return (Var)(num+VarPAR0+npar);
	}
 private:
	int npar;
	bool sysvar;
};

class NConstants
{
	public:
		
		NConstants( ) : cpMap( 0, false ), parxMap( 0, false ), varsMap( 0, true ),
			label(0), pttype(SolUser), cp(VarPAR0), branchsw(NOSwitch), neqns(0), nint(0), ndeg(0), nmul(0),
			stab(0), nmat(0), nint1(0), nint2(0), ndeg1(0), ndeg2(0), steps(0), cpMin(0.0),
			cpMax(0.0), ds(0.0), dsMin(0.0), dsMax(0.0), dsStart(0.0),
			epsC(0.0), epsR(0.0), epsS(0.0), nitC(0), nitR(0), nitS(0),
			nsym(0), npar(0), ndim(0) {  }
		
		NConstants( const NConstants& ct ) : inputFile(ct.inputFile), outputFile(ct.outputFile), sysname(ct.sysname),
			cpMap( ct.npar, false ), parxMap( ct.npar, false ), varsMap( ct.npar, true ),
			label(ct.label), pttype(ct.pttype), cp(ct.cp), branchsw(ct.branchsw), neqns(ct.neqns), nint(ct.nint), ndeg(ct.ndeg),
			nmul(ct.nmul), stab(ct.stab), nmat(ct.nmat),
			nint1(ct.nint1), nint2(ct.nint2), ndeg1(ct.ndeg1), ndeg2(ct.ndeg2), steps(ct.steps),
			cpMin(ct.cpMin), cpMax(ct.cpMax), ds(ct.ds), dsMin(ct.dsMin), dsMax(ct.dsMax), dsStart(ct.dsStart),
			epsC(ct.epsC), epsR(ct.epsR), epsS(ct.epsS), nitC(ct.nitC), nitR(ct.nitR), nitS(ct.nitS),
			nsym(ct.nsym), npar(ct.npar), ndim(ct.ndim),
			parxtype(ct.parxtype), parxnum(ct.parxnum), eqnstype(ct.eqnstype), eqnsnum(ct.eqnsnum), varstype(ct.varstype), varsnum(ct.varsnum),
			symRe(ct.symRe), symIm(ct.symIm) {  }

		int toEqnVar( System& sys,
		              Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
		              Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
		              Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN );
		inline bool inputAssert( std::istream& is );
		void loadFile(const std::string &fileName);
		void saveFile(const std::string &fileName);
		void printFile( std::ostream& file );

		const std::string& getInputFile( ) const { return inputFile; }
		const std::string& getOutputFile( ) const { return outputFile; }
		const std::string& getSysName( ) const { return sysname; }
		int  getLabel( ) const { return label; }
		PtType  getPointType( ) const { return pttype; }
		int  getPointTypeIdx( ) const { return pointMap.indexof( pttype ); }
		Var  getCp( ) const { return cp; }
		char getCpType( ) const { return cpMap.getType( cpMap.indexof( cp ) ); }
		int  getCpNum( ) const { return cpMap.getNum( cpMap.indexof( cp ) ); }
		int  getCpIdx( ) const { return cpMap.indexof(cp); }
		int  getBranchSW( ) const { return branchsw; }
		int  getBranchSWIdx( ) const { return branchswMap.indexof(branchsw); }
		int  getNEqns( ) const { return neqns; }
// 		void getParX( std::vector<Var>& p ) const { p = parx; }
		Var  getParX( int i ) const { return parxMap.fromTypeNum( parxtype[i], parxnum[i] ); }
		char getParXType( int i ) const { return parxtype[i]; }
		int  getParXNum( int i ) const { return parxnum[i]; }
// 		void getEqns( std::vector<Eqn>& e ) const { e = eqns; }
		Eqn  getEqns( int i ) const { if( eqnstype[i] == 'E' ) return (Eqn)eqnsnum[i]; else return EqnNone; }
		char getEqnsType( int i ) const { return eqnstype[i]; }
		int  getEqnsNum( int i ) const { return  eqnsnum[i]; }
// 		void getVars( std::vector<Var>& v ) const { v = vars; }
		Var  getVars( int i ) const { return varsMap.fromTypeNum( varstype[i], varsnum[i] ); }
		char getVarsType( int i ) const { return varstype[i]; }
		int  getVarsNum( int i ) const { return varsnum[i]; }
		int  getNInt( ) const { return nint; }
		int  getNDeg( ) const { return ndeg; }
		int  getNMul( ) const { return nmul; }
		bool getStab( ) const { return stab; }
		int  getNMat( ) const { return nmat; }
		int  getNInt1( ) const { return nint1; }
		int  getNInt2( ) const { return nint2; }
		int  getNDeg1( ) const { return ndeg1; }
		int  getNDeg2( ) const { return ndeg2; }
		int  getSteps( ) const { return steps; }
		double getCpMin( ) const { return cpMin; }
		double getCpMax( ) const { return cpMax; }
		double getDs( ) const { return ds; }
		double getDsMin( ) const { return dsMin; }
		double getDsMax( ) const { return dsMax; }
		double getDsStart( ) const { return dsStart; }
		double getEpsC( ) const { return epsC; }
		double getEpsR( ) const { return epsR; }
		double getEpsS( ) const { return epsS; }
		int getNItC( ) const { return nitC; }
		int getNItR( ) const { return nitR; }
		int getNItS( ) const { return nitS; }
		int getNSym( ) const { return nsym; }
		void getSymRe( std::vector<int>& re ) const { re = symRe; }
		int  getSymRe( int n ) const { return symRe[n]; }
		void getSymIm( std::vector<int>& im ) const { im = symIm; }
		int  getSymIm( int n ) const { return symIm[n]; }
		int getNPar( ) const { return npar; }
		int getNDim( ) const { return ndim; }

		const std::string& pointString( int i ) { return pointMap.string( i ); }
		const std::string& cpString( int i ) { return cpMap.string( i ); }
		const std::string& branchswString( int i ) { return branchswMap.string( i ); }
		const std::string& parxString( int i ) { return parxMap.string( i ); }
		const std::string& eqnsString( int i ) { return eqnsMap.string( i ); }
		const std::string& varsString( int i ) { return varsMap.string( i ); }

		const std::string& findPointString( PtType i ) { return pointMap.map( i ); }
		const std::string& findCpString( Var i ) { return cpMap.map( i ); }
		const std::string& findParXString( Var i ) { return parxMap.map( i ); }
		const std::string& findEqnsString( Eqn i ) { return eqnsMap.map( i ); }
		const std::string& findVarsString( Var i ) { return varsMap.map( i ); }

		int pointSize( ) { return pointMap.size(); }
		int cpSize( ) { return cpMap.size(); }
		int branchswSize( ) { return branchswMap.size(); }
		int parxSize( ) { return parxMap.size(); }
		int eqnsSize( ) { return eqnsMap.size(); }
		int varsSize( ) { return varsMap.size(); }
		
		// set routines
		virtual void setInputFileText( const std::string& str ) { inputFile = str; }
		virtual void setOutputFileText( const std::string& str ) { outputFile = str; }
		virtual void setSysNameText( const std::string& str )
		{
			sysname = str;
			System* sys = 0;
			try{ sys = new System( sysname ); }
			catch( pddeException ex )
			{
				npar = 0; ndim = 0;
				delete sys;
				throw( ex );
				return;
			}
			npar = sys->npar();
			ndim = sys->ndim();
			delete sys;
// 			std::cout<<"NDIM "<<ndim<<"\n";
			cpMap.setPar( npar );
			parxMap.setPar( npar );
			varsMap.setPar( npar );
		}
		virtual void setLabel( int i ) { label = i; }
		virtual void setPointType( PtType p ) { pttype = p; }
		virtual void setPointTypeIdx( int p ) { pttype = pointMap.getkey( p ); }
		virtual void setCp( Var v ) { cp = v; }
		virtual void setCp( char tp, int n ) { cp = cpMap.fromTypeNum( tp, n ); }
		virtual void setCpIdx( int i ) { cp = cpMap.getkey( i ); }
		virtual void setBranchSW( BranchSW i ) { branchsw = i; }
		virtual void setBranchSWIdx( int i ) { setBranchSW(branchswMap.getkey( i )); }
		virtual void setNEqns( int i )
		{
			parxnum.resize(i);
			eqnsnum.resize(i);
			varsnum.resize(i);
			parxtype.resize(i);
			eqnstype.resize(i);
			varstype.resize(i);
			for( int j = neqns; j < i; ++j )
			{
				parxtype[j] = 'S';
				parxnum[j] = VarNone;
				eqnstype[j] = 'E';
				eqnsnum[j] = EqnNone;
				varstype[j] = 'S';
				varsnum[j] = VarNone;
			}
			neqns = i;
		}
		virtual void setParX( int i, char tp, int n ) { parxtype[i] = tp; parxnum[i] = n; }
		virtual void setParXIdx( int i, int p ) { parxtype[i] = parxMap.getType(p); parxnum[i] = parxMap.getNum( p ); }
		virtual void setEqns( int i, char tp, int n ) { eqnstype[i] = tp; eqnsnum[i] = n; }
		virtual void setEqnsIdx( int i, int e ) { eqnstype[i] = eqnsMap.getType(e); eqnsnum[i] = eqnsMap.getNum( e ); }
		virtual void setVars( int i, char tp, int n ){ varstype[i] = tp; varsnum[i] = n; }
		virtual void setVarsIdx( int i, int v ) { varstype[i] = varsMap.getType(v); varsnum[i] = varsMap.getNum( v ); }
		virtual void setNInt( int i ) { nint = i; }
		virtual void setNDeg( int i ) { ndeg = i; }
		virtual void setNMul( int i ) { nmul = i; }
		virtual void setStab( bool s ) { stab = s; }
		virtual void setNMat( int i ) { nmat = i; }
		virtual void setNInt1( int i ) { nint1 = i; }
		virtual void setNInt2( int i ) { nint2 = i; }
		virtual void setNDeg1( int i ) { ndeg1 = i; }
		virtual void setNDeg2( int i ) { ndeg2 = i; }
		virtual void setSteps( int i ) { steps = i; }
		virtual void setCpMin( double d ) { cpMin = d; }
// 		virtual void setCpMin( const std::string& d ) { setCpMin(atof(d.c_str())); }
		virtual void setCpMax( double d ) { cpMax = d; }
// 		virtual void setCpMax( const std::string& d ) { setCpMax(atof(d.c_str())); }
		virtual void setDs( double d ) { ds = d; }
// 		virtual void setDs( const std::string& d ) { setDs(atof(d.c_str())); }
		virtual void setDsMin( double d ) { dsMin = d; }
// 		virtual void setDsMin( const std::string& d ) { setDsMin(atof(d.c_str())); }
		virtual void setDsMax( double d ) { dsMax = d; }
// 		virtual void setDsMax( const std::string& d ) { setDsMax(atof(d.c_str())); }
		virtual void setDsStart( double d ) { dsStart = d; }
// 		virtual void setDsStart( const std::string& d ) { setDsStart(atof(d.c_str())); }
		virtual void setEpsC( double d ) { epsC = d; }
// 		virtual void setEpsC( const std::string& d  ) { setEpsC(atof(d.c_str())); }
		virtual void setEpsR( double d ) { epsR = d; }
// 		virtual void setEpsR( const std::string& d ) { setEpsR(atof(d.c_str())); }
		virtual void setEpsS( double d ) { epsS = d; }
// 		virtual void setEpsS( const std::string& d ) { setEpsS(atof(d.c_str())); }
		virtual void setNItC( int i ) { nitC = i; }
		virtual void setNItR( int i ) { nitR = i; }
		virtual void setNItS( int i ) { nitS = i; }
		virtual void setNSym( int i )
		{
			symRe.resize(i);
			symIm.resize(i);
			for( int j = nsym; j < i; ++j )
			{
				symRe[j] = 0;
				symIm[j] = 0;
			}
			nsym = i;
		}
		virtual void setSymRe( int i, int v ) { symRe[i] = v; }
		virtual void setSymIm( int i, int v ) { symIm[i] = v; }
	
	private:

		std::string  inputFile;
		std::string  outputFile;
		std::string  sysname;
		int          label;
		PtType       pttype;
		PointType    pointMap;
		Var          cp;
		VarType      cpMap;
		BranchSW     branchsw;
		brswType     branchswMap;
		int          neqns;
// 		std::vector<Var> parx;
		std::vector<char> parxtype;
		std::vector<int>  parxnum;
		VarType      parxMap;
// 		std::vector<Eqn> eqns;
		std::vector<char> eqnstype;
		std::vector<int>  eqnsnum;
		EqnType      eqnsMap;
// 		std::vector<Var> vars;
		std::vector<char> varstype;
		std::vector<int>  varsnum;
		VarType      varsMap;
		int    nint;
		int    ndeg;
		int    nmul;
		bool   stab;
		int    nmat;
		int    nint1;
		int    nint2;
		int    ndeg1;
		int    ndeg2;
		int    steps;
		double cpMin;
		double cpMax;
		double ds;
		double dsMin;
		double dsMax;
		double dsStart;
		double epsC;
		double epsR;
		double epsS;
		int    nitC;
		int    nitR;
		int    nitS;
		int    nsym;
		std::vector<int> symRe;
		std::vector<int> symIm;
		// from sysname
		int    npar;
		int    ndim;
};

#endif
