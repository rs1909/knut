#include <fstream>
#include "constants.h"

#ifdef PDDE_GUI
  #include <QtXml>
#endif


class toNumber
{
 public:
	toNumber( int i ) { int n=snprintf(num,31,"%d", i); }
	toNumber( double d ) { int n=snprintf(num,31,"%lf", d); }
	const char* c_str() const { return num; }
 private:
	char num[32];
};

#ifdef PDDE_GUI

void NConstants::saveXmlFile(const std::string &fileName)
{
	QDomDocument doc("cfile");
	QDomElement root = doc.createElement("pdde");
	doc.appendChild(root);

	QDomElement __input_tag__ = doc.createElement("input");
	root.appendChild( __input_tag__ );
	QDomText __input__ = doc.createTextNode( getInputFile().c_str() );
	__input_tag__.appendChild( __input__ );
	
	QDomElement __output_tag__ = doc.createElement("output");
	root.appendChild( __output_tag__ );
	QDomText __output__ = doc.createTextNode( getOutputFile().c_str() );
	__output_tag__.appendChild( __output__ );
	
	QDomElement __sysname_tag__ = doc.createElement("sysname");
	root.appendChild( __sysname_tag__ );
	QDomText __sysname__ = doc.createTextNode( getSysName().c_str() );
	__sysname_tag__.appendChild( __sysname__ );
	
	QDomElement __ptlabel_tag__ = doc.createElement("label");
	root.appendChild( __ptlabel_tag__ );
	QDomText __ptlabel__ = doc.createTextNode( toNumber(getLabel()).c_str() );
	__ptlabel_tag__.appendChild( __ptlabel__ );
	
	QDomElement __pointtype_tag__ = doc.createElement("pointtype");
	root.appendChild( __pointtype_tag__ );
	QDomText __pointtype__ = doc.createTextNode( toNumber(getPointType()).c_str() );
	__pointtype_tag__.appendChild( __pointtype__ );
	
	QDomElement __cptype_tag__ = doc.createElement("cptype");
	root.appendChild( __cptype_tag__ );
	QDomText __cptype__ = doc.createTextNode( QChar(getCpType()) );
	__cptype_tag__.appendChild( __cptype__ );
	
	QDomElement __cpnum_tag__ = doc.createElement("cpnum");
	root.appendChild( __cpnum_tag__ );
	QDomText __cpnum__ = doc.createTextNode( toNumber(getCpNum()).c_str() );
	__cpnum_tag__.appendChild( __cpnum__ );

	QDomElement __brsw_tag__ = doc.createElement("switch");
	root.appendChild( __brsw_tag__ );
	QDomText __brsw__ = doc.createTextNode( toNumber(getBranchSW()).c_str() );
	__brsw_tag__.appendChild( __brsw__ );
	
	if( getPointType() != SolUser )
	{
		QDomElement __nparx_tag__ = doc.createElement("nparx");
		root.appendChild( __nparx_tag__ );
		QDomText __nparx__ = doc.createTextNode( toNumber(getNEqns()).c_str() );
		__nparx_tag__.appendChild( __nparx__ );
		
		if( getNEqns() != 0 )
		{
			QDomElement __parx_tag__ = doc.createElement("parx");
			root.appendChild( __parx_tag__ );
			for( int i = 0; i < getNEqns(); ++i )
			{
				QDomElement __parxelem_tag__ = doc.createElement( "par" );
				__parx_tag__.appendChild( __parxelem_tag__ );
				QDomElement __parxelemtype_tag__ = doc.createElement("type");
				__parxelem_tag__.appendChild( __parxelemtype_tag__ );
				QDomText __parxelemtype__ = doc.createTextNode( QChar(getParXType(i)) );
				__parxelemtype_tag__.appendChild( __parxelemtype__ );
				QDomElement __parxelemnum_tag__ = doc.createElement("num");
				__parxelem_tag__.appendChild( __parxelemnum_tag__ );
				QDomText __parxelemnum__ = doc.createTextNode( toNumber(getParXNum(i)).c_str() );
				__parxelemnum_tag__.appendChild( __parxelemnum__ );
			}
		}
	}else
	{
		QDomElement __neqns_tag__ = doc.createElement("neqns");
		root.appendChild( __neqns_tag__ );
		QDomText __neqns__ = doc.createTextNode( toNumber(getNEqns()).c_str() );
		__neqns_tag__.appendChild( __neqns__ );
		
		if( getNEqns() != 0 )
		{
			QDomElement __eqns_tag__ = doc.createElement("eqns");
			root.appendChild( __eqns_tag__ );
	
			QDomElement __vars_tag__ = doc.createElement("vars");
			root.appendChild( __vars_tag__ );
			
			for( int i = 0; i < getNEqns(); ++i )
			{
				QDomElement __eqnelem_tag__ = doc.createElement( "eqn" );
				__eqns_tag__.appendChild( __eqnelem_tag__ );
				QDomElement __eqnelemtype_tag__ = doc.createElement("type");
				__eqnelem_tag__.appendChild( __eqnelemtype_tag__ );
				QDomText __eqnelemtype__ = doc.createTextNode( QChar(getEqnsType(i)) );
				__eqnelemtype_tag__.appendChild( __eqnelemtype__ );
				QDomElement __eqnelemnum_tag__ = doc.createElement("num");
				__eqnelem_tag__.appendChild( __eqnelemnum_tag__ );
				QDomText __eqnelemnum__ = doc.createTextNode( toNumber(getEqnsNum(i)).c_str() );
				__eqnelemnum_tag__.appendChild( __eqnelemnum__ );
	
				QDomElement __varelem_tag__ = doc.createElement( "var" );
				__vars_tag__.appendChild( __varelem_tag__ );
				QDomElement __varelemtype_tag__ = doc.createElement("type");
				__varelem_tag__.appendChild( __varelemtype_tag__ );
				QDomText __varelemtype__ = doc.createTextNode( QChar(getVarsType(i)) );
				__varelemtype_tag__.appendChild( __varelemtype__ );
				QDomElement __varelemnum_tag__ = doc.createElement("num");
				__varelem_tag__.appendChild( __varelemnum_tag__ );
				QDomText __varelemnum__ = doc.createTextNode( toNumber(getVarsNum(i)).c_str() );
				__varelemnum_tag__.appendChild( __varelemnum__ );
			}
		}
	}
	
	QDomElement __nint_tag__ = doc.createElement("nint");
	root.appendChild( __nint_tag__ );
	QDomText __nint__ = doc.createTextNode( toNumber(getNInt()).c_str() );
	__nint_tag__.appendChild( __nint__ );

	QDomElement __ndeg_tag__ = doc.createElement("ndeg");
	root.appendChild( __ndeg_tag__ );
	QDomText __ndeg__ = doc.createTextNode( toNumber(getNDeg()).c_str() );
	__ndeg_tag__.appendChild( __ndeg__ );

	QDomElement __nmul_tag__ = doc.createElement("nmul");
	root.appendChild( __nmul_tag__ );
	QDomText __nmul__ = doc.createTextNode( toNumber(getNMul()).c_str() );
	__nmul_tag__.appendChild( __nmul__ );

	QDomElement __stab_tag__ = doc.createElement("stab");
	root.appendChild( __stab_tag__ );
	QDomText __stab__ = doc.createTextNode( toNumber(getStab()).c_str() );
	__stab_tag__.appendChild( __stab__ );
	
	QDomElement __nmat_tag__ = doc.createElement("nmat");
	root.appendChild( __nmat_tag__ );
	QDomText __nmat__ = doc.createTextNode( toNumber(getNMat()).c_str() );
	__nmat_tag__.appendChild( __nmat__ );
	
	QDomElement __nint1_tag__ = doc.createElement("nint1");
	root.appendChild( __nint1_tag__ );
	QDomText __nint1__ = doc.createTextNode( toNumber(getNInt1()).c_str() );
	__nint1_tag__.appendChild( __nint1__ );
	
	QDomElement __nint2_tag__ = doc.createElement("nint2");
	root.appendChild( __nint2_tag__ );
	QDomText __nint2__ = doc.createTextNode( toNumber(getNInt2()).c_str() );
	__nint2_tag__.appendChild( __nint2__ );
	
	QDomElement __ndeg1_tag__ = doc.createElement("ndeg1");
	root.appendChild( __ndeg1_tag__ );
	QDomText __ndeg1__ = doc.createTextNode( toNumber(getNDeg1()).c_str() );
	__ndeg1_tag__.appendChild( __ndeg1__ );
	
	QDomElement __ndeg2_tag__ = doc.createElement("ndeg2");
	root.appendChild( __ndeg2_tag__ );
	QDomText __ndeg2__ = doc.createTextNode( toNumber(getNDeg2()).c_str() );
	__ndeg2_tag__.appendChild( __ndeg2__ );
	
	QDomElement __steps_tag__ = doc.createElement("steps");
	root.appendChild( __steps_tag__ );
	QDomText __steps__ = doc.createTextNode( toNumber(getSteps()).c_str() );
	__steps_tag__.appendChild( __steps__ );
	
	QDomElement __cpmin_tag__ = doc.createElement("cpmin");
	root.appendChild( __cpmin_tag__ );
	QDomText __cpmin__ = doc.createTextNode( toNumber(getCpMin()).c_str() );
	__cpmin_tag__.appendChild( __cpmin__ );

	QDomElement __cpmax_tag__ = doc.createElement("cpmax");
	root.appendChild( __cpmax_tag__ );
	QDomText __cpmax__ = doc.createTextNode( toNumber(getCpMax()).c_str() );
	__cpmax_tag__.appendChild( __cpmax__ );

	QDomElement __ds_tag__ = doc.createElement("ds");
	root.appendChild( __ds_tag__ );
	QDomText __ds__ = doc.createTextNode( toNumber(getDs()).c_str() );
	__ds_tag__.appendChild( __ds__ );

	QDomElement __dsmin_tag__ = doc.createElement("dsmin");
	root.appendChild( __dsmin_tag__ );
	QDomText __dsmin__ = doc.createTextNode( toNumber(getDsMin()).c_str() );
	__dsmin_tag__.appendChild( __dsmin__ );

	QDomElement __dsmax_tag__ = doc.createElement("dsmax");
	root.appendChild( __dsmin_tag__ );
	QDomText __dsmax__ = doc.createTextNode( toNumber(getDsMax()).c_str() );
	__dsmax_tag__.appendChild( __dsmax__ );

	QDomElement __dsstart_tag__ = doc.createElement("dsstart");
	root.appendChild( __dsstart_tag__ );
	QDomText __dsstart__ = doc.createTextNode( toNumber(getDsStart()).c_str() );
	__dsstart_tag__.appendChild( __dsstart__ );
	
	QDomElement __epsc_tag__ = doc.createElement("epsc");
	root.appendChild( __epsc_tag__ );
	QDomText __epsc__ = doc.createTextNode( toNumber(getEpsC()).c_str() );
	__epsc_tag__.appendChild( __epsc__ );
	
	QDomElement __epsr_tag__ = doc.createElement("epsr");
	root.appendChild( __epsr_tag__ );
	QDomText __epsr__ = doc.createTextNode( toNumber(getEpsR()).c_str() );
	__epsr_tag__.appendChild( __epsr__ );
	
	QDomElement __epss_tag__ = doc.createElement("epss");
	root.appendChild( __epss_tag__ );
	QDomText __epss__ = doc.createTextNode( toNumber(getEpsS()).c_str() );
	__epss_tag__.appendChild( __epss__ );
	
	QDomElement __nitc_tag__ = doc.createElement("nitc");
	root.appendChild( __nitc_tag__ );
	QDomText __nitc__ = doc.createTextNode( toNumber(getNItC()).c_str() );
	__nitc_tag__.appendChild( __nitc__ );

	QDomElement __nitr_tag__ = doc.createElement("nitr");
	root.appendChild( __nitr_tag__ );
	QDomText __nitr__ = doc.createTextNode( toNumber(getNItR()).c_str() );
	__nitr_tag__.appendChild( __nitr__ );
	
	QDomElement __nits_tag__ = doc.createElement("nits");
	root.appendChild( __nits_tag__ );
	QDomText __nits__ = doc.createTextNode( toNumber(getNItS()).c_str() );
	__nits_tag__.appendChild( __nits__ );

	QDomElement __nsym_tag__ = doc.createElement("nsym");
	root.appendChild( __nsym_tag__ );
	QDomText __nsym__ = doc.createTextNode( toNumber(getNSym()).c_str() );
	__nsym_tag__.appendChild( __nsym__ );
	
	QDomElement __sym_tag__ = doc.createElement("sym");
	root.appendChild( __sym_tag__ );
	for( int i = 0; i < getNSym(); ++i )
	{
		QDomElement __symelem_tag__ = doc.createElement( "dim" );
		__sym_tag__.appendChild( __symelem_tag__ );
		QDomElement __symelemreal_tag__ = doc.createElement( "real" );
		__symelem_tag__.appendChild( __symelemreal_tag__ );
		QDomText __symelemreal__ = doc.createTextNode( toNumber(getSymRe(i)).c_str() );
		__symelemreal_tag__.appendChild( __symelemreal__ );
		QDomElement __symelemimag_tag__ = doc.createElement( "imag" );
		__symelem_tag__.appendChild( __symelemimag_tag__ );
		QDomText __symelemimag__ = doc.createTextNode( toNumber(getSymIm(i)).c_str() );
		__symelemimag_tag__.appendChild( __symelemimag__ );

	}
	std::ofstream file( fileName.c_str() );
	file<<doc.toString().toStdString().c_str();
}

void NConstants::loadXmlFile(const std::string &fileName)
{
	QDomDocument doc;
	
	QFile file( fileName.c_str() );
	
	if (!file.open(QIODevice::ReadOnly))
	{
		std::cout<<"Readonly\n"; std::cout.flush();
		return;
	}
	QString errorMsg; int errorLine; int errorColumn;
	if (!doc.setContent(&file, false, &errorMsg, &errorLine, &errorColumn ))
	{
		std::cout<<"Some other problem"<<errorMsg.toStdString()<<" "<<errorLine<<" "<<errorColumn<<"\n"; std::cout.flush();
		file.close();
		return;
	}
	
	QDomElement root = doc.documentElement();
	setInputFileText( root.firstChildElement("input").firstChild().toText().data().toStdString() );
	setOutputFileText( root.firstChildElement("output").firstChild().toText().data().toStdString() );
	setSysNameText( root.firstChildElement("sysname").firstChild().toText().data().toStdString() );
	setLabel( root.firstChildElement("label").firstChild().toText().data().toInt() );
	setPointType( (PtType)(root.firstChildElement("pointtype").firstChild().toText().data().toInt()) );
	setCp( root.firstChildElement("cptype").firstChild().toText().data()[0].toAscii(),
					root.firstChildElement("cpnum").firstChild().toText().data().toInt() );
	setBranchSW( (BranchSW)(root.firstChildElement("switch").firstChild().toText().data().toInt()) );

	if( getPointType() != SolUser )
	{
		setNEqns( root.firstChildElement("nparx").firstChild().toText().data().toInt() );
		if( getNEqns() != 0 )
		{
			QDomElement __par__ = root.firstChildElement("parx").firstChildElement("par");
			int it = 0;
			while( !__par__.isNull() )
			{
				setParX( it, __par__.firstChildElement("type").toText().data()[0].toAscii(),
										__par__.firstChildElement("num").toText().data().toInt() );
				__par__ = __par__.nextSiblingElement();
				++it;
			}
		}
	}else
	{
		setNEqns( root.firstChildElement("neqns").firstChild().toText().data().toInt() );
		if( getNEqns() != 0 )
		{
			QDomElement __eqn__ = root.firstChildElement("eqns").firstChildElement("eqn");
			int it = 0;
			while( !__eqn__.isNull() )
			{
// 				std::cout<<"setting up EQN\n";
				if( __eqn__.firstChildElement("type").firstChild().toText().data()[0].toAscii() != 'E' ) return;
				setEqns( it, 'E', __eqn__.firstChildElement("num").firstChild().toText().data().toInt() );
				__eqn__ = __eqn__.nextSiblingElement();
				++it;
			}
			QDomElement __var__ = root.firstChildElement("vars").firstChildElement("var");
			it = 0;
			while( !__var__.isNull() )
			{
// 				std::cout<<"setting up VAR\n";
				setVars( it, __var__.firstChildElement("type").firstChild().toText().data()[0].toAscii(),
										__var__.firstChildElement("num").firstChild().toText().data().toInt() );
				__var__ = __var__.nextSiblingElement();
				++it;
			}
		}
	}
	setNInt( root.firstChildElement("nint").firstChild().toText().data().toInt() );
	setNDeg( root.firstChildElement("ndeg").firstChild().toText().data().toInt() );
	setNMul( root.firstChildElement("nmul").firstChild().toText().data().toInt() );
	setStab( root.firstChildElement("stab").firstChild().toText().data().toInt() != 0 );
	setNMat( root.firstChildElement("nmat").firstChild().toText().data().toInt() );

	setNInt1( root.firstChildElement("nint1").firstChild().toText().data().toInt() );
	setNInt2( root.firstChildElement("nint2").firstChild().toText().data().toInt() );
	setNDeg1( root.firstChildElement("ndeg1").firstChild().toText().data().toInt() );
	setNDeg2( root.firstChildElement("ndeg2").firstChild().toText().data().toInt() );

	setSteps( root.firstChildElement("steps").firstChild().toText().data().toInt() );
	setCpMin( root.firstChildElement("cpmin").firstChild().toText().data().toDouble() );
	setCpMax( root.firstChildElement("cpmax").firstChild().toText().data().toDouble() );
	setDs( root.firstChildElement("ds").firstChild().toText().data().toDouble() );
	setDsMin( root.firstChildElement("dsmin").firstChild().toText().data().toDouble() );
	setDsMax( root.firstChildElement("dsmax").firstChild().toText().data().toDouble() );
	setDsStart( root.firstChildElement("dsstart").firstChild().toText().data().toDouble() );
	setEpsC( root.firstChildElement("epsc").firstChild().toText().data().toDouble() );
	setEpsR( root.firstChildElement("epsr").firstChild().toText().data().toDouble() );
	setEpsS( root.firstChildElement("epss").firstChild().toText().data().toDouble() );
	setNItC( root.firstChildElement("nitc").firstChild().toText().data().toInt() );
	setNItR( root.firstChildElement("nitr").firstChild().toText().data().toInt() );
	setNItS( root.firstChildElement("nits").firstChild().toText().data().toInt() );
	
	setNSym( root.firstChildElement("nsym").firstChild().toText().data().toInt() );
	
	QDomElement __sym__ = root.firstChildElement("sym").firstChildElement("dim");
	int it = 0;
	while( !__sym__.isNull() && it < getNSym() )
	{
		setSymRe( it, __sym__.firstChildElement("real").firstChild().toText().data().toInt() );
		setSymIm( it, __sym__.firstChildElement("imag").firstChild().toText().data().toInt() );
		__sym__ = __sym__.nextSiblingElement();
		++it;
	}
	
	file.close();
}

#endif // PDDE_GUI

inline bool NConstants::inputAssert( std::istream& is )
{
	if( is.rdstate() & (std::istream::eofbit | std::istream::failbit | std::istream::badbit) )
	{
		switch( is.rdstate() )
		{
			case std::istream::eofbit:
				std::cout<<"Unexpected end of file\n";
// 				QMessageBox::critical( this, std::string("IO Error"), std::string("Unexpected end of file"), QMessageBox::Ok, QMessageBox::NoButton );
				return true;
				break;
			case std::istream::failbit:
				std::cout<<"Input failed\n";
// 				QMessageBox::critical( this, std::string("IO Error"), std::string("Input failed"), QMessageBox::Ok, QMessageBox::NoButton );
				return true;
				break;
			case std::istream::badbit:
				std::cout<<"Bad input\n";
// 				QMessageBox::critical( this, std::string("IO Error"), std::string("Bad input"), QMessageBox::Ok, QMessageBox::NoButton );
				return true;
				break;
			default:
// 				QMessageBox::critical( this, std::string("IO Error"), std::string("Unexpected error"), QMessageBox::Ok, QMessageBox::NoButton );
				return true;
				break;
		}
	}
	return false;
}

void NConstants::loadFile(const std::string &fileName)
{
	std::ifstream file;
	
	file.open( fileName.c_str() );
	if(!file){ std::cout<<"Cannot open "<<fileName<<"\n"; return; }
	
	std::string __sysname__;
	file >> __sysname__; if( inputAssert( file ) ) return;
	while( file.get() != '\n' ); if( inputAssert( file ) ) return;

	if( __sysname__.find('/') == std::string::npos )
	{
		__sysname__ = QFileInfo( fileName.c_str() ).absolutePath( ).toStdString() + "/" + __sysname__;
	}
	setSysNameText( std::string( __sysname__.c_str() ) );

	int __ptlabel__;
	file >> __ptlabel__; if( inputAssert( file ) ) return;
	setLabel( __ptlabel__ );
	while( file.get() != '\n' ); if( inputAssert( file ) ) return;
	
	int __pttype__;
	char __cptype__;
	int __cp__;
	file >> __pttype__ >> __cptype__ >> __cp__; if( inputAssert( file ) ) return;
	if( ( __cptype__ != 'P' ) && ( __cptype__ != 'I' ) ) { std::cout<<"Error: CP: Bad parameter type."; }
	setPointType( (PtType)__pttype__ );
	setCp( __cptype__, __cp__ );
	
	bool loadsym = false;
	if( __pttype__ == SolUser )
	{
		int __brswitch__;
		int __neqn__;
		file >> __brswitch__; if( inputAssert( file ) ) return;
		file >> __neqn__; if( inputAssert( file ) ) return;
		setNEqns( __neqn__ );
		for( int i = 0; i < __neqn__; i++ )
		{
			char __eqntype__;
			int  __eqnnum__;
			file >> __eqntype__ >> __eqnnum__; if( inputAssert( file ) ) return;
			if( __eqntype__ != 'E' ) { std::cout<<"Error: EQN: Bad equation type."; }
			setEqns( i, 'E', __eqnnum__ );
			if( __eqnnum__ == EqnPhaseRot ) loadsym = true;
		}
		int __nvar__;
		file >> __nvar__; if( inputAssert( file ) ) return;
		if( __nvar__ != __neqn__ ) { std::cout << "Error: NVAR: Number of equations and variables are not the same."; }
		for( int i = 0; i < __neqn__; i++ )
		{
			char __vartype__;
			int  __varnum__;
			file >> __vartype__ >> __varnum__; if( inputAssert( file ) ) return;
			if( ( __vartype__ != 'S' ) && 
			    ( __vartype__ != 'P' ) && 
			    ( __vartype__ != 'I' ) )
				{ std::cout<<"Error: VAR: Bad variable/parameter type."; return; }
			setVars( i, __vartype__, __varnum__ );
		}
		while( file.get() != '\n' );
	}else
	{
		int __nvarx__;
		file >> __nvarx__; if( inputAssert( file ) ) return;
		setNEqns( __nvarx__ );
		for( int i = 0; i < __nvarx__; i++ )
		{
			char __vartype__;
			int  __varnum__;
			file >> __vartype__ >> __varnum__; if( inputAssert( file ) ) return;
			if( ( __vartype__ != 'P' ) && 
			    ( __vartype__ != 'I' ) )
				{ std::cout<<"Error: PARX: Bad parameter type."; }
			setParX( i, __vartype__, __varnum__ );
		}
		while( file.get() != '\n' );
	}
	
	int __nint__, __ndeg__, __nmul__, __stab__, __nmat__;
	file >> __nint__ >> __ndeg__ >> __nmul__ >> __stab__ >> __nmat__; if( inputAssert( file ) ) return;
	while( file.get() != '\n' );
	setNInt( __nint__ );
	setNDeg( __ndeg__ );
	setNMul( __nmul__ );
	
	setStab( __stab__ != 0 );
	setNMat( __nmat__ );
	
	int __nint1__, __nint2__, __ndeg1__, __ndeg2__; if( inputAssert( file ) ) return;
	file >> __nint1__ >> __nint2__ >> __ndeg1__ >> __ndeg2__;
	while( file.get() != '\n' );
	setNInt1( __nint1__ );
	setNInt2( __nint2__ );
	setNDeg1( __ndeg1__ );
	setNDeg2( __ndeg2__ );
	
	int __steps__;
	double __cpmin__, __cpmax__;
	file >> __steps__ >> __cpmin__ >> __cpmax__; if( inputAssert( file ) ) return;
	while( file.get() != '\n' );
	setSteps( __steps__ );
	setCpMin( __cpmin__ );
	setCpMax( __cpmax__ );
	
	double __ds__, __dsmin__, __dsmax__, __dsstart__;
	file >> __ds__ >> __dsmin__ >> __dsmax__ >> __dsstart__; if( inputAssert( file ) ) return;
	while( file.get() != '\n' );
	setDs( __ds__ );
	setDsMin( __dsmin__ );
	setDsMax( __dsmax__ );
	setDsStart( __dsstart__ );
	
	double __epsc__, __epsr__, __epss__;
	file >> __epsc__ >> __epsr__ >> __epss__; if( inputAssert( file ) ) return;
	while( file.get() != '\n' );
	setEpsC( __epsc__ );
	setEpsR( __epsr__ );
	setEpsS( __epss__ );
	
	int __nitc__, __nitr__, __nits__;
	file >> __nitc__ >> __nitr__ >> __nits__; if( inputAssert( file ) ) return;
	setNItC( __nitc__ );
	setNItR( __nitr__ );
	setNItS( __nits__ );

	if( loadsym )
	{
		while( file.get() != '\n' );
		int __nsym__;
		file >> __nsym__; if( inputAssert( file ) ) return;
		setNSym( __nsym__ );
		int __re__, __im__;
		for( int i = 0; i < __nsym__; ++i )
		{
			file >> __re__ >> __im__; if( inputAssert( file ) ) return;
			setSymRe( i, __re__ );
			setSymIm( i, __im__ );
		}
	}else
	{
		setNSym( 0 );
	}
}

void NConstants::saveFile(const std::string &fileName)
{
	std::ofstream file( fileName.c_str() );
	
	file << getSysName() << " \t\tSYSNAME\n";
	
	file << getLabel() << " \t\t\tLABEL\n";
	
	const PtType __pttype__ = getPointType();
	if( __pttype__ != SolUser )
	{
		file << (int)__pttype__ << " " << getCpType() << getCpNum() <<" "<< getNEqns() <<" ";
		for( int i = 0; i < getNEqns(); ++i ) file << getParXType(i) << getParXNum(i) <<" ";
		file <<"\t\tTYPE, CP, NPARX, PARX[NPARX]\n";
	}else
	{
		file << (int)__pttype__ << " " << getCpType() << getCpNum() <<" "
		     << 0 /*switch */ << " " << getNEqns() << " ";
		for( int i = 0; i < getNEqns(); i++ ) file << getEqnsType(i) << getEqnsNum(i) <<" ";
		file<< getNEqns() <<" ";
		for( int i = 0; i < getNEqns(); i++ ) file << getVarsType(i) << getVarsNum(i) <<" ";
		file<<"\tTYPE, CP, SWITCH, NEQN, EQN[NEQN], NVAR, VAR[NVAR]\n";
	}

	file << getNInt() << " "
	     << getNDeg() << " "
	     << getNMul() << " "
	     << getStab() << " "
	     << getNMat() << " \t\tNINT, NDEG, NMUL, STAB, NMAT\n";

	file << getNInt1() << " "
	     << getNInt2() << " "
	     << getNDeg2() << " "
	     << getNDeg2() << " \t\tNINT1, NINT2, NDEG1, NDEG2\n";
	
	file << getSteps() << " "
	     << getCpMin() << " "
	     << getCpMax() << " \t\tSTEPS, CPMIN, CPMAX\n";
	
	file << getDs() << " "
	     << getDsMin() << " "
	     << getDsMax() << " "
	     << getDsStart() << " \tDS, DSMIN, DSMAX, DSSTART \n";
	
	file << getEpsC() << " "
	     << getEpsR() << " "
	     << getEpsS() << " \tEPSC, EPSR, EPSS\n";
	
	file << getNItC() <<" "<< getNItR() <<" "<< getNItS() <<" \t\tNITC, NITR, NITS\n";

	file << getNSym() << " ";
	for( int i = 0; i < getNSym(); ++i ) file << getSymRe(i) <<" "<< getSymIm(i) <<" ";
	file<<"\t\tNSYM";
	for( int i = 0; i < getNSym(); ++i ) file << ", R" << i << ", I" << i;
	file << "\n";

}

int NConstants::toEqnVar( System& sys,
                Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN )
{
	// initializing the equations and variables
	if( getPointType() == SolUser )
	{
		eqn.Init(getNEqns());
		var.Init(getNEqns());
		for( int i = 0; i < getNEqns(); i++ )
		{
			eqn(i) = getEqns(i);
			var(i) = getVars(i);
		}
	}else
	{
		Array1D<Var> L_PARX(getNEqns());
		for( int i = 0; i < getNEqns(); i++ )
		{
			L_PARX(i) = getParX(i);
		}
		setBranchSW( PtToEqnVar( eqn, var, getPointType(), L_PARX, sys.npar() ) );
	}
	
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
			else P_MESSAGE("too many test functionals");
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
	
	if( getBranchSW() == TFHBSwitch )
	  { Array1D<Var> d(0); PtToEqnVar( eqn_refine, var_refine, SolTF, d, sys.npar() ); }// NPARX == 0
	
	// Here, we set up the branch switching.
	// We suppose that if there is a switch we use one parameter continuation afterwards
	// without using characteristic matrices. This means that we can switch on the characteristic matrix,
	// include the equation for the eigenvector norm before the other equations and
	// add CP to the variables as a normal parameter.
	Eqn eqn_temp;
	switch( getBranchSW() )
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
			var_start( var_refine.Size() ) = getCp();
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
			var_start( var_refine.Size() + 1 ) = getCp();
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
