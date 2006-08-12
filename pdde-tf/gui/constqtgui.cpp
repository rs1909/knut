
#include "constqtgui.h"
#include <QtXml>

void NConstantsQtGui::saveXmlFile(const std::string &fileName)
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
	root.appendChild( __dsmax_tag__ );
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

void NConstantsQtGui::loadXmlFile(const std::string &fileName)
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
				setParX( it, __par__.firstChildElement("type").firstChild().toText().data()[0].toAscii(),
										__par__.firstChildElement("num").firstChild().toText().data().toInt() );
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
