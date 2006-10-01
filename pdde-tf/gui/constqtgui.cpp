// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "constqtgui.h"
#include <fstream>
#include <QtXml/QtXml>

void NConstantsQtGui::saveXmlFile(const std::string &fileName)
{
	QDomDocument doc("cfile");
	QDomElement root = doc.createElement("pdde");
	doc.appendChild(root);

	QDomElement xml_input_tag = doc.createElement("input");
	root.appendChild( xml_input_tag );
	QDomText xml_input_var = doc.createTextNode( getInputFile().c_str() );
	xml_input_tag.appendChild( xml_input_var );
	
	QDomElement xml_output_tag = doc.createElement("output");
	root.appendChild( xml_output_tag );
	QDomText xml_output_var = doc.createTextNode( getOutputFile().c_str() );
	xml_output_tag.appendChild( xml_output_var );
	
	QDomElement xml_sysname_tag = doc.createElement("sysname");
	root.appendChild( xml_sysname_tag );
	QDomText xml_sysname_var = doc.createTextNode( getSysName().c_str() );
	xml_sysname_tag.appendChild( xml_sysname_var );
	
	QDomElement xml_ptlabel_tag = doc.createElement("label");
	root.appendChild( xml_ptlabel_tag );
	QDomText xml_ptlabel_var = doc.createTextNode( QString::number(getLabel()) );
	xml_ptlabel_tag.appendChild( xml_ptlabel_var );
	
	QDomElement xml_pointtype_tag = doc.createElement("pointtype");
	root.appendChild( xml_pointtype_tag );
	QDomText xml_pointtype_var = doc.createTextNode( QString::number(getPointType()) );
	xml_pointtype_tag.appendChild( xml_pointtype_var );
	
	QDomElement xml_cptype_tag = doc.createElement("cptype");
	root.appendChild( xml_cptype_tag );
	QDomText xml_cptype_var = doc.createTextNode( QChar(getCpType()) );
	xml_cptype_tag.appendChild( xml_cptype_var );
	
	QDomElement xml_cpnum_tag = doc.createElement("cpnum");
	root.appendChild( xml_cpnum_tag );
	QDomText xml_cpnum_var = doc.createTextNode( QString::number(getCpNum()) );
	xml_cpnum_tag.appendChild( xml_cpnum_var );

	QDomElement xml_brsw_tag = doc.createElement("switch");
	root.appendChild( xml_brsw_tag );
	QDomText xml_brsw_var = doc.createTextNode( QString::number(getBranchSW()) );
	xml_brsw_tag.appendChild( xml_brsw_var );
	
	if( getPointType() != SolUser )
	{
		QDomElement xml_nparx_tag = doc.createElement("nparx");
		root.appendChild( xml_nparx_tag );
		QDomText xml_nparx_var = doc.createTextNode( QString::number(getNEqns()) );
		xml_nparx_tag.appendChild( xml_nparx_var );
		
		if( getNEqns() != 0 )
		{
			QDomElement xml_parx_tag = doc.createElement("parx");
			root.appendChild( xml_parx_tag );
			for( int i = 0; i < getNEqns(); ++i )
			{
				QDomElement xml_parxelem_tag = doc.createElement( "par" );
				xml_parx_tag.appendChild( xml_parxelem_tag );
				QDomElement xml_parxelemtype_tag = doc.createElement("type");
				xml_parxelem_tag.appendChild( xml_parxelemtype_tag );
				QDomText xml_parxelemtype_var = doc.createTextNode( QChar(getParXType(i)) );
				xml_parxelemtype_tag.appendChild( xml_parxelemtype_var );
				QDomElement xml_parxelemnum_tag = doc.createElement("num");
				xml_parxelem_tag.appendChild( xml_parxelemnum_tag );
				QDomText xml_parxelemnum_var = doc.createTextNode( QString::number(getParXNum(i)) );
				xml_parxelemnum_tag.appendChild( xml_parxelemnum_var );
			}
		}
	}else
	{
		QDomElement xml_neqns_tag = doc.createElement("neqns");
		root.appendChild( xml_neqns_tag );
		QDomText xml_neqns_var = doc.createTextNode( QString::number(getNEqns()) );
		xml_neqns_tag.appendChild( xml_neqns_var );
		
		if( getNEqns() != 0 )
		{
			QDomElement xml_eqns_tag = doc.createElement("eqns");
			root.appendChild( xml_eqns_tag );
	
			QDomElement xml_vars_tag = doc.createElement("vars");
			root.appendChild( xml_vars_tag );
			
			for( int i = 0; i < getNEqns(); ++i )
			{
				QDomElement xml_eqnelem_tag = doc.createElement( "eqn" );
				xml_eqns_tag.appendChild( xml_eqnelem_tag );
				QDomElement xml_eqnelemtype_tag = doc.createElement("type");
				xml_eqnelem_tag.appendChild( xml_eqnelemtype_tag );
				QDomText xml_eqnelemtype_var = doc.createTextNode( QChar(getEqnsType(i)) );
				xml_eqnelemtype_tag.appendChild( xml_eqnelemtype_var );
				QDomElement xml_eqnelemnum_tag = doc.createElement("num");
				xml_eqnelem_tag.appendChild( xml_eqnelemnum_tag );
				QDomText xml_eqnelemnum_var = doc.createTextNode( QString::number(getEqnsNum(i)) );
				xml_eqnelemnum_tag.appendChild( xml_eqnelemnum_var );
	
				QDomElement xml_varelem_tag = doc.createElement( "var" );
				xml_vars_tag.appendChild( xml_varelem_tag );
				QDomElement xml_varelemtype_tag = doc.createElement("type");
				xml_varelem_tag.appendChild( xml_varelemtype_tag );
				QDomText xml_varelemtype_var = doc.createTextNode( QChar(getVarsType(i)) );
				xml_varelemtype_tag.appendChild( xml_varelemtype_var );
				QDomElement xml_varelemnum_tag = doc.createElement("num");
				xml_varelem_tag.appendChild( xml_varelemnum_tag );
				QDomText xml_varelemnum_var = doc.createTextNode( QString::number(getVarsNum(i)) );
				xml_varelemnum_tag.appendChild( xml_varelemnum_var );
			}
		}
	}
	
	QDomElement xml_nint_tag = doc.createElement("nint");
	root.appendChild( xml_nint_tag );
	QDomText xml_nint_var = doc.createTextNode( QString::number(getNInt()) );
	xml_nint_tag.appendChild( xml_nint_var );

	QDomElement xml_ndeg_tag = doc.createElement("ndeg");
	root.appendChild( xml_ndeg_tag );
	QDomText xml_ndeg_var = doc.createTextNode( QString::number(getNDeg()) );
	xml_ndeg_tag.appendChild( xml_ndeg_var );

	QDomElement xml_nmul_tag = doc.createElement("nmul");
	root.appendChild( xml_nmul_tag );
	QDomText xml_nmul_var = doc.createTextNode( QString::number(getNMul()) );
	xml_nmul_tag.appendChild( xml_nmul_var );

	QDomElement xml_stab_tag = doc.createElement("stab");
	root.appendChild( xml_stab_tag );
	QDomText xml_stab_var = doc.createTextNode( QString::number(getStab()) );
	xml_stab_tag.appendChild( xml_stab_var );
	
	QDomElement xml_nmat_tag = doc.createElement("nmat");
	root.appendChild( xml_nmat_tag );
	QDomText xml_nmat_var = doc.createTextNode( QString::number(getNMat()) );
	xml_nmat_tag.appendChild( xml_nmat_var );
	
	QDomElement __nint1_tag__ = doc.createElement("nint1");
	root.appendChild( __nint1_tag__ );
	QDomText __nint1__ = doc.createTextNode( QString::number(getNInt1()) );
	__nint1_tag__.appendChild( __nint1__ );
	
	QDomElement __nint2_tag__ = doc.createElement("nint2");
	root.appendChild( __nint2_tag__ );
	QDomText __nint2__ = doc.createTextNode( QString::number(getNInt2()) );
	__nint2_tag__.appendChild( __nint2__ );
	
	QDomElement __ndeg1_tag__ = doc.createElement("ndeg1");
	root.appendChild( __ndeg1_tag__ );
	QDomText __ndeg1__ = doc.createTextNode( QString::number(getNDeg1()) );
	__ndeg1_tag__.appendChild( __ndeg1__ );
	
	QDomElement __ndeg2_tag__ = doc.createElement("ndeg2");
	root.appendChild( __ndeg2_tag__ );
	QDomText __ndeg2__ = doc.createTextNode( QString::number(getNDeg2()) );
	__ndeg2_tag__.appendChild( __ndeg2__ );
	
	QDomElement xml_steps_tag = doc.createElement("steps");
	root.appendChild( xml_steps_tag );
	QDomText xml_steps_var = doc.createTextNode( QString::number(getSteps()) );
	xml_steps_tag.appendChild( xml_steps_var );
	
	QDomElement xml_cpmin_tag = doc.createElement("cpmin");
	root.appendChild( xml_cpmin_tag );
	QDomText xml_cpmin_var = doc.createTextNode( QString::number(getCpMin()) );
	xml_cpmin_tag.appendChild( xml_cpmin_var );

	QDomElement xml_cpmax_tag = doc.createElement("cpmax");
	root.appendChild( xml_cpmax_tag );
	QDomText xml_cpmax_var = doc.createTextNode( QString::number(getCpMax()) );
	xml_cpmax_tag.appendChild( xml_cpmax_var );

	QDomElement xml_ds_tag = doc.createElement("ds");
	root.appendChild( xml_ds_tag );
	QDomText xml_ds_var = doc.createTextNode( QString::number(getDs()) );
	xml_ds_tag.appendChild( xml_ds_var );

	QDomElement xml_dsmin_tag = doc.createElement("dsmin");
	root.appendChild( xml_dsmin_tag );
	QDomText xml_dsmin_var = doc.createTextNode( QString::number(getDsMin()) );
	xml_dsmin_tag.appendChild( xml_dsmin_var );

	QDomElement xml_dsmax_tag = doc.createElement("dsmax");
	root.appendChild( xml_dsmax_tag );
	QDomText xml_dsmax_var = doc.createTextNode( QString::number(getDsMax()) );
	xml_dsmax_tag.appendChild( xml_dsmax_var );

	QDomElement xml_dsstart_tag = doc.createElement("dsstart");
	root.appendChild( xml_dsstart_tag );
	QDomText xml_dsstart_var = doc.createTextNode( QString::number(getDsStart()) );
	xml_dsstart_tag.appendChild( xml_dsstart_var );
	
	QDomElement xml_epsc_tag = doc.createElement("epsc");
	root.appendChild( xml_epsc_tag );
	QDomText xml_epsc_var = doc.createTextNode( QString::number(getEpsC()) );
	xml_epsc_tag.appendChild( xml_epsc_var );
	
	QDomElement xml_epsr_tag = doc.createElement("epsr");
	root.appendChild( xml_epsr_tag );
	QDomText xml_epsr_var = doc.createTextNode( QString::number(getEpsR()) );
	xml_epsr_tag.appendChild( xml_epsr_var );
	
	QDomElement xml_epss_tag = doc.createElement("epss");
	root.appendChild( xml_epss_tag );
	QDomText xml_epss_var = doc.createTextNode( QString::number(getEpsS()) );
	xml_epss_tag.appendChild( xml_epss_var );
	
	QDomElement xml_nitc_tag = doc.createElement("nitc");
	root.appendChild( xml_nitc_tag );
	QDomText xml_nitc_var = doc.createTextNode( QString::number(getNItC()) );
	xml_nitc_tag.appendChild( xml_nitc_var );

	QDomElement xml_nitr_tag = doc.createElement("nitr");
	root.appendChild( xml_nitr_tag );
	QDomText xml_nitr_var = doc.createTextNode( QString::number(getNItR()) );
	xml_nitr_tag.appendChild( xml_nitr_var );
	
	QDomElement xml_nits_tag = doc.createElement("nits");
	root.appendChild( xml_nits_tag );
	QDomText xml_nits_var = doc.createTextNode( QString::number(getNItS()) );
	xml_nits_tag.appendChild( xml_nits_var );

	QDomElement xml_nsym_tag = doc.createElement("nsym");
	root.appendChild( xml_nsym_tag );
	QDomText xml_nsym_var = doc.createTextNode( QString::number(getNSym()) );
	xml_nsym_tag.appendChild( xml_nsym_var );
	
	QDomElement xml_sym_tag = doc.createElement("sym");
	root.appendChild( xml_sym_tag );
	for( int i = 0; i < getNSym(); ++i )
	{
		QDomElement xml_symelem_tag = doc.createElement( "dim" );
		xml_sym_tag.appendChild( xml_symelem_tag );
		QDomElement xml_symelemreal_tag = doc.createElement( "real" );
		xml_symelem_tag.appendChild( xml_symelemreal_tag );
		QDomText xml_symelemreal_var = doc.createTextNode( QString::number(getSymRe(i)) );
		xml_symelemreal_tag.appendChild( xml_symelemreal_var );
		QDomElement xml_symelemimag_tag = doc.createElement( "imag" );
		xml_symelem_tag.appendChild( xml_symelemimag_tag );
		QDomText xml_symelemimag_var = doc.createTextNode( QString::number(getSymIm(i)) );
		xml_symelemimag_tag.appendChild( xml_symelemimag_var );

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
			QDomElement xml_par_var = root.firstChildElement("parx").firstChildElement("par");
			int it = 0;
			while( !xml_par_var.isNull() )
			{
				setParX( it, xml_par_var.firstChildElement("type").firstChild().toText().data()[0].toAscii(),
										xml_par_var.firstChildElement("num").firstChild().toText().data().toInt() );
				xml_par_var = xml_par_var.nextSiblingElement();
				++it;
			}
		}
	}else
	{
		setNEqns( root.firstChildElement("neqns").firstChild().toText().data().toInt() );
		if( getNEqns() != 0 )
		{
			QDomElement xml_eqn_var = root.firstChildElement("eqns").firstChildElement("eqn");
			int it = 0;
			while( !xml_eqn_var.isNull() )
			{
// 				std::cout<<"setting up EQN\n";
				if( xml_eqn_var.firstChildElement("type").firstChild().toText().data()[0].toAscii() != 'E' ) return;
				setEqns( it, 'E', xml_eqn_var.firstChildElement("num").firstChild().toText().data().toInt() );
				xml_eqn_var = xml_eqn_var.nextSiblingElement();
				++it;
			}
			QDomElement xml_var_var = root.firstChildElement("vars").firstChildElement("var");
			it = 0;
			while( !xml_var_var.isNull() )
			{
// 				std::cout<<"setting up VAR\n";
				setVars( it, xml_var_var.firstChildElement("type").firstChild().toText().data()[0].toAscii(),
										xml_var_var.firstChildElement("num").firstChild().toText().data().toInt() );
				xml_var_var = xml_var_var.nextSiblingElement();
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
	
	QDomElement xml_sym_var = root.firstChildElement("sym").firstChildElement("dim");
	int it = 0;
	while( !xml_sym_var.isNull() && it < getNSym() )
	{
		setSymRe( it, xml_sym_var.firstChildElement("real").firstChild().toText().data().toInt() );
		setSymIm( it, xml_sym_var.firstChildElement("imag").firstChild().toText().data().toInt() );
		xml_sym_var = xml_sym_var.nextSiblingElement();
		++it;
	}
	
	file.close();
}
