#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "pderror.h"
#include "system.h"
#include "point.h"
#include "torpoint.h"
#include "parameters.h"
#include "compthread.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <QErrorMessage>

inline void parNamePrint( Vector& /*par*/, int npar, Array1D<Var>& var )
{
	for( int j = 1; j < var.Size(); j++ ) std::cout<<"\t"<<parType( npar, var(j) - VarPAR0 )<<parNum( npar, var(j) - VarPAR0 )<<"\t";
}

inline void parValuePrint( Vector& par, int /*npar*/, Array1D<Var>& var )
{
	for( int j = 1; j < var.Size(); j++ ) std::cout<<"\t"<<par( var(j) - VarPAR0 );
}

void MThread::run()
{
	try
	{
		System sys( params->getSysName() );
		if( sys.ndim() == 0 ) P_MESSAGE("zerodimensions");
	
		Vector par( params->getNPar()+ParEnd );
		mat4Data out ( params->getOutputFile(),
						params->getSteps(), sys.ndim(), sys.npar()+ParEnd,
						params->getNInt(), params->getNDeg(), params->getNMul() );
	// 	std::ofstream out( params->getOutputFile().c_str() );
	// 	out<<std::scientific;
	// 	out.precision(12);
		
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
		int trivial = params->toEqnVar( sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, testFN );
		const int npar = params->getNPar();
	
// 	for( int i=0; i<eqn.Size(); i++ ) std::cout<<EqnToStr( eqn(i) )<<", ";
// 	std::cout<<'\n';
// 	for( int i=0; i<var.Size(); i++ ) std::cout<<VarToStr( var(i) )<<", ";
// 	std::cout<<'\n';

		///!!!!!!!!!!!!!!!!!!
		//-----------------------------------------------------------------------------------------------------------
		//
		// END of initialization
		//
		//-----------------------------------------------------------------------------------------------------------
		
		Point* pt_ptr = new Point( sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul(), params->getNMat() );
		Point& pt = *pt_ptr;
		
		pt.setContIter( params->getNItC() );
		pt.setRefIter( params->getNItR() );
		pt.setStartIter( params->getNItS() );
		pt.setRefEps( params->getEpsR() );
		pt.setContEps( params->getEpsC() );
		pt.setStartEps( params->getEpsS() );
		pt.setCont( params->getCp()-VarPAR0 );
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Array1D<int> sre(params->getNSym()), sim(params->getNSym());
		for( int i=0; i<sre.Size(); ++i ) { sre(i) = params->getSymRe(i); sim(i) = params->getSymIm(i); }
		if( params->getNSym() != 0 ) pt.setSym( sre, sim );
		
		std::cout<<std::scientific;
		std::cout.precision(6);
		
		// load the initial guess
		if( params->getLabel() != 0 )
		{
			mat4Data istr( params->getInputFile() );
			pt.BinaryRead( istr, params->getLabel() );
		}
		
		pt.Refine();
		if( testFN != EqnNone )
		{
			std::cout<<"\n--- Finding the bifurcation point (TF) ---\n";
			pt.Reset( eqn_start, var_start );
			pt.setCont( params->getCp()-VarPAR0 );
			pt.StartTF( testFN ); // it only computes the characteristic multiplier refines the solution
		}
		
		// start the continuation!
		if( params->getBranchSW() != TFTRSwitch )
		{
			std::cout<<"\n--- Starting the continuation ---\n";
			
			for( int j = 0; j < par.Size(); j++ ) par(j) = pt.getPar()(j);
			//
			std::cout<<"\nLABEL\t"<<"   NORM\t\t"<<parType( npar, params->getCp()-VarPAR0 )<<parNum( npar, params->getCp()-VarPAR0 )<<"\t";
			parNamePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<parType( npar, (params->PARX)[j] )<<parNum( npar, (params->PARX)[j] )<<"\t";
			std::cout<<"\n";
			//
			std::cout<<"  "<<0<<"\t"<<pt.Norm()<<"\t"<<par(params->getCp()-VarPAR0);
			parValuePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<par((params->PARX)[j]);
			std::cout<<"\n";
			
			// making tangents
			if( params->getBranchSW() == TFPDSwitch )
			{
				std::cout<<"\nSwitching to the period two branch (TF).\n";
				pt.SwitchTFPD( params->getDsStart() );
				pt.setCont( params->getCp()-VarPAR0 );
			}
			else if( params->getBranchSW() == TFBRSwitch )
			{
				std::cout<<"\nSwitching to the other branch (TF).\n";
				pt.SwitchTFLP( params->getDsStart() );
				pt.setCont( params->getCp()-VarPAR0 );
			}
			else if( params->getBranchSW() == TFHBSwitch )
			{
				std::cout<<"\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
				pt.SwitchTFHB( params->getDsStart() );
				pt.setCont( params->getCp()-VarPAR0 );
			}
			else
			{
				std::cout<<"\nFinding the tangent.\n";
				pt.setCont( params->getCp()-VarPAR0 );
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
			double ds = params->getDs();
			for( int i = 0; i < params->getSteps(); i++ ) // 35
			{
				if( i % 24 == 0 )
				{
					std::cout<<"LABEL\t"<<"   NORM\t\t"<<(char)params->getCpType()<<params->getCpNum()<<"\t";
					parNamePrint( par, npar, var ); // for( int j = 0; j < params->NPARX; j++ ) std::cout<<"\t"<<parType( npar, (params->PARX)[j] )<<parNum( npar, (params->PARX)[j] )<<"\t";
					std::cout<<"\tUSTAB\tIT\n";
				}
				itpos = (itpos+1) % ithist;
				//
				it( itpos ) = pt.Continue( ds, (i == 0) && (params->getBranchSW() == TFHBSwitch) );
				//
				if( params->getStab() ) pt.Stability();
				ustabprev = ustab;
				if( trivial == 0 ) ustab = pt.UStab(); else if( trivial == 1 ) ustab = pt.UStabAUT(); else ustab = pt.UStabAUTRot();
				for( int j = 0; j < par.Size(); j++ ) par(j) = pt.getPar()(j);
				norm = pt.Norm();

				// console output
				std::cout<<"  "<<i+1<<"\t"<<norm<<"\t"<<par(params->getCp()-VarPAR0);
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
				pt.BinaryWrite( out, i );
				
				// branch output
// 				for( int j = 0; j < par.Size(); j++ ) ff<<par(j)<<"\t";
// 				ff<<"\t"<<norm<<"\t"<<pt.NormMX()<<"\t"<<ustab<<"\n";
// 				ff.flush();
				int itc = it( itpos );
				if( (itc > 3)&&(fabs(ds)/1.414 > params->getDsMin())&&(fabs(ds)/1.414 < params->getDsMax()) ) ds /= 1.414;
				if( (itc > 5)&&(fabs(ds)/2.0 > params->getDsMin())&&(fabs(ds)/2.0 < params->getDsMax()) ) ds /= 2.0;
				bool decr = true;
				for( int l=0; l < it.Size(); l++ ) if( it(l) > 3 ) decr = false;
				if( decr&&(fabs(ds)*1.414 > params->getDsMin())&&(fabs(ds)*1.414 < params->getDsMax()) ) ds *= 1.414;
				if( (itc >= params->getNItC())&&(fabs(ds)/2.0 < params->getDsMin()) )
				{
					P_MESSAGE("reached minimum stepsize (DSMIN)");
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
			Vector meshint(params->getNInt1() + 1), meshdeg( params->getNDeg1() + 1 );
			for( int i=0; i<meshint.Size(); i++ ) meshint(i) = (double)i/(params->getNInt1());
			for( int i=0; i<meshdeg.Size(); i++ ) meshdeg(i) = (double)i/(params->getNDeg1());

			// getting the sol and tangents
			pt.SwitchTRSol( Sol, meshint, meshdeg );
			pt.SwitchTFTRTan( TRe, TIm, alpha, meshint, meshdeg );

			// getting the parameters
			for( int j = 0; j < par.Size(); j++ ) par(j) = pt.getPar()(j);

			// destroy point, construct PointTR
			delete pt_ptr; pt_ptr = 0;
			PointTR pttr( sys, eqn, var, params->getNDeg1(),params->getNDeg2(), params->getNInt1(), params->getNInt2() );

			// construct the solution tangent from the eigenvectors
			// these next three functions could be only one
			pttr.ImportSol( Sol );
			pttr.ImportTan( TRe, TIm, alpha );
			pttr.Start( params->getDsStart() );
			pttr.setPar( par );
			
			pttr.setRho( alpha/(2.0*M_PI) );
			pttr.setCont( params->getCp()-VarPAR0 );

			double ds = params->getDs();
			std::ostringstream fdata, fidx;
			for( int i = 0; i < params->getSteps(); i++ )
			{
				pttr.Continue( ds, false );
				
				// write out the results
				for( int j = 0; j < npar; j++ ) std::cout<<par(j)<<"\t";
				std::cout<<std::endl;
				for( int j = 0; j < npar; j++ ) par(j) = pttr.getPar()(j);
// 				for( int j = 0; j < npar; j++ ) ff<<par(j)<<"\t";
// 				ff<<pttr.Norm()<<"\n";
// 				ff.flush();
				std::ostringstream fdata, fidx;
				fdata << "sol-" << i << ".dat";
				fidx << "sol-" << i << ".idx";
				pttr.SaveSol( fdata.str().c_str(), fidx.str().c_str() );
			}
		}
		// **********************************************************************************************************
		delete pt_ptr;
	}
	catch( pddeException ex )
	{
// 		QErrorMessage error;
// 		error.exec();
// 		error.showMessage( QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()) );
		std::cout<<ex.file<<":"<<ex.line<<" "<<ex.message.message;
		return;
	}
}
