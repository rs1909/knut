// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "pderror.h"
#include "system.h"
#include "point.h"
#include "torpoint.h"
//#include "parameters.h"
#include "constants.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

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

inline void parNamePrint( Vector& /*par*/, int npar, Array1D<Var>& var )
{
	for( int j = 1; j < var.Size(); j++ ) std::cout<<"\t"<<parType( npar, var(j) - VarPAR0 )<<parNum( npar, var(j) - VarPAR0 )<<"\t";
}

inline void parValuePrint( Vector& par, int /*npar*/, Array1D<Var>& var )
{
	for( int j = 1; j < var.Size(); j++ ) std::cout<<"\t"<<par( var(j) - VarPAR0 );
}


int main( int argc, const char** argv )
{
	// parameters
	NConstants*  params = 0;
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
					params = new NConstants;
					params->loadFile( constFile );
					params->printFile( std::cout );
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
				 #ifdef HAVE_CONFIG_H
					std::cout<<"This is "<<PACKAGE_NAME<<" version "<<PACKAGE_VERSION<<" ("<<PACKAGE_REVISION<<")\n";
				 #endif
					exit(0);
					break;
				default:
					P_MESSAGE("Unexpected command line argument.\n");
					break;
			}
		}else
		{
			P_MESSAGE("Unexpected command line argument.\n");
		}
	}
	

	if( params == 0 )
	{
		std::cout<<"Error: missing constants file.\n";
		exit(1);
	}
	if( (inFile == 0) && (params->getLabel() != 0) )
	{
		std::cout<<"Error: missing input file.\n";
		exit(1);
	}
	if( outFile == 0)
	{
		outFile ="out.pdde";
		std::cout<<"Warning: missing output file, using \""<<outFile<<"\" instead.\n";
	}
	if( branchFile == 0)
	{
		branchFile ="branch";
		std::cout<<"Warning: missing branch file, using \""<<branchFile<<"\" instead.\n";
	}
	
	// **********************************************************************************************************
	
	System sys( params->getSysName() );
	if( sys.ndim() == 0 ) P_MESSAGE("Number of dimensions are set to zero.");

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
	
	
	//-----------------------------------------------------------------------------------------------------------
	//
	// END of initialization
	//
	//-----------------------------------------------------------------------------------------------------------
	
	// just a block to contain pt, which eats too much memory
	try
	{
		int trivial = params->toEqnVar( sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, testFN );
		const int npar = sys.npar();
	
// 	for( int i=0; i<eqn.Size(); i++ ) std::cout<<EqnToStr( eqn(i) )<<", ";
// 	std::cout<<'\n';
// 	for( int i=0; i<var.Size(); i++ ) std::cout<<VarToStr( var(i) )<<", ";
// 	std::cout<<'\n';

		Point* pt_ptr = new Point( sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul(), params->getNMat() );
		Point& pt = *pt_ptr;
		
		pt.setContIter( params->getNItC() );
		pt.setRefIter( params->getNItR() );
		pt.setStartIter( params->getNItS() );
		pt.setRefEps( params->getEpsR() );
		pt.setContEps( params->getEpsC() );
		pt.setStartEps( params->getEpsS() );
		pt.setCont( params->getCp()-VarPAR0 );
		
		// setting the symmetric components
		Array1D<int> sre(params->getNSym()), sim(params->getNSym());
		for( int i=0; i<sre.Size(); ++i ) { sre(i) = params->getSymRe(i); sim(i) = params->getSymIm(i); }
		if( params->getNSym() != 0 ) pt.setSym( sre, sim );
		
		std::cout<<std::scientific;
		std::cout.precision(6);
		
		// load the initial guess
		if( params->getLabel() != 0 )
		{
			std::ifstream istr( inFile );
			for( int i=0; i<params->getLabel()-1; i++ )
			{
				pt.ReadNull( istr );
			}
			pt.Read( istr );
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
			parNamePrint( std::cout, npar, params->getCp(), var );
			std::cout<<"\n";
			parValuePrint( std::cout, par, params->getCp(), var, 0, pt.Norm(), 0, 0 );
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
				pt.SwitchTFLP( params->getDsStart());
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
					parNamePrint( std::cout, npar, params->getCp(), var );
					std::cout<<"\n";
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
				parValuePrint( std::cout, par, params->getCp(), var, i, norm, ustab, it( itpos ) );
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
				for( int j = 0; j < par.Size(); j++ ) ff<<par(j)<<"\t";
				ff<<"\t"<<norm<<"\t"<<pt.NormMX()<<"\t"<<ustab<<"\n";
				ff.flush();
				int itc = it( itpos );
				if( (itc > 3)&&(fabs(ds)/1.414 > params->getDsMin())&&(fabs(ds)/1.414 < params->getDsMax()) ) ds /= 1.414;
				if( (itc > 5)&&(fabs(ds)/2.0 > params->getDsMin())&&(fabs(ds)/2.0 < params->getDsMax()) ) ds /= 2.0;
				bool decr = true;
				for( int l=0; l < it.Size(); l++ ) if( it(l) > 3 ) decr = false;
				if( decr&&(fabs(ds)*1.414 > params->getDsMin())&&(fabs(ds)*1.414 < params->getDsMax()) ) ds *= 1.414;
				if( (itc >= params->getNItC())&&(fabs(ds)/2.0 < params->getDsMin()) )
				{
					P_MESSAGE("No convergence. The minimum arclength step size (DSMIN) has been reached.");
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
		delete pt_ptr;
	}
	catch( pddeException ex )
	{
		std::cout<<ex.file<<":"<<ex.line<<" "<<ex.message.message;
		exit(-1);
	}
	
	delete params;
}
