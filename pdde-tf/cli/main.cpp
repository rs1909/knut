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
#include "parameters.h"
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
	ConstFile*   params = 0;
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
					params = new ConstFile;
					params->getParams( constFile );
					params->printParams( );
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
					std::cout<<"This is "<<PACKAGE_NAME<<" version "<<PACKAGE_VERSION<<" ("<<PKG_REV<<")\n";
				 #endif
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
	
	int trivial = params->toEqnVar( sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, testFN );
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
			delete pt_ptr; pt_ptr = 0;
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
		delete pt_ptr;
	}

	delete params;
}
