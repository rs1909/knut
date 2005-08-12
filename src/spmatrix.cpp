// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages packages root directory
//
// ------------------------------------------------------------------------- //

#include <iostream>
#include "spmatrix.h" // itt mar meghivtuk <vector>-t
#include <cmath>

using namespace std;

extern "C" {

  /* ARPACK */
  int dnaupd_(integer *ido, char *bmat, integer *n, char *
      which, integer *nev, double *tol, double *resid, integer *ncv,
      double *v, integer *ldv, integer *iparam, integer *ipntr,
      double *workd, double *workl, integer *lworkl, integer *info,
      ftnlen bmat_len, ftnlen which_len);

  int dneupd_(logical *rvec, char *howmny, logical *select,
      double *dr, double *di, double *z__, integer *ldz,
      double *sigmar, double *sigmai, double *workev, 
      char *bmat, integer *n, char *which, integer *nev, double *tol,
      double *resid, integer *ncv, double *v, integer *ldv, 
      integer *iparam, integer *ipntr, double *workd, double *workl,
      integer *lworkl, integer *info, ftnlen howmny_len, ftnlen bmat_len,
      ftnlen which_len);
}

// **************************************************************************//
//                                                                           //
// **********************        SpMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

void SpMatrix::Init( char F, int n_, int nz )
{
	if( (F != 'R')&&(F != 'C') ) cout<<"SpMatrix::CONSTRUCTOR: invalid format specification.\n";
	format = F;
	n = 0;
	size = nz;
	Ap = new int[n_+1];
	Ai = new int[nz];
	Ax = new double[nz];
	for( int i = 0; i < nz; i++ ){ Ai[i] = 0; Ax[i] = 0.0; }
	for( int i = 0; i < n_+1; i++ ){ Ap[i] = 0; }
}

void SpMatrix::Init( const SpMatrix& M )
{
	format = M.format;
	n = M.n;
	size = M.size;
	
	Ap = new int[M.n+1];
	Ai = new int[size];
	Ax = new double[size];
	for( int i = 0; i < size; i++ ){ Ai[i] = M.Ai[i]; Ax[i] = M.Ax[i]; }
	for( int i = 0; i < n+1; i++ ){ Ap[i] = M.Ap[i]; }
}

SpMatrix::~SpMatrix()
{
	delete []Ap;
	delete []Ai;
	delete []Ax;
}

void SpMatrix::Check()
{
	for( int j = 0; j < n; j++ )
	{
		for( int k = Ap[j]+1; k < Ap[j+1]; k++ )
		{
			if( Ai[k-1] >= Ai[k] ) cout<<"SpMatrix::Check(): "<<j<<","<<Ai[k-1]<<"\n";
		}
	}
}

void SpMatrix::Clear()
{
	
	for( int i = 0; i < n+1; i++ ){ Ap[i] = 0; }
	n = 0;
	for( int i = 0; i < size; i++ ){ Ai[i] = 0; Ax[i] = 0.0; }
}

void SpMatrix::Clear(char F)
{
	
	for( int i = 0; i < n+1; i++ ){ Ap[i] = 0; }
	n = 0;
	format = F;
	for( int i = 0; i < size; i++ ){ Ai[i] = 0; Ax[i] = 0.0; }
}

// // multiplictation for row ordered OR transposed multiplication for column ordered
// void inline spmulR( double* out, const double* in, const int* Ap, const int* Ai, const int* Ax )
// {
// 	for( int j = 0; j < this->n; j++ )
// 	{
// 		out[j] = 0.0;
// 		for( int k = Ap[j]; k < Ap[j+1]; k++ )
// 		{
// #ifdef DEBUG
// 			if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AX bd\n"; throw(1); }
// #endif
// 			out[j] += Ax[k]*in[Ai[k]];
// 		}
// 	}
// }
// 
// // multiplictation for column ordered OR transposed multiplication for row ordered
// void inline spmulC( double* out, const double* in, const int* Ap, const int* Ai, const int* Ax )
// {
// 	for( int i = 0; i < this->n; i++ ) out[i] = 0.0;
// 	for( int j = 0; j < this->n; j++ )
// 	{
// 		for( int k = Ap[j]; k < Ap[j+1]; k++ )
// 		{
// #ifdef DEBUG
// 			if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AX bd\n"; throw(1); }
// #endif
// 			out[Ai[k]] += Ax[k]*in[j];
// 		}
// 	}
// }

// redefinition of purely virtual functions

void SpMatrix::AX( double* out, const double* in, double alpha, bool trans ) const
{
	// ez itt rossz, mert nincs megadva a vektorok merete, vagy a matrix masik dimenzioja
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int i = 0; i < this->n; i++ ) out[i] = 0.0;
			for( int j = 0; j < this->n; j++ )
			{
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AX bd\n"; throw(1); }
#endif
					out[Ai[k]] += Ax[k]*in[j];
				}
			}
		}else
		{
			// std::cout<<"SpAXt";
			for( int j = 0; j < this->n; j++ )
			{
				out[j] = 0.0;
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AX bd\n"; throw(1); }
#endif
					out[j] += Ax[k]*in[Ai[k]];
				}
			}
		}
	}else
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int i = 0; i < this->n; i++ ) out[i] = 0.0;
			for( int j = 0; j < this->n; j++ )
			{
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AX bd\n"; throw(1); }
#endif
					out[Ai[k]] += alpha*Ax[k]*in[j];
				}
			}
		}else
		{
			for( int j = 0; j < this->n; j++ )
			{
				out[j] = 0.0;
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AX bd\n"; throw(1); }
#endif
					out[j] += alpha*Ax[k]*in[Ai[k]];
				}
			}
		}
	}
}

void SpMatrix::AXpY( double* out, const double* in, const double* C, double alpha, double beta, bool trans ) const
{
	if( alpha == 1.0 ){
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int i = 0; i < this->n; i++ ) out[i] = beta*C[i];
			for( int j = 0; j < this->n; j++ )
			{
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AXpY bd\n"; throw(1); }
#endif
					out[Ai[k]] += Ax[k]*in[j];
				}
			}
		}else
		{
			for( int j = 0; j < this->n; j++ )
			{
				out[j] = beta*C[j];
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AXpY bd\n"; throw(1); }
#endif
					out[j] += Ax[k]*in[Ai[k]];
				}
			}
		}
	}else
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int i = 0; i < this->n; i++ ) out[i] = beta*C[i];
			for( int j = 0; j < this->n; j++ )
			{
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AXpY bd\n"; throw(1); }
#endif
					out[Ai[k]] += alpha*Ax[k]*in[j];
				}
			}
		}else
		{
			for( int j = 0; j < this->n; j++ )
			{
				out[j] = beta*C[j];
				for( int k = Ap[j]; k < Ap[j+1]; k++ )
				{
#ifdef DEBUG
					if( (Ai[k] < 0)||(Ai[k] >= this->n) ){ std::cout<<"Sp::AXpY bd\n"; throw(1); }
#endif
					out[j] += alpha*Ax[k]*in[Ai[k]];
				}
			}
		}
	}
}

void SpMatrix::AX( Matrix& out, const Matrix& in, double alpha, bool trans ) const
{
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int i = 0; i < out.Row(); i++ ) out(i,c) = 0.0;
				for( int j = 0; j < this->n; j++ )
				{
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(Ai[k],c) += Ax[k]*in(j,c);
					}
				}
			}
		}else
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int j = 0; j < this->n; j++ )
				{
					out(j,c) = 0.0;
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(j,c) += Ax[k]*in(Ai[k],c);
					}
				}
			}
		}
	}else
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int i = 0; i < out.Row(); i++ ) out(i,c) = 0.0;
				for( int j = 0; j < this->n; j++ )
				{
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(Ai[k],c) += alpha*Ax[k]*in(j,c);
					}
				}
			}
		}else
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int j = 0; j < this->n; j++ )
				{
					out(j,c) = 0.0;
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(j,c) += alpha*Ax[k]*in(Ai[k],c);
					}
				}
			}
		}
	}
}

void SpMatrix::AXpY( Matrix& out, const Matrix& in, const Matrix& C, double alpha, double beta, bool trans ) const
{
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int i = 0; i < in.Row(); i++ ) out(i,c) = beta*C(i,c);
				for( int j = 0; j < in.Row(); j++ )
				{
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(Ai[k],c) += Ax[k]*in(j,c);
					}
				}
			}
		}else
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int j = 0; j < this->n; j++ )
				{
					out(j,c) = beta*C(j,c);
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(j,c) += Ax[k]*in(Ai[k],c);
					}
				}
			}
		}
	}else
	{
		if( ((format == 'C')&&(trans == false))||((format == 'R')&&(trans == true)) )
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int i = 0; i < in.Row(); i++ ) out(i,c) = beta*C(i,c);
				for( int j = 0; j < in.Row(); j++ )
				{
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(Ai[k],c) += alpha*Ax[k]*in(j,c);
					}
				}
			}
		}else
		{
			for( int c = 0; c < in.Col(); c++ )
			{
				for( int j = 0; j < this->n; j++ )
				{
					out(j,c) = beta*C(j,c);
					for( int k = Ap[j]; k < Ap[j+1]; k++ )
					{
						out(j,c) += alpha*Ax[k]*in(Ai[k],c);
					}
				}
			}
		}
	}
}
// end of the purely virtual functions

void SpMatrix::Swap()
{
	int *Rp = new int[n+1];
	int *Ri = new int[size];
	double *Rx = new double[size];
	
	umfpack_di_transpose( n, n, Ap, Ai, Ax, 0, 0, Rp, Ri, Rx );
	
	delete []Ap; 
	delete []Ai; 
	delete []Ax;
	
	Ap = Rp; 
	Ai = Ri; 
	Ax = Rx;
	
	if( format == 'C' ) format = 'R'; else format = 'C';
}

void SpMatrix::StrPlot( GnuPlot& pl )
{
	pl.SetPointSize( 0.8 );
	pl.Plot(0, "with points");
	
	for( int i = 0; i < n; i++ )
	{
		for( int j = Ap[i]; j < Ap[i+1]; j++ )
		{
			if(Ax[j] != 0.0) pl.AddData( 0, (double)i, (double)Ai[j] ); 
		}
	}
	pl.Show();
}

void SpMatrix::Print()
{
	for( int i = 0; i < n; i++ )
	{
		for( int j = Ap[i]; j < Ap[i+1]; j++ )
		{
			cout<<"("<<i<<","<<Ai[j]<<")="<<Ax[j]<<'\t'; 
			cout.flush();
		}
		cout<<"\n";
	}
}

// **************************************************************************//
//                                                                           //
// **********************         SpFact Class         **********************//
//                                                                           //
// **************************************************************************//

SpFact::SpFact( char F, int nn_, int nz ) : SpMatrix( F, nn_, nz ), BaseFact()
{
	fact = false;

	Numeric = 0;

	Wi = new int[5*nn_+1];
	W  = new double[10*nn_+1];
  
	umfpack_di_defaults(Control);
	Control[UMFPACK_IRSTEP] = 0;
	Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
// 	Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
}

SpFact::SpFact( SpMatrix& M ) : SpMatrix( M ), BaseFact()
{
	fact = false;
  
	Numeric = 0;

	Wi = new int[5*n+1];
	W  = new double[10*n+1];
  
	umfpack_di_defaults(Control);
	Control[UMFPACK_IRSTEP] = 0;
	Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
//	Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
}

SpFact::~SpFact()
{
	if( Numeric != 0 ) umfpack_di_free_numeric( &Numeric );
	Numeric = 0;

	delete[] W;
	delete[] Wi;
}

void SpFact::Clear( )
{ 
	SpMatrix::Clear( ); 
	if( Numeric != 0 )
	{
		if( !fact ){ std::cout<<"GEBASZ\n"; throw(1); }
		umfpack_di_free_numeric( &Numeric ); 
		Numeric = 0;
	}
	fact = false;
}

void SpFact::Clear( char F )
{ 
	SpMatrix::Clear( F ); 
	if( Numeric != 0 )
	{
		if( !fact ){ std::cout<<"GEBASZ\n"; throw(1); }
		umfpack_di_free_numeric( &Numeric ); 
		Numeric = 0;
	}
	fact = false;
}


void SpFact::Fact()
{
	if( !fact )
	{
		void  *Symbolic = 0;
		if( Numeric != 0 )
		{
			cout<<"GEMASZ \n";
			exit(-1);
		}
		int   status;

		status = umfpack_di_symbolic(n, n, this->Ap, this->Ai, this->Ax, &Symbolic, Control, 0);
		if( status != 0 )
		{
			cout<<"SpFact::SpFact(Spmatrix& ): Symbolic "<<status<<"\n";
			throw(-1);
		}
		status = umfpack_di_numeric(this->Ap, this->Ai, this->Ax, Symbolic, &Numeric, Control, 0);
		if( status != 0 ) {
			cout<<"SpFact::SpFact(Spmatrix& ): Numeric "<<status<<"\n";
			throw(-1);
		}
		umfpack_di_free_symbolic(&Symbolic); Symbolic = 0;
		fact = true;
	}
}

void SpFact::Solve( double* x, double* b, bool trans )
{
	if( !fact ) Fact();
	int status,sys;
	
	if( format == 'C' )
	{ 
		if( trans ) sys = UMFPACK_At; else sys = UMFPACK_A; 
	}else
	{
		if( trans ) sys = UMFPACK_A; else sys = UMFPACK_At;
		if( format != 'R' )
		{
			cout<<"SpFact::Solve(Vector& , Vector& ): Wrong format";
			return; 
		}
	}
	status = umfpack_di_wsolve( sys, Ap, Ai, Ax, x, b, Numeric, Control, 0, Wi, W );
}

void SpFact::Solve( Vector& x, const Vector& b, bool trans )
{
	if( !fact ) Fact();
	int status,sys;
	
	if( format == 'C' )
	{ 
		if( trans ) sys = UMFPACK_At; else sys = UMFPACK_A; 
	}else
	{
		if( trans ) sys = UMFPACK_A; else sys = UMFPACK_At;
		if( format != 'R' )
		{
			cout<<"SpFact::Solve(Vector& , Vector& ): Wrong format";
			return; 
		}
	}
  
	status = umfpack_di_wsolve( sys, Ap, Ai, Ax, x.Pointer(), b.v, Numeric, Control, 0, Wi, W );
}

void SpFact::Scale( Vector& x, const Vector& b )
{
	if( !fact ) Fact();
	umfpack_di_scale( x.Pointer(), b.v, Numeric );
}

void SpFact::Solve( Matrix& x, const Matrix &b, bool trans )
{
	if( !fact ) Fact();
	int status,sys;
	if( format == 'C' )
	{ 
		if( trans ) sys = UMFPACK_At; else sys = UMFPACK_A; 
	}else
	{
		if( trans ) sys = UMFPACK_A; else sys = UMFPACK_At;
		if( format != 'R' )
		{
			cout<<"SpFact::Solve(Vector& , Vector& ): Wrong format";
			return;
		}
	}
  
	for( int i = 0; i < b.c; i++ )
	{
		status = umfpack_di_wsolve( sys, Ap, Ai, Ax, x.m+i*x.r, b.m+i*b.r, Numeric, Control, 0, Wi, W);
	}
}


void StabMatrix::Eigval( Vector& wr, Vector& wi )
{
	if( ( wr.Size() != wi.Size() )||( wr.Size() < 2 ) )
	{
		cout<<"Bad sizes of wr or wi\n";
		return;
	}
	integer IDO      = 0;
	char    BMAT     = 'I';
	integer N        = AI.Size() * A0.Size();
	char    WHICH[]  = {'L','M'};
	integer NEV      = wr.Size()-1;
	double  TOL      = 0.0;
//	 RESID was defined in the class
	integer NCV      = 2*NEV + 1;
	double* V        = new double[(N+1)*(NCV+1)];
	integer LDV      = N;
	integer IPARAM[] = {1,     // ISHIFT
	                    0,     // n.a.ff
	                    300,   // MXITER
	                    1,     // NB
	                    0,     // NCONV (out)
	                    0,     // IUPD n.a.
	                    1,     // MODE
	                    0,     // NP if there is shift ()
	                    0,     // NUMOP
	                    0,     // NUMOPB
	                    0,     // NUMREO
							  0, };  // padding
	integer IPNTR[14];
	double* WORKD     = new double[4*N+1]; //3*N
	integer LWORKL    = 3*NCV*NCV + 6*NCV;
	double* WORKL     = new double[LWORKL+1];
	integer INFO; // RESID is computed by the routine,
                         // but it would be better to use from a prev. run
	if( isINIT ) INFO = 1;
	else INFO = 0;

	double* tvec = new double[N+1];
	double* tvec2 = new double[N+1];
	
	// int it=0;
	do
	{
		dnaupd_( &IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, 
		         IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO, 1, 2 );
		if( ( IDO == 1 )||( IDO==-1 ) )
		{
			double* in = &(WORKD[IPNTR[0]-1]);
			double* out = &(WORKD[IPNTR[1]-1]);
			// these are the unit operators above the diagonal
			for( int i=0; i < AI.Size()-1; i++ )
			{
				for( int j=0; j<A0.Size(); j++ ) out[j + i*A0.Size()] = in[j + (i+1)*A0.Size()];
			}
			// the last row: multiplication and solution
			for( int j = 0; j < A0.Size(); j++ ) out[j+(AI.Size()-1)*A0.Size()] = 0.0;
			for( int i = 0; i < AI.Size(); i++ )
			{
				if( AI(AI.Size()-1-i).GetNZ() != 0 )
				{
					AI(AI.Size()-1-i).AX( tvec+i*A0.Size(), in+i*A0.Size(), 1.0, false );
					A0.Solve( tvec2+i*A0.Size(), tvec+i*A0.Size() );
					// summing up the last row
					for( int j=0; j<A0.Size(); j++ ) out[j+(AI.Size()-1)*A0.Size()] += tvec2[j+i*A0.Size()];
				}
			}
		}
	}
	while( ( IDO == 1 )||( IDO == -1 ) );
	if( IDO != 99 ){ std::cout<<"IDO = "<<IDO<<" is not expected\n"; throw(1); }
	delete[] tvec2;
	delete[] tvec;
	
	logical  RVEC = FALSE_;
	char     HOWMNY   = 'A';
	logical* SELECT   = new logical[NCV+1];
	double*  DR       = new double[NEV+2];
	double*  DI       = new double[NEV+2];
	double*  Z        = new double[N*(NEV+2)];
	integer  LDZ      = N;
	double   SIGMAR   = 0.0;
	double   SIGMAI   = 0.0;
	double*  WORKEV   = new double[4*NCV]; // 3*NCV

// 	std::cout<<"INFO1:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
// 	std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();
	
	if( true /*IPARAM[4] >= NEV*/ )
	{
		dneupd_( &RVEC, &HOWMNY, SELECT, DR, DI, Z, &LDZ,
		         &SIGMAR, &SIGMAI, WORKEV, 
		         &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, 
		         IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO, 1, 1, 2 );
// 		std::cout<<"INFO2:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
// 		std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();
		wr.Clear(); wi.Clear();
		if( IPARAM[4] > wr.Size() )
		{
			std::cout<<" more eigenvalues than expected: N="<<IPARAM[4]<<"\n"; std::cout.flush();
			IPARAM[4]=wr.Size();
		}
		for( int i = 0; i < IPARAM[4]; i++ )
		{
			wr(i) = DR[i];
			wi(i) = DI[i];
		}
	}
 
	// from the second
	delete[] WORKEV;
	delete[] Z;
	delete[] DI;
	delete[] DR;
	delete[] SELECT;
	// from the first
	delete[] WORKL;
	delete[] WORKD;
	delete[] V;
	// now we have a RESID
	isINIT = true;
}
