// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "plot.h"
 
using namespace std;
extern "C" {

  /* LAPACK */
  int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
      a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl,
      integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work,
      integer *lwork, integer *info, ftnlen jobvl_len, ftnlen jobvr_len);

  int dgesv_( integer *n, integer *nrhs, doublereal *a, integer *lda,
      integer *ipiv, doublereal *b, integer *ldb, integer *info);

  int dgesvx_(char *fact, char *trans, integer *n, integer *nrhs,
              doublereal *a, integer *lda, doublereal *af, integer *ldaf,
              integer *ipiv, char *equed, doublereal *r__, doublereal *c__,
              const doublereal *b, integer *ldb, doublereal *x, integer *ldx,
              doublereal *rcond, doublereal *ferr, doublereal *berr,
              doublereal *work, integer *iwork, integer *info,
              ftnlen fact_len, ftnlen trans_len, ftnlen equed_len);

}

// **************************************************************************//
//                                                                           //
// **********************          Matrix Class        **********************//
//                                                                           //
// **************************************************************************//


void Matrix::StrPlot( GnuPlot& pl )
{
	pl.Plot( 0, "with points ps 0.1" );
	
	for( int i = 0; i < this->r; i++)
	{
		for( int j = 0; j < this->c; j++)
		{
			if( (*this)(i,j) != 0.0 )
			{
				pl.AddData( 0, i, j );
			}
		}
	}
	pl.Show();
}

void Matrix::Eigval( Vector& wr, Vector& wi )
{
	if( (this->r != this->c)||(this->r != wr.Size())||(this->r != wi.Size()) )
	{
		cout << "eigval: wrong dimensions" << '\n';
		exit(-1);
	}
	const Matrix& AA = *this;

	if( this->r == 1 )
	{
		wr(0) = AA(0,0);
		wi(0) = 0.0;
		return;
	}
	if( this->r == 2 )
	{
		double tr = AA(0,0)+AA(1,1), det = AA(0,0)*AA(1,1) - AA(0,1)*AA(1,0);
		if( tr*tr - 4.0*det >= 0 )
		{
			wr(0) = 0.5*(tr-sqrt(tr*tr - 4.0*det));
			wi(0) = 0.0;
			wr(1) = 0.5*(tr+sqrt(tr*tr - 4.0*det));
			wi(1) = 0.0;
		}else
		{
			wr(0) = 0.5*tr;
			wi(0) = 0.5*sqrt(-tr*tr + 4.0*det);
			wr(1) = 0.5*tr;
			wi(1) = -0.5*sqrt(-tr*tr + 4.0*det);
		}
		return;
	}
	
	char    jobvl = 'N', jobvr = 'N';
	Matrix  vl(this->r,1), vr(this->r,1);
	integer n = this->r, lda = this->r;
	integer ldvl = 1, ldvr = 1;
	Matrix  work(this->r,4);
	integer lwork = 4*(this->r);
	integer info;

	dgeev_( &jobvl, &jobvr, &n, this->m, &lda, 
	        wr.Pointer(), wi.Pointer(), vl.m, &ldvl, vr.m, &ldvr,
	        work.m, &lwork, &info, 1, 1);
	if(info != 0) cout << "Error code" << info << '\n';
}

void Matrix::Eigval( Vector& wr, Vector& wi, Matrix& vl, Matrix& vr )
{
	if( (this->r != this->c)||(this->r != wr.Size())||(this->r != wi.Size()) )
	{
		cout << "eigval: wrong dimensions" << '\n';
		exit(-1);
	}
	
	if( this->r == 1 )
	{
		wr(0) = *(this->m);
		wi(0) = 0.0;
		vl(0,0) = 1.0;
		vr(0,0) = 1.0;
		return;
	}
	
	char    jobvl = 'V',jobvr = 'V';
	integer n = this->r, lda = this->r;
	integer ldvl = vl.Row(), ldvr = vr.Row();
	Matrix  work(this->r,10);
	integer lwork = 10*(this->r);
	integer info;

	dgeev_( &jobvl, &jobvr, &n, this->m, &lda, 
	        wr.Pointer(), wi.Pointer(), vl.m, &ldvl, vr.m, &ldvr,
           work.m, &lwork, &info, 1, 1);
	if(info != 0) cout << "Error code" << info << '\n';
}

// MatFact routines

void MatFact::Fact()
{
	if( (this->r != this->c) )
	{
		cout << "MatFact::Fact wrong dimensions" << '\n';
		throw(-1);
	}

	char FACT='N',trans='N',equed='N';
	integer n=this->r;
	integer nrhs=0;
  
	dgesvx_( &FACT, &trans, &n, &nrhs, this->m, &n, this->mf, &n,
	         ipiv, &equed, NULL, NULL, NULL, &n, NULL, &n, 
	         &rcond, NULL, NULL, work, iwork, &info, 1, 1, 1);
	
	// cout<<"work(1): "<<work[0]<<" rcond< "<<rcond<<"\n";
	if(info != 0) cout << "MatFact::Fact: Error code" << info << '\n';
	fact = true;
}

void MatFact::Solve( Vector& x, const Vector& b, bool trans )
{
	if( this->r != this->c ) { cout<<"MatFact::Solve not a square matrix\n"; throw(1); return; }
	if( (b.Size() != this->c)||(b.Size() != x.Size()) ) { cout<<"MatFact::Solve Vector sizes differ\n"; throw(1); return; }
	
	
	const Matrix& AA = *this;
	double det;
	switch( this->r )
	{
		case 1:
			x(0) = b(0)/AA(0,0);
			break;
		case 2:
			det = 1.0/( AA(0,0)*AA(1,1) - AA(0,1)*AA(1,0) );
			if( !trans )
			{
				x(0) = det*( AA(1,1)*b(0) - AA(0,1)*b(1) );
				x(1) = det*( -AA(1,0)*b(0) + AA(0,0)*b(1) );
			}else
			{
				x(0) = det*( AA(1,1)*b(0) - AA(1,0)*b(1) );
				x(1) = det*( -AA(0,1)*b(0) + AA(0,0)*b(1) );
			}
			break;
		default:
			if( !fact ) Fact();

			char FACT='F', equed='N', TRANS;
			if( trans ) TRANS = 'T'; else TRANS = 'N';
			integer n=this->r;
			integer nrhs=1;
			
			dgesvx_( &FACT, &TRANS, &n, &nrhs, this->m, &n, this->mf, &n,
			         ipiv, &equed, NULL, NULL, b.v, &n, x.Pointer(), &n, 
			         &rcond, &ferr, &berr, work, iwork, &info, 1, 1, 1);
			if(info != 0) cout << "MatFact::Solve: Error code" << info << '\n';
			break;
	}
// 	Vector err(b);
// 	AXpY( err, x, b, 1.0, -1.0, trans );
// 	cout<<"dim: "<<this->r<<" b: "<<sqrt(b*b)<<" e: "<<sqrt(err*err)<<"\n";
}

void MatFact::Solve( Matrix& x, const Matrix& b, bool trans )
{
	if( this->r != this->c ) { cout<<"MatFact::Solve not a square matrix\n"; throw(1); return; }
	if( b.Col() != x.Col() ) { cout<<"MatFact::Solve Matrix columns differ\n"; throw(1); return; }
	if( (b.Row() != this->c)||(b.Row() != x.Row()) ) { cout<<"MatFact::Solve Matrix rows differ\n"; throw(1); return; }
  
	const Matrix& AA = *this;
	double det;
	switch( this->r )
	{
		case 1:
			for( int i=0; i< b.Col(); i++ )
			{
				x(0,i) = b(0,i)/AA(0,0);
			}
			break;
		case 2:
			det = 1.0/( AA(0,0)*AA(1,1) - AA(0,1)*AA(1,0) );
			if( !trans )
			{
				for( int i=0; i< b.Col(); i++ )
				{
					x(0,i) = det*( AA(1,1)*b(0,i) - AA(0,1)*b(1,i) );
					x(1,i) = det*( -AA(1,0)*b(0,i) + AA(0,0)*b(1,i) );
				}
			}else
			{
				for( int i=0; i< b.Col(); i++ )
				{
					x(0,i) = det*( AA(1,1)*b(0,i) - AA(1,0)*b(1,i) );
					x(1,i) = det*( -AA(0,1)*b(0,i) + AA(0,0)*b(1,i) );
				}
			}
			break;
		default:
			if( !fact ) Fact();
			
			char FACT='F', equed='N', TRANS;
			if( trans ) TRANS = 'T'; else TRANS = 'N';
			double *fberr = new double[2*b.Col()+1]; 
			integer n=this->r;
			integer nrhs=b.Col();
			
			dgesvx_( &FACT, &TRANS, &n, &nrhs, this->m, &n, this->mf, &n,
			         ipiv, &equed, NULL, NULL, b.m, &n, x.m, &n, 
			         &rcond, fberr, fberr+b.Col(), work, iwork, &info, 1, 1, 1);
	
			if(info != 0) cout << "MatFact::Solve: Error code" << info << '\n';
			delete[] fberr;
			break;
	}
}
