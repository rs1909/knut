// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef SYSTEM_H
#define SYSTEM_H

#include "matrix.h"

#ifndef _WIN32

extern "C" {
#include <dlfcn.h>
}

typedef void* tdlhand;

#else

#include <windows.h>

typedef HMODULE tdlhand;

#endif

inline tdlhand tdlopen( const char* fname )
{
 #ifndef _WIN32
	return dlopen(fname, RTLD_NOW);
 #else
	return LoadLibrary( fname );
 #endif
}

inline void* tdlsym( tdlhand h, const char* sname )
{
 #ifndef _WIN32
	return dlsym( h, sname );
 #else
	return (void*) GetProcAddress( h, sname );
 #endif
}

inline int tdlclose( tdlhand h )
{
 #ifndef _WIN32
	return dlclose( h );
 #else
	return (int) FreeLibrary( h );
 #endif
}

inline const char* tdlerror( )
{
 #ifndef _WIN32
	return dlerror( );
 #else
	DWORD errcode = GetLastError( );
	SetLastError( 0 );
	if( errcode != 0 ) return "Windows System Error\n";
	else return 0;
 #endif
}

class System
{
	public:
	
		System( char* shobj );
		~System();
		
		int    ndim() { return (*v_ndim)( ); }
		int    npar() { return (*v_npar)( ); }
		int    ntau() { return (*v_ntau)( ); }
		void   tau( Vector& out, double t, const Vector& par ) 
		{
			(*v_tau)( out, t, par );
		}
		void   dtau( Vector& out, double t, const Vector& par, int vp )
		{
			(*v_dtau)( out, t, par, vp );
		}
		void   rhs( Vector& out, double t, const Matrix& x, const Vector& par )
		{
			(*v_rhs)( out, t, x, par );
		}
		void   deri( Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v )
		{
			if( nderi == 2 ) (*v_deri)( out, t, x, par, nx, vx, np, vp, v );
			else if( nderi == 0 ) discrderi( out, t, x, par, nx, vx, np, vp, v );
			else if( nderi == 1 && ( (nx == 1 && np == 0) || (nx == 0 && np == 1) ) ) (*v_deri)( out, t, x, par, nx, vx, np, vp, v );
			else discrderi( out, t, x, par, nx, vx, np, vp, v );
		}
		void   stpar( Vector& par )
		{
			(*v_stpar)( par );
		}
		void   stsol( Vector& out, double t )
		{
			(*v_stsol)( out, t );
		}
		
		void   discrderi( Matrix &out, double t, const Matrix& xx, const Vector& par,
		                  int nx, const int* vx, int np, const int* vp, const Matrix& vv );
		
	private:
	
		tdlhand handle;
		const char* error;
		
		int     nderi;
		
		Vector  f, f_eps;
		Vector  f2, f_eps2;
		Matrix  xx_eps;
		Matrix  xx_eps2;
		Vector  par_eps;
		Matrix  dxx2, dxx_eps2;
		Vector  vt;
		
		int     (*v_ndim)();
		int     (*v_npar)();
		int     (*v_ntau)();
		int     (*v_nderi)();
		void    (*v_tau)( Vector& out, double t, const Vector& par );
		void    (*v_dtau)( Vector& out, double t, const Vector& par, int vp );
		void    (*v_rhs)( Vector& out, double t, const Matrix& x, const Vector& par );
		void    (*v_deri)( Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v );
		void    (*v_stpar)( Vector& par );
		void    (*v_stsol)( Vector& out, double t );
};

#endif
