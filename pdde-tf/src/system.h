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
#else
#include <windows.h>
#endif

class System
{
	public:
	
		System( const std::string& shobj );
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
		
		typedef int  (*tp_sys_ndim)();
		typedef int  (*tp_sys_npar)();
		typedef int  (*tp_sys_ntau)();
		typedef int  (*tp_sys_nderi)();
		typedef void (*tp_sys_tau)( Vector& out, double t, const Vector& par );
		typedef void (*tp_sys_dtau)( Vector& out, double t, const Vector& par, int vp );
		typedef void (*tp_sys_rhs)( Vector& out, double t, const Matrix& x, const Vector& par );
		typedef void (*tp_sys_deri)( Matrix& out, double t, const Matrix& x, const Vector& par, int nx, const int* vx, int np, const int* vp, const Matrix& v );
		typedef void (*tp_sys_stpar)( Vector& par );
		typedef void (*tp_sys_stsol)( Vector& out, double t );
		
		typedef void (*FPTR)();
		union punned {
			void *obj;
			FPTR fun;
		};
		FPTR fptr( void * ptr );
	
	#ifndef _WIN32
		typedef void*   tdlhand;
	#else
		typedef HMODULE tdlhand;
	#endif
		
		tdlhand     tdlopen( const char* fname );
		void*       tdlsym( tdlhand h, const char* sname );
		int         tdlclose( tdlhand h );
		const char* tdlerror( );

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
		
		tp_sys_ndim   v_ndim;
		tp_sys_npar   v_npar;
		tp_sys_ntau   v_ntau;
		tp_sys_nderi  v_nderi;
		tp_sys_tau    v_tau;
		tp_sys_dtau   v_dtau;
		tp_sys_rhs    v_rhs;
		tp_sys_deri   v_deri;
		tp_sys_stpar  v_stpar;
		tp_sys_stsol  v_stsol;
};

inline System::FPTR System::fptr( void * ptr )
{
	punned tmp;
	tmp.obj = ptr;
	return tmp.fun;
}

inline System::tdlhand System::tdlopen( const char* fname )
{
 #ifndef _WIN32
	return dlopen(fname, RTLD_NOW);
 #else
	return LoadLibrary( fname );
 #endif
}

inline void* System::tdlsym( tdlhand h, const char* sname )
{
 #ifndef _WIN32
	return dlsym( h, sname );
 #else
	return (void*) GetProcAddress( h, sname );
 #endif
}

inline int System::tdlclose( System::tdlhand h )
{
 #ifndef _WIN32
	return dlclose( h );
 #else
	return (int) FreeLibrary( h );
 #endif
}

inline const char* System::tdlerror( )
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

#endif
