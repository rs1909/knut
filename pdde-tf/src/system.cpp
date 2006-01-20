// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cmath>
#include "system.h"

#include "matrix.h"
#include "pointtype.h"

extern "C" {
#include <stdlib.h>
#include <stdio.h>
}

#include <string>

System::System( char* shobj )
{
	std::string objname(shobj);
	if( objname.find('/') == std::string::npos ) objname.insert(0, "./");
	
	handle = tdlopen ( objname.c_str() );
	if (!handle)
	{
		std::cerr<<"Cannot open system definition file: "<<tdlerror()<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_ndim = (tp_sys_ndim) fptr( tdlsym(handle, "sys_ndim" ) );
	if ((error = tdlerror()) != NULL)
	{
		std::cerr<<"Cannot find sys_ndim(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_npar = (tp_sys_npar) fptr( tdlsym(handle, "sys_npar" ) );
	punned tmp; tmp.obj = tdlsym(handle, "sys_ndim" ); v_ndim = (int (*)())tmp.fun;
	if ((error = tdlerror()) != NULL)
	{
		std::cerr<<"Cannot find sys_npar(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_ntau = (tp_sys_ntau) fptr( tdlsym(handle, "sys_ntau" ) );
	if ((error = tdlerror()) != NULL)
	{
		std::cerr<<"Cannot find sys_ntau(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_nderi = (tp_sys_nderi) fptr( tdlsym(handle, "sys_nderi" ) );
	if ((error = tdlerror()) != NULL)
	{
		std::cerr<<"Cannot find sys_nderi(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_tau = (tp_sys_tau) fptr( tdlsym(handle, "sys_tau" ) );
	if ((error = tdlerror()) != 0 )
	{
		std::cerr<<"Cannot find sys_tau(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_dtau = (tp_sys_dtau) fptr( tdlsym(handle, "sys_dtau" ) );
	if ((error = tdlerror()) != 0 )
	{
		std::cerr<<"Cannot find sys_dtau(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_rhs = (tp_sys_rhs) fptr( tdlsym(handle, "sys_rhs" ) );
	if ((error = tdlerror()) != 0 )
	{
		std::cerr<<"Cannot find sys_rhs(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_deri = (tp_sys_deri) fptr( tdlsym(handle, "sys_deri" ) );
	if ((error = tdlerror()) != 0 )
	{
		std::cerr<<"Cannot find sys_deri(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_stpar = (tp_sys_stpar) fptr( tdlsym(handle, "sys_stpar" ) );
	if ((error = tdlerror()) != 0 )
	{
		std::cerr<<"Cannot find sys_stpar(): "<<error<<"\n";
		exit(1);
	}
	
	tdlerror();    /* Clear any existing error */
	v_stsol = (tp_sys_stsol) fptr( tdlsym(handle, "sys_stsol" ) );
	if ((error = tdlerror()) != 0 )
	{
		std::cerr<<"Cannot find sys_stsol(): "<<error<<"\n";
		exit(1);
	}
	
	nderi = (*v_nderi)( );
// 	std::cout<<"The order of supplied derivatives is "<<nderi<<".\n";
	f.Init( ndim() ), f_eps.Init( ndim() );
	f2.Init( ndim() ), f_eps2.Init( ndim() );
	xx_eps.Init( ndim(), 2 * ntau() + 1 );
	xx_eps2.Init( ndim(), 2 * ntau() + 1 );
	par_eps.Init( npar() + ParEnd );
	dxx2.Init( ndim(), ndim() ), dxx_eps2.Init( ndim(), ndim() );
	vt.Init( ndim() );
}

System::~System()
{
	tdlclose(handle);
}


static inline void AX( Vector & res, const Matrix& M, const Vector& v )
{
	for( int i=0; i<M.Row(); i++ )
	{
		for( int j=0; j<M.Col(); j++ )
		{
			res(i) = 0.0;
			res(i) += M(i,j)*v(j);
		}
	}
}

void System::discrderi( Matrix &out, double t, const Matrix& xx, const Vector& par, 
	       int nx, const int* vx, int np, const int* vp, const Matrix& vv )
{
	const double abs_eps_x1=1e-6;
	const double rel_eps_x1=1e-6;
	const double abs_eps_p1=1e-6;
	const double rel_eps_p1=1e-6;
	const double abs_eps_x2=2e-6;
	const double rel_eps_x2=2e-6;
	
	const int n = ndim();
	const int m = 2 * ntau() + 1;
	
	// f, f_eps, xx_eps
	// derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
	if( (nx == 1) && (np == 0) )
	{
		rhs( f, t, xx, par );
		for( int j = 0; j < n; j++ )
		{
			// xx_eps = xx;
			for( int q = 0; q < m; q++ ) for( int p = 0; p < n; p++ ) xx_eps(p,q) = xx(p,q);
			const double eps = abs_eps_x1 + rel_eps_x1 * fabs(xx(j,vx[0]));
			xx_eps(j,vx[0]) = xx(j,vx[0]) + eps;
			rhs( f_eps, t, xx_eps, par );
			for( int p = 0; p < n; p++ ) out(p,j) = ( f_eps(p) - f(p) ) / eps;
		}
	}
	// f, f_eps, par_eps
	// derivatives w.r.t. the parameters, purely, so this results a vector
	if( (nx == 0) && (np == 1) )
	{
		rhs( f, t, xx, par );
		par_eps = par;
		const double eps = abs_eps_p1 + rel_eps_p1*fabs(par(vp[0]));
		par_eps(vp[0]) = par(vp[0]) + eps;
		rhs( f_eps, t, xx, par_eps );
		for( int p = 0; p < n; p++ ) out(p) = ( f_eps(p) - f(p) ) / eps;
	}
	// f2, f_eps2, dxx2, dxx_eps2, xx_eps2, vt
	// second derivatives w.r.t. x
	if( (nx == 2) && (np == 0) )
	{
		for( int j = 0; j < n; j++ )
		{
			deri( dxx2, t, xx, par, 1, &vx[0], 0, vp, vv );
			xx_eps2 = xx;
			const double eps2 = abs_eps_x2 + rel_eps_x2*fabs(xx(j,vx[1]));
			xx_eps2(j,vx[1]) +=eps2;
			deri( dxx_eps2, t, xx_eps2, par, 1, &vx[0], 0, vp, vv );
			for( int p = 0; p < n; p++ )
			{
				out(p,j) = 0.0;
				for( int q = 0; q < n; q++ )
				{
					out(p,j) += ( dxx_eps2(p,q) - dxx2(p,q) )*vv(q,vx[0]);
				}
				out(p,j) /= eps2;
			}
		}
	}
	// mixed derivative w.r.t. x and the parameters
	if( (nx == 1) && (np == 1) )
	{
		deri( dxx2, t, xx, par, 1, vx, 0, vp, vv );
		par_eps = par;
		const double eps = abs_eps_p1 + rel_eps_p1*fabs(par(vp[0]));
		par_eps(vp[0]) = par(vp[0]) + eps;
		deri( dxx_eps2, t, xx, par_eps, 1, vx, 0, vp, vv );
		for( int p = 0; p < n; p++ )
		{
			for( int q = 0; q < n; q++ )
			{
				out(p,q) = ( dxx_eps2(p,q) - dxx2(p,q) ) / eps;
			}
		}
	}
}