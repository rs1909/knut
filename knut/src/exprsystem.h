// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 - 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef EXPRSYSTEM_H
#define EXPRSYSTEM_H

#include "expr.h"

#ifndef _WIN32
extern "C"
{
#include <dlfcn.h>
}
#else
#include <windows.h>
#endif

template <class T> class KNArray1D;
template <class T> class KNArray2D;
template <class T> class KNArray3D;
class KNVector;

class KNExprSystem
{
public:

  // sysName : the name of the vector field file
  KNExprSystem (); // this is to facilitate inheritance
  KNExprSystem (const std::string& vfexpr, bool trycompile = true);
  ~KNExprSystem ();

  size_t  ndim() const ;
  size_t  npar() const ;
  size_t  ntau() const ;
  // Vectorized versions
  void   p_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par );
  void   p_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp );
  void   mass(KNArray1D<double>& out);
  void   p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel );
  void   p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
                  size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
  // Setting the starting point
  void   stpar(KNVector& par) const;
  void   stsol(KNArray2D<double>& out, const KNArray1D<double>& time);
  void   parnames(const char *names[]) const;
  // create a C++ file
  void   toCxx (std::string& cxx) const;
  void   toString (std::string& vfexpr);
  // diagnostics
//   void   varprint( const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t pos );
protected:
  // vfexpr : the vector field expression
  // sysName : the file where vfexpr comes from (to track changes)
  // trycompile : whether to attempt to compile
  void constructor (const std::string& vfexpr, const std::string& sysName, bool trycompile = true);
    
private:
  ExpTree::Expression expr;
  std::vector<std::string> varName;
  std::vector<ExpTree::Expression> varDotExpr;
  std::vector<ExpTree::Expression> varInit;
  std::vector<double> varMass;
  std::vector<ExpTree::Expression> delayExpr;
  std::vector<std::string> parName;
  std::vector<double> parInit;
  
  // Delay W.R.T. parameters
  std::vector<ExpTree::Expression> delayExprDeri;
  // RHS W.R.T. parameters
  std::vector<ExpTree::Expression> varDot_p;
  // Jacobian
  std::vector<ExpTree::Expression> varDot_x;
  // mixed
  std::vector<ExpTree::Expression> varDot_x_p;
  // Hessian
  std::vector<ExpTree::Expression> varDot_xx;
  // Hessian (ORIG)
  std::vector<ExpTree::Expression> varDot_hess;

  // stack size
  size_t vectorSize;
  size_t stackSize;
  std::vector<ExpTree::Value> stack;
  
  void resizeStackVector (size_t vectorlen);
  
  std::vector<ExpTree::Value> varArray; 
  std::vector<double> parArray;
  
  void fillTime (const KNArray1D<double>& time);
  void fillVar (const KNArray3D<double>& x);
  void fillVar2 (const KNArray3D<double>& vv);
  void fillPar (const KNVector& par);
  void setupFunctions (const std::string& shobj);
  
  // functions
  std::function<size_t ()> fp_ndim;
  std::function<size_t ()> fp_npar;
  std::function<size_t ()> fp_ntau;
  std::function<void (KNArray2D<double>&, const KNArray1D<double>&, const KNVector&)> fp_p_tau;
  std::function<void (KNArray2D<double>&, const KNArray1D<double>&, const KNVector&, size_t)> fp_p_dtau;
  std::function<void (KNArray1D<double>&)> fp_mass;
  std::function<void (KNArray2D<double>&, const KNArray1D<double>&, const KNArray3D<double>&, const KNVector&, size_t)> fp_p_rhs;
  std::function<void (KNArray3D<double>&, const KNArray1D<double>&, const KNArray3D<double>&, const KNVector&,
                      size_t, size_t, const size_t*, size_t, const size_t*, const KNArray3D<double>&)> fp_p_deri;
  std::function<void (KNVector&)> fp_stpar;
  std::function<void (KNArray2D<double>& , const KNArray1D<double>&)> fp_p_stsol;
  std::function<void (const char *[])> fp_parnames;
  static bool workingCompiler;
  
#ifndef _WIN32
  typedef void*   tdlhand;
#else
  typedef HMODULE tdlhand;
#endif

  tdlhand     tdlopen(const char* fname);
  void*       tdlsym(tdlhand h, const char* sname);
  int         tdlclose(tdlhand h);
  const char* tdlerror();

  tdlhand handle;
  const char* error;
};

class KNSystem : public KNExprSystem
{
public:
  KNSystem ();
  KNSystem (const std::string& sysName, bool trycompile = true);
};

#endif
