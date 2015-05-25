// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#define KNUTSYS_H 

#include "expr.h"
#include "exprsystem.h"
#include "knerror.h"
#include "matrix.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <limits>
#include <iomanip>

// for calling the compiler

#include <cstring>
#include <cstdlib>
#include <cstdio>

extern "C"{
#include <sys/stat.h>
#ifndef _WIN32
#include <sys/wait.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <unistd.h>
#else
#include <windows.h>
#endif 
}

// for the compile time check results
#include "config.h"

using namespace ExpTree;

static void runCompiler(const std::string& cxxstring, const std::string& shobj, const std::string& executableDir);
static void getExecutableDir (std::string& executableDir);
static bool to_compile(const struct stat *sbuf_so, const struct stat *sbuf_src);

KNSystem::KNSystem (const std::string& sysName, bool trycompile)
{
  std::string contents;
  std::ifstream in (sysName);
  if (in)
  {
    in.seekg (0, std::ios::end);
    contents.resize (in.tellg());
    in.seekg (0, std::ios::beg);
    in.read (&contents[0], contents.size());
    in.close ();
  } else
  {
    P_MESSAGE2("Cannot open ", sysName);
  }
  
  std::string vfexpr;
  if (Expression::fromXML (vfexpr, contents)) {}
  else vfexpr = contents;

//   std::cout << "-BEGIN EX " << sysName << "\n";
//   std::cout << vfexpr << "\n";
//   std::cout << "-MID EX\n";
//   std::cout << contents << "\n";
//   std::cout << "-END EX\n";
//   std::cout.flush ();
  try
  {
    KNExprSystem::constructor (vfexpr, sysName, trycompile);
  }
  catch (KNException& ex)
  {
    ex.exprStr (vfexpr);
    throw ex;
  }
}

KNExprSystem::KNExprSystem () : handle(nullptr)
{
}

KNExprSystem::KNExprSystem (const std::string& vfexpr, bool trycompile) : handle(nullptr)
{
  constructor (vfexpr, std::string(), trycompile);
}
  
void KNExprSystem::constructor (const std::string& vfexpr, const std::string& sysName, bool trycompile)
{
  std::string vfName;
//   std::cout << "BEGIN EX\n";
//   std::cout << vfexpr << "\n";
//   std::cout << "MID EX\n";
//   std::cout.flush ();
  expr.fromString (vfexpr);
//   expr.print (std::cout);
//   std::cout << "END EX\n";
//   std::cout.flush ();
//  try
//  {
    expr.knutSplit (vfName, varName, varDotExpr, varInit, varMass, delayExpr, parName, parInit );
//  }
//  catch (KNException& ex)
//  {
//    std::cout << ex.str();
//  }
  
#ifdef DEBUG
  // print everything to see if there's an error
  std::cout << "NAME = "<< vfName << "\n";
  for (std::string& it : varName ) std::cout << "VAR = " << it << "\n";
  for (ExpTree::Expression& it : varDotExpr ) { it.optimize (); std::cout << "DOT = "; it.print (std::cout); std::cout << "\n"; }
  for (ExpTree::Expression& it : varInit ) { it.optimize (); std::cout << "INIT = "; it.print (std::cout); std::cout << "\n"; }
  for (double& it : varMass ) std::cout << "MASS = " << it << "\n";
  for (ExpTree::Expression& it : delayExpr ) { it.optimize ();std::cout << "DELAY = "; it.print (std::cout); std::cout << "\n"; }
  for (std::string& it : parName ) std::cout << "PAR = " << it << "\n";
  for (double& it : parInit ) std::cout << "PARINIT = " << it << "\n";
#endif
  
  // no name was provided
  if (vfName.empty ()) trycompile = false;
  
  // Delay W.R.T. parameters
  delayExprDeri.resize (delayExpr.size()*parName.size());
  for (size_t p = 0; p < parName.size(); p++)
  {
    NodePar pid(p, 0);
    for (size_t q = 0; q < delayExpr.size(); q++)
    {
      delayExpr[q].derivative(delayExprDeri[p + q*parName.size()], &pid);
    }
  }
  
  // W.R.T. parameters
  varDot_p.resize (varDotExpr.size()*parName.size());
  for (size_t q = 0; q < varDotExpr.size(); q++)
  {
    for (size_t p = 0; p < parName.size(); p++)
    {
      NodePar pid (p, 0);
      varDotExpr[q].derivative(varDot_p[p + q*parName.size()], &pid);
    }
  }

  // Jacobian
  varDot_x.resize (varDotExpr.size()*varDotExpr.size()*(delayExpr.size() + 1));
  for (size_t q = 0; q < varDotExpr.size(); q++)
  {
    for (size_t p = 0; p < varDotExpr.size()*(delayExpr.size() + 1); p++)
    {
      NodeVar vid (1 + p, 0);
      varDotExpr[q].derivative(varDot_x[p + q*varDotExpr.size()*(delayExpr.size() + 1)], &vid);
    }
  }
  // mixed
  varDot_x_p.resize (varDotExpr.size()*varDotExpr.size()*(delayExpr.size() + 1)*parName.size());
  for (size_t r = 0; r < parName.size(); r++)
  {
    NodePar pid (r, 0);
    for (size_t p = 0; p < varDotExpr.size()*(delayExpr.size() + 1); p++)
    {
      for (size_t q = 0; q < varDotExpr.size(); q++)
      {
        varDot_x[p + q*varDotExpr.size()*(delayExpr.size() + 1)].
          derivative(varDot_x_p[r + (p + q*varDotExpr.size()*(delayExpr.size() + 1))*parName.size()], &pid);
      }
    }
  }  
  // Hessian
//  varDot_xx.resize (varDotExpr.size()*varDotExpr.size()*(delayExpr.size() + 1));
//  // Matrix - Vector multiplication with the Jacobian & taking the derivative
//  std::vector<NodeVar*> newvar(varDotExpr.size()*(delayExpr.size() + 1));
//  for (size_t p = 0; p < varDotExpr.size()*(delayExpr.size() + 1); p++)
//  {
//    newvar[p] = new NodeVar(1 + p + varDotExpr.size()*(delayExpr.size() + 1));
//  }
//  
//  for (size_t q = 0; q < varDotExpr.size(); q++)
//  {
//    NodeAdd* ta = new NodeAdd(varDotExpr.size()*(delayExpr.size() + 1));
//    for (size_t p = 0; p < varDotExpr.size()*(delayExpr.size() + 1); p++)
//    {
//      NodeTimes* tm = new NodeTimes(2);
//      tm -> addArgument (0, varDot_x[p + q*varDotExpr.size()*(delayExpr.size() + 1)].copy(), 1.0, false);
//      tm -> addArgument (1, newvar[p] -> copy(), 1.0, false);
//      ta -> addArgument (p, tm, 1.0);
//    }
//    Node* opta = ta;
//    ta = 0;
//    opta -> optimize (&opta);
//    for (size_t k = 0; k < varDotExpr.size()*(delayExpr.size() + 1); k++)
//    {
//      NodeVar vid(1 + k);
//      Node* deri = opta -> derivative (&vid);
////      Node* deri = opta -> copy ();
//      deri -> optimize (&deri);
//      varDot_xx[k + q*varDotExpr.size()*(delayExpr.size() + 1)].fromNode (deri);
//    }
//    opta -> deleteTree ();
//    delete opta;
//  }
//  // cleaning up
//  for (size_t p = 0; p < varDotExpr.size()*(delayExpr.size() + 1); p++)
//  {
//    delete newvar[p];
//  }

  // hessian original interpretation (INEFFICIENT)
  varDot_hess.resize (varDotExpr.size()*varDotExpr.size()*(delayExpr.size() + 1)*(delayExpr.size() + 1));
  // k : dim1 = varDotExpr.size()*(delayExpr.size() + 1) deri w.r.t
  // l : dim2 = vx[0] which delay of vv
  // q : dim3 = rhs index
  // Matrix - Vector multiplication with the Jacobian & taking the derivative
  std::vector<NodeVar*> vvsub(varDotExpr.size()*(delayExpr.size() + 1));
  for (size_t p = 0; p < varDotExpr.size()*(delayExpr.size() + 1); p++)
  {
    vvsub[p] = new NodeVar (1 + p + varDotExpr.size()*(delayExpr.size() + 1), 0);
  }
  for (size_t q = 0; q < varDotExpr.size(); q++)
  {
    for (size_t l = 0; l < (delayExpr.size() + 1); l++)
    {
      NodeAdd* ta = new NodeAdd (varDotExpr.size(), 0);
      for (size_t p = 0; p < varDotExpr.size(); p++)
      {
        NodeTimes* tm = new NodeTimes (2, 0);
        tm -> addArgument (0, varDot_x[p + l*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1)].copy(), 1.0, false);
        tm -> addArgument (1, vvsub[p + l*varDotExpr.size()] -> copy(), 1.0, false);
        ta -> addArgument (p, tm, 1.0);
      }
      Node* opta = ta;
      ta = 0;
//      opta -> optimize (&opta);
      for (size_t k = 0; k < varDotExpr.size()*(delayExpr.size() + 1); k++)
      {
        NodeVar vid (1 + k, 0);
        Node* deri = opta -> derivative (&vid);
//        Node* deri = opta -> copy ();
        deri -> optimize (&deri);
        size_t idx = k + l*varDotExpr.size()*(delayExpr.size() + 1) + 
           q*varDotExpr.size()*(delayExpr.size() + 1)*(delayExpr.size() + 1);
        varDot_hess[idx].fromNode (deri);
      }
      opta -> deleteTree ();
      delete opta;
    }
  }
  // cleaning up
  for (size_t p = 0; p < vvsub.size(); p++)
  {
    delete vvsub[p];
  }
  
  stackSize = 0;
  vectorSize = 0;
  for (Expression& it : delayExpr) { it.stackCount (stackSize); }
  for (Expression& it : delayExprDeri) { it.stackCount (stackSize); }
  for (Expression& it : varDotExpr) { it.stackCount (stackSize); }
  for (Expression& it : varDot_p) { it.stackCount (stackSize); }
  for (Expression& it : varDot_x) { it.stackCount (stackSize); }
  for (Expression& it : varDot_x_p) { it.stackCount (stackSize); }
  for (Expression& it : varDot_xx) { it.stackCount (stackSize); }
  for (Expression& it : varDot_hess) { it.stackCount (stackSize); }
  stack.resize (stackSize+1);
  
  varArray.resize (1 + 2*varDotExpr.size()*(delayExpr.size() + 1));
  parArray.resize (parName.size());

  // try to compile
  if (workingCompiler && trycompile)
  {
    std::string shobj(vfName);
    shobj += ".so";
    
    struct stat sbuf_so;
    struct stat sbuf_vf;
    
    bool compile_cxx = true;
    if (!sysName.empty ())
    {
      int res_so = stat(shobj.c_str(), &sbuf_so);
      int res_vf = stat(sysName.c_str(), &sbuf_vf);
      P_ERROR_X3 (res_vf == 0, "Vector field definition '", sysName, "' is missing.");
      if (res_so != 0) compile_cxx = true;
      else compile_cxx = to_compile(&sbuf_so, &sbuf_vf);
    }
    
    if (compile_cxx)
    {
      std::string executableDir;
      getExecutableDir (executableDir);
      std::string cxxstring;
      toCxx (cxxstring);
      try {
        runCompiler (cxxstring, shobj, executableDir);
        setupFunctions (shobj);
      }
      catch (KNException& ex)
      {
        workingCompiler = false;
        std::cerr << "Cannot compile the C++ system definition.\n" << ex.str();
        std::cerr.flush ();
      }
    }
  }
}

KNExprSystem::~KNExprSystem()
{
  for (size_t k = 1; k < stack.size(); k++)
  {
    delete[] stack [k].data;
  }
  if (handle != 0) tdlclose(handle);
}

void KNExprSystem::toString (std::string& vfexpr)
{
  std::ostringstream ostr;
  expr.print (ostr);
  vfexpr = ostr.str ();
}

// private
void KNExprSystem::resizeStackVector (size_t vectorlen)
{
  if (vectorlen > vectorSize)
  {
    for (size_t k = 1; k < stack.size(); k++)
    {
      delete[] stack [k].data;
      stack [k].data = new double[vectorlen];
      stack [k].skip = 1;
    }
    vectorSize = vectorlen;
  }
}

// private
void KNExprSystem::fillTime (const KNArray1D<double>& time)
{
  varArray[0].data = const_cast<KNArray1D<double>&>(time).pointer ();
  varArray[0].skip = 1;
}

// private
void KNExprSystem::fillVar (const KNArray3D<double>& var)
{
  for (size_t p = 0; p < varDotExpr.size(); p++)
  {
    for (size_t q = 0; q < delayExpr.size() + 1; q++)
    {
      varArray[1 + p + q*varDotExpr.size()].data 
        = const_cast<KNArray3D<double>&>(var).pointer (p,q,0);
      varArray[1 + p + q*varDotExpr.size()].skip 
        = (((size_t)const_cast<KNArray3D<double>&>(var).pointer(p,q,1)) - ((size_t)const_cast<KNArray3D<double>&>(var).pointer(p,q,0)))/sizeof(double);
    }
  }
}

// private
void KNExprSystem::fillVar2 (const KNArray3D<double>& vv)
{
  const size_t skip = 1 + varDotExpr.size()*(delayExpr.size() + 1);
  for (size_t p = 0; p < varDotExpr.size(); p++)
  {
    for (size_t q = 0; q < delayExpr.size() + 1; q++)
    {
      varArray[skip + p + q*varDotExpr.size()].data 
        = const_cast<KNArray3D<double>&>(vv).pointer (p,q,0);
      varArray[skip + p + q*varDotExpr.size()].skip 
        = (((size_t)const_cast<KNArray3D<double>&>(vv).pointer(p,q,1)) - ((size_t)const_cast<KNArray3D<double>&>(vv).pointer(p,q,0)))/sizeof(double);
    }
  }
}

// private
void KNExprSystem::fillPar (const KNVector& par)
{
  for (size_t p = 0; p < parName.size(); p++)
  {
    parArray[p] = par(p);
  }
}

size_t KNExprSystem::ndim() const 
{
  if (fp_ndim) return fp_ndim ();
  else return varDotExpr.size();
}

/// WARNING Why do we need to add 3 to it?
size_t KNExprSystem::npar() const 
{
  if (fp_npar) return fp_npar ();
  else return parName.size() - 3;
}

size_t KNExprSystem::ntau() const 
{
  if (fp_ntau) return fp_ntau ();
  else return delayExpr.size() + 1;
}

// delay values
void KNExprSystem::p_tau ( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par )
{
  if (fp_p_tau) fp_p_tau (out, time, par);
  else
  {
    resizeStackVector (time.size ());
    stack[0].skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
    fillTime (time);
    fillPar (par);
    
    for (size_t k = 0; k < time.size (); k++) out (0,k) = 0.0;
    for (size_t k = 0; k < delayExpr.size(); k++)
    {
      stack[0].data = out.pointer(1 + k,0);
      delayExpr[k].evaluate (stack, varArray, parArray, time.size ());
    }
  }
}

void KNExprSystem::p_dtau ( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp )
{
  if (fp_p_dtau) fp_p_dtau (out, time, par, vp);
  else
  {
    P_ERROR(vp < parName.size() + 1);
    resizeStackVector (time.size ());
    stack[0].skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
    fillPar (par);
    fillTime (time);
    for (size_t k = 0; k < time.size (); k++) out (0,k) = 0.0;
    for (size_t k = 0; k < delayExpr.size(); k++)
    {
      stack[0].data = out.pointer(1 + k,0);
      delayExprDeri[vp + k*parName.size()].evaluate (stack, varArray, parArray, time.size ());
    }
  }
}

// no mass is implemented yet
void KNExprSystem::mass(KNArray1D<double>& out)
{
  for (size_t k = 0; k < varDotExpr.size(); k++)
  {
    out (k) = varMass[k];
  }
}

void KNExprSystem::p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel )
{
  if (fp_p_rhs) fp_p_rhs (out, time, x, par, sel);
  else
  {
    resizeStackVector (time.size ());
    stack[0].skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
    fillTime (time);
    fillPar (par);
    fillVar (x);

    for (size_t k = 0; k < varDotExpr.size(); k++)
    {
      stack[0].data = out.pointer (k,0);
      varDotExpr[k].evaluate (stack, varArray, parArray, time.size ());
    }
  }
}

void KNExprSystem::p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
  size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv )
{
  if (fp_p_deri) fp_p_deri (out, time, x, par, sel, nx, vx, np, vp, vv);
  else
  {
    resizeStackVector (time.size ());
    stack[0].skip = (((size_t)out.pointer(0,0,1)) - ((size_t)out.pointer(0,0,0)))/sizeof(double);
    fillTime (time);
    fillPar (par);
    fillVar (x);
    P_ERROR( par.size() >= parInit.size () );
    P_ERROR( x.dim3() == time.size () );
    P_ERROR( x.dim2() >= (delayExpr.size() + 1) );
    P_ERROR( x.dim1() == varDotExpr.size() );
    
    // parameters
    if (nx == 0 && np == 1)
    {
      P_ERROR( out.dim3() == time.size () );
      P_ERROR( out.dim2() == 1 );
      P_ERROR( out.dim1() == varDotExpr.size() );
      for (size_t k = 0; k < varDotExpr.size(); k++)
      {
        stack[0].data = out.pointer (k,0,0);
        stack[0].skip = (((size_t)out.pointer(k,0,1)) - ((size_t)out.pointer(k,0,0)))/sizeof(double);
        varDot_p[vp[0] + k*parName.size()].evaluate (stack, varArray, parArray, time.size ());
      }    
    }
    // jacobian
    if (nx == 1 && np == 0)
    {
      P_ERROR( out.dim3() == time.size () );
      P_ERROR( out.dim2() == varDotExpr.size() );
      P_ERROR( out.dim1() == varDotExpr.size() );
      for (size_t p = 0; p < varDotExpr.size(); p++)
      {
        for (size_t q = 0; q < varDotExpr.size(); q++)
        {
          stack[0].data = out.pointer (q,p,0);
          stack[0].skip = (((size_t)out.pointer(q,p,1)) - ((size_t)out.pointer(q,p,0)))/sizeof(double);
          varDot_x[p + vx[0]*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1)].evaluate (stack, varArray, parArray, time.size ());
        }
      }
    }
    // mixed derivatives
    if (nx == 1 && np == 1)
    {
      for (size_t p = 0; p < varDotExpr.size(); p++)
      {
        for (size_t q = 0; q < varDotExpr.size(); q++)
        {
          stack[0].data = out.pointer (q,p,0);
          stack[0].skip = (((size_t)out.pointer(q,p,1)) - ((size_t)out.pointer(q,p,0)))/sizeof(double);
          varDot_x_p[vp[0] + (p + vx[0]*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1))*parName.size()].evaluate (stack, varArray, parArray, time.size ());
        }
      }
    }
    // hessian NEW
  //   if (nx == 2 && np == 0)
  //   {
  //     if (vx[0] != 0) { out.clear (); return; }
  // //    std::cout << "waring: HESSIAN called\n"; // currently incompatible ...
  //     fillVar2 (vv);
  //     for (size_t q = 0; q < varDotExpr.size(); q++)
  //     {
  //       for (size_t p = 0; p < varDotExpr.size(); p++)
  //       {
  //         stack[0].data = out.pointer (q,p,0);
  //         varDot_xx[p + vx[1]*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1)].evaluate (stack, varArray, parArray, time.size ());
  // //        std::cout << "d_xx("<< q << "," << p << ") = "; std::cout.flush();
  // //        varDot_xx[p + vx[1]*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1)].print(std::cout);
  // //        std::cout << " = " << stack[0].data[0] << "\n";
  //       }
  //     }
  //   }
  // hessian OLD TYPE
    if (nx == 2 && np == 0)
    {
      fillVar2 (vv);
      P_ERROR( out.dim3() == time.size () );
      P_ERROR( out.dim2() == varDotExpr.size() );
      P_ERROR( out.dim1() == varDotExpr.size() );
      P_ERROR( vv.dim3() == time.size () );
      P_ERROR( vv.dim2() >= (delayExpr.size() + 1) );
      P_ERROR( vv.dim1() == varDotExpr.size() );

      out.clear ();
      for (size_t q = 0; q < varDotExpr.size(); q++)
      {
        for (size_t p = 0; p < varDotExpr.size(); p++)
        {
          stack[0].data = out.pointer (q,p,0);
          stack[0].skip = (((size_t)out.pointer(q,p,1)) - ((size_t)out.pointer(q,p,0)))/sizeof(double);
          const size_t idx = p + vx[1]*varDotExpr.size() + vx[0]*varDotExpr.size()*(delayExpr.size() + 1) + 
            q*varDotExpr.size()*(delayExpr.size() + 1)*(delayExpr.size() + 1);
          varDot_hess[idx].evaluate (stack, varArray, parArray, time.size ());
        }
      }
    }
  }
}

// for previous debug purposes
// void KNExprSystem::varprint( const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t pos )
// {
//   resizeStackVector (time.size ());
//   fillTime (time);
//   fillPar (par);
//   fillVar (x);
//   std::vector<double> locst (stack.size());
//   std::vector<double> locvar (varArray.size());
//   for (size_t k = 0; k < 1 + varDotExpr.size()*(delayExpr.size() + 1); k++) locvar[k] = varArray[k].data[pos];
// 
//   std::cout << "{";
//   for (size_t q = 0; q < parArray.size(); q++)
//   {
//     if (q != 0) std::cout << ",";
//     std::cout << "P[" << q << "]->" << parArray[q];
//   }
//   std::cout << "V[0]->" << time (pos) << ",";
//   for (size_t q = 0; q < (delayExpr.size() + 1); q++)
//   {
//     for (size_t p = 0; p < varDotExpr.size(); p++)
//     {
//       if (q != 0 || p != 0) std::cout << ",";
//       std::cout << "V[" << 1 + p + q*varDotExpr.size() << "]->" << x(p,q,pos);
//     }
//   }
//   std::cout << "}\n";
// }

// Setting the starting point
void KNExprSystem::stpar(KNVector& par) const
{
  if (fp_stpar) fp_stpar (par);
  else
  {
    for (size_t p = 0; p < parInit.size()-3; p++)
    {
      par(p) = parInit[p];
    }
  }
}

void KNExprSystem::stsol (KNArray2D<double>& out, const KNArray1D<double>& time)
{
  if (fp_p_stsol) fp_p_stsol (out, time);
  else
  {
    resizeStackVector (time.size ());
    stack[0].skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
    fillTime (time);
    for (size_t k = 0; k < varInit.size(); k++)
    {
      stack[0].data = out.pointer(k,0);
      varInit[k].evaluate (stack, varArray, parArray, time.size ());
    }
  }
}

void KNExprSystem::parnames (const char *names[]) const
{
  if (fp_parnames) fp_parnames (names);
  else
  {
    for (size_t p = 0; p < parName.size() - 3; p++)
    {
      names[p] = parName[p].c_str ();
    }
  }
}

static Node* indicesToCxx(const Expression& expr, 
                         const std::vector<NodePar*>& parlist,
                         const std::vector<NodeSymbol*>& parsymlist,
                         const std::vector<NodeVar*>& varlist,
                         const std::vector<NodeSymbol*>& varsymlist)
{
  Node* nd = expr.copy ();
  for (size_t k = 0; k < parlist.size(); k++)
  {
    nd -> replaceSymbol (*parlist[k], *parsymlist[k], &nd);
  }
  for (size_t k = 0; k < varlist.size(); k++)
  {
    nd -> replaceSymbol (*varlist[k], *varsymlist[k], &nd);
  }
  return nd;
}

void KNExprSystem::toCxx (std::string& cxx) const
{
  std::ostringstream out;
  std::vector<NodePar*> parlist (parName.size());
  std::vector<NodeSymbol*> parsymlist (parName.size());
  std::vector<NodeVar*> varlist (1+2*(delayExpr.size() + 1)*varDotExpr.size());
  std::vector<NodeSymbol*> varsymlist (1+2*(delayExpr.size() + 1)*varDotExpr.size());

  out.precision (std::numeric_limits<double>::digits10 + 1);
  out.setf (std::ios::scientific);
  
  std::ostringstream symtext;
  for (size_t q = 0; q < parName.size(); q++)
  {
    parlist[q] = new NodePar (q, 0);
    symtext << "par(" << q << ")";
    parsymlist[q] = new NodeSymbol (symtext.str (), 0);
    symtext.str (std::string());
  }

  varlist[0] = new NodeVar (0, 0);
  varsymlist[0] = new NodeSymbol ("time(idx)", 0);
  for (size_t q = 0; q < (delayExpr.size() + 1); q++)
  {
    for (size_t p = 0; p < varDotExpr.size(); p++)
    {
      varlist[1 + p + q*varDotExpr.size()] = new NodeVar (1 + p + q*varDotExpr.size(), 0);
      symtext << "x(" << p << "," << q << "," << "idx)";
      varsymlist[1 + p + q*varDotExpr.size()] = new NodeSymbol (symtext.str (), 0);
      symtext.str (std::string());
    }
  }
  const size_t skip = 1+(delayExpr.size() + 1)*varDotExpr.size();
  for (size_t q = 0; q < (delayExpr.size() + 1); q++)
  {
    for (size_t p = 0; p < varDotExpr.size(); p++)
    {
      varlist[skip + p + q*varDotExpr.size()] = new NodeVar (skip + p + q*varDotExpr.size(), 0);
      symtext << "v(" << p << "," << q << "," << "idx)";
      varsymlist[skip + p + q*varDotExpr.size()] = new NodeSymbol (symtext.str (), 0);
      symtext.str (std::string());
    }
  }
  // -- printing starts
  const size_t NPAR = parName.size()-3;
  out << "#include <cmath>\n"
         "#include \"knutsys.h\"\n"
         "\n"
         "static inline double heaviside(double x)\n"
         "{\n"
         "  if (x > 0) return 1.0;\n"
         "  else return 0.0;\n"
         "}\n"
         "\n"
         "extern \"C\" {\n"
         "\n"
         "size_t sys_ndim()  { return " << varDotExpr.size() << "; }  // KNExprSystem dimension\n"
         "size_t sys_npar()  { return " << NPAR << "; }  // Number of parameters, plus one (for the period)\n"
         "size_t sys_ntau()  { return " << delayExpr.size() + 1 << "; }  // Number of delays, plus one\n"
         "size_t sys_nderi() { return 2; }  // Order of derivatives computed here\n"
         "\n";
  out << "void sys_p_tau(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par)\n"
         "{\n"
         "  for (size_t idx = 0; idx < time.size(); idx++)\n"
         "  {\n"
         "    out(0,idx) = 0.0;\n";
  for (size_t q = 0; q < delayExpr.size(); q++)
  {
    Node* nd = indicesToCxx(delayExpr[q], parlist, parsymlist, varlist, varsymlist);
    out << "    out(" << 1+q << ",idx) = ";
    nd -> print (out);
    out << ";\n";
    nd -> deleteTree();
    delete nd;    
  }
  out << "  }\n"
         "}\n";
  out << "void sys_p_dtau(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par, size_t k)\n"
         "{\n";
  for (size_t p = 0; p < NPAR; p++)
  {
    out << "  if (k == " << p << ")\n"
           "  {\n"
           "    for (size_t idx = 0; idx < time.size(); idx++)\n"
           "    {\n"
           "      out(0,idx) = 0.0;\n";
    for (size_t q = 0; q < delayExpr.size(); q++)
    {
      Node* nd = indicesToCxx(delayExprDeri[p + q*parName.size()], parlist, parsymlist, varlist, varsymlist);
      out << "      out(" << 1+q << ",idx) = ";
      nd -> print (out);
      out << ";\n";
      nd -> deleteTree();
      delete nd;    
    }
    out << "    }\n"
           "  }\n";
  }
  out << "}\n";
  out << "void sys_p_rhs(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel)\n"
         "{\n"
         "  // Compute the vector field\n"
         "  for (size_t idx = 0; idx < time.size(); idx++)\n"
         "  {\n";
  for (size_t k = 0; k < varDotExpr.size(); k++)
  {
    Node* nd = indicesToCxx(varDotExpr[k], parlist, parsymlist, varlist, varsymlist);
    out << "    out(" << k << ",idx) = ";
    nd -> print (out);
    out << ";\n";
    nd -> deleteTree();
    delete nd;
  }
  out << "  }\n"
         "}\n";
  // --- derivatives
  out << "void sys_p_deri(KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par,\n" 
         "                size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& v)\n"
         "{\n"
         "  if (nx == 1 && np == 0)\n"
         "  {\n";
  for (size_t k = 0; k < delayExpr.size()+1; k++)
  {
    out << "    if (vx[0] == " << k << ")\n"
           "    {\n"
           "      for (size_t idx = 0; idx < time.size(); idx++)\n"
           "      {\n";
    for (size_t q = 0; q < varDotExpr.size(); q++)
    {
      for (size_t p = 0; p < varDotExpr.size(); p++)
      {
        Node* nd = indicesToCxx(varDot_x[p + k*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1)], parlist, parsymlist, varlist, varsymlist);
        out << "        out(" << q << "," << p << ",idx) = ";
        nd -> print (out);
        out << ";\n";
        nd -> deleteTree();
        delete nd;
      }
    }
    out << "      }\n"
           "    }\n";   
  }
  out << "  }\n"
         "  else if (nx == 0 && np == 1)\n"
         "  {\n";
  for (size_t p = 0; p < NPAR; p++)
  {
    out << "    if (vp[0] == " << p << ")\n"
           "    {\n"
           "      for (size_t idx = 0; idx < time.size(); idx++)\n"
           "      {\n";
    for (size_t k = 0; k < varDotExpr.size(); k++)
    {   
      Node* nd = indicesToCxx(varDot_p[p + k*parName.size()], parlist, parsymlist, varlist, varsymlist);
      out << "        out(" << k << "," << 0 << ",idx) = ";
      nd -> print (out);
      out << ";\n";
      nd -> deleteTree();
      delete nd;
    }
    out << "      }\n"
           "    }\n";
  }
  out << "  }\n"
         "  else if (nx == 1 && np == 1)\n"
         "  {\n";
  for (size_t p = 0; p < NPAR; p++)
  {
    out << "    if (vp[0] == " << p << ")\n"
           "    {\n";
    for (size_t k = 0; k < delayExpr.size()+1; k++)
    {
      out << "      if (vx[0] == " << k << ")\n"
             "      {\n"
             "        for (size_t idx = 0; idx < time.size(); idx++)\n"
             "        {\n";
      for (size_t q = 0; q < varDotExpr.size(); q++)
      {
        for (size_t r = 0; r < varDotExpr.size(); r++)
        {
          Node* nd = indicesToCxx(varDot_x_p[p + (r + k*varDotExpr.size() + q*varDotExpr.size()*(delayExpr.size() + 1))*parName.size()], 
                                  parlist, parsymlist, varlist, varsymlist);
          out << "          out(" << q << "," << r << ",idx) = ";
          nd -> print (out);
          out << ";\n";
          nd -> deleteTree();
          delete nd;
        }
      }
      out << "        }\n"
             "      }\n";   
    }
    out << "    }\n";
  }
   
  out << "  }\n"
         "  else if (nx == 2 && np == 0)\n"
         "  {\n";
  for (size_t r = 0; r < delayExpr.size()+1; r++)
  {
    out << "    if (vx[0] == " << r << ")\n"
           "    {\n";
    for (size_t k = 0; k < delayExpr.size()+1; k++)
    {
      out << "      if (vx[1] == " << k << ")\n"
             "      {\n"
             "        for (size_t idx = 0; idx < time.size(); idx++)\n"
             "        {\n";
      for (size_t q = 0; q < varDotExpr.size(); q++)
      {
        for (size_t p = 0; p < varDotExpr.size(); p++)
        {
          const size_t idx = p + k*varDotExpr.size() + r*varDotExpr.size()*(delayExpr.size() + 1) + 
            q*varDotExpr.size()*(delayExpr.size() + 1)*(delayExpr.size() + 1);
          Node* nd = indicesToCxx(varDot_hess[idx], parlist, parsymlist, varlist, varsymlist);
          out << "          out(" << q << "," << p << ",idx) = ";
          nd -> print (out);
          out << ";\n";
          nd -> deleteTree();
          delete nd;
        }
      }
      out << "        }\n"
             "      }\n";

    }
    out << "    }\n";
  }
  out << "  }\n"
         "}\n";
  // --- starting values
  out << "void sys_stpar(KNVector& par)\n"
         "{\n";
  for (size_t p = 0; p < NPAR; p++)
  {
    out << "  par(" << p << ") = " << parInit[p] << ";\n";
  } 
  out << "}\n"
         "\n"
         "void sys_p_stsol(KNArray2D<double>& out, const KNArray1D<double>& time)\n"
         "{\n"
         "  for (size_t idx = 0; idx < time.size(); idx++)\n"
         "  {\n";
  for (size_t k = 0; k < varDotExpr.size(); k++)
  {
    Node* nd = indicesToCxx(varInit[k], parlist, parsymlist, varlist, varsymlist);
    out << "    out(" << k << ",idx) = ";
    nd -> print (out);
    out << ";\n";
    nd -> deleteTree();
    delete nd;
  }
  out << "  }\n"
         "}\n"
         "\n"
         "void sys_parnames( const char *out[] )\n"
         "{\n";
  for (size_t p = 0; p < NPAR; p++)
  {
    out << "  out[" << p << "] = \"" << parName[p] << "\";\n";
  } 
  out << "}\n"
         "\n"
         "}  // extern \"C\"\n";
  
  
  cxx = out.str ();
  for (NodePar* nd: parlist) delete nd;
  for (NodeVar* nd: varlist) delete nd;
  for (NodeSymbol* nd: parsymlist) delete nd;
  for (NodeSymbol* nd: varsymlist) delete nd;  
}

#ifdef __APPLE__
#include <CoreFoundation/CFBase.h>
#include <CoreFoundation/CFString.h>
#include <CoreFoundation/CFBundle.h>
#endif

#ifdef _WIN32
#define DIRSEP '\\'
#else
#define DIRSEP '/'
#endif

bool KNExprSystem::workingCompiler = true;

static void getExecutableDir (std::string& executableDir)
{
  // Finding the executable directory
  std::string executableFile;
#ifdef __APPLE__
  CFURLRef bundleURL = CFBundleCopyExecutableURL(CFBundleGetMainBundle());
  if(bundleURL)
  {
    CFStringRef cfPath = CFURLCopyFileSystemPath(bundleURL, kCFURLPOSIXPathStyle);
    if(cfPath) 
    {
      CFIndex ssize = CFStringGetLength(cfPath)+1;
      char *buf = new char[ssize];
      CFStringGetCString(cfPath, buf, ssize, kCFStringEncodingMacRoman);
      executableFile = buf;
      delete[] buf; buf = 0;
      CFRelease(cfPath);
    }
    CFRelease(bundleURL);
  }
#elif __linux__
  std::string workingDir;
  std::ostringstream procf;
  procf << "/proc/" << getpid() << "/exe";
  char *buf = new char[512];
  const ssize_t bsfn = readlink(procf.str().c_str(), buf, 511);
  if ( bsfn != -1) { buf[bsfn] = '\0'; executableFile = buf; }
  delete[] buf; buf = 0;
#elif _WIN32
  char *buf = new char[MAX_PATH]; //always use MAX_PATH for filepaths
  GetModuleFileName(NULL, buf, MAX_PATH*sizeof(char));
  executableFile = buf;
  delete[] buf; buf = 0;
#endif
  // remove the filename from the end
  std::string::size_type sidx = executableFile.find_last_of(DIRSEP);
  if (sidx != std::string::npos) executableFile.erase(sidx,std::string::npos);
  executableDir = executableFile;
}

// Separates the string into arguments and add to the list of arguments.
// It does not work if the argument has a space inside.
static void addArgList(std::list<std::string>& argv, const std::string& str)
{
  size_t start = 0;
  size_t i=0;
  while (i < str.size())
  {    
    while (((str[i] == ' ')||(str[i] == '\t'))&&(i < str.size())) ++i;
    start = i;
    while (((str[i] != ' ')&&(str[i] != '\t'))&&(i < str.size())) ++i;
    argv.push_back(str.substr(start, i-start)); 
  }
}

// Constructs the command line and the argument list of the GNU C++ compiler (g++).
static inline void mkArgListCommandLine(std::list<std::string>& arglist, std::string& cmdline, const std::string& shobj, const std::string& executableDir)
{
#ifdef _WIN32
  std::string cxxcomp("g++");
  arglist.push_back(cxxcomp);
#else
  std::string cxxcomp(CMAKE_CXX_COMPILER);
  arglist.push_back(cxxcomp.substr(cxxcomp.find_last_of(DIRSEP)+1,std::string::npos));
#endif
  // constructing the command line
  addArgList(arglist, std::string( 
    CMAKE_CXX_FLAGS " "
    CMAKE_SHARED_LIBRARY_C_FLAGS " " 
    CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS " "
    KNUT_NO_EXCEPTIONS));
  std::string includeArg("-I");
  includeArg += executableDir.substr(0,executableDir.find_last_of(DIRSEP));
  includeArg += DIRSEP;
  includeArg += KNUT_INCLUDE_DIR;
  arglist.push_back(includeArg);
  arglist.push_back("-pipe");
  arglist.push_back("-x");
  arglist.push_back("c++");
  arglist.push_back("-o");
  arglist.push_back(shobj);
  arglist.push_back("-");

  cmdline.erase();
  for(std::list<std::string>::const_iterator it=arglist.begin(); it != arglist.end(); ++it)
  {
    if (it != arglist.begin()) cmdline += ' ';
#ifdef _WIN32
    if (it != arglist.begin()) cmdline += "\"";
    cmdline += *it;
    if (it != arglist.begin()) cmdline += "\"";
#else
    cmdline += *it;
#endif
  }
}

static bool to_compile(const struct stat *sbuf_so, const struct stat *sbuf_src)
{
#ifdef __APPLE__
  if (sbuf_so->st_mtimespec.tv_sec < sbuf_src->st_mtimespec.tv_sec) return true;
  else if (sbuf_so->st_mtimespec.tv_sec == sbuf_src->st_mtimespec.tv_sec) return (sbuf_so->st_mtimespec.tv_nsec <= sbuf_src->st_mtimespec.tv_nsec);
  else return false;
#else
  return (sbuf_so->st_mtime <= sbuf_src->st_mtime);
#endif
}

#ifdef _WIN32

// Run the compiler with cxxstring as the program to compile.
// Catch the standard input, output and error streams to interact.
// This is the Windows version.
static void runCompiler(const std::string& cxxstring, const std::string& shobj, const std::string& executableDir)
{
  // setting up the command line
  std::list<std::string> arglist;
  std::string cmdline;
  mkArgListCommandLine(arglist, cmdline, shobj, executableDir);
  
  HANDLE g_hChildStd_IN_Rd = NULL;
  HANDLE g_hChildStd_IN_Wr = NULL;
  HANDLE g_hChildStd_OUT_Rd = NULL;
  HANDLE g_hChildStd_OUT_Wr = NULL;
  HANDLE g_hChildStd_ERR_Rd = NULL;
  HANDLE g_hChildStd_ERR_Wr = NULL;

  HANDLE g_hInputFile = NULL;

  SECURITY_ATTRIBUTES saAttr;
 
  // Set the bInheritHandle flag so pipe handles are inherited. 
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
  saAttr.bInheritHandle = TRUE;
  saAttr.lpSecurityDescriptor = NULL;

  // Create a pipe for the child process's STDOUT. 
  if (!CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &saAttr, 0)) 
    P_MESSAGE2("StdoutRd CreatePipe", (int)GetLastError());
  // Ensure the read handle to the pipe for STDOUT is not inherited.
  if (!SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
    P_MESSAGE2("Stdout SetHandleInformation", (int)GetLastError());
  // Create a pipe for the child process's STDIN. 
  if (!CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &saAttr, 0))
    P_MESSAGE2("Stdin CreatePipe", (int)GetLastError());
  // Ensure the write handle to the pipe for STDIN is not inherited.
  if (!SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0))
    P_MESSAGE2("Stdin SetHandleInformation", (int)GetLastError());
  // Create a pipe for the child process's STDERR.
  if (!CreatePipe(&g_hChildStd_ERR_Rd, &g_hChildStd_ERR_Wr, &saAttr, 0)) 
    P_MESSAGE2("Stdin CreatePipe", (int)GetLastError());
  // Ensure the read handle to the pipe for STDERR is not inherited.
  if (!SetHandleInformation(g_hChildStd_ERR_Rd, HANDLE_FLAG_INHERIT, 0))
    P_MESSAGE2("Stdin SetHandleInformation", (int)GetLastError());
 
   // Create the child process. 
  CHAR *szCmdline = new CHAR[cmdline.size()+1];
  strcpy(szCmdline, cmdline.c_str());
  
  PROCESS_INFORMATION piProcInfo; 
  STARTUPINFO siStartInfo;
  BOOL bSuccess = FALSE;
 
  // Set up members of the PROCESS_INFORMATION structure. 
  ZeroMemory( &piProcInfo, sizeof(PROCESS_INFORMATION) );
 
  // Set up members of the STARTUPINFO structure. 
  // This structure specifies the STDIN and STDOUT handles for redirection.
  ZeroMemory( &siStartInfo, sizeof(STARTUPINFO) );
  siStartInfo.cb = sizeof(STARTUPINFO); 
  siStartInfo.hStdError = g_hChildStd_ERR_Wr;
  siStartInfo.hStdOutput = g_hChildStd_OUT_Wr;
  siStartInfo.hStdInput = g_hChildStd_IN_Rd;
  siStartInfo.dwFlags |= STARTF_USESTDHANDLES;
 
  // Create the child process. 
  bSuccess = CreateProcess(NULL, 
    szCmdline,     // command line 
    NULL,          // process security attributes 
    NULL,          // primary thread security attributes 
    TRUE,          // handles are inherited 
    0,             // creation flags 
    NULL,          // use parent's environment 
    NULL,          // use parent's current directory 
    &siStartInfo,  // STARTUPINFO pointer 
    &piProcInfo);  // receives PROCESS_INFORMATION 
   
  // If an error occurs, exit the application. 
  if ( ! bSuccess )
    P_MESSAGE2("CreateProcess", (int)GetLastError());
 
  // Write to the pipe that is the standard input for a child process. 
  // Data is written to the pipe's buffers, so it is not necessary to wait
  // until the child process is running before writing data.
  DWORD dwRead = cxxstring.size(), dwWritten;
  bSuccess = WriteFile(g_hChildStd_IN_Wr, cxxstring.c_str(), dwRead, &dwWritten, NULL);
  if ( ! bSuccess ) 
    P_MESSAGE2("WriteFile", (int)GetLastError());

  if (!CloseHandle(g_hChildStd_IN_Wr)) 
    P_MESSAGE2("StdInWr CloseHandle", (int)GetLastError());
  if (!CloseHandle(g_hChildStd_IN_Rd)) 
    P_MESSAGE2("StdInRd CloseHandle", (int)GetLastError());
 
  // Read from pipe that is the standard output for child process.
  const size_t bufsize = 2048;
  CHAR *out_buf = new char[bufsize];
  CHAR *err_buf = new char[bufsize];
  dwRead = 0;
  dwWritten = 0;
  bSuccess = FALSE;
  HANDLE hParentStdOut = GetStdHandle(STD_OUTPUT_HANDLE);

  // Close the write end of the pipe before reading from the 
  // read end of the pipe, to control child process execution.
  // The pipe is assumed to have enough buffer space to hold the
  // data the child process has already written to it.
  if (!CloseHandle(g_hChildStd_OUT_Wr)) 
    P_MESSAGE2("StdOutWr CloseHandle", (int)GetLastError());
   
  std::string outstr;
  std::string errstr;
  for (;;) 
  {
    bSuccess = ReadFile( g_hChildStd_OUT_Rd, out_buf, bufsize-1, &dwRead, NULL);
    err_buf[dwRead] = '\0';
    if( ! bSuccess || dwRead == 0 ) break;
    outstr.append(out_buf);
  }
  
  if (!CloseHandle(g_hChildStd_OUT_Rd)) 
    P_MESSAGE2("StdOutRd CloseHandle", (int)GetLastError());
  
  if (!CloseHandle(g_hChildStd_ERR_Wr)) 
    P_MESSAGE2("StdErrWr CloseHandle", (int)GetLastError());
   
  for (;;) 
  {
    bSuccess = ReadFile( g_hChildStd_ERR_Rd, err_buf, bufsize-1, &dwRead, NULL);
    err_buf[dwRead] = '\0';
    if( ! bSuccess || dwRead == 0 ) break;
    errstr.append(err_buf);
  }   
  
  if (!CloseHandle(g_hChildStd_ERR_Rd)) 
    P_MESSAGE2("StdErrRd CloseHandle", (int)GetLastError());
  
  WaitForSingleObject(piProcInfo.hProcess, INFINITE);
  DWORD exitStatus;
  GetExitCodeProcess(piProcInfo.hProcess, &exitStatus); 
  CloseHandle(piProcInfo.hProcess);
  CloseHandle(piProcInfo.hThread);

  if (exitStatus != 0) P_MESSAGE6("The error output of the compile command '",
      cmdline, "' is ", errstr.c_str(), " and the standard output is ", outstr.c_str());
  
  delete[] out_buf; out_buf = 0;
  delete[] err_buf; err_buf = 0;
  delete[] szCmdline; szCmdline = 0;
}

#else

// Opens the three pipes for interacting with a program.
// arglist is the argument list starting with the program name.
// the int* - s are the file numbers for the three pipes.
static pid_t pipeOpen(std::list<std::string>& arglist, int* input, int* output, int* error)
{
  if ( (input==0)&&(output==0)&&(error==0) ) return -1;
  char *argv[arglist.size()+1];
  int i = 0;
  for (std::list<std::string>::const_iterator it=arglist.begin(); it != arglist.end(); ++it)
  {
    argv[i] = new char[it->size()+1];
    strcpy(argv[i],it->c_str());
    ++i;
  }
  argv[i] = 0;

  int fds_output[2], fds_input[2], fds_error[2];
  int fdc_output[2], fdc_input[2], fdc_error[2];

  /* Create a pipe.  File descriptors for the two ends of the pipe are
     placed in fds.  */
  if (output) if (pipe (fds_output) == -1) P_MESSAGE2("pipe: ", strerror(errno));
  if (input)  if (pipe (fds_input) == -1) P_MESSAGE2("pipe: ", strerror(errno));
  if (error)  if (pipe (fds_error) == -1) P_MESSAGE2("pipe: ", strerror(errno));
  for (int i=0; i<2; ++i)
  {
    fdc_output[i] = fds_output[i]; 
    fdc_input[i] = fds_input[i];
    fdc_error[i] = fds_error[i];
  }
  /* Fork a child process.  */
  pid_t pid = fork ();
  if (pid == (pid_t) 0)
  {
//    std::cout << "Closing ends\n";
    /* This is the child process.  Close our copy of the read (write) end of
       the file descriptor.  */
    if (input)  if (close (fds_input[1]) == -1) P_MESSAGE2("close: ", strerror(errno));
    if (output) if (close (fds_output[0]) == -1) P_MESSAGE2("close: ", strerror(errno));
    if (error)  if (close (fds_error[0]) == -1) P_MESSAGE2("close ", strerror(errno));
//    std::cout << "Closed ends\n";
    /* Connect the read(write) end of the pipe to standard input.  */
    if (input)  if (dup2 (fds_input[0], STDIN_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
//    std::cout << "Dup input\n";
    if (output) if (dup2 (fds_output[1], STDOUT_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
//    std::cout << "Dup output\n";
    if (error)  if (dup2 (fds_error[1], STDERR_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
    /* Replace the child process with the "cat" program.  */
//    std::cout << "Starting process\n"; std::cout.flush();
    int st = execvp (argv[0], argv);
    // This wouldn't return unless there's a problem
    if (st == -1) P_MESSAGE4("Error executing command ", argv[0], ": ", strerror(errno));
//    std::cout << "Somehow returned\n"; std::cout.flush();
  }
  else
  {
    if (pid == -1) P_MESSAGE2("Failed to fork process: ", strerror(errno));
    /* Close our copy of the write (read) end of the file descriptor.  */
    if (input)
    {
      if (close (fdc_input[0]) == -1) P_MESSAGE2("close ", strerror(errno));
      int iflags = fcntl(fdc_input[1], F_GETFL, 0);
      iflags |= O_NONBLOCK;
      fcntl(fdc_input[1], F_SETFL, iflags);
      *input = fdc_input[1];
    }       
    if (output)
    {
      if (close (fdc_output[1]) == -1) P_MESSAGE2("close ", strerror(errno));
      int oflags = fcntl(fdc_output[0], F_GETFL, 0);
      oflags |= O_NONBLOCK;
      fcntl(fdc_output[0], F_SETFL, oflags);
      *output = fdc_output[0];
    }
    if (error)
    {
      if (close (fdc_error[1]) == -1) P_MESSAGE2("close ", strerror(errno));
      int eflags = fcntl(fdc_error[0], F_GETFL, 0);
      eflags |= O_NONBLOCK;
      fcntl(fdc_error[0], F_SETFL, eflags);
      *error = fdc_error[0];
    }
//    std::cout << "Somehow returned - other!!\n"; std::cout.flush();
  }
  for (size_t k=0; k<arglist.size(); ++k) { delete[] argv[k]; argv[k] = 0; }
  return pid;
}

// Run the compiler with cxxstring as the program to compile.
// Catch the standard input, output and error streams to interact.
// This is the Linux/Mac OS X version.
static void runCompiler(const std::string& cxxstring, const std::string& shobj, const std::string& executableDir)
{
  std::list<std::string> arglist;
  std::string cmdline;
  mkArgListCommandLine(arglist, cmdline, shobj, executableDir);
  
  //   std::cout << cmdline << std::endl;
  // running the command
  int input, output, error;
  // pid_t pid = pipeOpen ( ... could wait for this pid specifically
  pipeOpen(arglist, &input, &output, &error);
  // we need to write and read until all of them are finished
  bool outfin = false, infin = false, errfin = false;
  // input
  const char* cxxbuf = cxxstring.c_str();
  size_t cxxlen = cxxstring.size(), wbytes = 0;
  // output and err
  pollfd fds[3] = {{input, POLLOUT|POLLHUP, 0},{output, POLLIN|POLLHUP, 0},{error, POLLIN|POLLHUP, 0}};
  const size_t bufsize = 2048;
  char *out_buf = new char[bufsize];
  std::string out_str, err_str;
  int rct = 0, ect = 0;
  size_t iin = 0, iout = 1, ierr = 2, nfds = 3;
  do {
    int pact = poll(fds, nfds, 1000);
//     std::cerr << "Pipes ready: " << pact 
//      << " EV0 " <<  fds[iin].revents 
//      << " EV1 " <<  fds[iout].revents
//      << " EV2 " <<  fds[ierr].revents 
//      << " POLLOUT=" << POLLOUT << " POLLIN=" << POLLIN << " POLLPRI=" << POLLPRI 
//      << " POLLHUP=" << POLLHUP << " POLLERR=" << POLLERR 
//      << " POLLWRNORM=" << POLLWRNORM << " POLLWRBAND=" << POLLWRBAND 
//      << " POLLRDNORM=" << POLLRDNORM << " POLLRDBAND=" << POLLRDBAND << " POLLNVAL=" << POLLNVAL << "\n";
    P_ERROR_X2(pact != -1, "Error polling pipe status: ", strerror(errno));
    if (!infin)
    {
      P_ERROR_X2((fds[iin].revents & (POLLNVAL|POLLERR)) == 0, "Input pipe error -> ", err_str);
      if ((fds[iin].revents & POLLOUT) != 0)
      {
        ssize_t obytes = write (input, cxxbuf+wbytes, cxxlen-wbytes);
//        std::cerr << "runCompiler: Written bytes " << obytes << " out of " << cxxlen << "\n";
        if (obytes == -1) P_MESSAGE2("Error feeding the compiler: ", strerror(errno));
        else wbytes += static_cast<size_t>(obytes);
        if (wbytes == cxxlen) infin = true;
      }
      if ((fds[iin].revents & POLLHUP) != 0) infin = true;
      if (infin)
      {
        if (!outfin) { fds[iin] = fds[iout]; iout = iin; }
        if (!errfin) { fds[iout] = fds[ierr]; ierr = iout; }
        nfds -= 1;
        if (close (input) == -1) P_MESSAGE2("Error closing compiler input pipe: ", strerror(errno));
      }
    }
    if (!outfin)
    {
      P_ERROR_X2((fds[iout].revents & (POLLNVAL|POLLERR)) == 0, "Output pipe error -> ", err_str);
      if ((fds[iout].revents & POLLIN) != 0)
      {
        ssize_t bytes = read(output, out_buf, bufsize-1);
        ++rct;
        if (bytes > 0)
        {
//          std::cerr << "runCompiler: stdout bytes " << bytes << "\n";
          out_buf[bytes] = '\0'; 
          out_str.append(out_buf);
        } 
        else if ((bytes == -1) && (errno != EAGAIN)) P_MESSAGE2("Error reading standard output: ", strerror(errno));
      }
      if ((fds[iout].revents & POLLHUP) != 0)
      {
        outfin = true;
        if (!errfin) { fds[iout] = fds[ierr]; ierr = iout; }
        nfds -= 1;
        if (close (output) == -1) P_MESSAGE2("Error closing compiler output pipe: ", strerror(errno));
      }
    }
    if (!errfin)
    {
      P_ERROR_X2((fds[ierr].revents & (POLLNVAL|POLLERR)) == 0, "Error pipe error -> ", err_str);
      if ((fds[ierr].revents & POLLIN) != 0)
      {
        ssize_t bytes = read(error, out_buf, bufsize-1);
        ++ect;
        if (bytes > 0)
        {
//          std::cerr << "runCompiler: stderr bytes " << bytes << "\n";
          out_buf[bytes] = '\0'; 
          err_str.append(out_buf);
        } 
        else if ((bytes == -1) && (errno != EAGAIN)) P_MESSAGE2("Error reading standard error: ", strerror(errno));
      }
      if ((fds[ierr].revents & POLLHUP) != 0)
      {
        errfin = true;
        nfds -= 1;
        if (close (error) == -1) P_MESSAGE2("Error closing compiler error pipe: ", strerror(errno));
      }
    }
  } while (!outfin || !infin || !errfin);
  delete[] out_buf; out_buf = 0;
//  std::cerr << "Waited for read " << rct << " and error " << ect << " cycles."; 
  // checking status of compiler
  int status = 0;
  pid_t retval = waitpid(-1, &status, 0);
  if (retval == -1) P_MESSAGE2("Error while waiting for the compiler: ", strerror(errno));
  if (retval == 0) P_MESSAGE1("A mysterious error. (No child wanted to report status.)");
  // Check the exist status
  if (WIFEXITED(status))
  {
    if (WEXITSTATUS(status) != 0) P_MESSAGE8("Error code is ", WEXITSTATUS(status), " and the error output of the compile command '",
      cmdline, "' is\n", err_str.c_str(), "\nand the standard output is\n", out_str.c_str());
//    std::ofstream cxxfile("error.cxx");
//    cxxfile << cxxbuf;
  } else
  {
    P_MESSAGE1("The compiler was unexpectedly terminated.");
  }
}

#endif

KNExprSystem::tdlhand KNExprSystem::tdlopen(const char* fname)
{
#ifndef _WIN32
  return dlopen(fname, RTLD_NOW);
#else
  return LoadLibrary(fname);
#endif
}

void* KNExprSystem::tdlsym(tdlhand h, const char* sname)
{
#ifndef _WIN32
  void * res = dlsym(h, sname);
  return res;
#else
  return (void*) GetProcAddress(h, sname);
#endif
}

int KNExprSystem::tdlclose(KNExprSystem::tdlhand h)
{
#ifndef _WIN32
  return dlclose(h);
#else
  return (int) FreeLibrary(h);
#endif
}

const char* KNExprSystem::tdlerror()
{
#ifndef _WIN32
  return dlerror();
#else
  DWORD errcode = GetLastError();
  SetLastError(0);
  if (errcode != 0) return "Windows KNExprSystem Error\n";
  else return 0;
#endif
}

extern "C" {
  typedef size_t(*tp_sys_ndim)();
  typedef size_t(*tp_sys_npar)();
  typedef size_t(*tp_sys_ntau)();
  typedef void(*tp_sys_p_tau)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par );
  typedef void(*tp_sys_p_dtau)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par, size_t vp );
  typedef void(*tp_sys_p_rhs)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel );
  typedef void(*tp_sys_p_deri)( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
  typedef void(*tp_sys_stpar)(KNVector& par);
  typedef void(*tp_sys_p_stsol)(KNArray2D<double>& out, const KNArray1D<double>& time);
  typedef void(*tp_sys_parnames)(const char *names[]);
}

void KNExprSystem::setupFunctions (const std::string& shobj)
{
  handle = tdlopen(shobj.c_str());
  P_ERROR_X5(handle != nullptr, "Cannot open system definition file `", shobj, "'. Error code `", tdlerror(), "'.");
  
  tdlerror();    /* Clear any existing error */
  fp_ndim = reinterpret_cast<tp_sys_npar>(tdlsym(handle, "sys_ndim"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_ndim(): ", error, ".");
  
  tdlerror();    /* Clear any existing error */
  fp_npar = reinterpret_cast<tp_sys_npar>(tdlsym(handle, "sys_npar"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_npar(): ", error, ".");
  
  tdlerror();    /* Clear any existing error */
  fp_ntau = reinterpret_cast<tp_sys_ntau>(tdlsym(handle, "sys_ntau"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_ntau(): ", error, ".");
  
  /* Vectorized versions */
  tdlerror();    /* Clear any existing error */
  fp_p_tau = reinterpret_cast<tp_sys_p_tau>(tdlsym(handle, "sys_p_tau"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_p_tau(): ", error, ".");
  
  tdlerror();    /* Clear any existing error */
  fp_p_dtau = reinterpret_cast<tp_sys_p_dtau>(tdlsym(handle, "sys_p_dtau"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_p_dtau(): ", error, ".");
  
  tdlerror();    /* Clear any existing error */
  fp_p_rhs = reinterpret_cast<tp_sys_p_rhs>(tdlsym(handle, "sys_p_rhs"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_p_rhs(): ", error, ".");
  
  tdlerror();    /* Clear any existing error */
  fp_p_deri = reinterpret_cast<tp_sys_p_deri>(tdlsym(handle, "sys_p_deri"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_p_deri(): ", error, ".");
    
  tdlerror();    /* Clear any existing error */
  fp_stpar = reinterpret_cast<tp_sys_stpar>(tdlsym(handle, "sys_stpar"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_stpar(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  fp_p_stsol = reinterpret_cast<tp_sys_p_stsol>(tdlsym(handle, "sys_p_stsol"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_p_stsol(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  fp_parnames = reinterpret_cast<tp_sys_parnames>(tdlsym(handle, "sys_parnames"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_parnames(): ", error, ".");
}
