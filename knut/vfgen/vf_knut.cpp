
//  'vf_knut.cpp'
//  renamed from 'vf_pddecont.cpp'
//
//  This file defines the VectorField::PrintKnut method.
//
//  Copyright (C) 2008 Warren Weckesser
//  Copyright (C) 2009 Robert Szalai
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License, Version 2, as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program.  The file LICENSE-GPL2 in the VFGEN source directory
//  contains a copy. You may also obtain a copy by writing to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>
#include <ginac/ginac.h>

#include "vf.h"
#include "codegen_utils.h"

using namespace std;
using namespace GiNaC;


void VectorField::Knut_ConvertStateToZlags(ex& f)
{
  exset dlist;
  f.find(delay(wild(1), wild(2)), dlist);
  for (exset::const_iterator iter = dlist.begin(); iter != dlist.end(); ++iter)
  {
    ex delayfunc = *iter;
    ex delayexpr = delayfunc.op(0);
    lst vars = FindVarsInEx(delayexpr);
    ex del = delayfunc.op(1);
    int dindex = FindDelay(del);
    assert(dindex != -1);
    for (lst::const_iterator viter = vars.begin(); viter != vars.end(); ++viter)
    {
      int vindex = FindVar(ex_to<symbol>(*viter));
      delayexpr = delayexpr.subs(*viter == Zlags_(vindex, dindex + 1));
    }
    f = f.subs(delayfunc == delayexpr);
  }
  size_t nv = varname_list.nops();
  for (size_t i = 0; i < nv; ++i)
  {
    f = f.subs(varname_list[i] == Zlags_(i, 0));
  }
}

void VectorField::Knut_ConvertConstsPars(ex& f)
{
  for (size_t i = 0; i < conname_list.nops(); ++i)
  {
    f = f.subs(conname_list[i] == convalue_list[i]);
  }
  for (size_t i = 0; i < parname_list.nops(); ++i)
  {
    f = f.subs(parname_list[i] == par_(i));
  }
  for (size_t i = 0; i < internal_parname_list.nops(); ++i)
  {
    f = f.subs(internal_parname_list[i] == par_(parname_list.nops()+i));
  }
}

static void replace_Pi(std::vector<ex>& f)
{
  symbol sub_pi("M_PI");
  for (size_t i=0; i < f.size(); ++i)
  {
    f[i] = f[i].subs(Pi == sub_pi);
  }
}

//
// PrintKnut -- The Knut code generator.
//

void VectorField::PrintKnut(ostream& sys_out, map<string, string> /*options*/)
{
  size_t ex_ndim;
  size_t ex_npar;
  size_t ex_ntau;
  size_t ex_nevent;
  std::vector<GiNaC::ex> ex_tau;
  std::vector<GiNaC::ex> ex_tau_p;
  std::vector<GiNaC::ex> ex_rhs;
  std::vector<double> ex_mass;
  std::vector<GiNaC::ex> ex_rhs_p;
  std::vector<GiNaC::ex> ex_rhs_x;
  std::vector<GiNaC::ex> ex_rhs_xp;
  std::vector<GiNaC::ex> ex_rhs_xx;
  std::vector<double> ex_stpar;
  std::vector<GiNaC::ex> ex_stsol;
  std::vector<std::string> ex_parnames;

  ex_ndim = Knut_ndim();
  ex_npar = Knut_npar();
  ex_ntau = Knut_ntau();
  ex_nevent = Knut_nevent();
  Knut_tau(ex_tau);
  Knut_tau_p(ex_tau_p);
  Knut_RHS(ex_rhs);
  Knut_mass(ex_mass);
  Knut_RHS_p(ex_rhs_p);
  Knut_RHS_x(ex_rhs_x);
  Knut_RHS_xp(ex_rhs_xp);
  Knut_RHS_xx(ex_rhs_xx);
  Knut_stpar(ex_stpar);
  Knut_stsol(ex_stsol);
  Knut_parnames(ex_parnames);
  
  replace_Pi(ex_tau);
  replace_Pi(ex_tau_p);
  replace_Pi(ex_rhs);
  replace_Pi(ex_rhs_p);
  replace_Pi(ex_rhs_x);
  replace_Pi(ex_rhs_xp);
  replace_Pi(ex_rhs_xx);
  replace_Pi(ex_stsol);
  
  size_t nc = conname_list.nops();
  size_t nv = varname_list.nops();
  size_t np = parname_list.nops();
  size_t nf = funcname_list.nops();

  //
  //  Create the system definition file.
  //
  sys_out << csrc_double; // Make the output print C sorce formatting (double values)
  sys_out << "//" << endl;
  sys_out << "// Knut KNSystem Definition file for the VFGEN vector field: " << Name() << endl;
  sys_out << "//" << endl;
  PrintVFGENComment(sys_out, "// ");
  sys_out << "//" << endl;
  sys_out << endl;
  // sys_out << "#include <cstdlib>\n";
  sys_out << "#include <cmath>\n";
//  sys_out << "#include <iostream>\n";
  sys_out << "#include \"knutsys.h\"\n";
  sys_out << endl;
  sys_out << "static inline double heaviside(double x)\n";
  sys_out << "{\n";
  sys_out << "  if (x > 0) return 1.0;\n";
  sys_out << "  else return 0.0;\n";
  sys_out << "}\n";
  sys_out << endl;
  sys_out << "static inline double ramp(double x)\n";
  sys_out << "{\n";
  sys_out << "  if (x > 0) return x;\n";
  sys_out << "  else return 0.0;\n";
  sys_out << "}\n";
  sys_out << endl;
  sys_out << "extern \"C\" {\n";
  sys_out << endl;
  sys_out << "size_t sys_ndim()  { return " << ex_ndim << "; }  // KNSystem dimension\n";
  sys_out << "size_t sys_npar()  { return " << ex_npar << "; }  // Number of parameters, plus one (for the period)\n";
  sys_out << "size_t sys_ntau()  { return " << ex_ntau << "; }  // Number of delays, plus one\n";
  sys_out << "size_t sys_nderi() { return 2; }  // Order of derivatives computed here\n";
  sys_out << "size_t sys_nevent() { return " << ex_nevent << "; }  // Number of event functions\n";
  sys_out << endl;
    
  // ***************************************************************************
  //  Create the vector field function sys_rhs
  // ***************************************************************************
  sys_out << "//" << endl;
  sys_out << "// sys_p_rhs(...) computes the vector field.\n";
  sys_out << "//" << endl;
  sys_out << "void sys_p_rhs(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_, size_t sel)\n";
  sys_out << "{\n";
//  sys_out << "  std::cout << \" rhs\" << \"\\n\"; std::cout.flush();\n";
  sys_out   << "  // Compute the vector field\n";
  sys_out   << "  for (size_t idx=0; idx < time.size(); ++idx)\n";
  sys_out   << "  {\n";
  sys_out   << "    const double t = time(idx);\n";
  for (size_t i = 0; i < ex_rhs.size(); ++i)
  {
    sys_out << "    out(" << i << ",idx)" << " = " << ex_rhs[i] << ";" << endl;
  }
  sys_out   << "  }\n";
  sys_out   << "}\n";

  // ***************************************************************************
  // Mass vector
  // ***************************************************************************
  sys_out << "void sys_mass(KNArray1D<double>& out)\n";
  sys_out << "{\n";
  for (size_t j = 0; j < ex_rhs.size(); ++j)
  {
    sys_out << "  out(" << j << ") = " << ex_mass[j] << ";\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  
  // ***************************************************************************
  //  Derivative with respect to the delays
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// This function computes the Jacobian of the vector field\n";
  sys_out << "// with respect to the k-th delayed vector.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << "__jacx(KNArray3D<double>& out, const KNArray1D<double>& time, size_t k, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
//  sys_out << "  std::cout << \" del=\" << k << \"\\n\";\n";
  for (size_t del = 0; del < ex_ntau; ++del)
  {
    if (del == 0)
    {
      sys_out << "  if (k == 0)\n";
      sys_out << "  {\n";
      sys_out << "    // Derivatives wrt the state variables\n";
    }
    else
    {
      sys_out << "  else if (k == " << del << ")\n";
      sys_out << "  {\n";
      sys_out << "    // Derivatives wrt state variables with delay " << Delays[del-1] << endl;
    }
    sys_out << "    for (size_t idx=0; idx < time.size(); ++idx)\n";
    sys_out << "    {\n";
    sys_out << "      const double t = time(idx);\n";
    for (size_t i = 0; i < nv; ++i)
    {
      for (size_t j = 0; j < nv; ++j)
      {      
        sys_out << "      out(" << i << "," << j << ",idx)" << " = " << ex_rhs_x[i + (del + j*ex_ntau)*nv] << ";" << endl;
      }
    }	
    sys_out << "    }\n";
    sys_out << "  }\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  //  Derivative with respect to the parameters
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// This function computes the Jacobian of the vector field\n";
  sys_out << "// with respect to the j-th parameter.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << "__jacp(KNArray3D<double>& out, const KNArray1D<double>& time, size_t k, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  for (size_t alpha = 0; alpha < ex_npar; ++alpha)
  {
    if (alpha == 0)
    {
      sys_out << "  if (k == 0)\n";
      sys_out << "  {\n";
    }
    else
    {
      sys_out << "  else if (k == " << alpha << ")\n";
      sys_out << "  {\n";
    }
    sys_out << "    for (size_t idx=0; idx < time.size(); ++idx)\n";
    sys_out << "    {\n";
    sys_out << "      const double t = time(idx);\n";
    for (size_t i = 0; i < nv; ++i)
    {      
      sys_out << "      out(" << i << ",0,idx)" << " = " << ex_rhs_p[i + alpha*nv] << ";" << endl;
    }
    sys_out << "    }\n";
    sys_out << "  }\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  //  Derivative with respect to the delays and parameters (2nd order)
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// This function computes the Jacobian of the vector field\n";
  sys_out << "// with respect to the k-th delayed state vector and the j-th parameter.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << "__jacxp(KNArray3D<double>& out, const KNArray1D<double>& time, size_t k, size_t j, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  for (size_t alpha = 0; alpha < ex_npar; ++alpha)
  {
    if (alpha == 0)
    {
      sys_out << "  if (j == 0)\n";
      sys_out << "  {\n";
    }
    else
    {
      sys_out << "  else if (j == " << alpha << ")\n";
      sys_out << "  {\n";
    }
    for (size_t del = 0; del < ex_ntau; ++del)
    {
      if (del == 0)
      {
        sys_out << "    if (k == 0)\n";
        sys_out << "    {\n";
      }
      else
      {
        sys_out << "    else if (k == " << del << ")\n";
        sys_out << "    {\n";
      }
      sys_out << "      for (size_t idx=0; idx < time.size(); ++idx)\n";
      sys_out << "      {\n";
      sys_out << "        const double t = time(idx);\n";
      for (size_t i = 0; i < nv; ++i)
      {
        for (size_t j = 0; j < nv; ++j)
        {      
          sys_out << "        out(" << i << "," << j << ",idx)" << " = " << ex_rhs_xp[i + (del + (alpha + j*ex_npar)*ex_ntau)*nv] << ";" << endl;
        }
      }
      sys_out << "      }\n";
      sys_out << "    }\n";
    }
    sys_out << "  }\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  //  Second order derivative with respect to the delays 
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// This function computes the Hessian of the vector field\n";
  sys_out << "// with respect to the k1-th and k2-th delayed state vectors,\n";
  sys_out << "// then multiplies this by the m-th column of v to obtain a matrix.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << "__hess_times_v(KNArray3D<double>& out, const KNArray1D<double>& time, size_t k1, size_t k2, int m, const KNArray3D<double>& VZlags_, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  for (size_t del1 = 0; del1 < ex_ntau; ++del1)
  {
    if (del1 == 0)
    {
      sys_out << "  if (k1 == 0)\n";
      sys_out << "  {\n";
    }
    else
    {
      sys_out << "  else if (k1 == " << del1 << ")\n";
      sys_out << "  {\n";
    }
    for (size_t del2 = 0; del2 < ex_ntau; ++del2)
    {
      if (del2 == 0)
      {
        sys_out << "    if (k2 == 0)\n";
        sys_out << "    {\n";
      }
      else
      {
        sys_out << "    else if (k2 == " << del2 << ")\n";
        sys_out << "    {\n";
      }
      sys_out   << "      for (size_t idx=0; idx < time.size(); ++idx)\n";
      sys_out   << "      {\n";
      sys_out   << "        const double t = time(idx);\n";
      for (size_t i = 0; i < nv; ++i)
      {
        for (size_t j = 0; j < nv; ++j)
        {      
          sys_out << "        out(" << i << "," << j << ",idx)" << " = " << ex_rhs_xx[i + nv*(j + nv*(del1 + ex_ntau*del2))] << ";" << endl;
        }
      }
      sys_out << "      }\n";
      sys_out << "    }\n";
    }
    sys_out << "  }\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  //  Create the derivatives function sys_p_deri(...)
  // ***************************************************************************
  sys_out << "//" << endl;
  sys_out << "// sys_p_deri(...)\n";
  sys_out << "//" << endl;
  // Include a list of the lags in the comments.
  sys_out << "// The lags are: {";
  for (unsigned k = 0; k < Delays.size(); ++k)
  {
    sys_out << Delays[k];
    if (k < Delays.size() - 1)
      sys_out << ", ";
  }
  sys_out << "}" << endl;
  sys_out << "//\n";
  sys_out << "// If X(t) is the state vector at time t, then\n";
  sys_out << "//    Zlags_ = [ X(t) ";
  for (unsigned k = 0; k < Delays.size(); ++k)
  {
    sys_out << "X(t-(" << Delays[k] << "))";
    if (k < Delays.size() - 1)
      sys_out << " ";
  }
  sys_out << " ]\n";
  sys_out << "//\n";
  sys_out << "// The state vector:\n";
  GetFromVector2(sys_out, "// ", varname_list, "Zlags_", "(", ",0,idx)", 0, ";");
  sys_out << "//\n";
  // Function definition starts here.
  // ( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
  sys_out << "//" << endl;
  sys_out << "void sys_p_deri(KNArray3D<double>& jac_, const KNArray1D<double>& time, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_, size_t sel,\nsize_t nx_, const size_t* vx_, size_t np_, const size_t* vp_, const KNArray3D<double>& v_)\n";
  sys_out << "{\n";
//  sys_out << "  std::cout << \" nx=\" << nx_ << \" np=\" << np_ << \"\\n\"; std::cout.flush();\n";
  sys_out << "    if (nx_ == 1 && np_ == 0)\n";
  sys_out << "        " << "__jacx(jac_, time, vx_[0], Zlags_, par_);\n";
  sys_out << "    else if (nx_ == 0 && np_ == 1)\n";
  sys_out << "        " << "__jacp(jac_, time, vp_[0], Zlags_, par_);\n";
  sys_out << "    else if (nx_ == 1 && np_ == 1)\n";
  sys_out << "        " << "__jacxp(jac_,time, vx_[0], vp_[0], Zlags_, par_);\n";
  sys_out << "    else if (nx_ == 2 && np_ == 0)\n";
  sys_out << "        " << "__hess_times_v(jac_, time, vx_[0], vx_[1], vx_[0], v_, Zlags_, par_);\n";
  sys_out << "    else\n";
  sys_out << "    {\n";
  sys_out << "        // std::cerr << \"sys_deri: Requested derivative has not been implemented.\\n\";\n";
  sys_out << "        // exit(-1);\n";
  sys_out << "    }\n";
  sys_out << "}\n";

  // ***************************************************************************
  // Create sys_p_tau()
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//" << endl;
  sys_out << "// sys_p_tau(...) computes the vector of delays.\n";
  sys_out << "//" << endl;
  sys_out << endl;
  sys_out << "void sys_p_tau(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  sys_out << "  for (size_t idx=0; idx < time.size(); ++idx)\n";
  sys_out << "  {\n";
  for (size_t j = 0; j < ex_ntau; j++)
  {
    sys_out << "    out(" << j << ",idx) = " <<  ex_tau[j] << ";\n";
  }
  sys_out << "  }\n";
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  // Create sys_p_dtau()
  // ***************************************************************************
  sys_out << "//\n";
  sys_out << "// sys_p_dtau(...) computes the derivatives of the delays with respect to the parameters.\n";
  sys_out << "//\n";
  sys_out << "// The delays are: ";
  for (unsigned k = 0; k < Delays.size(); ++k)
  {
    sys_out << Delays[k];
    if (k < Delays.size() - 1)
      sys_out << ", ";
  }
  sys_out << endl;
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "void sys_p_dtau(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par_, size_t k)\n";
  sys_out << "{\n";
  for (size_t alpha = 0; alpha < ex_npar; ++alpha)
  {
    if (alpha == 0)
    {
      sys_out << "  if (k == 0)\n";
      sys_out << "  {\n";
    }
    else
    {
      sys_out << "  else if (k == " << alpha << ")\n";
      sys_out << "  {\n";
    }
    sys_out << "    for (size_t idx=0; idx < time.size(); ++idx)\n";
    sys_out << "    {\n";
    sys_out << "      const double t = time(idx);\n";
    for (size_t j = 0; j < ex_ntau; j++)
    {
      sys_out << "      out(" << j << ",idx) = " << ex_tau_p[j + alpha*ex_ntau] << ";\n";
    }
    sys_out << "    }\n";
    sys_out << "  }\n";
  }
  sys_out << "}\n";
  sys_out << endl;

  // ***************************************************************************
  // Starting parameters
  // ***************************************************************************
  sys_out << "void sys_stpar(KNVector& par_)\n";
  sys_out << "{\n";
  for (size_t j = 0; j < ex_npar; ++j)
  {
    sys_out << "  par_(" << j << ") = " << ex_stpar[j] << ";\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  
  // ***************************************************************************
  // Starting solution
  // ***************************************************************************
  sys_out << "void sys_stsol(KNVector& out, double t)\n";
  sys_out << "{\n";
  for (size_t j = 0; j < nv; ++j)
  {
    sys_out << "  out(" << j << ") = " << ex_stsol[j] << ";\n";
  }
  sys_out << "}\n";
  sys_out << endl;

  // ***************************************************************************
  // Parameter names
  // ***************************************************************************
  sys_out << "void sys_parnames( const char *out[] )\n";
  sys_out << "{\n";
  for (size_t k = 0; k < ex_npar; ++k)
  {
    sys_out << "  out[" << k << "] = \"" << ex_parnames[k] << "\";\n";
  }
  sys_out << "}\n";

  // ***************************************************************************
  // Event functions
  // ***************************************************************************
  sys_out << endl;
  // TODO: implement event functions          
  sys_out << "}  // extern \"C\"\n";

  return;
}

// These routines for use without a compiler
void VectorField::Knut_tau(std::vector<GiNaC::ex>& exls)
{
  const size_t nd = Delays.size()+1;
  exls.resize(nd);
  exls[0] = 0;
  for (size_t i = 1; i < nd; ++i)
  {
    ex f = Delays[i-1];
    Knut_ConvertConstsPars(f);
    exls[i] = f;
  }
}

void VectorField::Knut_tau_p(std::vector<GiNaC::ex>& exls)
{
  const size_t nd = Delays.size()+1;
  const size_t np = parname_list.nops();
  size_t par_shift = 1;
  if (HasPeriod) par_shift = 0;
  exls.resize(nd*(np+par_shift));
  for (size_t i = 1; i < nd; ++i) 
  {
    ex f = Delays[i-1];
    for (size_t j = 0; j < np; ++j)
    {
      symbol ps = ex_to<symbol>(parname_list[j]);
      ex df = f.diff(ps); 
      Knut_ConvertConstsPars(df);
      exls[i + (j+par_shift)*nd] = df;
    }
  }
}

void VectorField::Knut_RHS(std::vector<GiNaC::ex>& exls)
{
  exls.resize(varvecfield_list.nops());
  for (size_t i = 0; i < varvecfield_list.nops(); ++i)
  {
    ex f = varvecfield_list[i];
    Knut_ConvertStateToZlags(f);
    Knut_ConvertConstsPars(f);
    exls[i] = f;
  }
}

void VectorField::Knut_RHS_p(std::vector<GiNaC::ex>& exls)
{
  const size_t nv = varvecfield_list.nops();
  const size_t np = parname_list.nops();
  size_t par_shift = 1;
  if (HasPeriod) par_shift = 0;
  exls.resize(nv*(np+par_shift));
  for (size_t i = 0; i < nv; ++i)
  {
    ex f = varvecfield_list[i];
    Knut_ConvertStateToZlags(f);
    for (size_t j = 0; j < np; ++j)
    {
      symbol ps = ex_to<symbol>(parname_list[j]);
      ex df = f.diff(ps);
      Knut_ConvertConstsPars(df);
      exls[i + (j+par_shift)*nv] = df;
//      std::cout << "PD: " << i << "," << j << " -> " << df << "\n"; 
    }
  }
}

void VectorField::Knut_RHS_x(std::vector<GiNaC::ex>& exls)
{
  const size_t nv = varname_list.nops();
  const size_t nd = Delays.size()+1;
  exls.resize(nd*nv*nv);
  for (size_t i = 0; i < nv; ++i)
  {
    ex f = varvecfield_list[i];
    Knut_ConvertStateToZlags(f);
    Knut_ConvertConstsPars(f);
    for (size_t j = 0; j < nv; ++j)
    {
      for (size_t k = 0; k < nd; ++k)
      {
        symbol vtmp_("vtmp_");
        ex fj = f.subs(Zlags_(j, k) == vtmp_);
        ex df = fj.diff(vtmp_);
        df = df.subs(vtmp_ == Zlags_(j, k));
        exls[i + (k + j*nd)*nv] = df;
//        std::cout << "PX: " << i << "," << j << "," << k << " -> " << df << "\n";
      }
    }
  }
}

void VectorField::Knut_RHS_xp(std::vector<GiNaC::ex>& exls)
{
  const size_t nv = varname_list.nops();
  const size_t nd = Delays.size()+1;
  const size_t np = parname_list.nops();
  size_t par_shift = 1;
  if (HasPeriod) par_shift = 0;
  const size_t nps = np + par_shift;
  exls.resize(nd*nv*nv*nps);
  for (size_t i = 0; i < nv; ++i)
  {
    ex f = varvecfield_list[i];
    Knut_ConvertStateToZlags(f);
    for (size_t q = 0; q < np; ++q)
    {
      symbol ps = ex_to<symbol>(parname_list[q]);
      ex dfp = f.diff(ps);
      Knut_ConvertConstsPars(dfp);
      for (size_t j = 0; j < nv; ++j)
      {
        for (size_t k = 0; k < nd; ++k)
        {
          symbol vtmp_("vtmp_");
          ex fj = dfp.subs(Zlags_(j, k) == vtmp_);
          ex df = fj.diff(vtmp_);
          df = df.subs(vtmp_ == Zlags_(j, k));
          exls[i + (k + (q + par_shift + j*nps)*nd)*nv] = df;
//          std::cout << "PXp: " << nd << "|" << i << "," << j << "," << k << " -> " << df << "\n";
        }
      }
    }
  }
}

void VectorField::Knut_RHS_xx(std::vector<GiNaC::ex>& exls)
{
  const size_t nv = varname_list.nops();
  const size_t nd = Delays.size()+1;
  exls.resize(nd*nd*nv*nv);
  for (size_t i = 0; i < nv; ++i)
  {
    ex f = varvecfield_list[i];
    Knut_ConvertStateToZlags(f);
    Knut_ConvertConstsPars(f);
    for (size_t k1 = 0; k1 < nd; ++k1)
    {
      ex pdf1;
      for (size_t j = 0; j < nv; ++j)
      {
        symbol vtmp_("vtmp_");
        ex fj = f.subs(Zlags_(j, k1) == vtmp_);
        ex df = fj.diff(vtmp_);
        df = df.subs(vtmp_ == Zlags_(j, k1));
        pdf1 += df*VZlags_(j,k1); 
      }
      for (size_t j = 0; j < nv; ++j)
      {
        for (size_t k2 = 0; k2 < nd; ++k2)
        {
          symbol vtmp_("vtmp_");
          ex fj = pdf1.subs(Zlags_(j, k2) == vtmp_);
          ex df = fj.diff(vtmp_);
          df = df.subs(vtmp_ == Zlags_(j, k2));
          exls[i + nv*(j + nv*(k1 + nd*k2))] = df;
//          std::cout << "PXX: " << i << "," << j << ", (" << k1 << "," << k2 << ") -> " << df << "\n";
        }
      }
    }
  }
}

void VectorField::Knut_mass(std::vector<double>& out)
{  
  const size_t nc = conname_list.nops();
  const size_t nv = varmass_list.nops();
  out.resize(nv);
  for (size_t j = 0; j < nv; ++j)
  {
    // Make a replacement list with all constant values
    GiNaC::lst parsubs;
    for (size_t k = 0; k < nc; ++k)
    {
      parsubs.append(conname_list[k] == convalue_list[k]);
    }
    ex mass = varmass_list[j];
    mass = iterated_subs(mass, parsubs).evalf();
    if (is_exactly_a<numeric>(mass)) out[j] = ex_to<numeric>(mass).to_double();
    else P_MESSAGE1("Mass constant is not a number.");
  }
}

void VectorField::Knut_stpar(std::vector<double>& par)
{  
  const size_t nc = conname_list.nops();
  const size_t np = parname_list.nops();
  size_t par_shift = 1;
  if (HasPeriod) par_shift = 0;
  else par[0] = 1.0;
  par.resize(np+par_shift);
  for (size_t j = 0; j < np; ++j)
  {
    // First make a replacement list with all the other parameters
    GiNaC::lst parsubs;
    // add internal parameters
    for (size_t k = 0; k < internal_pardefval_list.nops(); ++k)
    {
      parsubs.append(internal_parname_list[k] == internal_pardefval_list[k]);
    }
    // add other parameters
    for (size_t k = 0; k < np; ++k)
    {
      if (k != j)
      {
        parsubs.append(parname_list[k] == pardefval_list[k]);
      }
    }
    // add constants
    for (size_t k = 0; k < nc; ++k)
    {
      parsubs.append(conname_list[k] == convalue_list[k]);
    }
    ex defval = pardefval_list[j];
    defval = iterated_subs(defval, parsubs).evalf();
    if (is_exactly_a<numeric>(defval)) par[j + par_shift] = ex_to<numeric>(defval).to_double();
    else P_MESSAGE1("Starting parameter is not a number.");
  }
}

void VectorField::Knut_stsol(std::vector<GiNaC::ex>& exls)
{
  const size_t nv = varname_list.nops();
  const size_t nc = conname_list.nops();
  const size_t np = parname_list.nops();
  exls.resize(nv);
  // parameter substitution
  GiNaC::lst parsubs;
  // add other parameters
  for (size_t k = 0; k < np; ++k)
  {
    parsubs.append(parname_list[k] == pardefval_list[k]);
  }
  // add constants
  for (size_t k = 0; k < nc; ++k)
  {
    parsubs.append(conname_list[k] == convalue_list[k]);
  }
  for (size_t j = 0; j < nv; ++j)
  {
    ex defsol = vardefic_list[j];
    defsol = iterated_subs(defsol, parsubs);
    exls[j] = defsol;
  }
}

void VectorField::Knut_parnames(std::vector<std::string>& pnames)
{
  const size_t np = parname_list.nops();  
  size_t par_shift = 1;
  if (HasPeriod) par_shift = 0;
  pnames.resize(np+par_shift);
  for (size_t k = 0; k < np; ++k)
  {
    if (is_exactly_a<symbol>(parname_list[k])) pnames[k+par_shift] = ex_to<symbol>(parname_list[k]).get_name();
    else P_MESSAGE1("Parameter is not a symbol.");
  }
}

static double loc_heaviside(double x)
{
  if (x > 0) return 1.0;
  else return 0.0;
}

static void expeval(const GiNaC::ex& expr, double* res, size_t skip, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray3D<double>& v, const KNVector& par)
{
  const size_t vsize = time.size();
  if (is_exactly_a<add>(expr))
  {
    double *res_loc = new double[vsize];
    for (size_t k = 0; k < vsize; ++k) res[k*skip] = 0; 
//    std::cout << " add -> (";
    for (const_iterator i = expr.begin(); i != expr.end(); ++i)
    {
      expeval(*i, res_loc, 1, time, x, v, par);
      for (size_t k = 0; k < vsize; ++k)
      {
        res[k*skip] += res_loc[k];
      }
    }
//    std::cout << ")";
    delete[] res_loc;
    return;
  } else if (is_exactly_a<mul>(expr))
  {
    double *res_loc = new double[vsize];
    for (size_t k = 0; k < vsize; ++k) res[k*skip] = 1.0; 
//    std::cout << " mul -> (";
    for (const_iterator i = expr.begin(); i != expr.end(); ++i)
    {
      expeval(*i, res_loc, 1, time, x, v, par);
      for (size_t k = 0; k < vsize; ++k)
      {
        res[k*skip] *= res_loc[k];
      }
    }
//    std::cout << ")";
    delete[] res_loc;
    return;
  } else if (is_exactly_a<power>(expr))
  {
    double *res_loc = new double[vsize];
    double *res_loc_p = new double[vsize];
//    std::cout << " power -> " << expr.nops() << "(" ;
    expeval(expr.op(0), res_loc, 1, time, x, v, par);
    expeval(expr.op(1), res_loc_p, 1, time, x, v, par);
    for (size_t k = 0; k < vsize; ++k)
    {
      res[k*skip] = std::pow(res_loc[k], res_loc_p[k]);
    }
//    std::cout << ")";
    delete[] res_loc_p;
    delete[] res_loc;
    return;
  } 
  else if (is_exactly_a<function>(expr))
  {
//    std::cout << " function -> " << ex_to<function>(expr).get_name() << "(" << expr.op(0) << ")";
    std::string name(ex_to<function>(expr).get_name());
    if (name == "par_")
    {
      size_t pid = static_cast<size_t>(ex_to<numeric>(expr.op(0)).to_int());
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = par(pid);
      return;
    }
    if (name == "Zlags_")
    {
      size_t d1 = static_cast<size_t>(ex_to<numeric>(expr.op(0)).to_int());
      size_t d2 = static_cast<size_t>(ex_to<numeric>(expr.op(1)).to_int());
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = x(d1,d2,k);
      return;
    }
    if (name == "VZlags_")
    {
      size_t d1 = static_cast<size_t>(ex_to<numeric>(expr.op(0)).to_int());
      size_t d2 = static_cast<size_t>(ex_to<numeric>(expr.op(1)).to_int());
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = v(d1,d2,k);
      return;
    }
    double (*fun_ptr)(double) = 0;
    if (name == "abs") fun_ptr = &fabs;
//    if (name == "step") return;
//    if (name == "csgn") return;
//    if (name == "eta") return;
//    if (name == "Li2") return;
//    if (name == "Li3") return;
//    if (name == "zetaderiv") return;
//    if (name == "factorial") return;
//    if (name == "binomial") return;
//    if (name == "Order") return;
    if (name == "exp") fun_ptr = &exp;
    if (name == "log") fun_ptr = &log;
    if (name == "sin") fun_ptr = &sin;
    if (name == "cos") fun_ptr = &cos;
    if (name == "tan") fun_ptr = &tan;
    if (name == "asin") fun_ptr = &asin;
    if (name == "acos") fun_ptr = &acos;
    if (name == "atan") fun_ptr = &atan;
//    if (name == "atan2") return;
    if (name == "sinh") fun_ptr = &sinh;
    if (name == "cosh") fun_ptr = &cosh;
    if (name == "tanh") fun_ptr = &tanh;
    if (name == "asinh") fun_ptr = &asinh;
    if (name == "acosh") fun_ptr = &acosh;
    if (name == "atanh") fun_ptr = &atanh;
    if (name == "lgamma") fun_ptr = &lgamma;
    if (name == "tgamma") fun_ptr = &tgamma;
//    if (name == "beta") return;
//    if (name == "psi") return;
//    if (name == "G") return;
//    if (name == "Li") return;
//    if (name == "S") return;
//    if (name == "H") return;
//    if (name == "zeta") return;
    if (name == "heaviside") fun_ptr = &loc_heaviside;
    if (fun_ptr)
    {
      double *res_loc = new double[vsize]; 
      expeval(expr.op(0), res_loc, 1, time, x, v, par); 
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = (*fun_ptr)(res_loc[k]); 
      delete[] res_loc; 
      return;
    } else
    {
      P_MESSAGE2("Not implemented function ", name);
    }
  }
  else if (is_exactly_a<numeric>(expr))
  {
//    std::cout << " numeric -> " << ex_to<numeric>(expr).to_double();
    for (size_t k = 0; k < vsize; ++k) res[k*skip] = ex_to<numeric>(expr).to_double();
    return;
  }
  else if (is_exactly_a<symbol>(expr))
  {
    std::string name = ex_to<symbol>(expr).get_name();
    if (name == "t")
    {
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = time(k);
      return;
    } else P_MESSAGE3("Unknown symbol ", name, ".");
  }
  else if (is_exactly_a<constant>(expr))
  {
    std::string name = ex_to<symbol>(expr).get_name();
    if ((name == "Pi") || (name == "\\pi"))
    {
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = M_PI;
      return;
    } else P_MESSAGE3("Unknown constant ", name, ".");
  } else
  {
    std::ostringstream os;
    os << expr;
    P_MESSAGE3("Unknown formula component: ", os.str(), ".");
  }
}

void VectorField::Knut_tau_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, const size_t, const size_t, const size_t nd )
{
  KNArray3D<double> dummy_x, dummy_v;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
  for (size_t i = 0; i < nd; ++i)
  {
    expeval(exls[i], out.pointer(i,0), skip, time, dummy_x, dummy_v, par);
  }
}

void VectorField::Knut_tau_p_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp, const size_t, const size_t, const size_t nd )
{
  KNArray3D<double> dummy_x, dummy_v;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
  for (size_t i = 0; i < nd; ++i)
  {
    expeval(exls[i + vp*nd], out.pointer(i,0), skip, time, dummy_x, dummy_v, par);
  }
}


void VectorField::Knut_RHS_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t /*sel*/ )
{
  KNArray3D<double> dummy_v;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
//  std::cout << "Knut_RHS_eval " << exls[0] << " " << skip << " \n";
  for (size_t n = 0; n < exls.size(); ++n)
  {
    expeval(exls[n], out.pointer(n,0), skip, time, x, dummy_v, par);
  }
}

void VectorField::Knut_RHS_p_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t /*sel*/, size_t alpha, const size_t nv, const size_t /*nps*/, const size_t /*nd*/)
{
  KNArray3D<double> dummy_v;
  const size_t skip = (((size_t)out.pointer(0,0,1)) - ((size_t)out.pointer(0,0,0)))/sizeof(double);
  for (size_t i = 0; i < nv; ++i)
  {
    expeval(exls[i + alpha*nv], out.pointer(i,0,0), skip, time, x, dummy_v, par);
  }	
}

void VectorField::Knut_RHS_x_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t /*sel*/, size_t del, const size_t nv, const size_t /*nps*/, const size_t nd)
{
  KNArray3D<double> dummy_v;
  for (size_t i = 0; i < nv; ++i)
  {
    for (size_t j = 0; j < nv; ++j)
    {
      const size_t skip = (((size_t)out.pointer(i,j,1)) - ((size_t)out.pointer(i,j,0)))/sizeof(double);
//      std::cout << "PiX: " << nd << "|" << i << "," << j << "," << del << " -> " << exls[i + (del + j*nd)*nv] << "\n";
      expeval(exls[i + (del + j*nd)*nv], out.pointer(i,j,0), skip, time, x, dummy_v, par);
    }
  }	
}

void VectorField::Knut_RHS_xp_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t /*sel*/, size_t del, size_t alpha, const size_t nv, const size_t nps, const size_t nd)
{
  KNArray3D<double> dummy_v;
  for (size_t i = 0; i < nv; ++i)
  {
    for (size_t j = 0; j < nv; ++j)
    {
      const size_t skip = (((size_t)out.pointer(i,j,1)) - ((size_t)out.pointer(i,j,0)))/sizeof(double);
//      std::cout << "PiXp: " << i << "," << j << "," << del << " -> " << exls[i + (del + (alpha + j*nps)*nd)*nv] << "\n";
      expeval(exls[i + (del + (alpha + j*nps)*nd)*nv], out.pointer(i,j,0), skip, time, x, dummy_v, par);
    }
  }	
}

void VectorField::Knut_RHS_xx_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray3D<double>& vv, const KNVector& par, size_t /*sel*/, size_t del1, size_t del2, const size_t nv, const size_t /*nps*/, const size_t nd)
{
  KNArray3D<double> dummy_v;
  for (size_t i = 0; i < nv; ++i)
  {
    for (size_t j = 0; j < nv; ++j)
    {
      const size_t skip = (((size_t)out.pointer(i,j,1)) - ((size_t)out.pointer(i,j,0)))/sizeof(double);
//      std::cout << "PiXX: " << nd << "|" << i << "," << j << "," << del << " -> " << exls[i + (del + j*nd)*nv] << "\n";
      expeval(exls[i + nv*(j + nv*(del1 + nd*del2))], out.pointer(i,j,0), skip, time, x, vv, par);
    }
  }	
}

void VectorField::Knut_RHS_stsol_eval( const std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time )
{
  KNArray3D<double> dummy_x, dummy_v;
  KNVector par;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
  for (size_t i = 0; i < exls.size(); ++i)
  {
    expeval(exls[i], out.pointer(i,0), skip, time, dummy_x, dummy_v, par);
  }
}
