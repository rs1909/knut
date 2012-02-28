
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

//
// Knut_ConvertDelaysToZlags(ex& f)
//
// This function converts each subexpression of the form delay(delayexpr,del)
// in f to an expression in terms of Zlags_(i,j).
//

void VectorField::Knut_ConvertDelaysToZlags(ex& f)
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
}

//
// Knut_ConvertStateToZlags(ex& f)
//
// This function converts all references to state variables--current time or
// delayed--to an expression in terms of Zlags_(i,j).
//

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
}

void VectorField::Knut_PrintParDerivs(ostream &dout, const vector<ex> &vf0)
{
  // int nc = conname_list.nops();
  size_t nv = varname_list.nops();
  size_t np = parname_list.nops();
  // int na = exprname_list.nops();
  // int nf = funcname_list.nops();
  int par_shift = 1;
  if (HasPeriod) par_shift = 0;

  for (size_t j = 0; j < np; ++j)
  {
    if (j == 0)
    {
      dout << "        if (j == " << j + par_shift << ")\n";
      dout << "        {\n";
    }
    else
    {
      dout << "        else if (j == " << j + par_shift << ")\n";
      dout << "        {\n";
    }
    dout << "            // Derivative wrt " << parname_list[j] << endl;
    for (size_t i = 0; i < nv; ++i)
    {
      symbol p = ex_to<symbol>(parname_list[j]);
      ex df = vf0[i].diff(p);
      dout << "            jac_(" << i << ",0,idx)" << " = " << df << ";" << endl;
    }
    dout << "        }\n";
  }
}

void VectorField::Knut_PrintJacobians(ostream &dout, const vector<ex> &vf0)
{
  // int nc = conname_list.nops();
  size_t nv = varname_list.nops();
  // int np = parname_list.nops();
  // int na = exprname_list.nops();
  // int nf = funcname_list.nops();
  size_t nd = Delays.size();

  std::vector<std::vector<size_t> > row(nd+1);
  std::vector<std::vector<size_t> > col(nd+1);

  for (size_t k = 0; k < nd + 1; ++k)
  {
    if (k == 0)
    {
      dout << "        if (k == 0)\n";
      dout << "        {\n";
      dout << "            // Derivatives wrt the state variables\n";
    }
    else
    {
      dout << "        else if (k == " << k << ")\n";
      dout << "        {\n";
      dout << "            // Derivatives wrt state variables with delay " << Delays[k-1] << endl;
    }
    for (size_t i = 0; i < nv; ++i)
    {
      ex f = vf0[i];
      for (size_t j = 0; j < nv; ++j)
      {
        symbol vtmp_("vtmp_");
        ex fj = f.subs(Zlags_(j, k) == vtmp_);
        ex df = fj.diff(vtmp_);
        df = df.subs(vtmp_ == Zlags_(j, k));
        if (!df.is_zero()) { row[k].push_back(i); col[k].push_back(j); }
        dout << "            jac_(" << i << "," << j << ",idx)" << " = " << df << ";" << endl;
      }
    }
    dout << "        }\n";
  }
  for (size_t i=0; i<row.size(); ++i)
  {
    dout << "        // Sparsity structure for tau = " << i << "\n";
    for (size_t j=0; j<row[i].size(); ++j)
    {
      dout << "        // (" << row[i][j] << ", " << col[i][j] << ")\n";
    }
  }
}

void VectorField::Knut_PrintXandParJacobians(ostream &dout, const vector<ex> &vf0)
{
  // int nc = conname_list.nops();
  size_t nv = varname_list.nops();
  size_t np = parname_list.nops();
  // int na = exprname_list.nops();
  // int nf = funcname_list.nops();
  size_t nd = Delays.size();
  
  int par_shift = 1;
  if (HasPeriod) par_shift = 0;

  std::vector<std::vector<size_t> > row(np*(nd+1));
  std::vector<std::vector<size_t> > col(np*(nd+1));

  for (size_t k = 0; k < nd + 1; ++k)
  {
    if (k == 0)
    {
      dout << "        if (k == 0)\n";
      dout << "        {\n";
      dout << "            // Derivatives wrt the state variables\n";
    }
    else
    {
      dout << "        }\n";
      dout << "        else if (k == " << k << ")\n";
      dout << "        {\n";
      dout << "            // Derivatives wrt state variables with delay " << Delays[k-1] << endl;
    }
    for (size_t m = 0; m < np; ++m)
    {
      if (m == 0)
      {
        dout << "            if (j == " << m + par_shift << ")\n";
        dout << "            {\n";
      }
      else
      {
        dout << "            else if (j == " << m + par_shift << ")\n";
        dout << "            {\n";
      }
      dout << "                // Derivative wrt " << parname_list[m] << "\n";
      for (size_t i = 0; i < nv; ++i)
      {
        ex f = vf0[i];
        for (size_t j = 0; j < nv; ++j)
        {
          symbol vtmp_("vtmp_");
          ex fj = f.subs(Zlags_(j, k) == vtmp_);
          ex df = fj.diff(vtmp_);
          df = df.subs(vtmp_ == Zlags_(j, k));
          symbol p = ex_to<symbol>(parname_list[m]);
          df = df.diff(p);
          // if (df != 0)
          if (!df.is_zero()) { row[k*np+m].push_back(i); col[k*np+m].push_back(j); }
          dout << "                jac_(" << i << "," << j << ",idx)" << " = " << df << ";" << endl;
        }
      }
      dout << "            }\n";
    }
  }
  dout << "        }\n";
  for (size_t i=0; i<row.size(); ++i)
  {
    dout << "        // Sparsity structure for tau = " << i/np << " and p = " << i%np + par_shift << "\n";
    for (size_t j=0; j<row[i].size(); ++j)
    {
      dout << "        // (" << row[i][j] << ", " << col[i][j] << ")\n";
    }
  }
}

ex pddec_second_deriv(const ex &f, size_t lag1, size_t var1, size_t lag2, size_t var2)
{
  symbol Z_("Z_");
  ex fj = f.subs(Zlags_(var1, lag1) == Z_);
  ex df = fj.diff(Z_);
  df = df.subs(Z_ == Zlags_(var1, lag1));
  df = df.subs(Zlags_(var2, lag2) == Z_);
  df = df.diff(Z_);
  df = df.subs(Z_ == Zlags_(var2, lag2));
  return df;
}

void VectorField::Knut_PrintHessiansTimesV(ostream &dout, const vector<ex> &vf0)
{
  // int nc = conname_list.nops();
  size_t nv = varname_list.nops();
  // int np = parname_list.nops();
  // int na = exprname_list.nops();
  // int nf = funcname_list.nops();
  size_t nd = Delays.size();
  
  std::vector<std::vector<size_t> > row((nd+1)*(nd+1));
  std::vector<std::vector<size_t> > col((nd+1)*(nd+1));
  
  for (size_t k1 = 0; k1 < nd + 1; ++k1)
  {
    if (k1 == 0)
    {
      dout << "        if (k1 == 0)\n";
      dout << "        {\n";
      dout << "            // Derivatives wrt the state variables\n";
    }
    else
    {
      dout << "        }\n";
      dout << "        else if (k1 == " << k1 << ")\n";
      dout << "        {\n";
      dout << "            // Derivatives wrt state variables with delay " << Delays[k1-1] << endl;
    }

    for (size_t k2 = 0; k2 < nd + 1; ++k2)
    {
      if (k2 == 0)
      {
        dout << "            if (k2 == 0)\n";
        dout << "            {\n";
        dout << "                // Derivatives wrt the state variables\n";
      }
      else
      {
        dout << "            else if (k2 == " << k2 << ")\n";
        dout << "            {\n";
        dout << "                // Derivatives wrt state variables with delay " << Delays[k2-1] << endl;
      }

      for (size_t i = 0; i < nv; ++i)
      {
        ex f = vf0[i];
        for (size_t j = 0; j < nv; ++j)
        {
          ostringstream os;
          os << csrc;
          for (size_t h = 0; h < nv; ++h)
          {
            ex d2f = pddec_second_deriv(f, k1, j, k2, h);
            if (d2f != 0)
            {
              if (os.str() != "")
                os << " + ";
              os << "(" << d2f << ")*v_(" << h << ",m,idx)";
            }
          }
          if (os.str() != "")
          {
            dout << "                jac_(" << i << "," << j << ",idx)" << " = " << os.str() << ";\n";
            row[k1*(nd+1)+k2].push_back(i); col[k1*(nd+1)+k2].push_back(j);
          } else
          {
            dout << "                jac_(" << i << "," << j << ",idx) = 0.0;\n";
          }
        }
      }
      dout << "            }\n";
    }
  }
  dout << "        }\n";
  for (size_t i=0; i<row.size(); ++i)
  {
    dout << "        // Sparsity structure for k1 = " << i/(nd+1) << " and k2 = " << i%(nd+1) << "\n";
    for (size_t j=0; j<row[i].size(); ++j)
    {
      dout << "        // (" << row[i][j] << ", " << col[i][j] << ")\n";
    }
  }

}

//
// PrintKnut -- The Knut code generator.
//

void VectorField::PrintKnut(ostream& sys_out, map<string, string> options)
{
  size_t nc = conname_list.nops();
  size_t nv = varname_list.nops();
  size_t np = parname_list.nops();
  size_t na = exprname_list.nops();
  size_t nf = funcname_list.nops();
  size_t par_shift = 1;
  if (HasPeriod) par_shift = 0;

  //
  //  Create the system definition file.
  //
  sys_out << csrc;
  sys_out << "//" << endl;
  sys_out << "// Knut KNSystem Definition file for the VFGEN vector field: " << Name() << endl;
  sys_out << "//" << endl;
  PrintVFGENComment(sys_out, "// ");
  sys_out << "//" << endl;
  sys_out << endl;
  // sys_out << "#include <cstdlib>\n";
  sys_out << "#include <cmath>\n";
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
  sys_out << "int sys_ndim()  { return " << nv << "; }  // KNSystem dimension\n";
  sys_out << "int sys_npar()  { return " << par_shift + np << "; }  // Number of parameters, plus one (for the period)\n";
  sys_out << "int sys_ntau()  { return " << 1 + Delays.size() << "; }  // Number of delays, plus one\n";
  sys_out << "int sys_nderi() { return 2; }  // Order of derivatives computed here\n";
  sys_out << "int sys_nevent() { return " << nf << "; }  // Number of event functions\n";
  sys_out << endl;
  
  bool hasMass = false;
  for (vector<StateVariable *>::iterator sv = StateVariables.begin(); sv != StateVariables.end(); ++sv)
  	if ((*sv)->Mass() != "") hasMass = true;
  if (hasMass)
  {
    sys_out << "void sys_mass(KNArray1D<double>& out)   // diagonal mass matrix\n";
    sys_out << "{\n";
    int k = 0;
  	for (vector<StateVariable *>::iterator sv = StateVariables.begin(); sv != StateVariables.end(); ++sv)
  	{
  	  if ((*sv)->Mass() != "") sys_out << "  out(" << k << ") = " << (*sv)->Mass() << ";\n";
  	  else sys_out << "  out(" << k << ") = 1.0;\n";
      ++k;
  	}
    sys_out << "}\n";
    sys_out << endl;
  }
  
  // ***************************************************************************
  //  Create the vector field function sys_rhs
  // ***************************************************************************
  sys_out << "//" << endl;
  sys_out << "// sys_p_rhs(...) computes the vector field.\n";
  sys_out << "//" << endl;
  sys_out << "void sys_p_rhs(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_, int sel)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  //
  // Parameters...
  //
  if (parname_list.nops() > 0)
  {
    sys_out << "    // Parameters (par_(0) is the period)\n";
    GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
    sys_out << endl;
  }
  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";

  // The expressions are already substituted into the equations
  //
  // StateVariables...
  //
  // sys_out << "    vf_ = zeros(" << nv << ",1);" << endl;
  sys_out << "        // Compute the vector field\n";
  for (size_t i = 0; i < nv; ++i)
  {
    ex f = varvecfield_list[i];
    Knut_ConvertStateToZlags(f);
//    Knut_ConvertConstsPars(f);
    sys_out << "        out(" << i << ",idx)" << " = " << f << ";" << endl;
  }
  sys_out << "    }\n";
  sys_out << "}\n";

  // ***************************************************************************
  //  Derivative with respect to the delays
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// " << Name() << "_jacx\n";
  sys_out << "//\n";
  sys_out << "// This function computes the Jacobian of the vector field\n";
  sys_out << "// with respect to the k-th delayed vector.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << Name() << "_jacx(KNArray3D<double>& jac_, const KNArray1D<double>& time, int k, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  //
  // Parameters...
  //
  if (parname_list.nops() > 0)
  {
    sys_out << "    // Parameters (par_(0) is the period)\n";
    GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
    // sys_out << endl;
  }

  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";
  //
  // After the following loop, vf0 will hold a version of the vector field
  // in which all occurrences of delay(expr,lag) will have been converted
  // to use Zlags_(i,j).  The i^th non-delayed state variable will be converted
  // to Zlags_(i,1).
  //
  vector<ex> vf0;
  for (size_t i = 0; i < nv; ++i)
  {
    ex f = varvecfield_list[i]; 
    // Convert the state variables and delay expressions to Zlags_
    Knut_ConvertStateToZlags(f);
    vf0.push_back(f);
  }
  Knut_PrintJacobians(sys_out, vf0);
  sys_out << "    }\n";
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  //  Derivative with respect to the parameters
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// " << Name() << "_jacp\n";
  sys_out << "//\n";
  sys_out << "// This function computes the Jacobian of the vector field\n";
  sys_out << "// with respect to the j-th parameter.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << Name() << "_jacp(KNArray3D<double>& jac_, const KNArray1D<double>& time, int j, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  //
  // Parameters...
  //
  if (parname_list.nops() > 0)
  {
    sys_out << "    // Parameters (par_(0) is the period)\n";
    GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
    // sys_out << endl;
  }
  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";
  Knut_PrintParDerivs(sys_out, vf0);
  sys_out << "    }\n";
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  //  Derivative with respect to the delays and parameters (2nd order)
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// " << Name() << "_jacxp\n";
  sys_out << "//\n";
  sys_out << "// This function computes the Jacobian of the vector field\n";
  sys_out << "// with respect to the k-th delayed state vector and the j-th parameter.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << Name() << "_jacxp(KNArray3D<double>& jac_, const KNArray1D<double>& time, int k, int j, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  // sys_out << endl;
  //
  // Parameters...
  //
  if (parname_list.nops() > 0)
  {
    sys_out << "    // Parameters (par_(0) is the period)\n";
    GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
    // sys_out << endl;
  }
  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";
  Knut_PrintXandParJacobians(sys_out, vf0);
  sys_out << "    }\n";
  sys_out << "}\n";
  sys_out << endl;

  // ***************************************************************************
  //  Second order derivative with respect to the delays 
  // ***************************************************************************
  sys_out << endl;
  sys_out << "//\n";
  sys_out << "// " << Name() << "_hess_times_v\n";
  sys_out << "//\n";
  sys_out << "// This function computes the Hessian of the vector field\n";
  sys_out << "// with respect to the k1-th and k2-th delayed state vectors,\n";
  sys_out << "// then multiplies this by the m-th column of v to obtain a matrix.\n";
  sys_out << "//\n";
  sys_out << endl;
  sys_out << "static inline void " << Name() << "_hess_times_v(KNArray3D<double>& jac_, const KNArray1D<double>& time, int k1, int k2, int m, const KNArray3D<double>& v_, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  sys_out << endl;

  //
  // Parameters...
  //
  if (parname_list.nops() > 0)
  {
    sys_out << "    // Parameters\n";
    GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
    sys_out << endl;
  }
  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";
  Knut_PrintHessiansTimesV(sys_out, vf0);
  sys_out << "    }\n";
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
  sys_out << "//" << endl;
  sys_out << "void sys_p_deri(KNArray3D<double>& jac_, const KNArray1D<double>& time, const KNArray3D<double>& Zlags_, const KNArray1D<double>& par_, int sel,\nint nx_, const int* vx_, int np_, const int* vp_, const KNArray3D<double>& v_)\n";
  sys_out << "{\n";
  sys_out << "    if (nx_ == 1 && np_ == 0)\n";
  sys_out << "        " << Name() << "_jacx(jac_, time, vx_[0], Zlags_, par_);\n";
  sys_out << "    else if (nx_ == 0 && np_ == 1)\n";
  sys_out << "        " << Name() << "_jacp(jac_, time, vp_[0], Zlags_, par_);\n";
  sys_out << "    else if (nx_ == 1 && np_ == 1)\n";
  sys_out << "        " << Name() << "_jacxp(jac_,time, vx_[0], vp_[0], Zlags_, par_);\n";
  sys_out << "    else if (nx_ == 2 && np_ == 0)\n";
  sys_out << "        " << Name() << "_hess_times_v(jac_,time,vx_[0],vx_[1],vx_[0],v_,Zlags_,par_);\n";
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
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  sys_out << "    // par_(0) is the period.\n";
  GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
  sys_out << endl;
  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";
  sys_out << "        out(0,idx) = 0.0;\n";
  int j = 1;
  for (vector<ex>::iterator p = Delays.begin(); p != Delays.end(); ++p)
  {
    // TODO: check that delays don't depend on states
    sys_out << "        out(" << j << ",idx) = " <<  *p << ";\n";
    ++j;
  }
  sys_out << "    }\n";
  sys_out << "}\n";
  sys_out << endl;
  // ***************************************************************************
  // Create sys_dtau()
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
  sys_out << "void sys_p_dtau(KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par_, int p_)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }

  sys_out << "    // par_(0) is the period.\n";
  GetFromVector(sys_out, "    const double ", parname_list, "par_", "()", par_shift, ";");
  sys_out << endl;
  sys_out << "    for (int idx=0; idx < time.size(); ++idx)\n";
  sys_out << "    {\n";
  sys_out << "        const double t = time(idx);\n";
  sys_out << "        out(0,idx) = 0.0;\n";
  sys_out << "        if (p_ == -1)\n";
  sys_out << "        {\n";
  for (unsigned k = 0; k < Delays.size(); ++k)
  {
    sys_out << "            out(" << k + 1 << ",idx) = 0.0;\n";
  }
  sys_out << "        }\n";
  for (size_t j = 0; j < np; ++j)
  {
    sys_out << "        else if (p_ == " << j + par_shift << ")\n";
    sys_out << "        {\n";
    sys_out << "            // Derivative wrt " << parname_list[j] << endl;
    for (unsigned k = 0; k < Delays.size(); ++k)
    {
      symbol p = ex_to<symbol>(parname_list[j]);
      ex df = Delays[k].diff(p);
      sys_out << "            out(" << k + 1 << ",idx) = " << df << ";\n";
    }
    sys_out << "        }\n";
  }

  sys_out << "    }\n";
  sys_out << "}\n";
  sys_out << endl;

  // ***************************************************************************
  // Starting parameters
  // ***************************************************************************
  sys_out << "void sys_stpar(KNVector& par_)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  sys_out << "    // VFGEN used the DefaultValues of the Parameters.\n";
  sys_out << "    // Change the following values to match your known solution.\n";
  if (!HasPeriod) sys_out << "    par_(0) = 1.0;\n";
  for (size_t j = 0; j < np; ++j)
  {
    // First make a replacement list with all the other parameters
    GiNaC::lst parsubs;
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
    defval = iterated_subs(defval, parsubs);
    sys_out << "    par_(" << j + par_shift << ") = " << defval << ";\n";
  }
  sys_out << "}\n";
  sys_out << endl;
  
  // ***************************************************************************
  // Starting solution
  // ***************************************************************************
  sys_out << "void sys_stsol(KNVector& out, double t)\n";
  sys_out << "{\n";
  if (HasPi)
  {
    sys_out << "    const double Pi = M_PI;\n";
  }
  //
  // Constants...
  //
  for (size_t i = 0; i < nc; ++i)
  {
    sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
  }
  //
  // Parameters...
  //
  for (size_t j = 0; j < np; ++j)
  {
    // First make a replacement list with all the other parameters
    GiNaC::lst parsubs;
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
    defval = iterated_subs(defval, parsubs);
    sys_out << "    const double " << parname_list[j] << " = " << defval << ";\n";
  }
  
  sys_out << "    // VFGEN used the DefaultInitialConditions of the StateVariables.\n";
  sys_out << "    // Change the following values to implement your known solution.\n";
  for (size_t j = 0; j < nv; ++j)
  {
    sys_out << "    out(" << j << ") = " << vardefic_list[j] << ";\n";
  }
  sys_out << "}\n";
  sys_out << endl;

  // ***************************************************************************
  // Parameter names
  // ***************************************************************************
  sys_out << "void sys_parnames( const char *out[] )\n";
  sys_out << "{\n";
  if (!HasPeriod) sys_out << "    out[0] = \"Period\";\n";
  for (unsigned int k = 0; k < parname_list.nops(); ++k)
  {
    sys_out << "    out[" << k + par_shift << "] = \"" << parname_list[k] << "\";\n";
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
  int par_shift = 1;
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
  int par_shift = 1;
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
    for (size_t q = 0; q < np; ++q)
    {
      symbol ps = ex_to<symbol>(parname_list[q]);
      ex dfp = f.diff(ps);
      Knut_ConvertStateToZlags(dfp);
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

void VectorField::Knut_stpar(std::vector<double>& par)
{  
  const size_t nc = conname_list.nops();
  const size_t np = parname_list.nops();
  int par_shift = 1;
  if (HasPeriod) par_shift = 0;
  else par[0] = 1.0;
  par.resize(np+par_shift);
  for (size_t j = 0; j < np; ++j)
  {
    // First make a replacement list with all the other parameters
    GiNaC::lst parsubs;
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
    defval = iterated_subs(defval, parsubs);
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
  int par_shift = 1;
  if (HasPeriod) par_shift = 0;
  pnames.resize(np+par_shift);
  for (size_t k = 0; k < np; ++k)
  {
    if (is_exactly_a<symbol>(parname_list[k])) pnames[k+par_shift] = ex_to<symbol>(parname_list[k]).get_name();
    else P_MESSAGE1("Parameter is not a symbol.");
  }
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
      int pid = ex_to<numeric>(expr.op(0)).to_int();
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = par(pid);
      return;
    }
    if (name == "Zlags_")
    {
      int d1 = ex_to<numeric>(expr.op(0)).to_int();
      int d2 = ex_to<numeric>(expr.op(1)).to_int();
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = x(d1,d2,k);
      return;
    }
    if (name == "VZlags_")
    {
      int d1 = ex_to<numeric>(expr.op(0)).to_int();
      int d2 = ex_to<numeric>(expr.op(1)).to_int();
      for (size_t k = 0; k < vsize; ++k) res[k*skip] = v(d1,d2,k);
      return;
    }
    double (*fun_ptr)(double) = 0;
    if (name == "abs") fun_ptr = &fabs;
    if (name == "step") return;
    if (name == "csgn") return;
    if (name == "eta") return;
    if (name == "Li2") return;
    if (name == "Li3") return;
    if (name == "zetaderiv") return;
    if (name == "factorial") return;
    if (name == "binomial") return;
    if (name == "Order") return;
    if (name == "exp") fun_ptr = &exp;
    if (name == "log") fun_ptr = &log;
    if (name == "sin") fun_ptr = &sin;
    if (name == "cos") fun_ptr = &cos;
    if (name == "tan") fun_ptr = &tan;
    if (name == "asin") fun_ptr = &asin;
    if (name == "acos") fun_ptr = &acos;
    if (name == "atan") fun_ptr = &atan;
    if (name == "atan2") return;
    if (name == "sinh") fun_ptr = &sinh;
    if (name == "cosh") fun_ptr = &cosh;
    if (name == "tanh") fun_ptr = &tanh;
    if (name == "asinh") fun_ptr = &asinh;
    if (name == "acosh") fun_ptr = &acosh;
    if (name == "atanh") fun_ptr = &atanh;
    if (name == "lgamma") fun_ptr = &lgamma;
    if (name == "tgamma") fun_ptr = &tgamma;
    if (name == "beta") return;
    if (name == "psi") return;
    if (name == "G") return;
    if (name == "Li") return;
    if (name == "S") return;
    if (name == "H") return;
    if (name == "zeta") return;
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
  } else
  {
    P_MESSAGE1("Unknown formula component.");
  }
}

void VectorField::Knut_tau_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, const size_t nv, const size_t nps, const size_t nd )
{
  KNArray3D<double> dummy_x, dummy_v;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
  for (size_t i = 0; i < nd; ++i)
  {
    expeval(exls[i], out.pointer(i,0), skip, time, dummy_x, dummy_v, par);
  }
}

void VectorField::Knut_tau_p_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp, const size_t nv, const size_t nps, const size_t nd )
{
  KNArray3D<double> dummy_x, dummy_v;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
  for (size_t i = 0; i < nd; ++i)
  {
    expeval(exls[i + vp*nd], out.pointer(i,0), skip, time, dummy_x, dummy_v, par);
  }
}


void VectorField::Knut_RHS_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel )
{
  KNArray3D<double> dummy_v;
  const size_t skip = (((size_t)out.pointer(0,1)) - ((size_t)out.pointer(0,0)))/sizeof(double);
//  std::cout << "Knut_RHS_eval " << exls[0] << " " << skip << " \n";
  for (size_t n = 0; n < exls.size(); ++n)
  {
    expeval(exls[n], out.pointer(n,0), skip, time, x, dummy_v, par);
  }
}

void VectorField::Knut_RHS_p_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel, int alpha, const size_t nv, const size_t nps, const size_t nd)
{
  KNArray3D<double> dummy_v;
  const size_t skip = (((size_t)out.pointer(0,0,1)) - ((size_t)out.pointer(0,0,0)))/sizeof(double);
  for (size_t i = 0; i < nv; ++i)
  {
    expeval(exls[i + alpha*nv], out.pointer(i,0,0), skip, time, x, dummy_v, par);
  }	
}

void VectorField::Knut_RHS_x_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel, int del, const size_t nv, const size_t nps, const size_t nd)
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

void VectorField::Knut_RHS_xp_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel, int del, int alpha, const size_t nv, const size_t nps, const size_t nd)
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

void VectorField::Knut_RHS_xx_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray3D<double>& vv, const KNVector& par, int sel, int del1, int del2, const size_t nv, const size_t nps, const size_t nd)
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
