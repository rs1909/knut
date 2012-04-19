
//
// vf.cpp
//
//
//  Copyright (C) 2008 Warren Weckesser
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

#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <ginac/ginac.h>

#include "vf.h"
#include "codegen_utils.h"
#include "knerror.h"

using namespace std;
using namespace GiNaC;

//
// This function overrides the default print method for the Zlags_ function.
// Ginac thinks Zlags_ is a function, but in fact, in the code generated in
// a couple of the VFGEN commands that handle delay equations, Zlags_ is a
// two-dimensional array.  So we want Zlags_(1,2) printed like that, not as
// Zlags_(1.0,2.0).
//

static void Zlags_print(const ex& arg1, const ex& arg2, const print_context& c)
{
  c.s << "Zlags_(";
  if (is_a<numeric>(arg1))
    c.s << ex_to<numeric>(arg1).to_int();
  else
    arg1.print(c);
  c.s << ",";
  if (is_a<numeric>(arg2))
    c.s << ex_to<numeric>(arg2).to_int();
  else
    arg2.print(c);
  c.s << ",idx)";
}

static void VZlags_print(const ex& arg1, const ex& arg2, const print_context& c)
{
  c.s << "VZlags_(";
  if (is_a<numeric>(arg1))
    c.s << ex_to<numeric>(arg1).to_int();
  else
    arg1.print(c);
  c.s << ",";
  if (is_a<numeric>(arg2))
    c.s << ex_to<numeric>(arg2).to_int();
  else
    arg2.print(c);
  c.s << ",idx)";
}

static void par_print(const ex& arg1, const print_context& c)
{
  c.s << "par_(";
  if (is_a<numeric>(arg1))
    c.s << ex_to<numeric>(arg1).to_int();
  else
    arg1.print(c);
  c.s << ")";
}

//
// Add a function called "delay" to the functions known by ginac.
// (The corresponding macro DECLARE_FUNCTION_2P(delay) is in vf.h.)
//

REGISTER_FUNCTION(delay, dummy())
REGISTER_FUNCTION(Zlags_, print_func<print_csrc_float>(Zlags_print).
                  print_func<print_csrc_double>(Zlags_print).
                  print_func<print_python>(Zlags_print))
REGISTER_FUNCTION(VZlags_, print_func<print_csrc_float>(VZlags_print).
                  print_func<print_csrc_double>(VZlags_print).
                  print_func<print_python>(VZlags_print))
REGISTER_FUNCTION(par_, print_func<print_csrc_float>(par_print).
                  print_func<print_csrc_double>(par_print).
                  print_func<print_python>(par_print))                  

ex delay_reader(const exvector& ev)
{
     return delay(ev[0],ev[1]);
}

ex Zlags__reader(const exvector& ev)
{
     return Zlags_(ev[0],ev[1]);
}

ex VZlags__reader(const exvector& ev)
{
     return VZlags_(ev[0],ev[1]);
}

ex par__reader(const exvector& ev)
{
     return par_(ev[0]);
}

static ex heaviside_evalf(const ex& x)
{
  if (is_exactly_a<numeric>(x))
  {
    if (ex_to<numeric>(x) > 0) return 1;
    else return 0;
  }
  else
    return heaviside(x).hold();
}

static ex heaviside_eval(const ex& x)
{
  return heaviside(x).hold();
}

static ex heaviside_deriv(const ex& /*x*/, unsigned /*deriv_param*/)
{
  return 0;
}

REGISTER_FUNCTION(heaviside, eval_func(heaviside_eval).
                  evalf_func(heaviside_evalf).
                  derivative_func(heaviside_deriv))

ex heaviside_reader(const exvector& ev)
{
     return heaviside(ev[0]);
}

static ex ramp_evalf(const ex& x)
{
  if (is_exactly_a<numeric>(x))
  {
    if (ex_to<numeric>(x) > 0) return x;
    else return 0;
  }
  else
    return ramp(x).hold();
}

static ex ramp_eval(const ex& x)
{
  return ramp(x).hold();
}

static ex ramp_deriv(const ex& x, unsigned /*deriv_param*/)
{
  return heaviside(x);
}

REGISTER_FUNCTION(ramp, eval_func(ramp_eval).
                  evalf_func(ramp_evalf).
                  derivative_func(ramp_deriv))

ex ramp_reader(const exvector& ev)
{
     return ramp(ev[0]);
}

//
// Symbol Methods
//

Symbol::Symbol() { }
Symbol::Symbol(const std::string& n) : name(n) { }
Symbol::Symbol(const std::string& n, const std::string& descr) : name(n), description(descr) { }

void   Symbol::Name(const std::string& n) { name = n; }
const std::string& Symbol::Name(void) const { return name; }
void   Symbol::Description(const std::string& descr) { description = descr; }
const std::string& Symbol::Description(void) const { return description; }
void   Symbol::Latex(const std::string& l) { latex = l; }
const std::string& Symbol::Latex(void) const { return latex; }

//
// FormulaSymbol Methods
//

FormulaSymbol::FormulaSymbol() : Symbol() { }
FormulaSymbol::FormulaSymbol(const std::string& name) : Symbol(name) { }
FormulaSymbol::FormulaSymbol(const std::string& name, const std::string& descr) : Symbol(name, descr) { }

void   FormulaSymbol::Formula(const std::string& f) { formula = f; }
const std::string& FormulaSymbol::Formula(void) const { return formula; }

//
// Constant Methods
//

Constant::Constant(const std::string& name) : Symbol(name), value("") { }
Constant::Constant(const std::string& name, const std::string& descr) : Symbol(name, descr), value("") { }

void   Constant::Value(const std::string& val) { value = val; }
const std::string& Constant::Value(void) const { return value; }


//
// Parameter Methods
//

Parameter::Parameter(const std::string& name) : Symbol(name), defaultvalue("") { }
Parameter::Parameter(const std::string& name, const std::string& descr) : Symbol(name, descr), defaultvalue("") { }

void   Parameter::DefaultValue(const std::string& val) { defaultvalue = val; }
const std::string& Parameter::DefaultValue(void) const { return defaultvalue; }

//
// Expression Methods
//
Expression::Expression(const std::string& name) : FormulaSymbol(name) { }
Expression::Expression(const std::string& name, const std::string& descr) : FormulaSymbol(name, descr) { }

//
// StateVariable Methods
//
StateVariable::StateVariable(const std::string& name) : FormulaSymbol(name), periodicfrom(""), periodicto(""), default_ic("") { }
StateVariable::StateVariable(const std::string& name, const std::string& descr) : FormulaSymbol(name, descr), periodicfrom(""), periodicto(""), default_ic("") { }

void   StateVariable::PeriodicFrom(const std::string& pfrom) { periodicfrom = pfrom; }
const std::string& StateVariable::PeriodicFrom(void) const { return periodicfrom; }
void   StateVariable::PeriodicTo(const std::string& pto) { periodicto = pto; }
const std::string& StateVariable::PeriodicTo(void) const { return periodicto; }
bool   StateVariable::IsPeriodic(void) { return periodicfrom != ""; }
void   StateVariable::DefaultInitialCondition(const std::string& ic) { default_ic = ic; }
const std::string& StateVariable::DefaultInitialCondition(void) const { return default_ic; }
void   StateVariable::DefaultHistory(const std::string& hist) { default_history = hist; }
const std::string& StateVariable::DefaultHistory() const { return default_history; }
void   StateVariable::Mass(const std::string& ms) { mass = ms; }
const std::string& StateVariable::Mass() const { return mass; }

//
// Function Methods
//
Function::Function(const std::string& name) : FormulaSymbol(name) { }
Function::Function(const std::string& name, const std::string& descr) : FormulaSymbol(name, descr) { }

//
// VectorField Methods
//

VectorField::VectorField(void) : IndependentVariable("t"), TimeDependentDelay(false), StateDependentDelay(false), IsAutonomous(true) { }
VectorField::VectorField(const std::string& name, const std::string& descr) : Symbol(name, descr), IndependentVariable("t"), TimeDependentDelay(false), StateDependentDelay(false), IsAutonomous(true) { }
VectorField::VectorField(const std::string& name, const std::string& descr, const std::string& /*indvar*/) : Symbol(name, descr), IndependentVariable("t"), TimeDependentDelay(false), StateDependentDelay(false), IsAutonomous(true) { }

VectorField::~VectorField()
{
  for (size_t i=0; i<Constants.size(); ++i)
  {
    delete Constants[i];
  }
  for (size_t i=0; i<Parameters.size(); ++i)
  {
    delete Parameters[i];
  }
  for (size_t i=0; i<Expressions.size(); ++i)
  {
    delete Expressions[i];
  }
  for (size_t i=0; i<StateVariables.size(); ++i)
  {
    delete StateVariables[i];
  }
  for (size_t i=0; i<Functions.size(); ++i)
  {
    delete Functions[i];
  }
}

void VectorField::AddConstant(Constant *p) { Constants.push_back(p); }
void VectorField::AddParameter(Parameter *p) { Parameters.push_back(p); }
void VectorField::AddExpression(Expression *e) { Expressions.push_back(e); }
void VectorField::AddStateVariable(StateVariable *sv) { StateVariables.push_back(sv); }

int VectorField::FindVar(const symbol &var)
{
  bool found = false;
  size_t k;
  for (k = 0; k < varname_list.nops(); ++k)
  {
    if (varname_list[k] == var)
    {
      // cerr << "Found the variable " << var << " at varname_list[" << k << "]\n";
      found = true;
      break;
    }
  }
  int vindex;
  if (found == false)
  {
    vindex = -1;
  }
  else
  {
    vindex = static_cast<int>(k);
  }
  return vindex;
}

GiNaC::lst VectorField::FindVarsInEx(const GiNaC::ex &e)
{
  GiNaC::lst vlist;
  for (size_t k = 0; k < varname_list.nops(); ++k)
  {
    if (e.has(varname_list[k])) vlist.append(varname_list[k]);
  }
  return vlist;
}

void VectorField::AddFunction(Function *f)
{
  Functions.push_back(f);
}

ex VectorField::SubsAllExpressions(const ex& e)
{
  const size_t na = exprname_list.nops();
  ex s = e;
  for (size_t k = 0; k < na; ++k)
  {
    s = s.subs(ex_to<symbol>(exprname_list[k]) == exprformula_list[na-k-1]);
  }
  return s;
}

int VectorField::FindDelay(ex &del)
{
  bool found = false;
  size_t k;
  for (k = 0; k < Delays.size(); ++k)
  {
//     cerr << "FindDelay: Delays[" << k << "] = " << Delays[k] << " " << del << endl;
    if (Delays[k] == del)
    {
//       cerr << "FindDelay: Found the delay " << del << " at Delays[" << k << "]\n";
      found = true;
      break;
    }
  }
  int dindex;
  if (found == false)
  {
    dindex = -1;
  }
  else
  {
    dindex = static_cast<int>(k);
  }
  return dindex;
}


int VectorField::AddDelay(ex &del)
{
  if (FindDelay(del) != -1)
    return 1;
  Delays.push_back(del);
  return 0;
}

//
// CheckForDelay(const ex& f)
//
// Look for expressions of the form delay(delayexpr,del) in f.
// For each of the these, call AddDelay.
// Also set the flags IsDelay and HasNonconstantDelay as appropriate.
//
// This function is called by ProcessSymbols for each Expression and
// StateVariable Formula.
//

void VectorField::CheckForDelay(const ex& f)
{
  exset occurrences;
  if (f.find(delay(wild(1), wild(2)), occurrences))
  {
    IsDelay = true;
    for (exset::const_iterator iter = occurrences.begin(); iter != occurrences.end(); ++iter)
    {
      ex del = iter->op(1);
      AddDelay(del);
      if (del.has(IndVar))
        TimeDependentDelay = true;  // time-dependent delay
      for (lst::const_iterator viter = varname_list.begin(); viter != varname_list.end(); ++viter)
      {
        if (del.has(*viter))
          StateDependentDelay = true; // state-dependent delay
      }
    }
  }
}

struct tup
{
  ex f;
  std::string eps;
};

static bool comp(const tup& a, const tup& b) { return a.eps < b.eps; }

int VectorField::ProcessSymbols(void)
{
  int rval = 0;

  IsDelay = false;
  HasPi   = false;
  HasPeriod = false;
  
  // Process the constants NAME
  for (vector<Constant *>::iterator c = Constants.begin(); c != Constants.end(); ++c)
  {
    symbol con(!(*c)->Latex().empty()
               ? symbol((*c)->Name(), (*c)->Latex())
               : symbol((*c)->Name()));
    conname_list.append(con);
    allsymbols.append(con);
  }

    // Process the constants VALUE
    for (vector<Constant *>::iterator c = Constants.begin(); c != Constants.end(); ++c)
    {
      string val = (*c)->Value();
      try
      {
        ex e(val,allsymbols);
        convalue_list.append(e);
        if (has(e, Pi)) HasPi = true;
      }
      catch (exception &p)
      {
        // cerr << "VectorField:ProcessSymbols: exception while processing constants\n";
        P_MESSAGE6("The Value \"", (*c)->Value(), "\" for the Constant ", (*c)->Name(), " has an error: ", p.what());
        rval = -1;
      }
    }

  // Process the parameters NAME
  for (vector<Parameter *>::iterator p = Parameters.begin(); p != Parameters.end(); ++p)
  {
    symbol par(!(*p)->Latex().empty()
               ? symbol((*p)->Name(), (*p)->Latex())
               : symbol((*p)->Name()));
    // symbol par((*p)->Name());
    bool isPeriod = !strcasecmp((*p)->Description().c_str(), "period");
    if (isPeriod) HasPeriod = true;
    if (isPeriod) parname_list.prepend(par); else parname_list.append(par);
    allsymbols.append(par);
  }
  // Process the parameters DEFAULTVALUE
  for (vector<Parameter *>::iterator p = Parameters.begin(); p != Parameters.end(); ++p)
  {
    bool isPeriod = !strcasecmp((*p)->Description().c_str(), "period");
    string defval = (*p)->DefaultValue();
    if (defval == "") defval = "0";
    try
    {
      ex e(defval, allsymbols);
      if (isPeriod) pardefval_list.prepend(e); else pardefval_list.append(e); 
      if (has(e, Pi)) HasPi = true;
    }
    catch (exception &ep)
    {
      P_MESSAGE6("The DefaultValue \"", (*p)->DefaultValue(), "\" for the Parameter ", (*p)->Name(), " has an error: ", ep.what());
      rval = -1;
    }
  }

  // At this point, allsymbols is a list of the ginac
  // symbols of the constants and parameters.
  // Now add IndVar to the list, because the expressions
  // for the default initial conditions, vector field
  // formulas, expressions, and functions can all be
  // functions of the independent variable.

  IndVar = symbol(IndependentVariable);
  allsymbols.append(IndVar);

  // Process the state variable names and default ICs
  //   (but not the formulas for the vector field)
  // NAME  __AND__ Initial condition
  size_t nv = 0;
  for (vector<StateVariable *>::iterator sv = StateVariables.begin(); sv != StateVariables.end(); ++sv, ++nv)
  {
    symbol var(!(*sv)->Latex().empty()
               ? symbol((*sv)->Name(), (*sv)->Latex())
               : symbol((*sv)->Name()));
    varname_list.append(var);
    string defic = (*sv)->DefaultInitialCondition();
    if (defic == "") defic = "0";
    try
    {
      ex e(defic, allsymbols);
      vardefic_list.append(e);
      if (has(e, Pi)) HasPi = true;
    }
    catch (exception &p)
    {
      P_MESSAGE6("The DefaultInitialCondition \"", (*sv)->DefaultInitialCondition(), "\" for the StateVariable ", (*sv)->Name(), " has an error: ", p.what());
      rval = -1;
    }
  }
  // History
  for (vector<StateVariable *>::iterator sv = StateVariables.begin(); sv != StateVariables.end(); ++sv)
  {
    string defhist = (*sv)->DefaultHistory();
    if (defhist == "")
      defhist = "0";  // The default DefaultHistory.  Shouldn't this be DefaultInitialCondition?
    // (And only be 0 if DefaultInitialCondition was also not given.)
    try
    {
      ex e(defhist, allsymbols);
      vardefhist_list.append(e);
      if (has(e, Pi))
        HasPi = true;
    }
    catch (exception &p)
    {
      P_MESSAGE6("The DefaultHistory \"", (*sv)->DefaultHistory(), "\" for the StateVariable ", (*sv)->Name(), " has an error: ", p.what());
      rval = -1;
    }
  }
  
  //
  // Now add the state variable symbols to allsymbols.
  // We didn't do this in the above loop because we don't want
  // to allow the default initial condition of one state variable
  // to be a function of another state variable. (In other words,
  // the default initial conditions can be functions of only the
  // constants, parameters, and the independent variable [for
  // delay equations].)
  //
  for (size_t j = 0; j < nv; ++j)
  {
    allsymbols.append(varname_list.op(j));
  }

  // Process the expressions NAME and FORMULA
  // Check for the delays later and replace the delay with the full expression at the end
  for (vector<Expression *>::iterator e = Expressions.begin(); e != Expressions.end(); ++e)
  {
    symbol auxe(!(*e)->Latex().empty()
                ? symbol((*e)->Name(), (*e)->Latex())
                : symbol((*e)->Name()));
    exprname_list.append(auxe);
    allsymbols.append(auxe);
    try
    {
      ex f((*e)->Formula(), allsymbols);
      exprformula_list.append(f);
      expreqn_list.append(auxe == f);
      if (has(f, Pi)) HasPi = true;
    }
    catch (exception &p)
    {
      P_MESSAGE6("The Formula \"", (*e)->Formula(), "\" for the Expression ", (*e)->Name(), " has an error: ", p.what());
      rval = -1;
    }
  }
  
  //
  // At this point, allsymbols holds the ginac symbols for
  // the independent variable, the constants, the parameters,
  // the state variables, and the auxiliary expressions.
  // No more symbols will be added to allsymbols.
  //

  // Process the vector field formulas
  for (vector<StateVariable *>::iterator sv = StateVariables.begin(); sv != StateVariables.end(); ++sv)
  {
    try
    {
      ex f((*sv)->Formula(), allsymbols);
      varvecfield_list.append(f);
      if (has(f, Pi)) HasPi = true;
      if (has(f, IndVar)) IsAutonomous = false;
    }
    catch (exception &p)
    {
      P_MESSAGE6("The Formula \"", (*sv)->Formula(), "\" for the StateVariable ", (*sv)->Name(), " has an error: ", p.what());
      rval = -1;
    }
  }

  // checking for delays with everything substituted in the vectorfields
  for (size_t i = 0; i < varvecfield_list.nops(); ++i)
  {
    ex f = iterated_subs(varvecfield_list[i], expreqn_list);
    varvecfield_list[i] = f;
    CheckForDelay(f);
  }
  
  // alphabetically ordering the delays. 
  // It is not necessary for functionality, was used to debug a wierd problem
  tup etup[Delays.size()];
  std::ostringstream os;
  for (size_t i = 0; i < Delays.size(); ++i)
  {
    os.str("");
    os << Delays[i];
    etup[i].eps = os.str();
    etup[i].f = Delays[i];
  }
  std::sort(etup, etup + Delays.size(), comp);
    
  // Functions
  for (vector<Function *>::iterator f = Functions.begin(); f != Functions.end(); ++f)
  {
    symbol funcname((*f)->Name());
    funcname_list.append(funcname);
    try
    {
      ex funcexpr((*f)->Formula(), allsymbols);
      funcformula_list.append(funcexpr);
      if (has(funcexpr, Pi)) HasPi = true;
    }
    catch (exception &p)
    {
      P_MESSAGE6("The Formula \"", (*f)->Formula(), "\" for the  Function ", (*f)->Name(), " has an error: ", p.what());
      rval = -1;
    }
  }
  return rval;
}

void VectorField::print(void)
{
  cout << "Name:       " << Name();
  cout << "   Independent Variable: " << IndependentVariable << endl;
  cout << "Constants:  " << conname_list;
  cout << "  Values: " << convalue_list << endl;

  cout << "Parameters: " << parname_list;
  cout << "  Default values: " << pardefval_list << endl;

  cout << "Variables:  " << varname_list << endl;
  cout << "  DefaultICs:     " << vardefic_list << endl;
  cout << "  DefaultHistory: " << vardefhist_list << endl;

  cout << "Expressions: " << endl;

  for (unsigned i = 0; i < exprname_list.nops(); ++i)
  {
    cout << "   " << exprname_list[i] << "=" << exprformula_list[i] << endl;
  }

  cout << "Expressions (equation list): " ;
  cout << expreqn_list << endl;

  cout << "KNVector field: " << endl;
  for (unsigned i = 0; i < varvecfield_list.nops(); ++i)
  {
    cout << "   " << varvecfield_list[i] << endl;
  }
  cout << "Functions: " << endl;
  for (unsigned i = 0; i < funcname_list.nops(); ++i)
  {
    cout << "   " << funcname_list[i] << "=" << funcformula_list[i] << endl;
  }
  string tf[2] = {"false", "true"};
  cout << "IsDelay: " << tf[IsDelay] << endl;
  cout << "TimeDependentDelay: " << tf[TimeDependentDelay] << endl;
  cout << "StateDependentDelay: " << tf[StateDependentDelay] << endl;
  cout << "IsAutonomous: " << tf[IsAutonomous] << endl;
  cout << "HasPi: " << tf[HasPi] << endl;
}
