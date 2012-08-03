
//
// vf.h
//
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


#ifndef VF_H_INCLUDED_
#define VF_H_INCLUDED_

#include <string>
#include <vector>
#include <map>
#include <ginac/ginac.h>

#include <matrix.h>

DECLARE_FUNCTION_2P(delay)
DECLARE_FUNCTION_2P(Zlags_)
DECLARE_FUNCTION_2P(VZlags_)
DECLARE_FUNCTION_1P(par_)

DECLARE_FUNCTION_1P(heaviside)
DECLARE_FUNCTION_1P(ramp)

//
// Symbol
//     FormulaSymbol
//         Expression
//         StateVariable
//         Function
//     Constant
//     Parameter
//     VectorField
//


//
// Symbol is the base class for the objects used to define a vector field.
// The data associated with a Symbol is the name of the symbol, a
// description, and a latex representation.
//

class Symbol
{
  private:

    std::string name;
    std::string description;
    std::string latex;

  public:

    // Constructors
    Symbol();
    Symbol(const std::string& name);
    Symbol(const std::string& name, const std::string& descr);

    // Get/Set methods
    void Name(const std::string& n);
    const std::string& Name(void) const;
    void Description(const std::string& descr);
    const std::string& Description(void) const;
    void Latex(const std::string& l);
    const std::string& Latex(void) const;
};


class FormulaSymbol : public Symbol
{
  private:

    std::string formula;

  public:

    // Constructors
    FormulaSymbol();
    FormulaSymbol(const std::string& name);
    FormulaSymbol(const std::string& name, const std::string& descr);

    // Get/Set methods
    void Formula(const std::string& f);
    const std::string& Formula(void) const;
};

class Constant : public Symbol
{
  private:

    std::string value;

  public:

    // Constructors
    Constant(const std::string& name);
    Constant(const std::string& name, const std::string& descr);

    // Get/Set methods
    void Value(const std::string& val);
    const std::string& Value(void) const;
};

class Parameter : public Symbol
{
  private:

    std::string defaultvalue;

  public:

    // Constructors
    Parameter(const std::string& name);
    Parameter(const std::string& name, const std::string& descr);

    // Get/Set methods
    void DefaultValue(const std::string& val);
    const std::string& DefaultValue(void) const;
};

class Expression : public FormulaSymbol
{
  public:

    // Constructors
    Expression(const std::string& name);
    Expression(const std::string& name, const std::string& descr);

};


class StateVariable : public FormulaSymbol
{
  private:

    std::string periodicfrom;
    std::string periodicto;
    std::string default_ic;
    std::string default_history;
    std::string mass;

  public:

    // Constructors
    StateVariable(const std::string& name);
    StateVariable(const std::string& name, const std::string& descr);

    // Get/Set methods
    void PeriodicFrom(const std::string& pfrom);
    const std::string& PeriodicFrom(void) const;
    void PeriodicTo(const std::string& pto);
    const std::string& PeriodicTo(void) const;
    bool IsPeriodic();
    void DefaultInitialCondition(const std::string& ic);
    const std::string& DefaultInitialCondition(void) const;
    void DefaultHistory(const std::string& hist);
    const std::string& DefaultHistory(void) const;
    void Mass(const std::string& hist);
    const std::string& Mass(void) const;
};


class Function : public FormulaSymbol
{
  public:

    // Constructors
    Function(const std::string& name);
    Function(const std::string& name, const std::string& descr);
};


class VectorField : public Symbol
{
  private:

    // There is no implementation of the copy constructor.
    VectorField(const VectorField& vf);
    // There is no implementation of the assignment operator.
    VectorField operator=(const VectorField& vf);

//  protected:

    GiNaC::lst conname_list;     // P
    GiNaC::lst convalue_list;    // P
    GiNaC::lst parname_list;     // P
    GiNaC::lst pardefval_list;   // P
    GiNaC::lst internal_parname_list;     // P
    GiNaC::lst internal_pardefval_list;   // P
    GiNaC::lst exprname_list;    // P
    GiNaC::lst exprformula_list; // P
    GiNaC::lst expreqn_list;     // P
    GiNaC::lst varname_list;     // P
    GiNaC::lst varvecfield_list; // P
    GiNaC::lst vardefic_list;    // P
    GiNaC::lst vardefhist_list;  // P
    GiNaC::lst varmass_list;     // P
    GiNaC::lst funcname_list;    // P
    GiNaC::lst funcformula_list; // P

    GiNaC::lst allsymbols;       // P
    GiNaC::symbol IndVar;        // P
    
    std::vector<std::vector<double> > mass_matrix;

    // Everything is public, for now.
//  public:
	// These are added by the ReadXML
    std::string IndependentVariable;             // X
    std::vector<Constant *>      Constants;      // X
    std::vector<Parameter *>     Parameters;     // X
    std::vector<Expression *>    Expressions;    // X
    std::vector<StateVariable *> StateVariables; // X
    std::vector<Function *>      Functions;      // X
    bool IsDelay;                                // P
    bool TimeDependentDelay;                     // P
    bool StateDependentDelay;                    // P
    bool HasPi;                                  // P
    bool IsAutonomous;                           // P
    bool HasPeriod;                              // P
    std::vector<GiNaC::ex> Delays;               // P

  public:
    // Constructors
    VectorField();
    VectorField(const std::string& name, const std::string& descr);
    VectorField(const std::string& name, const std::string& descr, const std::string& indvar);
    
    // Destructor
    ~VectorField();

    void PrintXML(const std::string& cmdstr);
    int  ReadXML(const std::string& xmlfilename);
    int  ProcessSymbols(void);
    void PrintKnut(std::ostream& sys_out, std::map<std::string, std::string> options);
    bool isTimeDependentDelay() { return TimeDependentDelay; }
    bool isStateDependentDelay() { return StateDependentDelay; }
    
    size_t Knut_ndim() { return varvecfield_list.nops(); }
    size_t Knut_npar() { if (HasPeriod) return parname_list.nops(); else return parname_list.nops()+1; }
    size_t Knut_ntau() { return Delays.size()+1; }
    size_t Knut_nevent() { return funcname_list.nops(); }
    void Knut_tau(std::vector<GiNaC::ex>& exls);
    void Knut_tau_p(std::vector<GiNaC::ex>& exls);
    // returns a vector of expressions for the right-hand-side
    // each expression is dependent on Zlags_(i,j) and par_(k) only
    // the expressions, constants and parameters are subsituted
    // This is used to evaluate the right-hand-side
    void Knut_RHS(std::vector<GiNaC::ex>&);
    void Knut_mass(std::vector<double>& out);
    void Knut_RHS_p(std::vector<GiNaC::ex>&);
    void Knut_RHS_x(std::vector<GiNaC::ex>&);
    void Knut_RHS_xp(std::vector<GiNaC::ex>&);
    void Knut_RHS_xx(std::vector<GiNaC::ex>&);
    void Knut_stpar(std::vector<double>&);
    void Knut_stsol(std::vector<GiNaC::ex>&);
    void Knut_parnames(std::vector<std::string>&);
    
    static void Knut_tau_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, const size_t nv, const size_t nps, const size_t nd );
    static void Knut_tau_p_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp, const size_t nv, const size_t nps, const size_t nd );
    static void Knut_RHS_eval( std::vector<GiNaC::ex>& exls, KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel );
    static void Knut_RHS_p_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel, size_t alpha, const size_t nv, const size_t nps, const size_t nd);
    static void Knut_RHS_x_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel, size_t del, const size_t nv, const size_t nps, const size_t nd);
    static void Knut_RHS_xp_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel, size_t del, size_t alpha, const size_t nv, const size_t nps, const size_t nd);
    static void Knut_RHS_xx_eval( std::vector<GiNaC::ex>& exls, KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray3D<double>& vv, const KNVector& par, size_t sel, size_t del1, size_t del2, const size_t nv, const size_t nps, const size_t nd);
    static void Knut_RHS_stsol_eval( const std::vector<GiNaC::ex>& exls, KNArray2D<double>& out_p, const KNArray1D<double>& time );

  private:
    void AddConstant(Constant *c);
    void AddParameter(Parameter *p);
    void AddExpression(Expression *e);
    void AddStateVariable(StateVariable *sv);
    int  FindVar(const GiNaC::symbol&);
    GiNaC::lst FindVarsInEx(const GiNaC::ex &e);
    GiNaC::ex SubsAllExpressions(const GiNaC::ex &e);

    void AddFunction(Function *f);

    int FindDelay(GiNaC::ex&);
    int AddDelay(GiNaC::ex&);

    void CheckForDelay(const GiNaC::ex& f);
    
    void print(void);
    void Knut_ConvertStateToZlags(GiNaC::ex& f);
    void Knut_ConvertConstsPars(GiNaC::ex& f);
};

#endif
