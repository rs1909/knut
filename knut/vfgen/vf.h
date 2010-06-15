
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

DECLARE_FUNCTION_2P(delay)
DECLARE_FUNCTION_2P(Zlags_)

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

    GiNaC::lst conname_list;
    GiNaC::lst convalue_list;
    GiNaC::lst parname_list;
    GiNaC::lst pardefval_list;
    GiNaC::lst exprname_list;
    GiNaC::lst exprformula_list;
    GiNaC::lst expreqn_list;
    GiNaC::lst varname_list;
    GiNaC::lst varvecfield_list;
    GiNaC::lst vardefic_list;
    GiNaC::lst vardefhist_list;
    GiNaC::lst funcname_list;
    GiNaC::lst funcformula_list;

    GiNaC::lst allsymbols;
    GiNaC::symbol IndVar;

    // Everything is public, for now.
//  public:

    std::string IndependentVariable;
    std::vector<Constant *>      Constants;
    std::vector<Parameter *>     Parameters;
    std::vector<Expression *>    Expressions;
    std::vector<StateVariable *> StateVariables;
    std::vector<Function *>      Functions;
    bool IsDelay;
    bool TimeDependentDelay;
    bool StateDependentDelay;
    bool HasPi;
    bool IsAutonomous;
    bool HasPeriod;

    std::vector<GiNaC::ex> Delays;

  public:
    // Constructors
    VectorField();
    VectorField(const std::string& name, const std::string& descr);
    VectorField(const std::string& name, const std::string& descr, const std::string& indvar);
    
    // Destructor
    ~VectorField();


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

    void PrintXML(const std::string& cmdstr);
    int  ReadXML(const std::string& xmlfilename);

    void CheckForDelay(const GiNaC::ex& f);
    int  ProcessSymbols(void);
    bool isTimeDependentDelay() { return TimeDependentDelay; }
    bool isStateDependentDelay() { return StateDependentDelay; }
    
    void print(void);
    void Knut_ConvertDelaysToZlags(GiNaC::ex& f);
    void Knut_ConvertStateToZlags(GiNaC::ex& f);
    void Knut_PrintParDerivs(std::ostream &dout, const std::vector<GiNaC::ex> &e);
    void Knut_PrintJacobians(std::ostream &dout, const std::vector<GiNaC::ex> &e);
    void Knut_PrintXandParJacobians(std::ostream &dout, const std::vector<GiNaC::ex> &e);
    void Knut_PrintHessiansTimesV(std::ostream &dout, const std::vector<GiNaC::ex> &e);
    void PrintKnut(std::ostream& sys_out, std::map<std::string, std::string> options);
};

#endif
