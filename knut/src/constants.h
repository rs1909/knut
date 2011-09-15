// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef CONSTANTS_H
#define CONSTANTS_H

// Includes local
#include "pointtype.h"
#include "system.h"

// Includes standard library
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cstdio>

#include <climits>
#include <cfloat>

static inline int toInt(long int x)
{
  if ((x <= INT_MAX)&&(x >= INT_MIN)) return (int)x;
  else { P_MESSAGE1("Integer is out of range."); return 0; }
}

static inline int toInt(size_t x)
{
  if (x < INT_MAX) return (int)x;
  else { P_MESSAGE1("Out of range."); return 0; }
}

#define KN_CONSTANT(type, name, capsName, defaultValue) private: \
    class var_##name { public: type value; \
      var_##name () : value(defaultValue) \
      { \
        KNConstantNames::addConstantName(#type, #name, &value, \
          reinterpret_cast<KNConstantNames::FPTR>(&KNConstants::snset##capsName)); \
      } \
      var_##name (const var_##name & src) { value = src.value; } \
    } name; \
  public: \
    const type & get##capsName() const { return name.value; } \
    void set##capsName(const type & d) { name.value = d; constantChanged(#name); } \
    void nset##capsName(const type & d) { name.value = d; } \
    void static snset##capsName(KNConstants& obj, const type & d) { obj.name.value = d; }

#define KN_ARRAY_CONSTANT( type, name, capsName, zero) private: \
    class var_##name { public: std::vector<type> value; \
      var_##name () { KNConstantNames::addConstantName("vector<" #type ">", #name, &value, (void(*)())&KNConstants::snset##capsName); } \
      var_##name (const var_##name & src) { value = src.value; } \
    } name; \
  public: \
    type get##capsName(int i) const \
      { if (i<get##capsName##Size() && i>=0) return name.value[static_cast<unsigned int>(i)]; else return zero; } \
    void set##capsName(int i, type d) { name.value[static_cast<unsigned int>(i)] = d; constantChanged(#name); } \
    void nset##capsName(int i, type d) { name.value[static_cast<unsigned int>(i)] = d; } \
    void static snset##capsName(KNConstants& obj, int i, type d) { obj.name.value[static_cast<unsigned int>(i)] = d; } \
    int get##capsName##Size() const { /*std::cout << name.value.size() << " " #name "MM\n";*/  return static_cast<int>(name.value.size()); } \
    void set##capsName##Size(int ns) \
    { \
      int ps = get##capsName##Size(); name.value.resize(static_cast<unsigned int>(ns)); \
      for (int i=ps; i<ns; ++i) name.value[static_cast<unsigned int>(i)] = zero; \
      constantChanged(#name); \
    }

class KNConstantNames
{
  public:
    typedef void(*FPTR)();
  private:
    class tuple_t {
      public:
        std::string type;
        void* var;
	FPTR nsetFun;
    };
    static std::map<std::string, tuple_t> nlist;
  public: 
    static void addConstantName(const char * type, const char * name, void * var, void(*nsetFun)());
    static const std::string& findType(const char * name);
    static const void* findValue(const char * name);
    static FPTR findFun(const char * name);
};

template<typename TP> struct TypeTuple
{
  int          index;
  TP           type;
  const char * code;
  const char * name;
};

template<typename TP> struct TypeTupleStr
{
  int          index;
  TP           type;
  std::string  code;
  std::string  name;
};

class KNConstants
{
  private:
    
    KN_CONSTANT(       std::string, inputFile,  InputFile,  "") 
    KN_CONSTANT(       std::string, outputFile, OutputFile, "")
    KN_CONSTANT(       std::string, sysname,    SysName,    "")
    KN_CONSTANT(       BifType,     fromType,   FromType,   BifNone)
    KN_CONSTANT(       int,         label,      Label,      0)
    KN_CONSTANT(       PtType,      pointType,  PointType,  SolUser)
    KN_CONSTANT(       Var,         cp,         Cp,         VarNone)
    KN_CONSTANT(       BranchSW,    branchSW,   BranchSW,   NOSwitch)
    KN_ARRAY_CONSTANT( Var,         parx,       Parx,       VarNone)
    KN_ARRAY_CONSTANT( Eqn,         eqns,       Eqns,       EqnNone)
    KN_ARRAY_CONSTANT( Var,         vars,       Vars,       VarNone)
    KN_CONSTANT(       int,         nInt,       NInt,       10)
    KN_CONSTANT(       int,         nDeg,       NDeg,       4)
    KN_CONSTANT(       int,         nMul,       NMul,       0)
    KN_CONSTANT(       bool,        stab,       Stab,       false)
    KN_CONSTANT(       int,         nInt1,      NInt1,      12)
    KN_CONSTANT(       int,         nInt2,      NInt2,      12)
    KN_CONSTANT(       int,         nDeg1,      NDeg1,      4)
    KN_CONSTANT(       int,         nDeg2,      NDeg2,      4)
    KN_CONSTANT(       int,         steps,      Steps,      10)
    KN_CONSTANT(       int,         iad,        Iad,        3)
    KN_CONSTANT(       int,         nPr,        NPr,        1)
    KN_CONSTANT(       double,      cpMin,      CpMin,      -DBL_MAX)
    KN_CONSTANT(       double,      cpMax,      CpMax,      DBL_MAX)
    KN_CONSTANT(       double,      ds,         Ds,         1e-3)
    KN_CONSTANT(       double,      dsMin,      DsMin,      1e-3)
    KN_CONSTANT(       double,      dsMax,      DsMax,      1e-3)
    KN_CONSTANT(       double,      dsStart,    DsStart,    1e-3)
    KN_CONSTANT(       double,      epsC,       EpsC,       1e-5)
    KN_CONSTANT(       double,      epsR,       EpsR,       1e-5)
    KN_CONSTANT(       double,      epsK,       EpsK,       1e-5)
    KN_CONSTANT(       int,         nItC,       NItC,       5)
    KN_CONSTANT(       int,         nItR,       NItR,       12)
    KN_CONSTANT(       int,         nItK,       NItK,       12)
    KN_ARRAY_CONSTANT( int,         symRe,      SymRe,      0)
    KN_ARRAY_CONSTANT( int,         symIm,      SymIm,      0)
    KN_CONSTANT(       int,         nDeri,      NDeri,      0)
    KN_CONSTANT(       double,      cAngle,     CAngle,     1.0/4) // Angle of the curvature
    // from sysname
    KN_CONSTANT(       int,         nPar,       NPar,       0)
    KN_CONSTANT(       int,         nDim,       NDim,       0)
    
  protected:
    std::vector<std::string> parNames;
    void addConstantName(const char* name) {}
 
  public:
    KNConstants() : PtTypeTable(this), EqnTable(this), BranchSWTable(this), BifTypeTable(this), VarTable(this) {}
    virtual ~KNConstants() { }
    // for setting up the continuation
    bool toEqnVar(KNSystem& sys,
                  KNArray1D<Eqn>& eqn, KNArray1D<Var>& var,                 // input
                  KNArray1D<Eqn>& eqn_refine, KNArray1D<Var>& var_refine,   // output
                  KNArray1D<Eqn>& eqn_start, KNArray1D<Var>& var_start, bool& findangle);
    void loadXmlFile(const std::string &fileName);
    void loadXmlFileOLD(const std::string &fileName);
    void saveXmlFile(const std::string &fileName);
    void printXmlFile(std::ostream& file);
    
    // converting equation numbers to what can be used in continuation
    Eqn eqnFromTypeNum(char type, int num) const
    {
      if (type == 'E') return (Eqn)num;
      else return EqnNone;	
    }
    // converting variable numbers to what can be used in continuation
    Var varFromTypeNum(char type, int num) const
    {
      Var tmp;
      if (type == 'S') tmp = (Var)num;
      else if (type == 'P') tmp = (Var)(num + VarPAR0);
      else if (type == 'I') tmp = (Var)(num + VarInternal);
      else tmp = VarNone;
//      std :: cout << type << num << (char)num << " -> " << tmp << " ";
      return tmp;
    }
    char eqnType(Eqn e)
    {
      return 'E';
    }
    int eqnNum(Eqn e)
    {
      return e;
    }
    char varType(Var v)
    {
      if (v < VarPAR0) return 'S';
      else if (v < VarInternal) return 'P';
      else return 'I';
    }
    int varNum(Var v)
    {
      if (v < VarPAR0) return v;
      else if (v < VarInternal) return v - VarPAR0;
      else return v - VarInternal;
    }
    
    void initParNames(int np)
    {
      parNames.resize(VarToIndex(VarEnd,np));
      for (int i=0; i<np; ++i)
      {
        char buf[16];
        snprintf(buf, 15, "%d", i);
        parNames[i] = "PAR(" + std::string(buf) + ")";
      }
      parNames[VarToIndex(VarPeriod,np)] = "Per. Mul.";
      parNames[VarToIndex(VarAngle,np)] = "Angle (NS)";
      parNames[VarToIndex(VarRot,np)] = "Rotation num.";      
    }
    void initDimensions(const KNSystem* sys);
    // This loads the shared object file
    virtual void setSysNameText(const std::string& str, bool testing = false);
    const std::vector<std::string>& getParNames() const { return parNames; }
    // this is to notify the system of changes
    virtual void constantChanged(const char* name) { }
    
    // type declarations
    template<typename TP> class TypeTupleTab
    {
    private:
      static const TypeTuple<TP> tabStatic[];
      std::vector<TypeTupleStr<TP> > tab;
      int sz;
      const KNConstants *parent;
    public:
      TypeTupleTab();
      TypeTupleTab(KNConstants *p) : tab(0), parent(p) { update(); }
      void         update(void); // needs to update 'tab'
      TP           CIndexToType(int idx) const;
      TP           NameToType(const char * name) const;
      TP           CodeToType(const char * name) const;
      int          TypeToCIndex(TP tp) const;
      std::string  CIndexToTypeName(int idx) const;
      std::string  TypeToName(TP tp) const;
      std::string  TypeToCode(TP tp) const;
      int          size() const;
    };
    
    const KNConstants::TypeTupleTab<PtType> PtTypeTable;
    const KNConstants::TypeTupleTab<Eqn> EqnTable;
    const KNConstants::TypeTupleTab<BranchSW> BranchSWTable;
    const KNConstants::TypeTupleTab<BifType> BifTypeTable;
    KNConstants::TypeTupleTab<Var> VarTable;
};

template<typename TP> void KNConstants::TypeTupleTab<TP>::update()
{
  int k = 0;
  while (tabStatic[k].index >= 0) ++k;
  tab.resize(k);
  for (size_t i=0; i<tab.size(); ++i)
  {
    tab[i].index = i; // tabStatic[i].index;
    tab[i].type = tabStatic[i].type;
    tab[i].name = tabStatic[i].name;
    tab[i].code = tabStatic[i].code;
  }
}

template<typename TP> TP KNConstants::TypeTupleTab<TP>::CIndexToType(int idx) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].index == idx) return tab[k].type;
  return tab[0].type;
}

template<typename TP> int KNConstants::TypeTupleTab<TP>::TypeToCIndex(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].index;
  return tab[0].index;
}

template<typename TP> std::string KNConstants::TypeTupleTab<TP>::TypeToCode(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].code;
  return tab[0].code;
}

template<typename TP> std::string KNConstants::TypeTupleTab<TP>::TypeToName(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].name;
  return tab[0].name;
}

template<typename TP> std::string KNConstants::TypeTupleTab<TP>::CIndexToTypeName(int idx) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].index == idx) return tab[k].name;
  return tab[0].name;
}

template<typename TP> TP KNConstants::TypeTupleTab<TP>::CodeToType(const char * name) const
{
  if (name)
  {
    for (size_t k=0; k<tab.size(); ++k)
      if (!tab[k].code.compare(name)) return tab[k].type;
  }
  return tab[0].type;
}

template<typename TP> TP KNConstants::TypeTupleTab<TP>::NameToType(const char * name) const
{
  if (name)
  {
    for (size_t k=0; k<tab.size(); ++k)
      if (!tab[k].name.compare(name)) return tab[k].type;
  }
  return tab[0].type;
}

template<typename TP> int KNConstants::TypeTupleTab<TP>::size() const
{
  return tab.size();
}

template<> void KNConstants::TypeTupleTab<Var>::update();

#endif
