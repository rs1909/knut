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
    class var_##name { public: type value; KNConstantsBase* base; \
      var_##name (); \
      var_##name (KNConstantsBase* b) : value(defaultValue), base(b) \
      { \
        b->cnames.addConstantName(#type, #name, &value, \
          reinterpret_cast<KNConstantNames::FPTR>(&KNConstantsBase::snset##capsName)); \
      } \
      var_##name (KNConstantsBase* b, const KNConstantsBase* src) : value(src->get##capsName()), base(b) \
      { \
        base->cnames.addConstantName(#type, #name, &value, \
          reinterpret_cast<KNConstantNames::FPTR>(&KNConstantsBase::snset##capsName)); \
      } \
    } name; \
  public: \
    const type & get##capsName() const { return name.value; } \
    void set##capsName(const type & d) { name.value = d; constantChanged(#name); } \
    void nset##capsName(const type & d) { name.value = d; } \
    void static snset##capsName(KNConstantsBase& obj, const type & d) { obj.name.value = d; }

#define KN_ARRAY_CONSTANT( type, name, capsName, zero) private: \
    class var_##name { public: std::vector<type> value; KNConstantsBase* base; \
      var_##name (); \
      var_##name (KNConstantsBase* b) : base(b) \
        { b->cnames.addConstantName("vector<" #type ">", #name, &value, (void(*)())&KNConstantsBase::snset##capsName); } \
      var_##name (KNConstantsBase* b, const KNConstantsBase* src) : value(src->get##capsName##Vector()), base(b) \
        { b->cnames.addConstantName("vector<" #type ">", #name, &value, (void(*)())&KNConstantsBase::snset##capsName); } \
    } name; \
  public: \
    const std::vector<type> & get##capsName##Vector() const { return name.value; } \
    type get##capsName(int i) const \
      { if (i<get##capsName##Size() && i>=0) return name.value[static_cast<unsigned int>(i)]; else return zero; } \
    void set##capsName(int i, type d) { name.value[static_cast<unsigned int>(i)] = d; constantChanged(#name); } \
    void nset##capsName(int i, type d) { name.value[static_cast<unsigned int>(i)] = d; } \
    void static snset##capsName(KNConstantsBase& obj, int i, type d) { obj.name.value[static_cast<unsigned int>(i)] = d; } \
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
    std::map<std::string, tuple_t> nlist;
  public: 
    void addConstantName(const char * type, const char * name, void * var, void(*nsetFun)());
    const std::string& findType(const char * name);
    const void* findValue(const char * name);
    FPTR findFun(const char * name);
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

class KNConstantsBase
{
  protected:
    KNConstantNames cnames;
  private:   
    KN_CONSTANT(       std::string, inputFile,  InputFile,  "") 
    KN_CONSTANT(       std::string, outputFile, OutputFile, "")
    KN_CONSTANT(       std::string, sysname,    SysName,    "")
    KN_CONSTANT(       std::string, sysType,    SysType,    "")
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
    KNConstantsBase() : 
      inputFile(this), outputFile(this), sysname(this), sysType(this), fromType(this), label(this), pointType(this),
      cp(this), branchSW(this), parx(this), eqns(this), vars(this), nInt(this), nDeg(this),
      nMul(this), stab(this), nInt1(this), nInt2(this), nDeg1(this), nDeg2(this), 
      steps(this), iad(this), nPr(this), cpMin(this), cpMax(this),
      ds(this), dsMin(this), dsMax(this), dsStart(this), epsC(this), epsR(this), epsK(this),
      nItC(this), nItR(this), nItK(this), symRe(this), symIm(this), nDeri(this), cAngle(this),
      nPar(this), nDim(this),
      PtTypeTable(this), EqnTable(this), BranchSWTable(this), BifTypeTable(this), VarTable(this) {}
    KNConstantsBase(const KNConstantsBase* src) : 
      inputFile(this,src), outputFile(this,src), sysname(this,src), sysType(this,src), fromType(this,src), label(this,src), pointType(this,src),
      cp(this,src), branchSW(this,src), parx(this,src), eqns(this,src), vars(this,src), nInt(this,src), nDeg(this,src),
      nMul(this,src), stab(this,src), nInt1(this,src), nInt2(this,src), nDeg1(this,src), nDeg2(this,src), 
      steps(this,src), iad(this,src), nPr(this,src), cpMin(this,src), cpMax(this,src),
      ds(this,src), dsMin(this,src), dsMax(this,src), dsStart(this,src), epsC(this,src), epsR(this,src), epsK(this,src),
      nItC(this,src), nItR(this,src), nItK(this,src), symRe(this,src), symIm(this,src), nDeri(this,src), cAngle(this,src),
      nPar(this,src), nDim(this,src),
      PtTypeTable(this), EqnTable(this), BranchSWTable(this), BifTypeTable(this), VarTable(this) {}
    virtual ~KNConstantsBase() { }
    // for setting up the continuation
    bool toEqnVar(KNSystem& sys,
                  KNArray1D<Eqn>& eqn, KNArray1D<Var>& var,                 // input
                  KNArray1D<Eqn>& eqn_refine, KNArray1D<Var>& var_refine,   // output
                  KNArray1D<Eqn>& eqn_start, KNArray1D<Var>& var_start, bool& findangle);
    void loadXmlFile(const std::string &fileName);
    void loadXmlFileOld(const std::string &fileName);
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
      const KNConstantsBase *parent;
    public:
      TypeTupleTab();
      TypeTupleTab(const KNConstantsBase *p) : tab(0), parent(p) { update(); }
      void         update(void); // needs to update 'tab'
      void         update(const KNConstantsBase *p) { parent = p; update(); }
      TP           CIndexToType(int idx) const;
      TP           NameToType(const char * name) const;
      TP           CodeToType(const char * name) const;
      int          TypeToCIndex(TP tp) const;
      std::string  CIndexToTypeName(int idx) const;
      std::string  TypeToName(TP tp) const;
      std::string  TypeToCode(TP tp) const;
      int          size() const;
    };  
    KNConstantsBase::TypeTupleTab<PtType> PtTypeTable;
    KNConstantsBase::TypeTupleTab<Eqn> EqnTable;
    KNConstantsBase::TypeTupleTab<BranchSW> BranchSWTable;
    KNConstantsBase::TypeTupleTab<BifType> BifTypeTable;
    KNConstantsBase::TypeTupleTab<Var> VarTable;
};

class KNConstants : public KNConstantsBase
{
  public:
    KNConstants() {}
    KNConstants( const KNConstants& cs ) :  KNConstantsBase(cs)
    {
      PtTypeTable.update(this);
      EqnTable.update(this);
      BranchSWTable.update(this);
      BifTypeTable.update(this);
      VarTable.update(this);
    }
    virtual ~KNConstants() { }
};

template<typename TP> void KNConstantsBase::TypeTupleTab<TP>::update()
{
  size_t k = 0;
  while (tabStatic[k].index >= 0) ++k;
  tab.resize(k);
  for (size_t i=0; i<tab.size(); ++i)
  {
    tab[i].index = static_cast<int>(i); // tabStatic[i].index;
    tab[i].type = tabStatic[i].type;
    tab[i].name = tabStatic[i].name;
    tab[i].code = tabStatic[i].code;
  }
}

template<typename TP> TP KNConstantsBase::TypeTupleTab<TP>::CIndexToType(int idx) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].index == idx) return tab[k].type;
  return tab[0].type;
}

template<typename TP> int KNConstantsBase::TypeTupleTab<TP>::TypeToCIndex(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].index;
  return tab[0].index;
}

template<typename TP> std::string KNConstantsBase::TypeTupleTab<TP>::TypeToCode(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].code;
  return tab[0].code;
}

template<typename TP> std::string KNConstantsBase::TypeTupleTab<TP>::TypeToName(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].name;
  return tab[0].name;
}

template<typename TP> std::string KNConstantsBase::TypeTupleTab<TP>::CIndexToTypeName(int idx) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].index == idx) return tab[k].name;
  return tab[0].name;
}

template<typename TP> TP KNConstantsBase::TypeTupleTab<TP>::CodeToType(const char * name) const
{
  if (name)
  {
    for (size_t k=0; k<tab.size(); ++k)
      if (!tab[k].code.compare(name)) return tab[k].type;
  }
  return tab[0].type;
}

template<typename TP> TP KNConstantsBase::TypeTupleTab<TP>::NameToType(const char * name) const
{
  if (name)
  {
    for (size_t k=0; k<tab.size(); ++k)
      if (!tab[k].name.compare(name)) return tab[k].type;
  }
  return tab[0].type;
}

template<typename TP> int KNConstantsBase::TypeTupleTab<TP>::size() const
{
  return tab.size();
}

template<> void KNConstantsBase::TypeTupleTab<Var>::update();

#endif
