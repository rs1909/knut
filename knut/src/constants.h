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

class KNConstants
{
  private:
    
    KN_CONSTANT(       std::string, inputFile,  InputFile,  "") 
    KN_CONSTANT(       std::string, outputFile, OutputFile, "")
    KN_CONSTANT(       std::string, sysname,    SysName,    "")
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
    // for setting up the continuation
    bool toEqnVar(KNSystem& sys,
                  KNArray1D<Eqn>& eqn, KNArray1D<Var>& var,                 // input
                  KNArray1D<Eqn>& eqn_refine, KNArray1D<Var>& var_refine,   // output
                  KNArray1D<Eqn>& eqn_start, KNArray1D<Var>& var_start, bool& findangle);
    void loadXmlFile(const std::string &fileName);
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
        parNames[i] = "Par. " + std::string(buf);
      }
      parNames[VarToIndex(VarPeriod,np)] = "Per. Mul.";
      parNames[VarToIndex(VarAngle,np)] = "Angle (NS)";
      parNames[VarToIndex(VarRot,np)] = "Rotation num.";      
    }
    void initDimensions(const KNSystem* sys)
    {
      initParNames(sys->npar());
      const char **names = new const char *[sys->npar()+1];
      for (int i=0; i<sys->npar()+1; ++i) names[i] = 0;
      sys->parnames(names);
      for (int i=0; i<sys->npar(); ++i) if (names[i] != 0) parNames[i] = names[i];
      delete[] names;
      // this will send a signal so parNames must be ready
      setNPar(sys->npar());
      setNDim(sys->ndim());
    }
    // This loads the shared object file
    virtual void setSysNameText(const std::string& str, bool testing = false)
    {
      setSysName(str);
      KNSystem* sys = 0;
      try
      {
        sys = new KNSystem(getSysName(), getNDeri());
        initDimensions(sys);
        delete sys;
      }
      catch (KNException ex)
      {
        delete sys;
        setNPar(0);
        setNDim(0);
        throw(ex);
      }
    }
    const std::vector<std::string>& getParNames() const { return parNames; }
    // this is to notify the system of changes
    virtual void constantChanged(const char* name) { }
    
    template<typename TP> TP           CIndexToType(const TypeTuple<TP>* tab, int idx);
    template<typename TP> int          TypeToCIndex(const TypeTuple<TP>* tab, TP tp);
    template<typename TP> const char * CIndexToTypeName(const TypeTuple<TP>* tab, int idx);
    template<typename TP> const char * TypeToName(const TypeTuple<TP>* tab, TP tp);
    template<typename TP> int          TypeTableSize(const TypeTuple<TP>* tab);

    static const TypeTuple<PtType> PtTypeTable[];
    PtType       CIndexToPtType(int idx) { return CIndexToType<PtType>(PtTypeTable, idx); }
    int          PtTypeToCIndex(PtType pt) { return TypeToCIndex<PtType>(PtTypeTable, pt); }
    const char * CIndexToPtTypeName(int idx) { return CIndexToTypeName<PtType>(PtTypeTable, idx); }
//    const char * PtTypeToName(PtType tp) { return TypeToName<PtType>(PtTypeTable, tp); }
    int          PtTypeTableSize() { return TypeTableSize<PtType>(PtTypeTable); }
  
    static const TypeTuple<Eqn> EqnTable[];
    Eqn          CIndexToEqn(int idx) { return CIndexToType<Eqn>(EqnTable, idx); }
    int          EqnToCIndex(Eqn pt) { return TypeToCIndex<Eqn>(EqnTable, pt); }
    const char * CIndexToEqnName(int idx) { return CIndexToTypeName<Eqn>(EqnTable, idx); }
//    const char * EqnToName(Eqn tp) { return TypeToName<Eqn>(EqnTable, tp); }
    int          EqnTableSize() { return TypeTableSize<Eqn>(EqnTable); }
    
    static const TypeTuple<BranchSW> BranchSWTable[];
    BranchSW     CIndexToBranchSW(int idx) { return CIndexToType<BranchSW>(BranchSWTable, idx); }
    int          BranchSWToCIndex(BranchSW pt) { return TypeToCIndex<BranchSW>(BranchSWTable, pt); }
    const char * CIndexToBranchSWName(int idx) { return CIndexToTypeName<BranchSW>(BranchSWTable, idx); }
//    const char * BranchSWToName(BranchSW tp) { return TypeToName<BranchSW>(BranchSWTable, tp); }
    int          BranchSWTableSize() { return TypeTableSize<BranchSW>(BranchSWTable); }

    static const TypeTuple<Var> VarTableS[];
    static const TypeTuple<Var> VarTableI[];
    Var          CIndexToVar(int idx);
    int          VarToCIndex(Var pt);
    const char * CIndexToVarName(int idx);
    int          VarTableSize();
};

template<typename TP> TP           KNConstants::CIndexToType(const TypeTuple<TP>* tab, int idx)
{
  for (int k=0; tab[k].index >= 0; ++k)
    if (tab[k].index == idx) return tab[k].type;
  return tab[0].type;
}

template<typename TP> int          KNConstants::TypeToCIndex(const TypeTuple<TP>* tab, TP tp)
{
  for (int k=0; tab[k].index >= 0; ++k)
    if (tab[k].type == tp) return tab[k].index;
  return tab[0].index;
}

template<typename TP> const char * KNConstants::CIndexToTypeName(const TypeTuple<TP>* tab, int idx)
{
  for (int k=0; tab[k].index >= 0; ++k)
    if (tab[k].index == idx) return tab[k].name;
  return tab[0].name;
}

template<typename TP> const char * KNConstants::TypeToName(const TypeTuple<TP>* tab, TP tp)
{
  for (int k=0; tab[k].index >= 0; ++k)
    if (tab[k].type == tp) return tab[k].name;
  return tab[0].name;
}

template<typename TP> int          KNConstants::TypeTableSize(const TypeTuple<TP>* tab)
{
  int k = 0;
  while (tab[k].index >= 0) ++k;
  return k;
}

#endif
