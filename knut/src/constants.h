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
#include "exprsystem.h"
#include "knerror.h"

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

static inline size_t toSizeT(long int x)
{
  if ((x <= INT_MAX)&&(x >= 0)) return static_cast<size_t>(x);
  else { P_MESSAGE1("Unsigned integer is out of range."); return 0; }
}

static inline int toInt(size_t x)
{
  if (x < INT_MAX) return (int)x;
  else { P_MESSAGE1("Out of range."); return 0; }
}

class KNAbstractConstantsReader
{
  public:
    virtual bool getSystem(std::string& strout) = 0;
    virtual bool getTextField(const char* field, std::string& strout) = 0;
    virtual bool getIndexField(const char* field, size_t& out) = 0;
    virtual bool getDoubleField(const char* field, double& out) = 0;
    virtual bool getStringListField(const char* field, std::vector<std::string>& out) = 0;
    virtual bool getIndexListField(const char* field, std::vector<size_t>& out) = 0;
};

class KNAbstractConstantsWriter
{
  public:
    virtual void setSystem(const std::string& str) = 0;
    virtual void setTextField(const char* fieldname, const std::string& str) = 0;
    virtual void setDoubleField(const char* fieldname, double val) = 0;
    virtual void setIndexField(const char* fieldname, size_t val) = 0;
    virtual void setStringListField(const char* field, const std::vector<std::string>& vec) = 0;
    virtual void setIndexListField(const char* field, const std::vector<size_t>& vec) = 0;
};

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
    void set##capsName##Vector(const std::vector<type> &v) { name.value = v; constantChanged(#name); } \
    type get##capsName(size_t i) const \
      { if (i<get##capsName##Size()) return name.value[i]; else return zero; } \
    void set##capsName(size_t i, type d) { name.value[i] = d; constantChanged(#name); } \
    void nset##capsName(size_t i, type d) { name.value[i] = d; } \
    void static snset##capsName(KNConstantsBase& obj, size_t i, type d) { obj.name.value[i] = d; } \
    size_t get##capsName##Size() const { /*std::cout << name.value.size() << " " #name "MM\n";*/  return name.value.size(); } \
    void set##capsName##Size(size_t ns) \
    { \
      size_t ps = get##capsName##Size(); name.value.resize(ns); \
      for (size_t i=ps; i<ns; ++i) name.value[i] = zero; \
      constantChanged(#name); \
    }

class KNConstantNames
{
  public:
    using FPTR = void (*)();
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

// forward declearation
class KNConstantsBase;

// type declarations
template<typename TP> class TypeTupleTab  : public TypeTupleTabBase<TP>
{
  private:
    std::vector<TypeTupleStr<TP> > tab;
    size_t sz;
    const KNConstantsBase *parent;
  public:
    TypeTupleTab();
    TypeTupleTab(const KNConstantsBase *p) : tab(0), parent(p) { update(); }
    void         update(void); // needs to update 'tab'
    void         update(const KNConstantsBase *p) { parent = p; update(); }
    TP           CIndexToType(size_t idx) const;
    TP           NameToType(const char * name) const;
    TP           CodeToType(const char * code) const;
    void         CodeToTypeVector (std::vector<TP>& tp, const std::vector<std::string>& code);
    size_t       TypeToCIndex(TP tp) const;
    std::string  CIndexToTypeName(size_t idx) const;
    std::string  TypeToName(TP tp) const;
    std::string  TypeToCode(TP tp) const;
    void         TypeToCodeVector (std::vector<std::string>& code, const std::vector<TP>& tp);
    size_t       size() const;
};  

class KNConstantsBase
{
  friend class TypeTupleTab<Var>;
  protected:
    KNConstantNames cnames;
  private:   
    KN_CONSTANT(       std::string, system,     System,     "")
    KN_CONSTANT(       std::string, inputFile,  InputFile,  "") 
    KN_CONSTANT(       std::string, outputFile, OutputFile, "")
    KN_CONSTANT(       std::string, sysname,    SysName,    "")
    KN_CONSTANT(       std::string, sysType,    SysType,    "")
    KN_CONSTANT(       BifType,     fromType,   FromType,   BifNone)
    KN_CONSTANT(       size_t,      label,      Label,      0)
    KN_CONSTANT(       PtType,      pointType,  PointType,  SolUser)
    KN_CONSTANT(       Var,         cp,         Cp,         VarNone)
    KN_CONSTANT(       BranchSW,    branchSW,   BranchSW,   NOSwitch)
    KN_ARRAY_CONSTANT( Var,         parx,       Parx,       VarNone)
    KN_ARRAY_CONSTANT( Eqn,         eqns,       Eqns,       EqnNone)
    KN_ARRAY_CONSTANT( Var,         vars,       Vars,       VarNone)
    KN_CONSTANT(       size_t,      nInt,       NInt,       10)
    KN_CONSTANT(       size_t,      nDeg,       NDeg,       4)
    KN_CONSTANT(       size_t,      nMul,       NMul,       0)
    KN_CONSTANT(       bool,        stab,       Stab,       false)
    KN_CONSTANT(       size_t,      nInt1,      NInt1,      12)
    KN_CONSTANT(       size_t,      nInt2,      NInt2,      12)
    KN_CONSTANT(       size_t,      nDeg1,      NDeg1,      4)
    KN_CONSTANT(       size_t,      nDeg2,      NDeg2,      4)
    KN_CONSTANT(       size_t,      steps,      Steps,      10)
    KN_CONSTANT(       size_t,      iad,        Iad,        3)
    KN_CONSTANT(       size_t,      nPr,        NPr,        1)
    KN_CONSTANT(       size_t,      SSM,        SSM,        0)
    KN_CONSTANT(       double,      cpMin,      CpMin,      -DBL_MAX)
    KN_CONSTANT(       double,      cpMax,      CpMax,      DBL_MAX)
    KN_CONSTANT(       double,      ds,         Ds,         1e-3)
    KN_CONSTANT(       double,      dsMin,      DsMin,      1e-3)
    KN_CONSTANT(       double,      dsMax,      DsMax,      1e-3)
    KN_CONSTANT(       double,      dsStart,    DsStart,    1e-3)
    KN_CONSTANT(       double,      epsC,       EpsC,       1e-5)
    KN_CONSTANT(       double,      epsR,       EpsR,       1e-5)
    KN_CONSTANT(       double,      epsK,       EpsK,       1e-5)
    KN_CONSTANT(       size_t,      nItC,       NItC,       5)
    KN_CONSTANT(       size_t,      nItR,       NItR,       12)
    KN_CONSTANT(       size_t,      nItK,       NItK,       12)
    KN_ARRAY_CONSTANT( size_t,      symRe,      SymRe,      0)
    KN_ARRAY_CONSTANT( size_t,      symIm,      SymIm,      0)
    KN_CONSTANT(       size_t,      nDeri,      NDeri,      0)
    KN_CONSTANT(       double,      cAngle,     CAngle,     1.0/4) // Angle of the curvature
    // from sysname
    KN_CONSTANT(       size_t,      nPar,       NPar,       0)
    KN_CONSTANT(       size_t,      nDim,       NDim,       0)
    
  protected:
    std::vector<std::string> parNames;
    void addConstantName(const char* /*name*/) {}
 
  public:
    KNConstantsBase() : 
      system(this), 
      inputFile(this), outputFile(this), sysname(this), sysType(this), fromType(this), label(this), pointType(this),
      cp(this), branchSW(this), parx(this), eqns(this), vars(this), nInt(this), nDeg(this),
      nMul(this), stab(this), nInt1(this), nInt2(this), nDeg1(this), nDeg2(this), 
      steps(this), iad(this), nPr(this), SSM(this), cpMin(this), cpMax(this),
      ds(this), dsMin(this), dsMax(this), dsStart(this), epsC(this), epsR(this), epsK(this),
      nItC(this), nItR(this), nItK(this), symRe(this), symIm(this), nDeri(this), cAngle(this),
      nPar(this), nDim(this),
      PtTypeTable(this), EqnTable(this), BranchSWTable(this), BifTypeTable(this), VarTable(this) {}
    KNConstantsBase(const KNConstantsBase* src) : 
      system(this, src), 
      inputFile(this,src), outputFile(this,src), sysname(this,src), sysType(this,src), fromType(this,src), label(this,src), pointType(this,src),
      cp(this,src), branchSW(this,src), parx(this,src), eqns(this,src), vars(this,src), nInt(this,src), nDeg(this,src),
      nMul(this,src), stab(this,src), nInt1(this,src), nInt2(this,src), nDeg1(this,src), nDeg2(this,src), 
      steps(this,src), iad(this,src), nPr(this,src), SSM(this,src), cpMin(this,src), cpMax(this,src),
      ds(this,src), dsMin(this,src), dsMax(this,src), dsStart(this,src), epsC(this,src), epsR(this,src), epsK(this,src),
      nItC(this,src), nItR(this,src), nItK(this,src), symRe(this,src), symIm(this,src), nDeri(this,src), cAngle(this,src),
      nPar(this,src), nDim(this,src),
      PtTypeTable(this), EqnTable(this), BranchSWTable(this), BifTypeTable(this), VarTable(this) {}
    virtual ~KNConstantsBase() { }
    // for setting up the continuation
    bool toEqnVar(KNExprSystem& sys,
                  KNArray1D<Eqn>& eqn, KNArray1D<Var>& var,                 // input
                  KNArray1D<Eqn>& eqn_refine, KNArray1D<Var>& var_refine,   // output
                  KNArray1D<Eqn>& eqn_start, KNArray1D<Var>& var_start, bool& findangle);
    void loadXmlFileV5(const std::string &fileName);
    void printXmlFileV5(std::ostream& file);
    void saveXmlFileV5(const std::string &fileName);
    void read(KNAbstractConstantsReader& reader);
    void write(KNAbstractConstantsWriter& writer);

    // historical    
    void loadXmlFileV4(const std::string &fileName);
    void saveXmlFileV4(const std::string &fileName);
    void printXmlFileV4(std::ostream& file);
    void loadXmlFileV2(const std::string &fileName);

    
    // converting equation numbers to what can be used in continuation
    Eqn eqnFromTypeNum(char type, size_t num) const
    {
      if (type == 'E') return (Eqn)num;
      else return EqnNone;	
    }
    // converting variable numbers to what can be used in continuation
    Var varFromTypeNum(char type, size_t num) const
    {
      Var tmp;
      if (type == 'S') tmp = (Var)num;
      else if (type == 'P') tmp = (Var)(num + VarPAR0);
      else if (type == 'I') tmp = (Var)(num + VarInternal);
      else tmp = VarNone;
//      std :: cout << type << num << (char)num << " -> " << tmp << " ";
      return tmp;
    }
    char eqnType(Eqn )
    {
      return 'E';
    }
    size_t eqnNum(Eqn e)
    {
      return e;
    }
    char varType(Var v)
    {
      if (v < VarPAR0) return 'S';
      else if (v < VarInternal) return 'P';
      else return 'I';
    }
    size_t varNum(Var v)
    {
      if (v < VarPAR0) return v;
      else if (v < VarInternal) return v - VarPAR0;
      else return v - VarInternal;
    }
    
    void initParNames(size_t np)
    {
      parNames.resize(VarToIndex(VarEnd,np));
      for (size_t i=0; i<np; ++i)
      {
        char buf[16];
        snprintf(buf, 15, "%zd", i);
        parNames[i] = "PAR(" + std::string(buf) + ")";
      }
      parNames[VarToIndex(VarPeriod,np)] = "Per. Mul.";
      parNames[VarToIndex(VarAngle,np)] = "Angle (NS)";
      parNames[VarToIndex(VarRot,np)] = "Rotation num.";      
    }
    void initDimensions (const KNExprSystem* sys);
    // This loads the shared object file
    virtual void setSystemText (const std::string& str);
    virtual void setSysNameText (const std::string& str, bool testing = false);
    const std::vector<std::string>& getParNames () const { return parNames; }
    // this is to notify the system of changes
    virtual void constantChanged(const char* /*name*/) { }
    
    TypeTupleTab<PtType> PtTypeTable;
    TypeTupleTab<Eqn> EqnTable;
    TypeTupleTab<BranchSW> BranchSWTable;
    TypeTupleTab<BifType> BifTypeTable;
    TypeTupleTab<Var> VarTable;
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
    ~KNConstants() override { }
};

template<typename TP> void TypeTupleTab<TP>::update()
{
  size_t k = 0;
  while (TypeTupleTabBase<TP>::tabStatic[k].index < ~(size_t)0) ++k;
  tab.resize(k);
  for (size_t i=0; i<tab.size(); ++i)
  {
    tab[i].index = i; // tabStatic[i].index;
    tab[i].type = TypeTupleTabBase<TP>::tabStatic[i].type;
    tab[i].name = TypeTupleTabBase<TP>::tabStatic[i].name;
    tab[i].code = TypeTupleTabBase<TP>::tabStatic[i].code;
  }
}

template<typename TP> TP TypeTupleTab<TP>::CIndexToType(size_t idx) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].index == idx) return tab[k].type;
  return tab[0].type;
}

template<typename TP> size_t TypeTupleTab<TP>::TypeToCIndex(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].index;
  return tab[0].index;
}

template<typename TP> std::string TypeTupleTab<TP>::TypeToCode(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].code;
  return tab[0].code;
}

template<typename TP> void TypeTupleTab<TP>::TypeToCodeVector (std::vector<std::string>& code, const std::vector<TP>& tp)
{
  code.resize(tp.size());
  for (size_t k=0; k < tp.size(); ++k) code[k] = TypeToCode(tp[k]);
}

template<typename TP> std::string TypeTupleTab<TP>::TypeToName(TP tp) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].type == tp) return tab[k].name;
  return tab[0].name;
}

template<typename TP> std::string TypeTupleTab<TP>::CIndexToTypeName(size_t idx) const
{
  for (size_t k=0; k<tab.size(); ++k)
    if (tab[k].index == idx) return tab[k].name;
  return tab[0].name;
}

template<typename TP> TP TypeTupleTab<TP>::CodeToType(const char * code) const
{
  if (code)
  {
    for (size_t k=0; k<tab.size(); ++k)
      if (!tab[k].code.compare(code)) return tab[k].type;
  }
  return tab[0].type;
}

template<typename TP> void TypeTupleTab<TP>::CodeToTypeVector (std::vector<TP>& tp, const std::vector<std::string>& code)
{
  tp.resize(code.size());
  for (size_t k=0; k < code.size(); ++k) tp[k] = CodeToType(code[k].c_str());
}

template<typename TP> TP TypeTupleTab<TP>::NameToType(const char * name) const
{
  if (name)
  {
    for (size_t k=0; k<tab.size(); ++k)
      if (!tab[k].name.compare(name)) return tab[k].type;
  }
  return tab[0].type;
}

template<typename TP> size_t TypeTupleTab<TP>::size() const
{
  return tab.size();
}

template<> void TypeTupleTab<Var>::update();

#endif
