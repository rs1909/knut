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

#ifndef EXTRACT_NAMES

  // Includes local
  #include "pointtype.h"
  #include "system.h"
  
  // Includes standard library
  #include <iostream>
  #include <vector>
  #include <string>
  #include <cstdio>

  #include <climits>

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

#endif

#ifdef REMOVE_TYPES
# define SELF(a)
#else 
# define SELF(a) a
#endif

#define KN_CONSTANT(name, type, getName, setName, nsetName) private: \
    SELF(type) name; \
  public: \
    type getName() const { return name; } \
    virtual void setName(type d) { name = d; constantChanged(#name); } \
    virtual void nsetName(type d) { name = d; }

#define KN_STRING_CONSTANT(name, getName, setName, nsetName) private: \
    SELF(std::string) name; \
  public: \
    const std::string& getName() const { return name; } \
    virtual void setName(const std::string& d) { name = d; constantChanged(#name); } \
    virtual void nsetName(const std::string& d) { name = d; }

#define KN_ARRAY_CONSTANT( name, type, getName, setName, getSize, setSize, zero) private: \
    SELF(std::vector<type>) name; \
  public: \
    type getName(int i) const \
      { if (i<getSize() && i>=0) return name[static_cast<unsigned int>(i)]; else return zero; } \
    virtual void setName(int i, type d) { name[static_cast<unsigned int>(i)] = d; constantChanged(#name); } \
    int getSize() const { return static_cast<int>(name.size()); } \
    virtual void setSize(int ns) \
    { \
      int ps = getSize(); name.resize(static_cast<unsigned int>(ns)); \
      for (int i=ps; i<ns; ++i) name[static_cast<unsigned int>(i)] = zero; \
      constantChanged(#name); \
    }

class NConstants
{
  private:
    
    KN_STRING_CONSTANT(inputFile,            getInputFile,  setInputFile,  nsetInputFile) 
    KN_STRING_CONSTANT(outputFile,           getOutputFile, setOutputFile, nsetOutputFile)
    KN_STRING_CONSTANT(sysname,              getSysName,    setSysName,    nsetSysName)
    KN_CONSTANT(       label,      int,      getLabel,      setLabel,      nsetLabel)
    KN_CONSTANT(       pointType,  PtType,   getPointType,  setPointType,  nsetPointType)
    KN_CONSTANT(       cpType,     char,     getCpType,     setCpType,     nsetCpType)
    KN_CONSTANT(       cpNum,      int,      getCpNum,      setCpNum,      nsetCpNum)
    KN_CONSTANT(       branchSW,   BranchSW, getBranchSW,   setBranchSW,   nsetBranchSW)
    KN_ARRAY_CONSTANT( parxType,   char,     getParxType,   setParxType,   getParxTypeSize, setParxTypeSize, 'P')
    KN_ARRAY_CONSTANT( parxNum,    int,      getParxNum,    setParxNum,    getParxNumSize,  setParxNumSize, 0)
    KN_ARRAY_CONSTANT( eqnsType,   char,     getEqnsType,   setEqnsType,   getEqnsTypeSize, setEqnsTypeSize, 'E')
    KN_ARRAY_CONSTANT( eqnsNum,    int,      getEqnsNum,    setEqnsNum,    getEqnsNumSize,  setEqnsNumSize, 0)
    KN_ARRAY_CONSTANT( varsType,   char,     getVarsType,   setVarsType,   getVarsTypeSize, setVarsTypeSize, 'P')
    KN_ARRAY_CONSTANT( varsNum,    int,      getVarsNum,    setVarsNum,    getVarsNumSize,  setVarsNumSize, 0)
    KN_CONSTANT(       nInt,       int,      getNInt,       setNInt,       nsetNInt)
    KN_CONSTANT(       nDeg,       int,      getNDeg,       setNDeg,       nsetNDeg)
    KN_CONSTANT(       nMul,       int,      getNMul,       setNMul,       nsetNMul)
    KN_CONSTANT(       stab,       bool,     getStab,       setStab,       nsetStab)
    KN_CONSTANT(       nMat,       int,      getNMat,       setNMat,       nsetNMat)
    KN_CONSTANT(       nInt1,      int,      getNInt1,      setNInt1,      nsetNInt1)
    KN_CONSTANT(       nInt2,      int,      getNInt2,      setNInt2,      nsetNInt2)
    KN_CONSTANT(       nDeg1,      int,      getNDeg1,      setNDeg1,      nsetNDeg1)
    KN_CONSTANT(       nDeg2,      int,      getNDeg2,      setNDeg2,      nsetNDeg2)
    KN_CONSTANT(       steps,      int,      getSteps,      setSteps,      nsetSteps)
    KN_CONSTANT(       iad,        int,      getIad,        setIad,        nsetIad)
    KN_CONSTANT(       nPr,        int,      getNPr,        setNPr,        nsetNPr)
    KN_CONSTANT(       cpMin,      double,   getCpMin,      setCpMin,      nsetCpMin)
    KN_CONSTANT(       cpMax,      double,   getCpMax,      setCpMax,      nsetCpMax)
    KN_CONSTANT(       ds,         double,   getDs,         setDs,         nsetDs)
    KN_CONSTANT(       dsMin,      double,   getDsMin,      setDsMin,      nsetDsMin)
    KN_CONSTANT(       dsMax,      double,   getDsMax,      setDsMax,      nsetDsMax)
    KN_CONSTANT(       dsStart,    double,   getDsStart,    setDsStart,    nsetDsStart)
    KN_CONSTANT(       epsC,       double,   getEpsC,       setEpsC,       nsetEpsC)
    KN_CONSTANT(       epsR,       double,   getEpsR,       setEpsR,       nsetEpsR)
    KN_CONSTANT(       epsK,       double,   getEpsK,       setEpsK,       nsetEpsK)
    KN_CONSTANT(       nItC,       int,      getNItC,       setNItC,       nsetNItC)
    KN_CONSTANT(       nItR,       int,      getNItR,       setNItR,       nsetNItR)
    KN_CONSTANT(       nItK,       int,      getNItK,       setNItK,       nsetNItK)
    KN_ARRAY_CONSTANT( symRe,      int,      getSymRe,      setSymRe,      getSymReSize, setSymReSize, 0)
    KN_ARRAY_CONSTANT( symIm,      int,      getSymIm,      setSymIm,      getSymImSize, setSymImSize, 0)
    // from sysname
    KN_CONSTANT(       nPar,       int,      getNPar,       setNPar,       nsetNPar)
    KN_CONSTANT(       nDim,       int,      getNDim,       setNDim,       nsetNDim)
    
  protected:
    std::vector<std::string> parNames;
 
  public:
    // for double, int, char
    // for i in `g++ -E -DEXTRACT_NAMES constants.h | grep -e "private:\ *double" -e "private:\ *int" -e "private:\ *char"| sed -e s/\;.*//g -e s/\ *private:\ *double//g -e s/\ *private:\ *int//g -e s/\ *private:\ *char//g`; do echo -n $i\(0\),\ ; done;
    NConstants() : label(0), cpType(0), cpNum(0), nInt(0), nDeg(0), nMul(0), nMat(0), nInt1(0), nInt2(0), nDeg1(0), nDeg2(0), steps(0), iad(0), nPr(0), cpMin(0), cpMax(0), ds(0), dsMin(0), dsMax(0), dsStart(0), epsC(0), epsR(0), epsK(0), nItC(0), nItR(0), nItK(0), nPar(0), nDim(0)
    {}
    // to extract the dependencies run:
    // for i in `g++ -E -DEXTRACT_NAMES -DREMOVE_TYPES constants.h | grep private | sed -e s/\;.*//g -e s/.*private:\ //g`; do echo -n $i\(src.$i\),\ ; done;
    NConstants(const NConstants& src) : inputFile(src.inputFile), outputFile(src.outputFile), sysname(src.sysname), label(src.label), pointType(src.pointType), cpType(src.cpType), cpNum(src.cpNum), branchSW(src.branchSW), parxType(src.parxType), parxNum(src.parxNum), eqnsType(src.eqnsType), eqnsNum(src.eqnsNum), varsType(src.varsType), varsNum(src.varsNum), nInt(src.nInt), nDeg(src.nDeg), nMul(src.nMul), stab(src.stab), nMat(src.nMat), nInt1(src.nInt1), nInt2(src.nInt2), nDeg1(src.nDeg1), nDeg2(src.nDeg2), steps(src.steps), iad(src.iad), nPr(src.nPr), cpMin(src.cpMin), cpMax(src.cpMax), ds(src.ds), dsMin(src.dsMin), dsMax(src.dsMax), dsStart(src.dsStart), epsC(src.epsC), epsR(src.epsR), epsK(src.epsK), nItC(src.nItC), nItR(src.nItR), nItK(src.nItK), symRe(src.symRe), symIm(src.symIm), nPar(src.nPar), nDim(src.nDim) {}
    ~NConstants() { initParNames(); }
    // for setting up the continuation
    bool toEqnVar(System& sys,
                  Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                  Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                  Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, bool& findangle);
    void loadXmlFile(const std::string &fileName);
    void saveXmlFile(const std::string &fileName);
    void printXmlFile(std::ostream& file);
    // this is to notify the system of changes
    virtual void constantChanged(const char* name) { }
	
	// converting equation numbers to what can be used in continuation
    Eqn eqnFromTypeNum(char type, int num) const
    {
      if (type == 'E') return (Eqn)num;
      else return EqnNone;	
    }
    // converting variable numbers to what can be used in continuation
    Var varFromTypeNum(char type, int num) const
    {
      if (type == 'S') return(Var)num;
      else if (type == 'P') return(Var)(num + VarPAR0);
      else if (type == 'I') return(Var)(num + VarPAR0 + getNPar());
      return VarNone;
    }
    // Convert type, number into enum that can be used in the computation
    // OR MAYBE MAKE THE COMPUTATION WORK WITH THESE DIRECTLY?
    Var  getCp() const { return varFromTypeNum(getCpType(),getCpNum()); }
    Var  getParx(int i) const 
    {
      if (i<getParxNumSize()) return varFromTypeNum(getParxType(i), getParxNum(i)); 
      else return VarNone;
    }
    Var  getVars(int i) const 
    {
      if (i<getVarsNumSize()) return varFromTypeNum(getVarsType(i), getVarsNum(i));
      else return VarNone;
    }
    Eqn  getEqns(int i) const 
    {
      if (i<getEqnsNumSize()) return eqnFromTypeNum(getEqnsType(i), getEqnsNum(i)); 
      else return EqnNone;
    }
    
    void initParNames()
    {
      parNames.resize(getNPar()+ParEnd);
      for (int i=0; i<getNPar(); ++i)
      {
        char buf[16];
        snprintf(buf, 15, "%d", i);
        parNames[i] = "Par. " + std::string(buf);
      }
      parNames[getNPar()+ParPeriod] = "Per. Mul.";
      parNames[getNPar()+ParAngle] = "Angle (NS)";
      parNames[getNPar()+ParRot] = "Rotation num.";      
    }
    void initDimensions(const System* sys)
    {
      setNPar(sys->npar());
      setNDim(sys->ndim());
      initParNames();
      const char **names = new const char *[sys->npar()];
      sys->parnames(names);
      for (int i=0; i<sys->npar(); ++i) parNames[i] = names[i];
      delete[] names;
    }
    // This loads the shared object file
    virtual void setSysNameText(const std::string& str)
    {
      setSysName(str);
      System* sys = 0;
      try
      {
        sys = new System(sysname);
        initDimensions(sys);
        delete sys;
      }
      catch (knutException ex)
      {
      	delete sys;
        setNPar(0);
        setNDim(0);
        throw(ex);
      }
    }
    const std::vector<std::string>& getParNames() const { return parNames; }
};

#endif
