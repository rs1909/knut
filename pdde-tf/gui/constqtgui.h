// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef CONSTQTGUI_H
#define CONSTQTGUI_H

#include <QObject>
#include "constants.h"

template<class KEY> class genMap
{
  public:
    genMap() : invalid("Invalid")
    { }
    // get the index of a key
    unsigned int indexof(const KEY& v) const
    {
      for (unsigned int i = 0; i < key.size(); ++i) if (key[i] == v) return i;
      return 0;
    }
    // get a key at a certain position
    const KEY& getkey(unsigned int i) const
    {
      if (i < key.size()) return key[i];
      else return key.at(0);
    }
    // a string from index
    const std::string& string(unsigned int i) const
    {
      if (i < key.size()) return dsc[i];
      else return invalid;
    }
    // a string from key
    const std::string& map(const KEY& v) const
    {
      return string(indexof(v));
    }
    unsigned int size() const
    {
      return key.size();
    }
  protected:
    std::vector<KEY>         key;
    std::vector<std::string> dsc;
    const std::string        invalid;
};

// for branch switching
class brswType : public genMap<BranchSW>
{
  public:
    brswType()
    {
      key.push_back(NOSwitch);
      dsc.push_back("No switch");
      key.push_back(TFBRSwitch);
      dsc.push_back("Branch");
      key.push_back(TFPDSwitch);
      dsc.push_back("Period doubling");
      key.push_back(TFHBSwitch);
      dsc.push_back("Hopf (aut)");
      key.push_back(TFTRSwitch);
      dsc.push_back("Torus");
      key.push_back(TFBRAUTSwitch);
      dsc.push_back("Branch (aut)");
      key.push_back(TFBRAUTROTSwitch);
      dsc.push_back("Branch (sym)");
    }
    BranchSW getNum(unsigned int i) const
    {
      return getkey(i);
    }
};

// for selecting the solution type
class PointType : public genMap<PtType>
{
  public:

    PointType()
    {
      dsc.push_back("User");
      key.push_back(SolUser);
      dsc.push_back("ODE limit cycle");
      key.push_back(SolODE);
      dsc.push_back("ODE limit cycle (aut)");
      key.push_back(SolAUTODE);
      dsc.push_back("Limit cycle");
      key.push_back(SolTF);
      dsc.push_back("Limit point");
      key.push_back(BifTFLP);
      dsc.push_back("Period doubling");
      key.push_back(BifTFPD);
      dsc.push_back("Neimark-Sacker");
      key.push_back(BifTFNS);
      dsc.push_back("Branch switch");
      key.push_back(SolTFBRSW);
      dsc.push_back("Period doubling switch");
      key.push_back(SolTFPDSW);
      dsc.push_back("Limit cycle (aut)");
      key.push_back(SolTFAUT);
      dsc.push_back("Limit point (aut)");
      key.push_back(BifTFAUTLP);
      dsc.push_back("Period doubling (aut)");
      key.push_back(BifTFAUTPD);
      dsc.push_back("Neimark-Sacker (aut)");
      key.push_back(BifTFAUTNS);
      dsc.push_back("Branch switch (aut)");
      key.push_back(SolTFAUTBRSW);
      dsc.push_back("Period doubling switch (aut)");
      key.push_back(SolTFAUTPDSW);
      dsc.push_back("Hopf switch (aut)");
      key.push_back(SolTFAUTHBSW);
      dsc.push_back("Torus");
      key.push_back(SolTor);
      dsc.push_back("Torus from NS");
      key.push_back(SolTorNS);
      dsc.push_back("Torus (aut)");
      key.push_back(SolAUTTor);
      dsc.push_back("Torus from NS (aut)");
      key.push_back(SolAUTTorNS);
    }
    PtType getNum(unsigned int i) const
    {
      return getkey(i);
    }
};

class EqnType : public genMap<Eqn>
{
  public:
    EqnType()
    {
      dsc.push_back("None");
      key.push_back(EqnNone);
      dsc.push_back("Limit cycle");
      key.push_back(EqnSol);
      dsc.push_back("ODE limit cycle");
      key.push_back(EqnODESol);
      dsc.push_back("PHCND (transl)");
      key.push_back(EqnPhase);
      dsc.push_back("PHCND (rot)");
      key.push_back(EqnPhaseRot);
      dsc.push_back("Torus");
      key.push_back(EqnTORSol);
      dsc.push_back("PHCND 1 (torus)");
      key.push_back(EqnTORPhase0);
      dsc.push_back("PHCND 2 (torus)");
      key.push_back(EqnTORPhase1);
      dsc.push_back("LP TF");
      key.push_back(EqnTFLP);
      dsc.push_back("PD TF");
      key.push_back(EqnTFPD);
      dsc.push_back("LP TF (aut)");
      key.push_back(EqnTFLPAUT);
      dsc.push_back("LP TF (sym)");
      key.push_back(EqnTFLPAUTROT);
      dsc.push_back("NS TF (real)");
      key.push_back(EqnTFCPLX_RE);
      dsc.push_back("NS TF (imag)");
      key.push_back(EqnTFCPLX_IM);
    }
    // these translate to Type and Num. Here it is trivial
    char getType(unsigned int) const
    {
      return 'E';
    }
    Eqn getNum(unsigned int i) const
    {
      return getkey(i);
    }
    Eqn fromTypeNum(char type, int num) const
    {
      if (type == 'E') return (Eqn)num;
      else return EqnNone;	
    }
};

class VarType : public genMap<Var>
{
  public:

    VarType(int npar_ = 0, bool full_ = true) : npar(npar_), sysvar(full_)
    {
      setPar(npar_);
    }

    void setPar(int npar_)
    {
      npar = npar_;
      key.clear();
      dsc.clear();
      dsc.push_back("None");
      key.push_back(VarNone);
      if (sysvar)
      {
        dsc.push_back("Limit cycle");
        key.push_back(VarSol);
        dsc.push_back("ODE limit cycle");
        key.push_back(VarODESol);
        dsc.push_back("Torus");
        key.push_back(VarTORSol);
      }
      for (int i = 0; i < npar; ++i)
      {
        char _buf[7+1];
        _buf[7] = '\0';
        snprintf(_buf, 7, "%d", i);
        dsc.push_back("P " + std::string(_buf));
        key.push_back((Var)(VarPAR0 + i));
      }
      dsc.push_back("Angle (NS)");
      key.push_back((Var)(VarPAR0 + npar + ParAngle));
      dsc.push_back("Rotation num.");
      key.push_back((Var)(VarPAR0 + npar + ParRot));
    }

    void setSysVar(bool b)
    {
      sysvar = b;
      setPar(npar);
    }

    char getType(unsigned int i) const
    {
      if (key[i] < VarPAR0) return 'S';
      else if (key[i] < VarPAR0 + npar) return 'P';
      else return 'I';
    }

    int  getNum(unsigned int i) const
    {
      if (key[i] < VarPAR0) return key[i];
      else if (key[i] < VarPAR0 + npar) return key[i] - VarPAR0;
      else return key[i] - VarPAR0 - npar;
    }

    Var fromTypeNum(char type, int num) const
    {
      if (type == 'S') return(Var)num;
      else if (type == 'P') return(Var)(num + VarPAR0);
      else if (type == 'I') return(Var)(num + VarPAR0 + npar);
      return VarNone;
    }
  private:
    int npar;
    bool sysvar;
};

// IMPLEMENT THE SIGNALS AND SLOTS IN MAINWINDOW.[H,CPP]
class NConstantsQtGui : public QObject, public NConstants
{
  Q_OBJECT
  private:
	// These are not necessary by the command line
    PointType          pointTypeMap;
    VarType            cpMap;
    brswType           branchSWMap;
    VarType            parxMap;
    EqnType            eqnsMap;
    VarType            varsMap;
  
  public:
    // for the continuation parameter
    unsigned int cpSize() { return cpMap.size(); }
    const std::string& cpString(unsigned int i) { return cpMap.string(i); }
    unsigned int getCpIdx() { return cpMap.indexof(getCp()); }
    void setCpIdx(unsigned int i)
    {
      setCpType(cpMap.getType(i));
      setCpNum(cpMap.getNum(i));
    }
    
    // for pointType
    unsigned int pointTypeSize() { return pointTypeMap.size(); }
    const std::string& pointTypeString(unsigned int i) { return pointTypeMap.string(i); }
    unsigned int getPointTypeIdx() { return pointTypeMap.indexof(getPointType()); }
    void setPointTypeIdx(unsigned int i) { setPointType(pointTypeMap.getNum(i)); }
    
    unsigned int branchSWSize() { return branchSWMap.size(); }
    const std::string& branchSWString(unsigned int i) { return branchSWMap.string(i); }
    unsigned int getBranchSWIdx() { return branchSWMap.indexof(getBranchSW()); }
    void setBranchSWIdx(unsigned int i) { setBranchSW(branchSWMap.getNum(i)); }

    unsigned int parxSize() { return parxMap.size(); }
    unsigned int eqnsSize() { return eqnsMap.size(); }
    unsigned int varsSize() { return varsMap.size(); }    
    const std::string& findParxString(Var i) { return parxMap.map(i); }
    const std::string& findEqnsString(Eqn i) { return eqnsMap.map(i); }
    const std::string& findVarsString(Var i) { return varsMap.map(i); }
    const std::string& parxString(unsigned int i) { return parxMap.string(i); }
    const std::string& eqnsString(unsigned int i) { return eqnsMap.string(i); }
    const std::string& varsString(unsigned int i) { return varsMap.string(i); } 
    
    virtual void setParxIdx(int i, unsigned int p)
    {
      setParxType(i, parxMap.getType(p));
      setParxNum(i, parxMap.getNum(p));
    }
    virtual void setEqnsIdx(int i, unsigned int e)
    {
      setEqnsType(i, eqnsMap.getType(e));
      setEqnsNum(i, eqnsMap.getNum(e));
    }
    virtual void setVarsIdx(int i, unsigned int v)
    {
      setVarsType(i, varsMap.getType(v));
      setVarsNum(i, varsMap.getNum(v));
    }
    virtual void setSysNameText(const std::string& str)
    {
      try
      {
        NConstants::setSysNameText(str);
      }
      catch(knutException ex)
      {
        emit exceptionOccured(ex);
      }
      cpMap.setPar(getNPar());
      parxMap.setPar(getNPar());
      varsMap.setPar(getNPar());
      emit constantChangedSignal("translationMaps");
    }
    virtual void constantChanged(const char* name) { emit constantChangedSignal(name); }
  signals:
    void constantChangedSignal(const char* name);
    void exceptionOccured(const knutException&);
};

#endif
