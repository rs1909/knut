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
    char getType(int) const
    {
      return 'E';
    }
    int  getNum(unsigned int i) const
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

    char getType(int i) const
    {
      if (key[i] < VarPAR0) return 'S';
      else if (key[i] < VarPAR0 + npar) return 'P';
      else return 'I';
    }

    int  getNum(int i) const
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
    int  cpSize() { return cpMap.size(); }
    const std::string& cpString(int i) { return cpMap.string(i); }
    unsigned int getCpIdx() { return cpMap.indexof(getCp()); }
    void setCpIdx(int i)
    {
      setCpType(cpMap.getType(i));
      setCpNum(cpMap.getNum(i));
    }
    
    // for pointType
    int  pointTypeSize() { return pointTypeMap.size(); }
    const std::string& pointTypeString(int i) { return pointTypeMap.string(i); }
    unsigned int getPointTypeIdx() { return pointTypeMap.indexof(getPointType()); }
    void setPointTypeIdx(int i) { setPointType(pointTypeMap.getNum(i)); }
    
    int  branchSWSize() { return branchSWMap.size(); }
    const std::string& branchSWString(int i) { return branchSWMap.string(i); }
    unsigned int getBranchSWIdx() { return branchSWMap.indexof(getBranchSW()); }
    void setBranchSWIdx(int i) { setBranchSW(branchSWMap.getNum(i)); }

    int parxSize() { return parxMap.size(); }
    int eqnsSize() { return eqnsMap.size(); }
    int varsSize() { return varsMap.size(); }    
    const std::string& findParXString(Var i) { return parxMap.map(i); }
    const std::string& findEqnsString(Eqn i) { return eqnsMap.map(i); }
    const std::string& findVarsString(Var i) { return varsMap.map(i); }
    const std::string& parxString(int i) { return parxMap.string(i); }
    const std::string& eqnsString(int i) { return eqnsMap.string(i); }
    const std::string& varsString(int i) { return varsMap.string(i); } 
    
    virtual void setParxIdx(int i, int p)
    {
      setParxType(i, parxMap.getType(p));
      setParxNum(i, parxMap.getNum(p));
    }
    virtual void setEqnsIdx(int i, int e)
    {
      setEqnsType(i, eqnsMap.getType(e));
      setEqnsNum(i, eqnsMap.getNum(e));
    }
    virtual void setVarsIdx(int i, int v)
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

//class NConstantsQtGui : public QObject, public NConstants
//{
//    Q_OBJECT
//  public:
//    NConstantsQtGui(QObject *parent = 0) : QObject(parent)
//    { }
//    NConstantsQtGui(const NConstants& ct, QObject *parent = 0) : QObject(parent), NConstants(ct)
//    { }
//
//  public slots:
//    void setInputFileText(const std::string& str)
//    {
//      NConstants::setInputFileText(str);
//      emit inputFileChanged(str);
//    }
//    void setOutputFileText(const std::string& str)
//    {
//      NConstants::setOutputFileText(str);
//      emit outputFileChanged(str);
//    }
//    void setSysNameText(const std::string& str)
//    {
//      try
//      {
//        NConstants::setSysNameText(str);
//      }
//      catch (knutException ex)
//      {
//        emit exceptionOccured(ex);
//        return;
//      }
//      emit sysnameChanged(str);
//    }
//    void setLabel(int i)
//    {
//      NConstants::setLabel(i) ;
//      emit labelChanged(NConstants::getLabel());
//    }
//    void setPointType(PtType p)
//    {
//      NConstants::setPointType(p);
//      emit pointTypeChangedIdx(getPointTypeIdx());
//    }
//    void setPointTypeIdx(int p)
//    {
//      NConstants::setPointTypeIdx(p);
//      emit pointTypeChangedIdx(getPointTypeIdx());
//    }
//    void setCp(char tp, int n)
//    {
//      NConstants::setCp(tp, n);
//      emit cpChangedIdx(getCpIdx());
//    }
//    void setCpIdx(int i)
//    {
//      NConstants::setCpIdx(i);
//      emit cpChangedIdx(getCpIdx());
//    }
//    void setBranchSW(BranchSW i)
//    {
//      NConstants::setBranchSW(i);
//      emit branchswChangedIdx(getBranchSWIdx());
//    }
//    void setBranchSWIdx(int i)
//    {
//      NConstants::setBranchSWIdx(i);
//      emit branchswChangedIdx(getBranchSWIdx());
//    }
//    void setNEqns(int i)
//    {
//      NConstants::setNEqns(i);
//      emit neqnsChanged(getNEqns());
//    }
//    void setNInt(int i)
//    {
//      NConstants::setNInt(i);
//      emit nintChanged(getNInt());
//    }
//    void setNDeg(int i)
//    {
//      NConstants::setNDeg(i);
//      emit ndegChanged(getNDeg());
//    }
//    void setNMul(int i)
//    {
//      NConstants::setNMul(i);
//      emit nmulChanged(getNMul());
//    }
//    void setStab(bool s)
//    {
//      NConstants::setStab(s);
//      emit stabChanged(getStab());
//    }
//    void setStab(int state)
//    {
//      setStab(state == Qt::Checked);
//    }
//    void setNMat(int i)
//    {
//      NConstants::setNMat(i);
//      emit nmatChanged(getNMat());
//    }
//    void setNInt1(int i)
//    {
//      NConstants::setNInt1(i);
//      emit nint1Changed(getNInt1());
//    }
//    void setNInt2(int i)
//    {
//      NConstants::setNInt2(i);
//      emit nint2Changed(getNInt2());
//    }
//    void setNDeg1(int i)
//    {
//      NConstants::setNDeg1(i);
//      emit ndeg1Changed(getNDeg1());
//    }
//    void setNDeg2(int i)
//    {
//      NConstants::setNDeg2(i);
//      emit ndeg2Changed(getNDeg2());
//    }
//    void setSteps(int i)
//    {
//      NConstants::setSteps(i);
//      emit stepsChanged(getSteps());
//    }
//    void setIad(int i)
//    {
//      NConstants::setIad(i);
//      emit iadChanged(getIad());
//    }
//    void setCpMin(double d)
//    {
//      NConstants::setCpMin(d);
//      emit cpMinChanged(QString::number(getCpMin(), 'g', 12));
//    }
//    void setCpMax(double d)
//    {
//      NConstants::setCpMax(d);
//      emit cpMaxChanged(QString::number(getCpMax(), 'g', 12));
//    }
//    void setDs(double d)
//    {
//      NConstants::setDs(d);
//      emit dsChanged(QString::number(getDs(), 'g', 12));
//    }
//    void setDsMin(double d)
//    {
//      NConstants::setDsMin(d);
//      emit dsMinChanged(QString::number(getDsMin(), 'g', 12));
//    }
//    void setDsMax(double d)
//    {
//      NConstants::setDsMax(d);
//      emit dsMaxChanged(QString::number(getDsMax(), 'g', 12));
//    }
//    void setDsStart(double d)
//    {
//      NConstants::setDsStart(d);
//      emit dsStartChanged(QString::number(getDsStart(), 'g', 12));
//    }
//    void setEpsC(double d)
//    {
//      NConstants::setEpsC(d);
//      emit epsCChanged(QString::number(getEpsC(), 'g', 12));
//    }
//    void setEpsR(double d)
//    {
//      NConstants::setEpsR(d);
//      emit epsRChanged(QString::number(getEpsR(), 'g', 12));
//    }
//    void setEpsK(double d)
//    {
//      NConstants::setEpsK(d);
//      emit epsKChanged(QString::number(getEpsK(), 'g', 12));
//    }
//    void setNItC(int i)
//    {
//      NConstants::setNItC(i);
//      emit nitCChanged(getNItC());
//    }
//    void setNItR(int i)
//    {
//      NConstants::setNItR(i);
//      emit nitRChanged(getNItR());
//    }
//    void setNItK(int i)
//    {
//      NConstants::setNItK(i);
//      emit nitKChanged(getNItK());
//    }
//    void setNSym(int i)
//    {
//      NConstants::setNSym(i);
//      emit nsymChanged(getNSym());
//    }
//
//    void setInputFileText(const QString& str)
//    {
//      setInputFileText(str.toStdString());
//    }
//    void setOutputFileText(const QString& str)
//    {
//      setOutputFileText(str.toStdString());
//    }
//    void setSysNameText(const QString& str)
//    {
//      setSysNameText(str.toStdString());
//    }
//    void setCpMin(const QString& d)
//    {
//      textCpMin = d;
//    }
//    void setCpMax(const QString& d)
//    {
//      textCpMax = d;
//    }
//    void setDs(const QString& d)
//    {
//      textDs = d;
//    }
//    void setDsMin(const QString& d)
//    {
//      textDsMin = d;
//    }
//    void setDsMax(const QString& d)
//    {
//      textDsMax = d;
//    }
//    void setDsStart(const QString& d)
//    {
//      textDsStart = d;
//    }
//    void setEpsC(const QString& d)
//    {
//      textEpsC = d;
//    }
//    void setEpsR(const QString& d)
//    {
//      textEpsR = d;
//    }
//    void setEpsK(const QString& d)
//    {
//      textEpsK = d;
//    }
//
//    void editedCpMin()
//    {
//      setCpMin(textCpMin.toDouble());
//    }
//    void editedCpMax()
//    {
//      setCpMax(textCpMax.toDouble());
//    }
//    void editedDs()
//    {
//      setDs(textDs.toDouble());
//    }
//    void editedDsMin()
//    {
//      setDsMin(textDsMin.toDouble());
//    }
//    void editedDsMax()
//    {
//      setDsMax(textDsMax.toDouble());
//    }
//    void editedDsStart()
//    {
//      setDsStart(textDsStart.toDouble());
//    }
//    void editedEpsC()
//    {
//      setEpsC(textEpsC.toDouble());
//    }
//    void editedEpsR()
//    {
//      setEpsR(textEpsR.toDouble());
//    }
//    void editedEpsK()
//    {
//      setEpsK(textEpsK.toDouble());
//    }
//
//  signals:
//    void inputFileChanged(const std::string& str);
//    void outputFileChanged(const std::string& str);
//    void sysnameChanged(const std::string& str);
//    void labelChanged(int i);
//    void pointTypeChangedIdx(int i);
//    void cpChangedIdx(int i);
//    void branchswChangedIdx(int i);
//    void neqnsChanged(int neqns_);
//    void nintChanged(int i);
//    void ndegChanged(int i);
//    void nmulChanged(int i);
//    void stabChanged(bool b);
//    void nmatChanged(int i);
//    void nint1Changed(int i);
//    void nint2Changed(int i);
//    void ndeg1Changed(int i);
//    void ndeg2Changed(int i);
//    void stepsChanged(int i);
//    void iadChanged(int i);
//    void cpMinChanged(const QString& d);
//    void cpMaxChanged(const QString& d);
//    void dsChanged(const QString& d);
//    void dsMinChanged(const QString& d);
//    void dsMaxChanged(const QString& d);
//    void dsStartChanged(const QString& d);
//    void epsCChanged(const QString& d);
//    void epsRChanged(const QString& d);
//    void epsKChanged(const QString& d);
//    void nitCChanged(int i);
//    void nitRChanged(int i);
//    void nitKChanged(int i);
//    void nsymChanged(int nsym_);
//    void exceptionOccured(const knutException& ex);
//
//  private:
//    QString textCpMin;
//    QString textCpMax;
//    QString textDs;
//    QString textDsMin;
//    QString textDsMax;
//    QString textDsStart;
//    QString textEpsC;
//    QString textEpsR;
//    QString textEpsK;
//};

#endif
