// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef CONSTQTGUI_H
#define CONSTQTGUI_H

#include "constants.h"
#include <string>
#include <QObject>

class NConstantsQtGui : public QObject, public NConstants
{
    Q_OBJECT
  public:
    NConstantsQtGui(QObject *parent = 0) : QObject(parent)
    { }
    NConstantsQtGui(const NConstants& ct, QObject *parent = 0) : QObject(parent), NConstants(ct)
    { }
    void loadXmlFile(const std::string &fileName);
    void saveXmlFile(const std::string &fileName);

  public slots:
    void setInputFileText(const std::string& str)
    {
      NConstants::setInputFileText(str);
      emit inputFileChanged(str);
    }
    void setOutputFileText(const std::string& str)
    {
      NConstants::setOutputFileText(str);
      emit outputFileChanged(str);
    }
    void setSysNameText(const std::string& str)
    {
      try
      {
        NConstants::setSysNameText(str);
      }
      catch (pddeException ex)
      {
        emit exceptionOccured(ex);
        return;
      }
      emit sysnameChanged(str);
    }
    void setLabel(int i)
    {
      NConstants::setLabel(i) ;
      emit labelChanged(NConstants::getLabel());
    }
    void setPointType(PtType p)
    {
      NConstants::setPointType(p);
      emit pointTypeChangedIdx(getPointTypeIdx());
    }
    void setPointTypeIdx(int p)
    {
      NConstants::setPointTypeIdx(p);
      emit pointTypeChangedIdx(getPointTypeIdx());
    }
    void setCp(Var v)
    {
      NConstants::setCp(v);
      emit cpChangedIdx(getCpIdx());
    }
    void setCp(char tp, int n)
    {
      NConstants::setCp(tp, n);
      emit cpChangedIdx(getCpIdx());
    }
    void setCpIdx(int i)
    {
      NConstants::setCpIdx(i);
      emit cpChangedIdx(getCpIdx());
    }
    void setBranchSW(BranchSW i)
    {
      NConstants::setBranchSW(i);
      emit branchswChangedIdx(getBranchSWIdx());
    }
    void setBranchSWIdx(int i)
    {
      NConstants::setBranchSWIdx(i);
      emit branchswChangedIdx(getBranchSWIdx());
    }
    void setNEqns(int i)
    {
      NConstants::setNEqns(i);
      emit neqnsChanged(getNEqns());
    }
    void setNInt(int i)
    {
      NConstants::setNInt(i);
      emit nintChanged(getNInt());
    }
    void setNDeg(int i)
    {
      NConstants::setNDeg(i);
      emit ndegChanged(getNDeg());
    }
    void setNMul(int i)
    {
      NConstants::setNMul(i);
      emit nmulChanged(getNMul());
    }
    void setStab(bool s)
    {
      NConstants::setStab(s);
      emit stabChanged(getStab());
    }
    void setStab(int state)
    {
      setStab(state == Qt::Checked);
    }
    void setNMat(int i)
    {
      NConstants::setNMat(i);
      emit nmatChanged(getNMat());
    }
    void setNInt1(int i)
    {
      NConstants::setNInt1(i);
      emit nint1Changed(getNInt1());
    }
    void setNInt2(int i)
    {
      NConstants::setNInt2(i);
      emit nint2Changed(getNInt2());
    }
    void setNDeg1(int i)
    {
      NConstants::setNDeg1(i);
      emit ndeg1Changed(getNDeg1());
    }
    void setNDeg2(int i)
    {
      NConstants::setNDeg2(i);
      emit ndeg2Changed(getNDeg2());
    }
    void setSteps(int i)
    {
      NConstants::setSteps(i);
      emit stepsChanged(getSteps());
    }
    void setIad(int i)
    {
      NConstants::setIad(i);
      emit iadChanged(getIad());
    }
    void setCpMin(double d)
    {
      NConstants::setCpMin(d);
      emit cpMinChanged(QString::number(getCpMin(), 'g', 12));
    }
    void setCpMax(double d)
    {
      NConstants::setCpMax(d);
      emit cpMaxChanged(QString::number(getCpMax(), 'g', 12));
    }
    void setDs(double d)
    {
      NConstants::setDs(d);
      emit dsChanged(QString::number(getDs(), 'g', 12));
    }
    void setDsMin(double d)
    {
      NConstants::setDsMin(d);
      emit dsMinChanged(QString::number(getDsMin(), 'g', 12));
    }
    void setDsMax(double d)
    {
      NConstants::setDsMax(d);
      emit dsMaxChanged(QString::number(getDsMax(), 'g', 12));
    }
    void setDsStart(double d)
    {
      NConstants::setDsStart(d);
      emit dsStartChanged(QString::number(getDsStart(), 'g', 12));
    }
    void setEpsC(double d)
    {
      NConstants::setEpsC(d);
      emit epsCChanged(QString::number(getEpsC(), 'g', 12));
    }
    void setEpsR(double d)
    {
      NConstants::setEpsR(d);
      emit epsRChanged(QString::number(getEpsR(), 'g', 12));
    }
    void setEpsK(double d)
    {
      NConstants::setEpsK(d);
      emit epsKChanged(QString::number(getEpsK(), 'g', 12));
    }
    void setNItC(int i)
    {
      NConstants::setNItC(i);
      emit nitCChanged(getNItC());
    }
    void setNItR(int i)
    {
      NConstants::setNItR(i);
      emit nitRChanged(getNItR());
    }
    void setNItK(int i)
    {
      NConstants::setNItK(i);
      emit nitKChanged(getNItK());
    }
    void setNSym(int i)
    {
      NConstants::setNSym(i);
      emit nsymChanged(getNSym());
    }

    void setInputFileText(const QString& str)
    {
      setInputFileText(str.toStdString());
    }
    void setOutputFileText(const QString& str)
    {
      setOutputFileText(str.toStdString());
    }
    void setSysNameText(const QString& str)
    {
      setSysNameText(str.toStdString());
    }
    void setCpMin(const QString& d)
    {
      textCpMin = d;
    }
    void setCpMax(const QString& d)
    {
      textCpMax = d;
    }
    void setDs(const QString& d)
    {
      textDs = d;
    }
    void setDsMin(const QString& d)
    {
      textDsMin = d;
    }
    void setDsMax(const QString& d)
    {
      textDsMax = d;
    }
    void setDsStart(const QString& d)
    {
      textDsStart = d;
    }
    void setEpsC(const QString& d)
    {
      textEpsC = d;
    }
    void setEpsR(const QString& d)
    {
      textEpsR = d;
    }
    void setEpsK(const QString& d)
    {
      textEpsK = d;
    }

    void editedCpMin()
    {
      setCpMin(textCpMin.toDouble());
    }
    void editedCpMax()
    {
      setCpMax(textCpMax.toDouble());
    }
    void editedDs()
    {
      setDs(textDs.toDouble());
    }
    void editedDsMin()
    {
      setDsMin(textDsMin.toDouble());
    }
    void editedDsMax()
    {
      setDsMax(textDsMax.toDouble());
    }
    void editedDsStart()
    {
      setDsStart(textDsStart.toDouble());
    }
    void editedEpsC()
    {
      setEpsC(textEpsC.toDouble());
    }
    void editedEpsR()
    {
      setEpsR(textEpsR.toDouble());
    }
    void editedEpsK()
    {
      setEpsK(textEpsK.toDouble());
    }

  signals:
    void inputFileChanged(const std::string& str);
    void outputFileChanged(const std::string& str);
    void sysnameChanged(const std::string& str);
    void labelChanged(int i);
    void pointTypeChangedIdx(int i);
    void cpChangedIdx(int i);
    void branchswChangedIdx(int i);
    void neqnsChanged(int neqns_);
    void nintChanged(int i);
    void ndegChanged(int i);
    void nmulChanged(int i);
    void stabChanged(bool b);
    void nmatChanged(int i);
    void nint1Changed(int i);
    void nint2Changed(int i);
    void ndeg1Changed(int i);
    void ndeg2Changed(int i);
    void stepsChanged(int i);
    void iadChanged(int i);
    void cpMinChanged(const QString& d);
    void cpMaxChanged(const QString& d);
    void dsChanged(const QString& d);
    void dsMinChanged(const QString& d);
    void dsMaxChanged(const QString& d);
    void dsStartChanged(const QString& d);
    void epsCChanged(const QString& d);
    void epsRChanged(const QString& d);
    void epsKChanged(const QString& d);
    void nitCChanged(int i);
    void nitRChanged(int i);
    void nitKChanged(int i);
    void nsymChanged(int nsym_);
    void exceptionOccured(const pddeException& ex);


  private:
    QString textCpMin;
    QString textCpMax;
    QString textDs;
    QString textDsMin;
    QString textDsMax;
    QString textDsStart;
    QString textEpsC;
    QString textEpsR;
    QString textEpsK;
};

#endif
