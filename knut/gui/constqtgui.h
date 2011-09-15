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

// IMPLEMENT THE SIGNALS AND SLOTS IN MAINWINDOW.[H,CPP]
class NConstantsQtGui : public QObject, public KNConstants
{
  Q_OBJECT
  private:
    class tuple_t {
      public:
        const char* type;
        const char* name;
        QObject*    object;
        const char* signalName;
        const char* slotName;
    };
    std::vector<tuple_t> connectionList;
  
  public:
    virtual void setSysNameText(const std::string& str, bool testing = false)
    {
      try
      {
        KNConstants::setSysNameText(str, testing);
      }
      catch(KNException ex)
      {
        if (testing) emit sendMessage(QString::fromStdString(ex.getMessage().str()));
        else emit exceptionOccured(ex);
        return;
      }
    }
    // connecting to signals and slots
    void registerCallback(const char * type, const char * name, QObject* object, const char * signal, const char * slot);
    virtual void constantChanged(const char* name);
    
  private:
    template<typename A> class callbackPtr
    {
      public:
        callbackPtr(KNConstantNames::FPTR v) : value((void (*)(KNConstants&, const A&))v) {}
        void (*value)(KNConstants&, const A&);
    };
    template<typename B> callbackPtr<B> findConnection(QObject *object, B val);
    
  private slots:
    void slotToInt(int val);
    void slotToBool(bool val);
    void slotToDouble(const QString& val);
    void slotToPtType(int val);
    void slotToEqn(int val);
    void slotToBranchSW(int val);
    void slotToBifType(int val);
    void slotToVar(int val);
    void slotToString(const QString& val);
    
  signals:
    void constantChangedSignal(const char* name);
    void exceptionOccured(const KNException&);
    void sendMessage(const QString & message);
};

#endif
