// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef COMPTHREAD_H
#define COMPTHREAD_H

#include "constqtgui.h"
#include "basecomp.h"
#include <string>
#include <sstream>
#include <QThread>
#include <QMutex>
#include <QEventLoop>
 
class KNDataFile;

class MThread : public QObject, public KNAbstractContinuation
{
    Q_OBJECT
  public:
    MThread(const KNConstants& constants, QObject* parent = 0);
    ~MThread();
    void printStream();
    void clearLastLine() { emit printClearLastLine(); }
    void storeCursor() { emit printStoreCursor(); }
    void raiseException(const KNException& ex);
    void createData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms);
//    void createDataLC (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms);
//    void createDataTR (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms);
    KNDataFile& data();
    const KNDataFile* dataPointer();
    void deleteData();
    void dataUpdated();
  public slots:
    void stopReq(bool flag);
    void process();
    void consolePrint(const std::string& str);
    void dataDeleteAck();
    void dataChangedAck();
    void dataCreated(KNDataFile* dataFile);

  signals:
    void finished();
    void exceptionOccured(const KNException& ex);
    void printStoreCursor();
    void printClearLastLine();
    void printToScreen(const std::string& str) const;
    void dataChanged(const KNDataFile* dataFile);
    // called by deleteData(). This notifies the system that the continuation has finished.
    void dataDeleteReq();
    // these have to be connected and in return call dataCreated()
    void createDataRequest (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms);
//    void createDataRequestLC (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms);
//    void createDataRequestTR (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms);

  private:
    KNDataFile* output;
    bool changeQueued;
    static const QEventLoop::ProcessEventsFlags waitFlag;
};

#endif
