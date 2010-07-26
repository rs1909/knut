// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "constqtgui.h"
#include "basecomp.h"
#include <string>
#include <sstream>
#include <QThread>
#include <QMutex>

class KNDataFile;

class MThread : public QThread, public KNAbstractContinuation
{
    Q_OBJECT
  public:
    MThread(const KNConstants& constants, QObject* parent = 0);
    ~MThread();
    void run();
    void printStream();
    void clearLastLine() { emit printClearLastLine(); }
    void raiseException(const KNException& ex);
    void setData(KNDataFile* data);
    KNDataFile& data();
    const KNDataFile* dataPointer();
    void deleteData();
    void dataUpdated();
  public slots:
    void consolePrint(const std::string& str);
    void dataDeleteAck();
    void dataChangedAck();

  signals:
    void exceptionOccured(const KNException& ex);
    void printClearLastLine();
    void printToScreen(const std::string& str) const;
    void dataChanged(const KNDataFile* dataFile);
    void dataDeleteReq();
    void dataSet(const KNDataFile* dataFile);

  private:
    KNDataFile* output;
    bool changeQueued;
    QMutex changeLock;
};
