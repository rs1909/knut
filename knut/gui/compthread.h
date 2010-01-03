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
#include <QSharedPointer>

class KNDataFile;

class MThread : public QThread, public KNAbstractContinuation
{
    Q_OBJECT
  public:
    MThread(const KNConstants& constants, QObject* parent = 0);
    ~MThread();
    void run();
    void print(std::ostringstream& str);
    void raiseException(const KNException& ex);
    void setData(KNDataFile* data);
    KNDataFile& data();
    const QSharedPointer<KNDataFile>& dataPointer();
    void deleteData();
    void dataUpdated();
  public slots:
    void consolePrint(const std::string& str);

  signals:
    void exceptionOccured(const KNException& ex);
    void printToScreen(const std::string& str);
    void dataChanged(const QSharedPointer<const KNDataFile>& mat);

  private:
    QSharedPointer<KNDataFile> output;
};
