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

class mat4Data;

class MThread : public QThread, public BaseComp
{
    Q_OBJECT
  public:
    MThread(const NConstants& constants, QObject* parent = 0);
    ~MThread();
    void run();
    void print(std::ostringstream& str);
    void raiseException(const knutException& ex);
    void setData(mat4Data* data);
    mat4Data& data();
    const QSharedPointer<mat4Data>& dataPointer();
    void deleteData();
    void dataUpdated();
  public slots:
    void consolePrint(const std::string& str);

  signals:
    void exceptionOccured(const knutException& ex);
    void printToScreen(const std::string& str);
    void dataChanged(const QSharedPointer<const mat4Data>& mat);

  private:
    QSharedPointer<mat4Data> output;
};
