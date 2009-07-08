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

class MThread : public QThread, public BaseComp
{
    Q_OBJECT
  public:
    MThread(const NConstants& constants, QObject* parent = 0) 
    : QThread(parent), BaseComp(constants)
    { }
    ~MThread() { }
    void run() { BaseComp::run(); }
    void print(std::ostringstream& str) { emit printToScreen(str.str()); str.str(""); }
    void raiseException(const knutException& ex) { emit exceptionOccured(ex); }
  public slots:
    void consolePrint(const std::string& str) { std::cout << str; }

  signals:
    void exceptionOccured(const knutException& ex);
    void printToScreen(const std::string& str);
};
