// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "mainwindow.h"
#include <string>
#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);
  app.setWindowIcon(QIcon(":/res/images/icon-pdde-cont.png"));
  MainWindow mainWin(app.applicationDirPath());
  qRegisterMetaType<pddeException>("pddeException");
  qRegisterMetaType<std::string>("std::string");
  mainWin.show();
  return app.exec();
}

