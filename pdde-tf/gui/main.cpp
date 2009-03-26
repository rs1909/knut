// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "config.h"
#include "mainwindow.h"
#include <string>
#include <QApplication>
#include <QDir>
#include <QFileInfo>

#ifdef Q_WS_MAC
#include "macopenevent.h"
#endif

// the command line can contain one xml file name.
// if -r is given then it will run the continuation and exit afterwards

QString constFile;

int main(int argc, char *argv[])
{
  // searching for the run option
  bool run =  false;
  for (int acnt = 1; acnt < argc;  acnt++)
  {
    if (argv[acnt][0] == '-')
    {
      switch (argv[acnt][1])
      {
        case 'r':
          if (argv[acnt][2] == 0) run = true;
          break;
#ifdef HAVE_CONFIG_H
        case 'v':
          if (argv[acnt][2] == 0)
          {
            std::cout << "This is " << PACKAGE_NAME 
                      << " version " << PACKAGE_VERSION 
                      << " (" << PACKAGE_REVISION << ")\n";
            exit(0);
          }
          break;
#endif
      }
    } else
    {
      constFile = argv[acnt];
    }
  }

  if (run)
  {
    // don't make it a gui application if we run it...
    QCoreApplication app(argc, argv);
    NConstantsQtGui params;
    QFileInfo cfInfo(constFile);
    QDir::setCurrent(cfInfo.absolutePath());
    
    params.loadXmlFile(cfInfo.fileName().toStdString());
    CLIComp comp(params);
    comp.run();
  } else
  {
    QApplication app(argc, argv);
    app.setWindowIcon(QIcon(":/res/images/icon-pdde-cont.png"));
    MainWindow mainWin(app.applicationDirPath());

#ifdef Q_WS_MAC
    FileOpenEvent* openEvent = new FileOpenEvent(&app);
    app.installEventFilter(openEvent);
    QObject::connect(openEvent, SIGNAL(open(const QString&)), &mainWin, SLOT(loadFile(const QString&)));
#endif

    qRegisterMetaType<pddeException>("pddeException");
    qRegisterMetaType<std::string>("std::string");
    mainWin.show();
    return app.exec();
  }
}

