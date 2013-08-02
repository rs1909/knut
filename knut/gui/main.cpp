// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "config.h"
#include "mainwindow.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <QApplication>
#include <QDir>
#include <QFileInfo>

#ifdef Q_WS_MAC
#include <cstdio>
#include <mach-o/dyld.h>
#include "macopenevent.h"
#endif

// the command line can contain one xml file name.
// if -r is given then it will run the continuation and exit afterwards

int main(int argc, char *argv[])
{
  QString constFile;
  // searching for the run option
  bool run = false;
#ifdef Q_WS_MAC
  bool printDyld = false;
#endif
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
#ifdef Q_WS_MAC
        case 'L':
          if (argv[acnt][2] == 0) printDyld = true;
          break;
#endif
      }
    } else
    {
      constFile = argv[acnt];
    }
  }
  
#ifdef Q_WS_MAC
  // gettig shared library names from 'Mac OS X interals' by A. Singh
  if (printDyld)
  {
    const char *s = 0;
    uint32_t i, image_max;
    image_max = _dyld_image_count();
    // start from 1 so that the program itself is not included
    for (i = 1; i < image_max; i++)
    {
      if ((s = _dyld_get_image_name(i)))
      {
        QFileInfo lib(s);
        // remove system libraries
        if (!lib.path().contains("/usr/lib") && !lib.path().contains("/KNSystem"))
        {
          std::cout << s << '\n';
        }
      } else
      {
        // No name
        std::cerr << "Library has no name\n";
      }
    }
    return 0;
  }
#endif

  if (run)
  {
    // don't make it a gui application if we run it...
    QCoreApplication app(argc, argv);
    NConstantsQtGui params;
    QFileInfo cfInfo(constFile);
    QDir::setCurrent(cfInfo.absolutePath());
    
    params.loadXmlFileV5(cfInfo.fileName().toStdString());
    KNCliContinuation comp(params);
//    params.printXmlFileV5(std::cout);
    comp.run();
  } else
  {
    qRegisterMetaType<KNException>("KNException");
    qRegisterMetaType<std::string>("std::string");
    qRegisterMetaType<KNDataFile*>("KNDataFile*");
    qRegisterMetaType<size_t>("size_t");
    qRegisterMetaType<KNConstants*>("KNConstants*");
    qRegisterMetaType<KNConstants*>("DataType");

    QApplication app(argc, argv);
    app.setWindowIcon(QIcon(":/res/images/icon-knut.png"));
#ifdef Q_WS_WIN
    app.setFont(QFont("Helvetica", 9));
#endif
    MainWindow mainWin(app.applicationDirPath());
    if (!constFile.isEmpty()) { mainWin.loadFile(constFile); /*mainWin.loadFile(constFile);*/ }
#ifdef Q_WS_MAC
    FileOpenEvent* openEvent = new FileOpenEvent(&app);
    app.installEventFilter(openEvent);
    QObject::connect(openEvent, SIGNAL(open(const QString&)), &mainWin, SLOT(loadFile(const QString&)));
#endif
    
    mainWin.show();
    return app.exec();
  }
}

