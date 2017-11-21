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

#ifdef Q_OS_OSX
#include <cstdio>
#include <mach-o/dyld.h>
#endif

#include <cfenv>

class KnutApplication : public QApplication
{
  public:
    KnutApplication(int &argc, char **argv, const QString& file)
        : QApplication(argc, argv)
    {
        mainWin = new MainWindow(this->applicationDirPath());
        mainWin -> show();
        if (!file.isEmpty()) mainWin->loadFile(file);
    }
#ifdef Q_OS_OSX
    bool event(QEvent *event) override
    {
        if (event->type() == QEvent::FileOpen) {
            QFileOpenEvent *openEvent = static_cast<QFileOpenEvent *>(event);
            mainWin->loadFile (openEvent->file());
        }

        return QApplication::event(event);
    }
    ~KnutApplication() override
    {
      delete mainWin;
    }
  signals:
    void open(const QString& file);
#endif
  MainWindow *mainWin;
};

// the command line can contain one xml file name.
// if -r is given then it will run the continuation and exit afterwards

int main(int argc, char *argv[])
{
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  QString constFile;
  // searching for the run option
  bool run = false;
#ifdef Q_OS_OSX
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
#ifdef Q_OS_OSX
        case 'L':
          if (argv[acnt][2] == 0) printDyld = true;
          break;
#endif
      }
    } else
    {
      constFile = argv[acnt];
      std::cout << constFile.toStdString() << "\n";
    }
  }
  
#ifdef Q_OS_OSX
  // gettig shared library names from 'Mac OS X interals' by A. Singh
  if (printDyld)
  {
    const char *s = nullptr;
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
  try {    
  if (run)
  {
    // don't make it a gui application if we run it...
    QCoreApplication app(argc, argv);
    NConstantsQtGui params;
    QFileInfo cfInfo(constFile);
    QDir::setCurrent(cfInfo.absolutePath());
    
    params.loadXmlFileV5(cfInfo.fileName().toStdString());
    KNCliContinuation comp;
    KNSystem *sys = new KNSystem (params.getSysName ());
//    params.printXmlFileV5(std::cout);
    KNDataFile *inputData = nullptr;
    if (params.getLabel() != 0)
    {
      inputData = new KNDataFile (params.getInputFile());
    }
    comp.run (sys, &params, inputData);
    delete inputData;
    delete sys;
  } else
  {
//    qRegisterMetaType<interprocess_exception>("interprocess_exception");
    qRegisterMetaType<KNException>("KNException");
    qRegisterMetaType<std::string>("std::string");
    qRegisterMetaType<KNDataFile*>("KNDataFile*");
    qRegisterMetaType<size_t>("size_t");
    qRegisterMetaType<KNConstants*>("KNConstants*");
    qRegisterMetaType<KNConstants*>("DataType");

    KnutApplication app(argc, argv, constFile);
    app.setWindowIcon(QIcon(":/res/images/icon-knut.png"));
#ifdef Q_WS_WIN
    app.setFont(QFont("Helvetica", 9));
#endif
    return app.exec();
  }
  } // try
  // catching the possible uncaught exceptions. There should not be any.
  catch (KNException& ex)
  {
    std::cerr << ex.exprStr ("");
  }
}

