// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "compthread.h"
#include "paramview.h"
#include "screendialog.h"
#include "system.h"

#include <QMainWindow>
#include <QCheckBox>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QMessageBox>

#include <QObject>

class QAction;
class QMenu;
class QTextEdit;
class QVBoxLayout;
class QGridLayout;
class QPushButton;
class QLabel;
class QProcess;
class KNDataFile;
class plotWindow;

class EqnVarTableView;

class MainWindow : public QMainWindow
{
    Q_OBJECT

  public:
    MainWindow(const QString& appDir);
        
    static void showException(QWidget* parent, const KNException& ex)
    {
      QMessageBox::critical(parent, "Critical error", 
        QString("%1\nThis has occurred in file '%2' at line %3.")
          .arg(ex.getMessage().str().c_str())
          .arg(ex.getFile().c_str()).arg(ex.getLine()), 
        QMessageBox::Ok, 0, 0);
    }

  public slots:
    // We need to be able to load a file from the outside
    // and also from an event on MacOS
    void loadFile(const QString& fileName);
    // this loads the file, not just sets the parameter
    // otherwise it would be `setSysName()'
    void setSysNameParameter();
    void setNEqnsParameter(int d)
    {
      parameters.setParxSize(d);
      parameters.setEqnsSize(d);
      parameters.setVarsSize(d);
    }
    void setNSymParameter(int d) { parameters.setSymReSize(d); parameters.setSymImSize(d); }
    
    void externalException(const KNException& ex)
    {
      showException(this, ex);
    }

  protected:
    void closeEvent(QCloseEvent *event);

  // We need to be able to run from the outside
  public slots:
    void run();

  private slots:
    void stopped();
    void threadDataDelete();
    void threadDataSet(const KNDataFile* dataFile);
    void stop();
    void newFile();
    void open();
    bool save();
    bool saveAs();
    void about();

    void setSysName();
    void setInputFile();
    void setOutputFile();
    void setPointType();
    void inputPlot();
    void inputPlotDestroyed();
    void outputPlot();
    void outputPlotDestroyed();
    void terminalView();
    void terminalViewDestroyed();
    void terminalTextAppend(const std::string& str)
    {
      terminalText.append(QString(str.c_str()));
    }
    void terminalClearLastLine()
    {
      int sz = terminalText.size();
      int i = 0;
      do { --sz; ++i; } while(terminalText[sz] != QChar('\n'));
      terminalText.truncate(sz+1);
    }
    void compileSystem();
    void generateSystem();
    void openInputMatFile(const QString& fileName);
    void closeInputMatFile(const KNDataFile* dataFile);
    void openOutputMatFile(const QString& fileName);
    void closeOutputMatFile(const KNDataFile* dataFile);
    void nParChanged(int n);
    
  private:
    inline bool inputAssert(std::istream& is);

    void createActions();
    void createMenus();
    void createToolBars();
    void createStatusBar();
    void readSettings();
    void writeSettings();
    bool maybeSave();
    bool saveFile(const QString &fileName);
    void setCurrentFile(const QString &fileName);
    QString strippedName(const QString &fullFileName);

    // path of the executable
    QString     executableDir;

    // all the parameters
    NConstantsQtGui parameters;
    MThread      compThread;
    QWidget     *paramsWidget;
    QGridLayout *paramsGrid;
    QLabel      *eqnsLabel;

    //the textual output
    QString      terminalText;

    // these contain the constants
    QLineEdit *inputFile;
    plotWindow *inputPlotWindow;
    QLineEdit *outputFile;
    plotWindow *outputPlotWindow;
    screenDialog* terminalDialog;
    QLineEdit *sysname;
    QSpinBox  *label;
    QComboBox *pttype;
    QComboBox *cp;
    QSpinBox  *eqns;
    QComboBox *branchsw;

    QSpinBox  *nint;
    QSpinBox  *ndeg;
    QSpinBox  *nmul;
    QCheckBox *stab;
    QLineEdit *curveAngle;
    QSpinBox  *nint1;
    QSpinBox  *nint2;
    QSpinBox  *ndeg1;
    QSpinBox  *ndeg2;
    QSpinBox  *steps;
    QSpinBox  *iad;
    QSpinBox  *nPr;
    QLineEdit *cpMin;
    QLineEdit *cpMax;
    QLineEdit *ds;
    QLineEdit *dsMin;
    QLineEdit *dsMax;
    QLineEdit *dsStart;
    QLineEdit *epsC;
    QLineEdit *epsR;
    QLineEdit *epsK;
    QSpinBox  *nitC;
    QSpinBox  *nitR;
    QSpinBox  *nitK;
    SYMTableView *sym;
    QSpinBox  *nsym;

    EqnVarTableView* table;

    QString  curFile;
    QProcess *compilerProcess;

    QMenu    *fileMenu;
    QMenu    *helpMenu;
    QToolBar *fileToolBar;
    QAction *runAct;
    QAction *stopAct;
    QAction *terminalAct;
    QAction *openAct;
    QAction *saveAct;
    QAction *saveAsAct;
    QAction *exitAct;
    QAction *aboutAct;
    QAction *aboutQtAct;
    
    const KNDataFile* inputMatFile;
    const KNDataFile* outputMatFile;
    const KNDataFile* inThreadMatFile;
};

#endif
