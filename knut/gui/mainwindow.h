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
#include "exprsystem.h"

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
      QString msg (ex.getMessage().str().c_str());
      msg.replace ('\n', "<br>");
      QMessageBox::critical(parent, "Critical error", 
        QString("<tt>%1</tt>\nThis has occurred in file '%2' at line %3.")
          .arg(msg)
          .arg(ex.getFile()).arg(ex.getLine()), 
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

  signals:
    void threadDataCreated(KNDataFile* dataFile);
    void threadDataDeleteAck();

  protected:
    void closeEvent(QCloseEvent *event) override;

  // We need to be able to run from the outside
  public slots:
    void run();

  private slots:
    void stopped();
    void threadDataDelete();
    void createThreadData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms);
//    void createThreadDataTR (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms);
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
    void outputPlot();
    void plotReq(const QString& fileName);
    void plotOpenFile(const QString& fileName);
    void plotDestroyed();
    void terminalView();
    void terminalViewDestroyed();
    void terminalTextAppend(const std::string& str)
    {
      terminalText.append(QString(str.c_str()));
    }
    void terminalStoreCursor()
    {
      terminalTextSize = terminalText.size();
    }
    void terminalClearLastLine()
    {
      terminalText.truncate(terminalTextSize);
    }
//     void compileSystem();
//     void generateSystem();
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
    
    void openMatFile (const KNDataFile** matFile, const QString& fileName);
    void raisePlot (plotWindow **window, const KNDataFile** matFile, const QString& fileName);

    // path of the executable
    QString     executableDir;

    // all the parameters
    NConstantsQtGui parameters;
    MThread     *compThread;
    QThread     *workThread;
    QWidget     *paramsWidget;
    QGridLayout *paramsGrid;
    QLabel      *eqnsLabel;

    //the textual output
    QString      terminalText;
    int          terminalTextSize;

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
//     QProcess *compilerProcess;

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
    KNDataFile* inThreadMatFile;
};

#endif
