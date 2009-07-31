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
class mat4Data;
class plotWindow;

class EqnVarTableView;

class MainWindow : public QMainWindow
{
    Q_OBJECT

  public:
    MainWindow(const QString& appDir, const QString& fileName);
        
    static void showException(QWidget* parent, const knutException& ex)
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

    void setConstant(const char* name)
    {
      if (!strcmp(name,"inputFile"))
      {
        inputFile->blockSignals(true);
        int cpos = inputFile->cursorPosition();
        inputFile->setText(parameters.getInputFile().c_str());
        inputFile->setCursorPosition(cpos);
        inputFile->blockSignals(false);
      }
      else if (!strcmp(name,"outputFile"))
      {
        outputFile->blockSignals(true);
        int cpos = outputFile->cursorPosition();
        outputFile->setText(parameters.getOutputFile().c_str());
        outputFile->setCursorPosition(cpos);
        outputFile->blockSignals(false);
      }
      else if (!strcmp(name,"sysname"))
      {
        sysname->blockSignals(true);
        int cpos = sysname->cursorPosition();
        sysname->setText(parameters.getSysName().c_str());
        sysname->setCursorPosition(cpos);
        sysname->blockSignals(false);
      }
      else if (!strcmp(name,"label"))
      {
        label->blockSignals(true);
        label->setValue(parameters.getLabel());
        label->blockSignals(false);
      }
      else if (!strcmp(name,"pointType"))
      {
        pttype->blockSignals(true);
        pttype->setCurrentIndex(static_cast<int>(parameters.getPointTypeIdx()));
        pttype->blockSignals(false);
      }
      else if (!strcmp(name,"cpType"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"cpNum"))
      {
        cp->blockSignals(true);
        cp->setCurrentIndex(static_cast<int>(parameters.getCpIdx()));
        cp->blockSignals(false);
      }
      else if (!strcmp(name,"branchSW"))
      {
        branchsw->blockSignals(true);
        branchsw->setCurrentIndex(static_cast<int>(parameters.getBranchSWIdx()));
        branchsw->blockSignals(false);
      }
      else if (!strcmp(name,"parxType"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"parxNum"))
      {
        eqns->blockSignals(true);
        eqns->setValue(parameters.getParxNumSize());
        eqns->blockSignals(false);
      }
      else if (!strcmp(name,"eqnsType"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"eqnsNum"))
      {
        eqns->blockSignals(true);
        eqns->setValue(parameters.getEqnsNumSize());
        eqns->blockSignals(false);
      }
      else if (!strcmp(name,"varsType"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"varsNum"))
      {
        eqns->blockSignals(true);
        eqns->setValue(parameters.getVarsNumSize());
        eqns->blockSignals(false);
      }
      else if (!strcmp(name,"nInt"))
      {
        nint->blockSignals(true);
        nint->setValue(parameters.getNInt());
        nint->blockSignals(false);
      }
      else if (!strcmp(name,"nDeg"))
      {
        ndeg->blockSignals(true);
        ndeg->setValue(parameters.getNDeg());
        ndeg->blockSignals(false);
      }
      else if (!strcmp(name,"nMul"))
      {
        nmul->blockSignals(true);
        nmul->setValue(parameters.getNMul());
        nmul->blockSignals(false);
      }
      else if (!strcmp(name,"stab"))
      {
        stab->blockSignals(true);
        if (parameters.getStab()) stab->setCheckState(Qt::Checked);
        else stab->setCheckState(Qt::Unchecked);
        stab->blockSignals(false);
      }
      else if (!strcmp(name,"nMat"))
      {
        nmat->blockSignals(true);
        nmat->setValue(parameters.getNMat());
        nmat->blockSignals(false);
      }
      else if (!strcmp(name,"nInt1"))
      {
        nint1->blockSignals(true);
        nint1->setValue(parameters.getNInt1());
        nint1->blockSignals(false);
      }
      else if (!strcmp(name,"nInt2"))
      {
        nint2->blockSignals(true);
        nint2->setValue(parameters.getNInt2());
        nint2->blockSignals(false);
      }
      else if (!strcmp(name,"nDeg1"))
      {
        ndeg1->blockSignals(true);
        ndeg1->setValue(parameters.getNDeg1());
        ndeg1->blockSignals(false);
      }
      else if (!strcmp(name,"nDeg2"))
      {
        ndeg2->blockSignals(true);
        ndeg2->setValue(parameters.getNDeg2());
        ndeg2->blockSignals(false);
      }
      else if (!strcmp(name,"steps"))
      {
        steps->blockSignals(true);
        steps->setValue(parameters.getSteps());
        steps->blockSignals(false);
      }
      else if (!strcmp(name,"iad"))
      {
        iad->blockSignals(true);
        iad->setValue(parameters.getIad());
        iad->blockSignals(false);
      }
      else if (!strcmp(name,"nPr"))
      {
        nPr->blockSignals(true);
        nPr->setValue(parameters.getNPr());
        nPr->blockSignals(false);
      }
      else if (!strcmp(name,"cpMin"))
      {
        cpMin->blockSignals(true);
        cpMin->setText(QString::number(parameters.getCpMin(), 'g', 12));
        cpMin->blockSignals(false);
      }
      else if (!strcmp(name,"cpMax"))
      {
        cpMax->blockSignals(true);
        cpMax->setText(QString::number(parameters.getCpMax(), 'g', 12));
        cpMax->blockSignals(false);
      }
      else if (!strcmp(name,"ds"))
      {
        ds->blockSignals(true);
        ds->setText(QString::number(parameters.getDs(), 'g', 12));
        ds->blockSignals(false);
      }
      else if (!strcmp(name,"dsMin"))
      {
        dsMin->blockSignals(true);
        dsMin->setText(QString::number(parameters.getDsMin(), 'g', 12));
        dsMin->blockSignals(false);
      }
      else if (!strcmp(name,"dsMax"))
      {
        dsMax->blockSignals(true);
        dsMax->setText(QString::number(parameters.getDsMax(), 'g', 12));
        dsMax->blockSignals(false);
      }
      else if (!strcmp(name,"dsStart"))
      {
        dsStart->blockSignals(true);
        dsStart->setText(QString::number(parameters.getDsStart(), 'g', 12));
        dsStart->blockSignals(false);
      }
      else if (!strcmp(name,"epsC"))
      {
        epsC->blockSignals(true);
        epsC->setText(QString::number(parameters.getEpsC(), 'g', 12));
        epsC->blockSignals(false);
      }
      else if (!strcmp(name,"epsR"))
      {
        epsR->blockSignals(true);
        epsR->setText(QString::number(parameters.getEpsR(), 'g', 12));
        epsR->blockSignals(false);
      }
      else if (!strcmp(name,"epsK"))
      {
        epsK->blockSignals(true);
        epsK->setText(QString::number(parameters.getEpsK(), 'g', 12));
        epsK->blockSignals(false);
      }
      else if (!strcmp(name,"nItC"))
      {
        nitC->blockSignals(true);
        nitC->setValue(parameters.getNItC());
        nitC->blockSignals(false);
      }
      else if (!strcmp(name,"nItR"))
      {
        nitR->blockSignals(true);
        nitR->setValue(parameters.getNItR());
        nitR->blockSignals(false);
      }
      else if (!strcmp(name,"nItK"))
      {
        nitK->blockSignals(true);
        nitK->setValue(parameters.getNItK());
        nitK->blockSignals(false);
      }
      else if (!strcmp(name,"symRe"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"symIm"))
      {
        nsym->blockSignals(true);
        nsym->setValue(parameters.getSymImSize());
        nsym->blockSignals(false);
      }
      else if (!strcmp(name,"nPar"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"nDim"))
      {
        // Don't react!
      }
      else if (!strcmp(name,"translationMaps"))
      {
        cp->blockSignals(true);
        int idx = cp->currentIndex();
        cp->clear();
        for (unsigned int i = 0; i < parameters.cpSize(); ++i) cp->addItem(parameters.cpString(i).c_str());
        if (idx < cp->count()) cp->setCurrentIndex(idx);
        cp->blockSignals(false);
      }
      else
      {
        P_MESSAGE3("No constant `", name, "' is defined in mainwindow.h.");
      }
    }
    
    // this loads the file, not just sets the parameter
    // otherwise it would be `setSysName()'
    void setSysNameParameter() { parameters.setSysNameText(sysname->text().toStdString()); }
    void setInputFileParameter() { parameters.setInputFile(inputFile->text().toStdString()); }
    void setOutputFileParameter() { parameters.setOutputFile(outputFile->text().toStdString()); }
    void setLabelParameter(int d) { parameters.setLabel(d); }
    void setPointTypeIdxParameter(int d)
      { parameters.setPointTypeIdx(static_cast<unsigned int>(d)); }
    void setCpIdxParameter(int d)
      { parameters.setCpIdx(static_cast<unsigned int>(d)); }
    void setBranchSWIdxParameter(int d)
      { parameters.setBranchSWIdx(static_cast<unsigned int>(d)); }
    void setNEqnsParameter(int d)
    {
      parameters.setParxTypeSize(d);
      parameters.setParxNumSize(d);
      parameters.setEqnsTypeSize(d);
      parameters.setEqnsNumSize(d);
      parameters.setVarsTypeSize(d);
      parameters.setVarsNumSize(d);
    }
    void setNIntParameter(int d) { parameters.setNInt(d); }
    void setNDegParameter(int d) { parameters.setNDeg(d); }
    void setNMulParameter(int d) { parameters.setNMul(d); }
    void setStabParameter(int d) { parameters.setStab(d); }
    void setNMatParameter(int d) { parameters.setNMat(d); }
    void setNInt1Parameter(int d) { parameters.setNInt1(d); }
    void setNInt2Parameter(int d) { parameters.setNInt2(d); }
    void setNDeg1Parameter(int d) { parameters.setNDeg1(d); }
    void setNDeg2Parameter(int d) { parameters.setNDeg2(d); }
    void setStepsParameter(int d) { parameters.setSteps(d); }
    void setDsParameter() { parameters.setDs(ds->text().toDouble()); }
    void setDsMinParameter() { parameters.setDsMin(dsMin->text().toDouble()); }
    void setDsMaxParameter() { parameters.setDsMax(dsMax->text().toDouble()); }
    void setDsStartParameter() { parameters.setDsStart(dsStart->text().toDouble()); }
    void setEpsCParameter() { parameters.setEpsC(epsC->text().toDouble()); }
    void setEpsRParameter() { parameters.setEpsR(epsR->text().toDouble()); }
    void setEpsKParameter() { parameters.setEpsK(epsK->text().toDouble()); }
    void setCpMinParameter() { parameters.setCpMin(cpMin->text().toDouble()); }
    void setCpMaxParameter() { parameters.setCpMax(cpMax->text().toDouble()); }
    void setNSymParameter(int d) { parameters.setSymReSize(d); parameters.setSymImSize(d); }
    void setNItCParameter(int d) { parameters.setNItC(d); }
    void setNItRParameter(int d) { parameters.setNItR(d); }
    void setNItKParameter(int d) { parameters.setNItK(d); }
    void setIadParameter(int d) { parameters.setIad(d); }
    void setNPrParameter(int d) { parameters.setNPr(d); }    
    
    void externalException(const knutException& ex)
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
    void compileSystem();

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
    QSpinBox  *nmat;
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
};

#endif
