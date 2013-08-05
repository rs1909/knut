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
#include "plotdata.h"
#include "plotwindow.h"
#include "paramview.h"
#include "constqtgui.h"

#include <fstream>
#include <cfloat>

#include <QtGui/QtGui>
#include <QLabel>
#include <QAction>
#include <QToolButton>
#include <QHBoxLayout>
#include <QToolBar>
#include <QMenu>
#include <QFileDialog>
#include <QMenuBar>
#include <QStatusBar>

MainWindow::MainWindow(const QString& appDir) :
    executableDir(appDir),
    terminalTextSize(0),
    inputPlotWindow(0), outputPlotWindow(0),
    terminalDialog(0), // compilerProcess(0),
    inputMatFile(0), outputMatFile(0), inThreadMatFile(0)
{
#if WIN32
  executableDir.replace('/','\\');
#endif
  // QTabWidget
  // a) files + equations b) numerics c) symmetry d) torus
  QTabWidget* tabWidget = new QTabWidget();
  
  // the container widgets
  QWidget* systemWidget = new QWidget();
  QWidget* numericsWidget = new QWidget();
  QWidget* symmetryWidget = new QWidget();
  QWidget* torusWidget = new QWidget();
  
  // the layout widgets
  QGridLayout* systemGrid = new QGridLayout();
  QGridLayout* numericsGrid = new QGridLayout();
  QGridLayout* symmetryGrid = new QGridLayout();
  QGridLayout* torusGrid = new QGridLayout();

  tabWidget->addTab(systemWidget, QString("System"));
  tabWidget->addTab(numericsWidget, QString("Numerics"));
  tabWidget->addTab(symmetryWidget, QString("Symmetry"));
  tabWidget->addTab(torusWidget, QString("Torus"));

  systemWidget->setLayout(systemGrid);
  numericsWidget->setLayout(numericsGrid);
  symmetryWidget->setLayout(symmetryGrid);
  torusWidget->setLayout(torusGrid);

  setCentralWidget(tabWidget);

  // Lable for the input file
  QHBoxLayout *getInputFileLayout = new QHBoxLayout;
  QLabel* inputFileLabel = new QLabel("INPUT");
  inputFileLabel->setToolTip(QString("Input file which contains the starting point"));
  QAction* inputFilePlotAct = new QAction(QIcon(":/res/images/cr16-action-pencil.png"), tr("&Plot"), this);
  inputFile = new QLineEdit();
  QAction* inputFileAct = new QAction(QIcon(":/res/images/cr16-action-fileopen.png"), tr("&Browse..."), this);
  QToolButton* getInputFile = new QToolButton();
  QToolButton* getInputFilePlot = new QToolButton();
  getInputFile->setDefaultAction(inputFileAct);
  getInputFilePlot->setDefaultAction(inputFilePlotAct);
  getInputFileLayout->addWidget(getInputFile);
  getInputFileLayout->addWidget(getInputFilePlot);
  getInputFileLayout->addStretch();
  systemGrid->addWidget(inputFileLabel, 0, 0, Qt::AlignLeft | Qt::AlignVCenter);
  systemGrid->addWidget(inputFile, 0, 1, 1, 3);
  systemGrid->addLayout(getInputFileLayout, 0, 4);
  connect(inputFileAct, SIGNAL(triggered()), this, SLOT(setInputFile())); // opening an input file
  connect(inputFilePlotAct, SIGNAL(triggered()), this, SLOT(inputPlot())); // plotting the input file
  parameters.registerCallback("std::string", "inputFile", inputFile, SIGNAL(textEdited(const QString &)), "setText");

  QHBoxLayout *getOutputFileLayout = new QHBoxLayout;
  QLabel* outputFileLabel = new QLabel("OUTPUT");
  outputFileLabel->setToolTip(QString("Output file, which contains the result of computation"));
  outputFile = new QLineEdit();
  QAction* outputFileAct = new QAction(QIcon(":/res/images/cr16-action-fileopen.png"), tr("&Browse..."), this);
  QAction* outputFilePlotAct = new QAction(QIcon(":/res/images/cr16-action-pencil.png"), tr("&Plot"), this);
  QToolButton* getOutputFile = new QToolButton();
  QToolButton* getOutputFilePlot = new QToolButton();
  getOutputFile->setDefaultAction(outputFileAct);
  getOutputFilePlot->setDefaultAction(outputFilePlotAct);
  getOutputFileLayout->addWidget(getOutputFile);
  getOutputFileLayout->addWidget(getOutputFilePlot);
  getOutputFileLayout->addStretch();
  systemGrid->addWidget(outputFileLabel, 1, 0, Qt::AlignLeft | Qt::AlignVCenter);
  systemGrid->addWidget(outputFile, 1, 1, 1, 3);
  systemGrid->addLayout(getOutputFileLayout, 1, 4);
  connect(outputFileAct, SIGNAL(triggered()), this, SLOT(setOutputFile()));
  connect(outputFilePlotAct, SIGNAL(triggered()), this, SLOT(outputPlot()));
  parameters.registerCallback("std::string", "outputFile", outputFile, SIGNAL(textEdited(const QString &)), "setText");

  // this only for setting SYSNAME
  QHBoxLayout *sysnameLayout = new QHBoxLayout;
  QLabel* sysnameLabel = new QLabel("SYSNAME");
  sysnameLabel->setToolTip(QString("Vector field file, e.g., \"sys-problem.vf\""));
  sysname = new QLineEdit();
  QAction* sysdefAct = new QAction(QIcon(":/res/images/cr16-action-fileopen.png"), tr("&Browse..."), this);
//   QAction* compileAct = new QAction(QIcon(":/res/images/cr22-action-compile.png"), tr("&Compile..."), this);
//   QAction* generateAct = new QAction(QIcon(":/res/images/cr22-action-generate.png"), tr("&Generate..."), this);  
  QToolButton* getSysdef = new QToolButton();
//   QToolButton* compile = new QToolButton();
//   QToolButton* generate = new QToolButton();
  getSysdef->setDefaultAction(sysdefAct);
//   compile->setDefaultAction(compileAct);
//   generate->setDefaultAction(generateAct);
  sysnameLayout->addWidget(getSysdef);
//   sysnameLayout->addWidget(compile);
//   sysnameLayout->addWidget(generate);
  sysnameLayout->addStretch();
  systemGrid->addWidget(sysnameLabel, 2, 0, Qt::AlignLeft | Qt::AlignVCenter);
  systemGrid->addWidget(sysname, 2, 1, 1, 3);
  systemGrid->addLayout(sysnameLayout, 2, 4);
  connect(sysdefAct, SIGNAL(triggered()), this, SLOT(setSysName()));
//   connect(compileAct, SIGNAL(triggered()), this, SLOT(compileSystem()));
//   connect(generateAct, SIGNAL(triggered()), this, SLOT(generateSystem()));
  // sets up a bidirectional connection
  connect(sysname, SIGNAL(editingFinished()), this, SLOT(setSysNameParameter()));
  parameters.registerCallback("std::string", "sysname", sysname, SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("size_t","nPar", this, 0, "nParChanged");

  // setting LABEL
  QLabel* labelLabel = new QLabel("LABEL");
  labelLabel->setToolTip(QString("The label of the solution in the input file to be used as a starting point."));
  label = new QSpinBox();
  label->setRange(0, 0xffff);
  parameters.registerCallback("size_t", "label", label, SIGNAL(valueChanged(int)), "setValue");

  QLabel* fromTypeLabel = new QLabel("FROM");
  QComboBox* fromType = new QComboBox();
  for (size_t i = 0; i < parameters.BifTypeTable.size(); ++i)
  {
    fromType->insertItem(i, parameters.BifTypeTable.CIndexToTypeName(i).c_str());
  }
  parameters.registerCallback("BifType", "fromType", fromType, SIGNAL(currentIndexChanged(int)), "setCurrentIndex");
  
  QLabel* pttypeLabel = new QLabel("POINT TYPE");
  pttypeLabel->setToolTip(QString("The type of a solution to be continued."));
  pttype = new QComboBox();
  for (size_t i = 0; i < parameters.PtTypeTable.size(); ++i)
  {
    pttype->insertItem(i, parameters.PtTypeTable.CIndexToTypeName(i).c_str());
  }
  parameters.registerCallback("PtType", "pointType", pttype, SIGNAL(currentIndexChanged(int)), "setCurrentIndex");

  QLabel* cpLabel = new QLabel("CP");
  cpLabel->setToolTip(QString("The continuation parameter."));
  cp = new QComboBox();

  QLabel* branchswLabel = new QLabel("SWITCH");
  branchswLabel->setToolTip("Switches to another branch at the bifurcation point.");
  branchsw = new QComboBox();
  for (size_t i = 0; i < parameters.BranchSWTable.size(); ++i)
  {
    branchsw->insertItem(i, parameters.BranchSWTable.CIndexToTypeName(i).c_str());
  }
  systemGrid->addWidget(branchswLabel, 3, 4, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(branchsw, 4, 4);
  parameters.registerCallback("BranchSW", "branchSW", branchsw, SIGNAL(currentIndexChanged(int)), "setCurrentIndex");

  eqnsLabel = new QLabel("NEQNS");
  eqnsLabel->setToolTip(QString("NPARX: Number of additional parameters to be used in the continuation.\n"
                                "NEQNS: Number of equations and variables that define the system."));
  eqns = new QSpinBox();
  eqns->setRange(0, 0xffff);
  systemGrid->addWidget(fromTypeLabel, 3, 0, 1, 1, Qt::AlignLeft | Qt::AlignBottom);
  systemGrid->addWidget(fromType, 4, 0);
  systemGrid->addWidget(labelLabel, 3, 1, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(label, 4, 1);
  systemGrid->addWidget(cpLabel, 3, 3, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(pttypeLabel, 3, 2, 1, 1, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(pttype, 4, 2);
  systemGrid->addWidget(cp, 4, 3);
  systemGrid->addWidget(eqnsLabel, 5, 0, Qt::AlignLeft | Qt::AlignBottom);
  systemGrid->addWidget(eqns, 6, 0, Qt::AlignTop);

  table = new EqnVarTableView(&parameters);
  systemGrid->addWidget(table, 5, 1, 2, 4, Qt::AlignVCenter);
  // this has to reconfigure the table
  connect(eqns, SIGNAL(valueChanged(int)), this, SLOT(setNEqnsParameter(int)));
  connect(table, SIGNAL(sizeChanged(int)), eqns, SLOT(setValue(int)));
  parameters.registerCallback("Var", "cp", cp, SIGNAL(currentIndexChanged(int)), "setCurrentIndex");
  parameters.registerCallback("vector<Var>", "parx", table, 0, "dataUpdated");
  parameters.registerCallback("vector<Var>", "vars", table, 0, "dataUpdated");
  parameters.registerCallback("vector<Eqn>", "eqns", table, 0, "dataUpdated");
  parameters.registerCallback("vector<Var>", "parx", eqns, 0, "setValue");
  parameters.registerCallback("vector<Var>", "vars", eqns, 0, "setValue");
  parameters.registerCallback("vector<Eqn>", "eqns", eqns, 0, "setValue");
  connect(pttype, SIGNAL(currentIndexChanged(int)), table, SLOT(dataUpdated(int)));
  
  QDoubleValidator* dbValid = new QDoubleValidator(-DBL_MAX, DBL_MAX, 16, this);

  // setting NINT, NDEG, NMUL, STAB
  QLabel* nintLabel = new QLabel("NINT");
  QLabel* ndegLabel = new QLabel("NDEG");
  QLabel* nmulLabel = new QLabel("NMUL");
  QLabel* stabLabel = new QLabel("STAB");
  QLabel* curveAngleLabel = new QLabel("CURVATURE");
  nintLabel->setToolTip(QString("Number of collocation intervals."));
  ndegLabel->setToolTip(QString("The degree of the piecewise polynomial."));
  nmulLabel->setToolTip(QString("Number of Floquet multipliers to be computed."));
  stabLabel->setToolTip(QString("Checked when stability information (the Floquet multipliers) is computed."));
  nint = new QSpinBox();
  nint->setRange(0, 0xffff);
  ndeg = new QSpinBox();
  ndeg->setRange(2, 24);
  nmul = new QSpinBox();
  nmul->setRange(0, 0xffff);
  stab = new QCheckBox();
  curveAngle = new QLineEdit();
  curveAngle->setValidator(dbValid);

  numericsGrid->addWidget(nintLabel, 0, 0, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(ndegLabel, 0, 1, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nmulLabel, 0, 2, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(stabLabel, 0, 3, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(curveAngleLabel, 0, 4, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nint, 1, 0);
  numericsGrid->addWidget(ndeg, 1, 1);
  numericsGrid->addWidget(nmul, 1, 2);
  numericsGrid->addWidget(stab, 1, 3, Qt::AlignHCenter);
  numericsGrid->addWidget(curveAngle, 1, 4);
  
  parameters.registerCallback("size_t",    "nInt",   nint,       SIGNAL(valueChanged(int)),           "setValue");
  parameters.registerCallback("size_t",    "nDeg",   ndeg,       SIGNAL(valueChanged(int)),           "setValue");
  parameters.registerCallback("size_t",    "nMul",   nmul,       SIGNAL(valueChanged(int)),           "setValue");
  parameters.registerCallback("double", "cAngle", curveAngle, SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("bool",   "stab",   stab,       SIGNAL(clicked(bool)),               "setChecked");

  QLabel* nint1Label = new QLabel("NINT1");
  QLabel* nint2Label = new QLabel("NINT2");
  QLabel* ndeg1Label = new QLabel("NDEG1");
  QLabel* ndeg2Label = new QLabel("NDEG2");
  nint1 = new QSpinBox();
  nint1->setRange(0, 0xffff);
  nint2 = new QSpinBox();
  nint1->setRange(0, 0xffff);
  ndeg1 = new QSpinBox();
  ndeg1->setRange(0, 0xffff);
  ndeg2 = new QSpinBox();
  ndeg2->setRange(0, 0xffff);
  torusGrid->addWidget(nint1Label, 0, 0, Qt::AlignHCenter | Qt::AlignBottom);
  torusGrid->addWidget(nint2Label, 0, 1, Qt::AlignHCenter | Qt::AlignBottom);
  torusGrid->addWidget(ndeg1Label, 0, 2, Qt::AlignHCenter | Qt::AlignBottom);
  torusGrid->addWidget(ndeg2Label, 0, 3, Qt::AlignHCenter | Qt::AlignBottom);
  torusGrid->addWidget(nint1, 1, 0);
  torusGrid->addWidget(nint2, 1, 1);
  torusGrid->addWidget(ndeg1, 1, 2);
  torusGrid->addWidget(ndeg2, 1, 3);
  parameters.registerCallback("size_t", "nInt1", nint1, SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "nInt2", nint2, SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "nDeg1", ndeg1, SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "nDeg2", ndeg2, SIGNAL(valueChanged(int)), "setValue");

  QLabel* stepsLabel = new QLabel("STEPS");
  QLabel* dsLabel = new QLabel("DS");
  QLabel* dsMinLabel = new QLabel("DSMIN");
  QLabel* dsMaxLabel = new QLabel("DSMAX");
  QLabel* dsStartLabel = new QLabel("DSSTART");
  stepsLabel->setToolTip(QString("Number of continuation steps to be taken."));
  dsLabel->setToolTip(QString("The default arc-length step size."));
  dsMinLabel->setToolTip(QString("The minimal arc-length step size."));
  dsMaxLabel->setToolTip(QString("The maximal arc-length step size."));
  dsStartLabel->setToolTip(QString("The amplitude of the initial guess along the critical eigenvector at branch switching."));
  steps = new QSpinBox();
  steps->setRange(0, 0xffff);
  ds = new QLineEdit();
  ds->setValidator(dbValid);
  dsMin = new QLineEdit();
  dsMin->setValidator(dbValid);
  dsMax = new QLineEdit();
  dsMax->setValidator(dbValid);
  dsStart = new QLineEdit();
  dsStart->setValidator(dbValid);
  numericsGrid->addWidget(stepsLabel, 2, 0, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(dsLabel, 2, 1, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(dsMinLabel, 2, 2, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(dsMaxLabel, 2, 3, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(dsStartLabel, 2, 4, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(steps, 3, 0);
  numericsGrid->addWidget(ds, 3, 1);
  numericsGrid->addWidget(dsMin, 3, 2);
  numericsGrid->addWidget(dsMax, 3, 3);
  numericsGrid->addWidget(dsStart, 3, 4);
  parameters.registerCallback("size_t", "steps", steps,   SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("double", "ds", ds,      SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "dsMin", dsMin,   SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "dsMax", dsMax,   SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "dsStart", dsStart, SIGNAL(textEdited(const QString &)), "setText");

  QLabel* epsCLabel = new QLabel("EPSC");
  QLabel* epsRLabel = new QLabel("EPSR");
  QLabel* epsKLabel = new QLabel("EPSK");
  QLabel* cpMinLabel = new QLabel("CPMIN");
  QLabel* cpMaxLabel = new QLabel("CPMAX");
  epsCLabel->setToolTip(QString("Convergence tolerance in continuation."));
  epsRLabel->setToolTip(QString("Convergence tolerance when converging to a solution."));
  epsKLabel->setToolTip(QString("Convergence tolerance when converging to a singular (e.g. tangent) vector."));
  cpMinLabel->setToolTip(QString("The minimal value of the continuation parameter."));
  cpMaxLabel->setToolTip(QString("The maximal value of the continuation parameter."));
  epsC = new QLineEdit();
  epsC->setValidator(dbValid);
  epsR = new QLineEdit();
  epsR->setValidator(dbValid);
  epsK = new QLineEdit();
  epsK->setValidator(dbValid);
  cpMin = new QLineEdit();
  cpMin->setValidator(dbValid);
  cpMax = new QLineEdit();
  cpMax->setValidator(dbValid);
  numericsGrid->addWidget(epsCLabel, 4, 0, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(epsRLabel, 4, 1, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(epsKLabel, 4, 2, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(cpMinLabel, 4, 3, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(cpMaxLabel, 4, 4, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(epsC, 5, 0);
  numericsGrid->addWidget(epsR, 5, 1);
  numericsGrid->addWidget(epsK, 5, 2);
  numericsGrid->addWidget(cpMin, 5, 3);
  numericsGrid->addWidget(cpMax, 5, 4);
  parameters.registerCallback("double", "epsC", epsC, SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "epsR", epsR, SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "epsK", epsK, SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "cpMin", cpMin, SIGNAL(textEdited(const QString &)), "setText");
  parameters.registerCallback("double", "cpMax", cpMax, SIGNAL(textEdited(const QString &)), "setText");

  QLabel* nitCLabel = new QLabel("NITC");
  QLabel* nitRLabel = new QLabel("NITR");
  QLabel* nitKLabel = new QLabel("NITK");
  QLabel* iadLabel = new QLabel("IAD");
  QLabel* nPrLabel = new QLabel("NPR");
  QLabel* nsymLabel = new QLabel("NSYM");
  nitCLabel->setToolTip(QString("Maximal number of iteration steps in continuation."));
  nitRLabel->setToolTip(QString("Maximal number of iteration when converging to a solution."));
  nitKLabel->setToolTip(QString("Maximal number of iteration steps when computing a singular (e.g. tangent) vector."));
  iadLabel->setToolTip(QString("Mesh adaptation in every IAD-th continuation steps."));
  nPrLabel->setToolTip(QString("Print every NPR-th continuation step."));
  nsymLabel->setToolTip(QString("Number of rotational symmetric complex dimensions."));

  nitC = new QSpinBox();
  nitC->setRange(0, 0xffff);
  nitR = new QSpinBox();
  nitR->setRange(0, 0xffff);
  nitK = new QSpinBox();
  nitK->setRange(0, 0xffff);
  iad = new QSpinBox();
  iad->setRange(0, 0xffff);
  nPr = new QSpinBox();
  nPr->setRange(1, 0xffff);
  nsym = new QSpinBox();
  nsym->setRange(0, 0xffff);
  numericsGrid->addWidget(nitCLabel, 6, 0, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nitRLabel, 6, 1, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nitKLabel, 6, 2, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(iadLabel, 6, 3, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nPrLabel, 6, 4, Qt::AlignHCenter | Qt::AlignBottom);
  symmetryGrid->addWidget(nsymLabel, 0, 0, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nitC, 7, 0);
  numericsGrid->addWidget(nitR, 7, 1);
  numericsGrid->addWidget(nitK, 7, 2);
  numericsGrid->addWidget(iad, 7, 3);
  numericsGrid->addWidget(nPr, 7, 4);
  symmetryGrid->addWidget(nsym, 2, 0);

  sym = new SYMTableView(&parameters);

  symmetryGrid->addWidget(sym, 3, 0, 2, 5);
  connect(nsym, SIGNAL(valueChanged(int)), this, SLOT(setNSymParameter(int)));
  connect(sym, SIGNAL(sizeChanged(int)), nsym, SLOT(setValue(int)));
  parameters.registerCallback("size_t", "nItC", nitC,   SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "nItR", nitR,   SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "nItK", nitK,   SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "iad", iad,   SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("size_t", "nPr", nPr,   SIGNAL(valueChanged(int)), "setValue");
  parameters.registerCallback("vector<size_t>", "symRe", sym, 0, "dataUpdated");
  parameters.registerCallback("vector<size_t>", "symIm", sym, 0, "dataUpdated");
  parameters.registerCallback("vector<size_t>", "symRe", nsym, 0, "setValue");
  parameters.registerCallback("vector<size_t>", "symIm", nsym, 0, "setValue");

  // connecting exceptions
  connect(&compThread, SIGNAL(exceptionOccured(const KNException&)), this, SLOT(externalException(const KNException&)), Qt::QueuedConnection);
  connect(&parameters, SIGNAL(exceptionOccured(const KNException&)), this, SLOT(externalException(const KNException&)), Qt::QueuedConnection);
  connect(&parameters, SIGNAL(sendMessage(const QString &)), statusBar(), SLOT(showMessage(const QString &)), Qt::QueuedConnection);
  // text output
  connect(&compThread, SIGNAL(printToScreen(const std::string&)), this, SLOT(terminalTextAppend(const std::string&)), Qt::QueuedConnection);
  connect(&compThread, SIGNAL(printStoreCursor()), this, SLOT(terminalStoreCursor()), Qt::QueuedConnection);
  connect(&compThread, SIGNAL(printClearLastLine()), this, SLOT(terminalClearLastLine()), Qt::QueuedConnection);

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  readSettings();
  setCurrentFile("");

  // compThread is a permanent object
  // should it be dynamic?
  workThread = new QThread;
  compThread.moveToThread(workThread);
  connect(workThread, SIGNAL(started()), &compThread, SLOT(process()), Qt::QueuedConnection);
  connect(&compThread, SIGNAL(finished()), workThread, SLOT(quit()), Qt::QueuedConnection);
  connect(workThread, SIGNAL(finished()), this, SLOT(stopped()), Qt::QueuedConnection);
  // opening the file
  connect(&compThread, SIGNAL(createDataRequest (const std::string&, DataType, size_t, size_t, KNConstants*)), 
          this,           SLOT(createThreadData (const std::string&, DataType, size_t, size_t, KNConstants*)), Qt::QueuedConnection);
//  connect(&compThread, SIGNAL(createDataRequestTR (const std::string&, size_t, size_t, KNConstants*)), 
//          this,           SLOT(createThreadDataTR (const std::string&, size_t, size_t, KNConstants*)), Qt::QueuedConnection);
  connect(this, SIGNAL(threadDataCreated(KNDataFile*)), &compThread, SLOT(dataCreated(KNDataFile*)), Qt::QueuedConnection);
  // closing the file
  connect(&compThread, SIGNAL(dataDeleteReq()), this, SLOT(threadDataDelete()), Qt::QueuedConnection);
  connect(this, SIGNAL(threadDataDeleteAck()), &compThread, SLOT(dataDeleteAck()), Qt::QueuedConnection);
}

void MainWindow::run()
{
  if (!workThread->isRunning())
  {
    terminalText.clear();
    setSysNameParameter(); // this is probably not updated by editingFinished()
    compThread.setConstants(parameters);
//    compThread.setStopFlag(false);
    const bool fflag = false;
    QMetaObject::invokeMethod (&compThread, "stopReq", Qt::QueuedConnection, QGenericArgument("bool", &fflag));
    workThread->start();
    stopAct->setEnabled(true);
  } else
  {
//     std::cout << "MainWindow::run running!\n";
  }
}

void MainWindow::stopped()
{
  stopAct->setEnabled(false);
//   std::cout << "MainWindow::stopped\n";
}

void MainWindow::threadDataDelete()
{
  if (inThreadMatFile == inputMatFile) inputMatFile = 0;
  if (inThreadMatFile == outputMatFile) outputMatFile = 0;
  delete inThreadMatFile;
  inThreadMatFile = 0;
  emit threadDataDeleteAck();
}

void MainWindow::createThreadData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms)
{
  if (t == DataType::LC) {
    inThreadMatFile = new KNDataFile(prms->getOutputFile(), prms->getParNames(),
                          prms->getSteps(), ndim, VarToIndex(VarEnd,npar),
                          prms->getNInt(), prms->getNDeg(), prms->getNMul());
  } else if(t == DataType::TR) {
    inThreadMatFile = new KNDataFile(prms->getOutputFile(), prms->getParNames(),
                          prms->getSteps(), ndim, VarToIndex(VarEnd,npar),
                          prms->getNInt1(), prms->getNInt2(), prms->getNDeg1(), prms->getNDeg2());
  } else if(t == DataType::ST) {
    inThreadMatFile = new KNDataFile(prms->getOutputFile(), prms->getParNames(),
                          prms->getSteps(), ndim, VarToIndex(VarEnd,npar),
                          0, 0, prms->getNMul());
  } else return;
  emit threadDataCreated (inThreadMatFile);
}

void MainWindow::stop()
{
  const bool fflag = true;
  QMetaObject::invokeMethod (&compThread, "stopReq", Qt::QueuedConnection, QGenericArgument("bool", &fflag));
//  compThread.setStopFlag(true);
}

void MainWindow::setSysName()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open system definition file",
                                                  QString(), "Vector field (*.vf);;All files (*)");
  if (!fileName.isEmpty())
  {
    sysname->setText(QDir::current().relativeFilePath(fileName));
    setSysNameParameter();
  }
}

void MainWindow::setInputFile()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open input file", QString(), "v4 MAT files (*.mat);;All files (*)");
  if (!fileName.isEmpty())
  {
    inputFile->setText(QDir::current().relativeFilePath(fileName));
    parameters.setInputFile(QDir::current().relativeFilePath(fileName).toStdString());
  }
}

void MainWindow::setOutputFile()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save output file as", QString(), "v4 MAT files (*.mat);;All files (*)");
  if (!fileName.isEmpty())
  {
    outputFile->setText(QDir::current().relativeFilePath(fileName));
    parameters.setOutputFile(QDir::current().relativeFilePath(fileName).toStdString());
  }
}

void MainWindow::setPointType()
{
  if (parameters.getPointType() == SolUser)
  {
    eqnsLabel->setText(QString("NEQNS"));
    branchsw->setEnabled(true);
  }
  else
  {
    eqnsLabel->setText(QString("NPARX"));
    branchsw->setEnabled(false);
  }
}

// --------------------------------------------------------------------------------------
void MainWindow::raisePlot (plotWindow **window, const KNDataFile** matFile, const QString& fileName)
{
//   std::cout << "raiseplot 0\n";
  if (*window == 0)
  {
//     std::cout << "raiseplot 1\n";
    // opening the data file
    openMatFile (matFile, fileName);
    if (*matFile)
    {
      *window = new plotWindow(fileName);
      connect (&compThread, SIGNAL(dataChanged(const KNDataFile*)), 
                *window, SLOT(updatePlot(const KNDataFile*)), Qt::QueuedConnection);
      connect (*window, SIGNAL(updated()),
                &compThread, SLOT(dataChangedAck()), Qt::QueuedConnection);
      connect (*window, SIGNAL(requestPlot(const QString&)), this, SLOT(plotReq(const QString&)));
      connect (*window, SIGNAL(openFile(const QString&)), this, SLOT(plotOpenFile(const QString&)));
      connect (*window, SIGNAL(windowClosed()), this, SLOT(plotDestroyed()));

      // calling the slot in a thread safe way (and not directly)
      QMetaObject::invokeMethod (&compThread, "dataChangedAck", Qt::QueuedConnection);
      (*window)->init (parameters.getCp());
      (*window)->setWindowTitle("Plot data");
      (*window)->show();
    }
  }
}

// a slot that is called back from plotWindow to open a file
// if the inThreadMatFile is the same as the input then we would need to use that
void MainWindow::openMatFile (const KNDataFile** matFile, const QString& fileName)
{
  if (fileName.isEmpty())
  {
    return;
  }
  if (*matFile != 0)
  {
    if (QFileInfo(fileName) == QFileInfo(QString::fromStdString((*matFile)->getFileName())))
    {
      return; // nothing to open
    }
    if (*matFile != inThreadMatFile)
    {
      delete *matFile;
      *matFile = 0;
    }
  }
  try {
    *matFile = new KNDataFile(fileName.toStdString());
  }
  catch (KNException ex) 
  {
    delete *matFile;
    *matFile = 0;
    externalException(ex);
  }
}

void MainWindow::inputPlot ()
{
  raisePlot (&inputPlotWindow, &inputMatFile, inputFile->text());
}

void MainWindow::outputPlot ()
{
  raisePlot (&outputPlotWindow, &outputMatFile, outputFile->text());
}

void MainWindow::plotOpenFile (const QString& fileName)
{
  if (sender() == inputPlotWindow)
  {
//     std::cout << "MainWindow::plotOpenFile inputPlotWindow\n";
    openMatFile(&inputMatFile, fileName);
    if (inputMatFile) inputPlotWindow->initPlot (inputMatFile);
  }
  if (sender() == outputPlotWindow)
  {
//     std::cout << "MainWindow::plotOpenFile outputPlotWindow\n";
    openMatFile(&outputMatFile, fileName);
    if (outputMatFile) outputPlotWindow->initPlot (outputMatFile);
  }
}

void MainWindow::plotReq (const QString& fileName)
{
  if (sender() == inputPlotWindow) 
  {
//     std::cout << "MainWindow::plotReq inputPlotWindow\n";
    if (inputMatFile == 0)
    {
      openMatFile(&inputMatFile, fileName);
      if (inputMatFile) inputPlotWindow->initPlot (inputMatFile);
    }
    inputPlotWindow->addPlot (inputMatFile);
  }
  if (sender() == outputPlotWindow)
  {
//     std::cout << "MainWindow::plotReq outputPlotWindow\n";
    if (outputMatFile == 0)
    {
      openMatFile(&outputMatFile, fileName);
      if (outputMatFile) outputPlotWindow->initPlot (outputMatFile);
    }
    outputPlotWindow->addPlot (outputMatFile);
  }
}

void MainWindow::plotDestroyed ()
{
  if (sender() == inputPlotWindow) 
  {
//     std::cout << "MainWindow::plotReq inputPlotWindow\n";
    if (inputMatFile != inThreadMatFile)
    {
      delete inputMatFile;
    }
    delete inputPlotWindow;
    inputPlotWindow = 0;
    inputMatFile = 0;
  }
  if (sender() == outputPlotWindow)
  {
//     std::cout << "MainWindow::plotReq outputPlotWindow\n";
    if (outputMatFile != inThreadMatFile)
    {
      delete outputMatFile;
    }
    delete outputPlotWindow;
    outputPlotWindow = 0;
    outputMatFile = 0;
  }
}

// --------------------------------------------------------------------------------------

void MainWindow::terminalViewDestroyed()
{
//  std::cout<<"TERM DESTROYED\n";
  delete terminalDialog;
  terminalDialog = 0;
}

void MainWindow::terminalView()
{
  if (terminalDialog == 0)
  {
    terminalDialog = new screenDialog(0);
    terminalDialog->setWindowTitle("terminal view");
    connect(&compThread, SIGNAL(printToScreen(const std::string&)), terminalDialog, SLOT(append(const std::string&)), Qt::QueuedConnection);
    connect(&compThread, SIGNAL(printClearLastLine()), terminalDialog, SLOT(clearLastLine()), Qt::QueuedConnection);
    connect(&compThread, SIGNAL(printStoreCursor()), terminalDialog, SLOT(storeCursor()), Qt::QueuedConnection);
    connect(terminalDialog, SIGNAL(windowClosed()), this, SLOT(terminalViewDestroyed()));
    terminalDialog->show();
    terminalDialog->setText(terminalText);
    terminalDialog->storeCursor();
  } else
  {
    delete terminalDialog;
    terminalDialog = 0;
  }
}

void MainWindow::setSysNameParameter()
{ 
  parameters.setSysNameText(sysname->text().toStdString(), true);
}

// void MainWindow::compileSystem()
// {
//   QString fileName = QFileDialog::getOpenFileName(this, "Open system definition", QString(), "C++ source (*.cpp)");
//   if (!fileName.isEmpty() && compilerProcess == 0)
//   {
//     QString newfile(fileName);
//     newfile.replace(QString(".cpp"), QString(".so"));
//     try {
//       KNSystem::compileSystem(fileName.toStdString(), newfile.toStdString(), executableDir.toStdString());
//       sysname->setText(QDir::current().relativeFilePath(fileName));
//       setSysNameParameter();
//     }
//     catch (KNException ex)
//     {
//       externalException(ex);
//     }
//   }
// }
// 
// void MainWindow::generateSystem()
// {
//   QString fileName = QFileDialog::getOpenFileName(this, "Open system definition", QString(), "Vector field definition (*.vf)");
//   if (!fileName.isEmpty() && compilerProcess == 0)
//   {
//     try {
//       KNSystem::generateSystem(fileName.toStdString(), SYS_TYPE_VFC, executableDir.toStdString());
//       sysname->setText(QDir::current().relativeFilePath(fileName));
//       setSysNameParameter();
//     }
//     catch (KNException ex)
//     {
//       externalException(ex);
//     }
//   }
// }

void MainWindow::closeEvent(QCloseEvent *event)
{
  if (workThread->isRunning())
  {
  	const bool fflag = true;
  	QMetaObject::invokeMethod (&compThread, "stopReq", Qt::QueuedConnection, QGenericArgument("bool", &fflag));
//    compThread.setStopFlag(true);
    workThread->wait();
  }
  delete workThread;
  workThread = 0;
  writeSettings();
  qApp->quit();
  event->accept();
}

void MainWindow::newFile()
{
  setCurrentFile("");
}

void MainWindow::open()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open constants file", QString(), "Constants files (*.knut);;All files (*)");
  if (!fileName.isEmpty())
  {
    loadFile(fileName);
  }
}

bool MainWindow::save()
{
  if (curFile.isEmpty())
  {
    return saveAs();
  }
  else
  {
    return saveFile(curFile);
  }
}

bool MainWindow::saveAs()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save constants file as", curFile, "Constants files (*.knut);;All files (*)");
  if (fileName.isEmpty())
    return false;

  return saveFile(fileName);
}

void MainWindow::about()
{
  QMessageBox::about(this, tr("About Knut"),
                     QString(tr(
                               "<h3>%1: A continuation software for delay-differential equations</h3>"
                               "<p>Version %2 (%3), %4</p>"
                               "<p>Available at <a href=%5> %5</a></p>"
                               "<p>This program is free software; you can redistribute it and/or "
                               "modify it under the terms of the "
                               "<a href=http://www.gnu.org/copyleft/gpl.html>GNU General Public License</a> "
                               "as published by the Free Software Foundation; either version 2 "
                               "of the License, or (at your option) any later version.</p>"
                               "<p>This program is distributed in the hope that it will be useful, "
                               "but WITHOUT ANY WARRANTY; without even the implied warranty of "
                               "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the "
                               "GNU General Public License for more details.</p>"
                               "<p>You should have received a copy of the GNU General Public License "
                               "along with this program; if not, write to the Free Software "
                               "Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.</p>"
                             )).arg(PACKAGE_NAME)
                             .arg(PACKAGE_VERSION)
                             .arg(PACKAGE_REVISION)
                             .arg(PACKAGE_COPYRIGHT)
                             .arg(PACKAGE_URL));
}

void MainWindow::createActions()
{
  runAct = new QAction(QIcon(":/res/images/cr22-action-launch.png"), tr("&Run"), this);
  runAct->setShortcut(tr("Ctrl+R"));
  runAct->setStatusTip(tr("Start the computation"));
  connect(runAct, SIGNAL(triggered()), this, SLOT(run()));

  stopAct = new QAction(QIcon(":/res/images/cr22-action-stop.png"), tr("&Stop"), this);
  stopAct->setShortcut(tr("Ctrl+T"));
  stopAct->setStatusTip(tr("Stop the computation"));
  stopAct->setEnabled(false);
  connect(stopAct, SIGNAL(triggered()), this, SLOT(stop()));

  terminalAct = new QAction(QIcon(":/res/images/cr22-action-view_text.png"), tr("&Text"), this);
  terminalAct->setShortcut(tr("Ctrl+V"));
  terminalAct->setStatusTip(tr("View textual output"));
  connect(terminalAct, SIGNAL(triggered()), this, SLOT(terminalView()));

  openAct = new QAction(QIcon(":/res/images/cr22-action-fileopen.png"), tr("&Open..."), this);
  openAct->setShortcut(tr("Ctrl+O"));
  openAct->setStatusTip(tr("Open an existing file"));
  connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

  saveAct = new QAction(QIcon(":/res/images/cr22-action-filesave.png"), tr("&Save"), this);
  saveAct->setShortcut(tr("Ctrl+S"));
  saveAct->setStatusTip(tr("Save the document to disk"));
  connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

  saveAsAct = new QAction(QIcon(":/res/images/cr22-action-filesaveas.png"), tr("Save &As..."), this);
  saveAsAct->setShortcut(tr("Ctrl+Shift+S"));
  saveAsAct->setStatusTip(tr("Save the document under a new name"));
  connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));

  exitAct = new QAction(QIcon(":/res/images/cr22-action-exit.png"), tr("E&xit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  exitAct->setStatusTip(tr("Exit the application"));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

  aboutAct = new QAction(QIcon(":/res/images/icon-knut.png"), tr("&About"), this);
  aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

  aboutQtAct = new QAction(tr("About &Qt"), this);
  aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
  connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(openAct);
  fileMenu->addAction(saveAct);
  fileMenu->addAction(saveAsAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);

  menuBar()->addSeparator();

  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutAct);
  helpMenu->addAction(aboutQtAct);
}

void MainWindow::createToolBars()
{
  fileToolBar = addToolBar(tr("File"));
  fileToolBar->addAction(openAct);
  fileToolBar->addAction(saveAct);
  fileToolBar->addAction(saveAsAct);
  fileToolBar->addAction(runAct);
  fileToolBar->addAction(stopAct);
  fileToolBar->addAction(terminalAct);
}

void MainWindow::createStatusBar()
{
  statusBar()->showMessage(tr("Ready"));
}

void MainWindow::readSettings()
{
  QSettings settings("Knut", "main window");
  QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
  QSize size = settings.value("size", QSize(400, 400)).toSize();
  resize(size);
  move(pos);
}

void MainWindow::writeSettings()
{
  QSettings settings("Knut", "main window");
  settings.setValue("pos", pos());
  settings.setValue("size", size());
}

void MainWindow::loadFile(const QString &fileName)
{
//   bool success = true;
  try
  {
    QDir::setCurrent(QFileInfo(fileName).absolutePath());
    parameters.loadXmlFileV5(fileName.toStdString());
  }
  catch (KNException ex)
  {
    externalException(ex);
  }
  setCurrentFile(fileName);
  statusBar()->showMessage(tr("File loaded"), 2000);
}

bool MainWindow::saveFile(const QString &fileName)
{
  parameters.saveXmlFileV5(fileName.toStdString());
  setCurrentFile(fileName);
  statusBar()->showMessage(tr("File saved"), 2000);

  return true;
}

void MainWindow::setCurrentFile(const QString &fileName)
{
  curFile = fileName;
  setWindowModified(false);

  QString shownName;
  if (curFile.isEmpty())
    shownName = "untitled";
  else
    shownName = strippedName(curFile);

  setWindowTitle(tr("%1[*] - %2").arg(shownName).arg(tr("Knut")));
}

QString MainWindow::strippedName(const QString &fullFileName)
{
  return QFileInfo(fullFileName).fileName();
}

void MainWindow::nParChanged(int n)
{
  cp->blockSignals(true);
  cp->clear();
  for (size_t i = 0; i < parameters.VarTable.size(); ++i)
  {
    cp->insertItem(i, parameters.VarTable.CIndexToTypeName(i).c_str());
  }
  cp->setCurrentIndex(parameters.VarTable.TypeToCIndex(parameters.getCp()));
  cp->blockSignals(false);
}
