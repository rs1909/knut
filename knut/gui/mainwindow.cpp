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

MainWindow::MainWindow(const QString& appDir, const QString& fileName) :
    executableDir(appDir), compThread(parameters),
    inputPlotWindow(0), outputPlotWindow(0),
    terminalDialog(0), compilerProcess(0),
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
  connect(inputFile, SIGNAL(textEdited(const QString &)), this, SLOT(setInputFileParameter(const QString &)));

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
  connect(outputFile, SIGNAL(textEdited(const QString &)), this, SLOT(setOutputFileParameter(const QString &)));

  // this only for setting SYSNAME
  QHBoxLayout *sysnameLayout = new QHBoxLayout;
  QLabel* sysnameLabel = new QLabel("SYSNAME");
  sysnameLabel->setToolTip(QString("The compiled system definition, e.g., \"sys-problem.so\""));
  sysname = new QLineEdit();
  QAction* sysdefAct = new QAction(QIcon(":/res/images/cr16-action-fileopen.png"), tr("&Browse..."), this);
  QAction* compileAct = new QAction(QIcon(":/res/images/cr22-action-compile.png"), tr("&Compile..."), this);
  QAction* generateAct = new QAction(QIcon(":/res/images/cr22-action-generate.png"), tr("&Generate..."), this);  
  QToolButton* getSysdef = new QToolButton();
  QToolButton* compile = new QToolButton();
  QToolButton* generate = new QToolButton();
  getSysdef->setDefaultAction(sysdefAct);
  compile->setDefaultAction(compileAct);
  generate->setDefaultAction(generateAct);
  sysnameLayout->addWidget(getSysdef);
  sysnameLayout->addWidget(compile);
  sysnameLayout->addWidget(generate);
  sysnameLayout->addStretch();
  systemGrid->addWidget(sysnameLabel, 2, 0, Qt::AlignLeft | Qt::AlignVCenter);
  systemGrid->addWidget(sysname, 2, 1, 1, 3);
  systemGrid->addLayout(sysnameLayout, 2, 4);
  connect(sysdefAct, SIGNAL(triggered()), this, SLOT(setSysName()));
  connect(compileAct, SIGNAL(triggered()), this, SLOT(compileSystem()));
  connect(generateAct, SIGNAL(triggered()), this, SLOT(generateSystem()));
  // sets up a bidirectional connection
  connect(sysname, SIGNAL(editingFinished()), this, SLOT(setSysNameParameter()));

  // setting LABEL
  QLabel* labelLabel = new QLabel("LABEL");
  labelLabel->setToolTip(QString("The label of the solution in the input file to be used as a starting point."));
  label = new QSpinBox();
  label->setRange(0, 0xffff);
  systemGrid->addWidget(labelLabel, 3, 0, Qt::AlignLeft | Qt::AlignBottom);
  systemGrid->addWidget(label, 4, 0);
  connect(label, SIGNAL(valueChanged(int)), this, SLOT(setLabelParameter(int)));

  QLabel* pttypeLabel = new QLabel("POINT TYPE");
  pttypeLabel->setToolTip(QString("The type of a solution to be continued."));
  pttype = new QComboBox();
  for (unsigned int i = 0; i < parameters.pointTypeSize(); ++i)
  {
    pttype->addItem(parameters.pointTypeString(i).c_str());
  }
  // thes set the values
  connect(pttype, SIGNAL(currentIndexChanged(int)), this, SLOT(setPointTypeIdxParameter(int)));

  QLabel* cpLabel = new QLabel("CP");
  cpLabel->setToolTip(QString("The continuation parameter."));
  cp = new QComboBox();
  connect(cp, SIGNAL(currentIndexChanged(int)), this, SLOT(setCpIdxParameter(int)));

  QLabel* branchswLabel = new QLabel("SWITCH");
  branchswLabel->setToolTip("Switches to another branch at the bifurcation point.");
  branchsw = new QComboBox();
  for (unsigned int i = 0; i < parameters.branchSWSize(); ++i)
  {
    branchsw->addItem(parameters.branchSWString(i).c_str());
  }
  systemGrid->addWidget(branchswLabel, 3, 4, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(branchsw, 4, 4);
  connect(branchsw, SIGNAL(currentIndexChanged(int)), this, SLOT(setBranchSWIdxParameter(int)));

  eqnsLabel = new QLabel("NEQNS");
  eqnsLabel->setToolTip(QString("NPARX: Number of additional parameters to be used in the continuation.\n"
                                "NEQNS: Number of equations and variables that define the system."));
  eqns = new QSpinBox();
  eqns->setRange(0, 0xffff);
  systemGrid->addWidget(pttypeLabel, 3, 1, 1, 2, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(cpLabel, 3, 3, Qt::AlignHCenter | Qt::AlignBottom);
  systemGrid->addWidget(pttype, 4, 1, 1, 2);
  systemGrid->addWidget(cp, 4, 3);
  systemGrid->addWidget(eqnsLabel, 5, 0, Qt::AlignLeft | Qt::AlignBottom);
  systemGrid->addWidget(eqns, 6, 0, Qt::AlignTop);

  table = new EqnVarTableView(&parameters);
  systemGrid->addWidget(table, 5, 1, 2, 4, Qt::AlignVCenter);
  // this has to reconfigure the table
  connect(eqns, SIGNAL(valueChanged(int)), this, SLOT(setNEqnsParameter(int)));

  // setting NINT, NDEG, NMUL, STAB
  QLabel* nintLabel = new QLabel("NINT");
  QLabel* ndegLabel = new QLabel("NDEG");
  QLabel* nmulLabel = new QLabel("NMUL");
  QLabel* stabLabel = new QLabel("STAB");
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

  numericsGrid->addWidget(nintLabel, 0, 0, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(ndegLabel, 0, 1, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nmulLabel, 0, 2, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(stabLabel, 0, 4, Qt::AlignHCenter | Qt::AlignBottom);
  numericsGrid->addWidget(nint, 1, 0);
  numericsGrid->addWidget(ndeg, 1, 1);
  numericsGrid->addWidget(nmul, 1, 2);
  numericsGrid->addWidget(stab, 1, 4, Qt::AlignHCenter);
  connect(nint, SIGNAL(valueChanged(int)), this, SLOT(setNIntParameter(int)));
  connect(ndeg, SIGNAL(valueChanged(int)), this, SLOT(setNDegParameter(int)));
  connect(nmul, SIGNAL(valueChanged(int)), this, SLOT(setNMulParameter(int)));
  connect(stab, SIGNAL(stateChanged(int)), this, SLOT(setStabParameter(int)));

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
  connect(nint1, SIGNAL(valueChanged(int)), this, SLOT(setNInt1Parameter(int)));
  connect(nint2, SIGNAL(valueChanged(int)), this, SLOT(setNInt2Parameter(int)));
  connect(ndeg1, SIGNAL(valueChanged(int)), this, SLOT(setNDeg1Parameter(int)));
  connect(ndeg2, SIGNAL(valueChanged(int)), this, SLOT(setNDeg2Parameter(int)));

  QDoubleValidator* dbValid = new QDoubleValidator(-DBL_MAX, DBL_MAX, 16, this);

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
  connect(steps, SIGNAL(valueChanged(int)), this, SLOT(setStepsParameter(int)));
  connect(ds, SIGNAL(textEdited(const QString &)), this, SLOT(setDsParameter(const QString &)));
  connect(dsMin, SIGNAL(textEdited(const QString &)), this, SLOT(setDsMinParameter(const QString &)));
  connect(dsMax, SIGNAL(textEdited(const QString &)), this, SLOT(setDsMaxParameter(const QString &)));
  connect(dsStart, SIGNAL(textEdited(const QString &)), this, SLOT(setDsStartParameter(const QString &)));

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
  connect(epsC, SIGNAL(textEdited(const QString &)), this, SLOT(setEpsCParameter(const QString &)));
  connect(epsR, SIGNAL(textEdited(const QString &)), this, SLOT(setEpsRParameter(const QString &)));
  connect(epsK, SIGNAL(textEdited(const QString &)), this, SLOT(setEpsKParameter(const QString &)));
  connect(cpMin, SIGNAL(textEdited(const QString &)), this, SLOT(setCpMinParameter(const QString &)));
  connect(cpMax, SIGNAL(textEdited(const QString &)), this, SLOT(setCpMaxParameter(const QString &)));

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
  connect(nitC, SIGNAL(valueChanged(int)), this, SLOT(setNItCParameter(int)));
  connect(nitR, SIGNAL(valueChanged(int)), this, SLOT(setNItRParameter(int)));
  connect(nitK, SIGNAL(valueChanged(int)), this, SLOT(setNItKParameter(int)));
  connect(iad, SIGNAL(valueChanged(int)), this, SLOT(setIadParameter(int)));
  connect(nPr, SIGNAL(valueChanged(int)), this, SLOT(setNPrParameter(int)));

  // update parameters in the GUI, when `parameters' have changed.
  connect(&parameters, SIGNAL(constantChangedSignal(const char*)), this, SLOT(setConstant(const char*)));

  // connecting exceptions
  connect(&compThread, SIGNAL(exceptionOccured(const KNException&)), this, SLOT(externalException(const KNException&)));
  connect(&parameters, SIGNAL(exceptionOccured(const KNException&)), this, SLOT(externalException(const KNException&)));
  connect(&parameters, SIGNAL(sendMessage(const QString &)), statusBar(), SLOT(showMessage(const QString &)));
  // text output
  connect(&compThread, SIGNAL(printToScreen(const std::string&)), this, SLOT(terminalTextAppend(const std::string&)));
  connect(&compThread, SIGNAL(printClearLastLine()), this, SLOT(terminalClearLastLine()));

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  readSettings();
  if (fileName.isEmpty()) setCurrentFile("");
  else loadFile(fileName);

  // compThread is a permanent object
  // should it be dynamic?
  connect(&compThread, SIGNAL(finished()), this, SLOT(stopped()));
  connect(&compThread, SIGNAL(dataDeleteReq()), this, SLOT(threadDataDelete()));
  connect(&compThread, SIGNAL(dataSet(const KNDataFile*)), this, SLOT(threadDataSet(const KNDataFile*)));
}

void MainWindow::run()
{
  if (!compThread.isRunning())
  {
    terminalText.clear();
    setSysNameParameter(); // this is probably not updated by editingFinished()
    compThread.setConstants(parameters);
    compThread.setStopFlag(false);
    compThread.start();
    stopAct->setEnabled(true);
  }
}

void MainWindow::stopped()
{
  stopAct->setEnabled(false);
}

void MainWindow::threadDataDelete()
{
  // don't let it delete if the plotter is using it,
  // take away ownership in any case
  bool inToDelete = true, outToDelete = true; 
  if (inputMatFile)
    if(inputMatFile == inThreadMatFile)
      inToDelete = false;
  if (outputMatFile)
    if(outputMatFile == inThreadMatFile)
      outToDelete = false;
  if (inToDelete && outToDelete)
  {
    compThread.dataDeleteAck();
//    std::cout << " --- Deleting\n";
  } else
  {
    const std::string fn(inThreadMatFile->getFileName());
//    std::cout << " --- reopening as read only " << fn << "\n";
    compThread.dataDeleteAck();
    const KNDataFile* df = new KNDataFile(fn);
    if (!inToDelete) 
    {
      inputMatFile = df;
      if (inputPlotWindow) inputPlotWindow->setData(df);
    }
    if (!outToDelete)
    {
      outputMatFile = df;
      if (outputPlotWindow) outputPlotWindow->setData(df);
    }
//    std::cout << " --- reopened as read only " << df->getFileName() << "\n";
  }
  inThreadMatFile = 0;
}

void MainWindow::threadDataSet(const KNDataFile* dataFile)
{
  if (inThreadMatFile)
  {
//    std::cout << "One is already in thread\n";
  } else
  {
//    std::cout << "setting inThreadMatFile\n";
    inThreadMatFile = dataFile;
    // is the opened data used somewhere?
    // if yes, set it to the new, the old will be deleted
    if (inputPlotWindow)
    {
//      std::cout << "intest " << inputMatFile << "\n";
      if (inputPlotWindow->isThisData(dataFile))
      {
//        std::cout << "Same as Inp. 1\n";
        if (inputMatFile) delete inputMatFile;
        inputMatFile = dataFile;
        inputPlotWindow->setData(dataFile);
//        std::cout << "Same as Inp. 2\n";
      }
    }
    if (outputPlotWindow)
    {
//      std::cout << "outtest " << outputMatFile << "\n";
      if (outputPlotWindow->isThisData(dataFile))
      {
//        std::cout << "Same as Outp. 1\n";
        if (outputMatFile) delete outputMatFile;
        outputMatFile = dataFile;
        outputPlotWindow->setData(dataFile);
//        std::cout << "Same as Outp. 2\n";
      }
    }
  }
}

void MainWindow::stop()
{
  compThread.setStopFlag(true);
}

void MainWindow::setSysName()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open system definition file",
                                                  QString(), "Shared object (*.so);;All files (*)");
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
    setInputFileParameter(QDir::current().relativeFilePath(fileName));
  }
}

void MainWindow::setOutputFile()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save output file as", QString(), "v4 MAT files (*.mat);;All files (*)");
  if (!fileName.isEmpty())
  {
    outputFile->setText(QDir::current().relativeFilePath(fileName));
    setOutputFileParameter(QDir::current().relativeFilePath(fileName));
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

void MainWindow::inputPlot()
{
  if (inputPlotWindow == 0)
  {
    inputPlotWindow = new plotWindow(inputFile->text());
    connect(inputPlotWindow, SIGNAL(openFile(const QString&)), this, SLOT(openInputMatFile(const QString&)));
    connect(inputPlotWindow, SIGNAL(closeFile(const KNDataFile*)), this, SLOT(closeInputMatFile(const KNDataFile*)));
    connect(&compThread, SIGNAL(dataChanged(const KNDataFile*)), 
            inputPlotWindow, SLOT(updatePlot(const KNDataFile*)), Qt::QueuedConnection);
    connect(inputPlotWindow, SIGNAL(updated()),
            &compThread, SLOT(dataChangedAck()), Qt::QueuedConnection);
    compThread.dataChangedAck();
    inputPlotWindow->init();
    if (inputPlotWindow->isDataSet())
    {
      connect(inputPlotWindow, SIGNAL(windowClosed()), this, SLOT(inputPlotDestroyed()));
      inputPlotWindow->setWindowTitle("Plot - input data");
      inputPlotWindow->show();
    }
    else
    {
      inputPlotDestroyed();
    }
  }
}

// a slot that is called back from plotWindow to open a file
void MainWindow::openInputMatFile(const QString& fileName)
{
  if (fileName.isEmpty())
  {
//    std::cout << "Empty fileName\n";
    return;
  }
  if (inThreadMatFile)
  {
    if (QFileInfo(fileName) == QFileInfo(QString::fromStdString(inThreadMatFile->getFileName())))
    {
      if (inputMatFile != inThreadMatFile && inputMatFile != 0) delete inputMatFile;
      inputMatFile = inThreadMatFile;
      if (inputPlotWindow) inputPlotWindow->setData(inputMatFile);
//      std::cout << "Already open in compThread " << fileName.toStdString() << "\n";
      return;
    }
  }
  if (inputMatFile)
  {
    if (QFileInfo(fileName) != QFileInfo(QString::fromStdString(inputMatFile->getFileName())))
    {
      // if it is not used in computation -> safe to delete;
      delete inputMatFile;
      try {
        inputMatFile = new KNDataFile(fileName.toStdString());
//        std::cout << "Openining (1)" << fileName.toStdString() << "\n";
      }
      catch (KNException ex) 
      {
        externalException(ex);
      }
      if (inputPlotWindow) inputPlotWindow->setData(inputMatFile);
      return;
    } else
    {
//      std::cout << "Already open in MainWindow" << fileName.toStdString() << "\n";
    }
  } else
  {
    try {
      inputMatFile = new KNDataFile(fileName.toStdString());
//      std::cout << "Openining (2)" << fileName.toStdString() << "\n";
      if (inputPlotWindow) inputPlotWindow->setData(inputMatFile);
      return;
    }
    catch (KNException ex) 
    {
      externalException(ex);
    }
  }
//  std::cout << "should not have ended up here.\n";
}

void MainWindow::closeInputMatFile(const KNDataFile* dataFile)
{
  if (dataFile == inputMatFile)
  {
    if (inputMatFile != inThreadMatFile) delete dataFile;
    inputMatFile = 0;
  } else
  {
//    std::cerr << "Mismatch of datafiles in plotWindow and MainWindow.\n";
  }
}

void MainWindow::inputPlotDestroyed()
{
//  std::cout<<"plot destroyed\n";
  delete inputPlotWindow;
  inputPlotWindow = 0;
  inputMatFile = 0;
}

// --------------------------------------------------------------------------------------

void MainWindow::outputPlot()
{
  if (outputPlotWindow == 0)
  {
    outputPlotWindow = new plotWindow(outputFile->text());
    connect(outputPlotWindow, SIGNAL(openFile(const QString&)), this, SLOT(openOutputMatFile(const QString&)));
    connect(outputPlotWindow, SIGNAL(closeFile(const KNDataFile*)), this, SLOT(closeOutputMatFile(const KNDataFile*)));
    connect(&compThread, SIGNAL(dataChanged(const KNDataFile*)), 
            outputPlotWindow, SLOT(updatePlot(const KNDataFile*)), Qt::QueuedConnection);
    connect(outputPlotWindow, SIGNAL(updated()),
            &compThread, SLOT(dataChangedAck()), Qt::QueuedConnection);
    compThread.dataChangedAck();
    outputPlotWindow->init();
    if (outputPlotWindow->isDataSet())
    {
      connect(outputPlotWindow, SIGNAL(windowClosed()), this, SLOT(outputPlotDestroyed()));
      outputPlotWindow->setWindowTitle("Plot - output data");
      outputPlotWindow->setWindowFlags(Qt::Window);
      outputPlotWindow->show();
    }
    else
    {
      outputPlotDestroyed();
    }
  }
}

// a slot that is called back from plotWindow to open a file
void MainWindow::openOutputMatFile(const QString& fileName)
{
  if (fileName.isEmpty())
  {
//    std::cout << "Empty fileName\n";
    return;
  }
  if (inThreadMatFile)
  {
    if (QFileInfo(fileName) == QFileInfo(QString::fromStdString(inThreadMatFile->getFileName())))
    {
      if (outputMatFile != inThreadMatFile && outputMatFile != 0) delete outputMatFile;
      outputMatFile = inThreadMatFile;
      if (outputPlotWindow) outputPlotWindow->setData(outputMatFile);
//      std::cout << "Already open in compThread " << fileName.toStdString() << "\n";
      return;
    }
  }
  if (outputMatFile)
  {
    if (QFileInfo(fileName) != QFileInfo(QString::fromStdString(outputMatFile->getFileName())))
    {
      // if it is not used in computation -> safe to delete;
      delete outputMatFile;
      try {
        outputMatFile = new KNDataFile(fileName.toStdString());
//        std::cout << "Openining (1)" << fileName.toStdString() << "\n";
      }
      catch (KNException ex) 
      {
        externalException(ex);
      }
      if (outputPlotWindow) outputPlotWindow->setData(outputMatFile);
      return;
    } else
    {
//      std::cout << "Already open in MainWindow" << fileName.toStdString() << "\n";
    }
  } else
  {
    try {
      outputMatFile = new KNDataFile(fileName.toStdString());
//      std::cout << "Openining (2)" << fileName.toStdString() << "\n";
      if (outputPlotWindow) outputPlotWindow->setData(outputMatFile);
      return;
    }
    catch (KNException ex) 
    {
      externalException(ex);
    }
  }
//  std::cout << "should not have ended up here.\n";
}

void MainWindow::closeOutputMatFile(const KNDataFile* dataFile)
{
  if (dataFile == outputMatFile)
  {
    if (outputMatFile != inThreadMatFile) delete dataFile;
    outputMatFile = 0;
  } else
  {
//    std::cerr << "Mismatch of datafiles in plotWindow and MainWindow.\n";
  }
}

void MainWindow::outputPlotDestroyed()
{
//  std::cout<<"plot destroyed\n";
  delete outputPlotWindow;
  outputPlotWindow = 0;
  outputMatFile = 0;
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
    connect(&compThread, SIGNAL(printToScreen(const std::string&)), terminalDialog, SLOT(append(const std::string&)));
    connect(&compThread, SIGNAL(printClearLastLine()), terminalDialog, SLOT(clearLastLine()));
    connect(terminalDialog, SIGNAL(windowClosed()), this, SLOT(terminalViewDestroyed()));
    terminalDialog->setWindowTitle("terminal view");
    terminalDialog->setText(terminalText);
    terminalDialog->show();
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

void MainWindow::compileSystem()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open system definition", QString(), "C++ source (*.cpp)");
  if (!fileName.isEmpty() && compilerProcess == 0)
  {
    QString newfile(fileName);
    newfile.replace(QString(".cpp"), QString(".so"));
    try {
      KNSystem::compileSystem(fileName.toStdString(), newfile.toStdString(), executableDir.toStdString());
      sysname->setText(QDir::current().relativeFilePath(newfile));
      setSysNameParameter();
    }
    catch (KNException ex)
    {
      externalException(ex);
    }
  }
}

void MainWindow::generateSystem()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open system definition", QString(), "Vector field definition (*.vf)");
  if (!fileName.isEmpty() && compilerProcess == 0)
  {
    try {
      KNSystem::generateSystem(fileName.toStdString(), executableDir.toStdString());
      sysname->setText(QDir::current().relativeFilePath(fileName +".so"));
      setSysNameParameter();
    }
    catch (KNException ex)
    {
      externalException(ex);
    }
  }
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  if (compThread.isRunning())
  {
    compThread.setStopFlag(true);
    compThread.wait();
  }
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
  bool success = true;
  try
  {
    QDir::setCurrent(QFileInfo(fileName).absolutePath());
    parameters.loadXmlFile(fileName.toStdString());
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
  parameters.saveXmlFile(fileName.toStdString());
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

void MainWindow::setConstant(const char* name)
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
  else if (!strcmp(name,"nDeri"))
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
