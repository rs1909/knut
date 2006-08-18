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
#include "plotdata.h"
#include "plotwindow.h"
#include "paramview.h"

#include <fstream>
#include <cfloat>

#include <QtGui>

MainWindow::MainWindow( const QString& appDir ) :
	executableDir(appDir), compThread(parameters),
	inputData(0), inputPlotWindow(0),
	outputData(0), outputPlotWindow(0),
	terminalDialog(0), compilerProcess(0)
{
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
	
	tabWidget->addTab( systemWidget, QString("System") );
	tabWidget->addTab( numericsWidget, QString("Numerics") );
	tabWidget->addTab( symmetryWidget, QString("Symmetry") );
	tabWidget->addTab( torusWidget, QString("Torus") );
	
	systemWidget->setLayout( systemGrid );
	numericsWidget->setLayout( numericsGrid );
	symmetryWidget->setLayout( symmetryGrid );
	torusWidget->setLayout( torusGrid );
	
	setCentralWidget( tabWidget );

	/// QTabWidget
	/// a) files + equations b) numerics c) symmetry d) torus

	QHBoxLayout *getInputFileLayout = new QHBoxLayout;
	QLabel* inputFileLabel = new QLabel("INPUT");
	inputFileLabel->setToolTip( QString("Input file which contains the starting point") );
	QAction* inputFilePlotAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Plot"), this);
	inputFile = new QLineEdit();
	QAction* inputFileAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Browse..."), this);
	QToolButton* getInputFile = new QToolButton( );
	QToolButton* getInputFilePlot = new QToolButton( );
	getInputFile->setDefaultAction( inputFileAct );
	getInputFilePlot->setDefaultAction( inputFilePlotAct );
	getInputFileLayout->addWidget( getInputFile );
	getInputFileLayout->addWidget( getInputFilePlot );
	getInputFileLayout->addStretch();
	systemGrid->addWidget( inputFileLabel, 0, 0, Qt::AlignLeft | Qt::AlignVCenter);
	systemGrid->addWidget( inputFile, 0, 1, 1, 3 );
	systemGrid->addLayout( getInputFileLayout, 0, 4 );
	connect( inputFileAct, SIGNAL(triggered()), this, SLOT(setInputFile()));
	connect( inputFilePlotAct, SIGNAL(triggered()), this, SLOT(inputPlot()));
	connect( inputFile, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setInputFileText(const QString&)) );
	connect( &parameters, SIGNAL(inputFileChanged(const std::string&)), this, SLOT(setInputFileText(const std::string&)) );

	QHBoxLayout *getOutputFileLayout = new QHBoxLayout;
	QLabel* outputFileLabel = new QLabel("OUTPUT");
	outputFileLabel->setToolTip( QString("Output file, which contains the result of computation") );
	outputFile = new QLineEdit();
	QAction* outputFileAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Browse..."), this);
	QAction* outputFilePlotAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Plot"), this);
	QToolButton* getOutputFile = new QToolButton( );
	QToolButton* getOutputFilePlot = new QToolButton( );
	getOutputFile->setDefaultAction( outputFileAct );
	getOutputFilePlot->setDefaultAction( outputFilePlotAct );
	getOutputFileLayout->addWidget( getOutputFile );
	getOutputFileLayout->addWidget( getOutputFilePlot );
	getOutputFileLayout->addStretch();
	systemGrid->addWidget( outputFileLabel, 1, 0, Qt::AlignLeft | Qt::AlignVCenter);
	systemGrid->addWidget( outputFile, 1, 1, 1, 3 );
	systemGrid->addLayout( getOutputFileLayout, 1, 4 );
	connect( outputFileAct, SIGNAL(triggered()), this, SLOT(setOutputFile()));
	connect( outputFilePlotAct, SIGNAL(triggered()), this, SLOT(outputPlot()));
	connect( outputFile, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setOutputFileText(const QString&)) );
	connect( &parameters, SIGNAL(outputFileChanged(const std::string&)), this, SLOT(setOutputFileText(const std::string&)) );

	// this only for setting SYSNAME
	QHBoxLayout *sysnameLayout = new QHBoxLayout;
	QLabel* sysnameLabel = new QLabel("SYSNAME");
	sysnameLabel->setToolTip( QString("The compiled system definition, e.g., \"sys-problem.so\"") );
	sysname = new QLineEdit();
	QAction* sysdefAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Browse..."), this);
	QAction* compileAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Compile..."), this);
	QToolButton* getSysdef = new QToolButton( );
	QToolButton* compile = new QToolButton( );
	getSysdef->setDefaultAction( sysdefAct );
	compile->setDefaultAction( compileAct );
	sysnameLayout->addWidget( getSysdef );
	sysnameLayout->addWidget( compile );
	sysnameLayout->addStretch();
	systemGrid->addWidget( sysnameLabel, 2, 0, Qt::AlignLeft | Qt::AlignVCenter);
	systemGrid->addWidget( sysname, 2, 1, 1, 3 );
	systemGrid->addLayout( sysnameLayout, 2, 4 );
	connect( sysdefAct, SIGNAL(triggered()), this, SLOT(setSysName()));
	connect( compileAct, SIGNAL(triggered()), this, SLOT(compileSystem()));
	// sets up a bidirectional connection
	connect( sysname, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setSysNameText(const QString&)) );
	connect( &parameters, SIGNAL(sysnameChanged(const std::string&)), this, SLOT(setSysNameText(const std::string&)) );

	// setting LABEL
	QLabel* labelLabel = new QLabel("LABEL");
	labelLabel->setToolTip( QString("The label of the solution in the input file to be used as a starting point.") );
	label = new QSpinBox();
	label->setRange( 0, 0xffff );
	systemGrid->addWidget( labelLabel, 3, 0, Qt::AlignLeft | Qt::AlignBottom );
	systemGrid->addWidget( label, 4, 0 );
	connect( label, SIGNAL(valueChanged(int)), &parameters, SLOT(setLabel(int)) );
	connect( &parameters, SIGNAL(labelChanged(int)), this, SLOT(setLabel(int)) );

	QLabel* pttypeLabel = new QLabel("POINT TYPE");
	pttypeLabel->setToolTip( QString("The type of a solution to be continued.") );
	pttype = new QComboBox();
	for( int i = 0; i < parameters.pointSize(); ++i )
	{
		pttype->addItem( parameters.pointString( i ).c_str() );
	}
	// thes set the values
	connect( pttype, SIGNAL(currentIndexChanged(int)), &parameters, SLOT(setPointTypeIdx(int)));
	connect( &parameters, SIGNAL(pointTypeChangedIdx(int)), this, SLOT(setPointTypeIdx(int)) );
	// arranges some widget changes
	connect( &parameters, SIGNAL(pointTypeChangedIdx(int)), this, SLOT(setPointType()));
	
	QLabel* cpLabel = new QLabel("CP");
	cpLabel->setToolTip( QString("The continuation parameter.") );
	cp = new QComboBox( );
	setupCp();
	// this sets up the indices
	connect( &parameters, SIGNAL(sysnameChanged(const std::string&)), this, SLOT(setupCp()) );
	// this is the bidirectional connection
	connect( cp, SIGNAL(currentIndexChanged(int)), &parameters, SLOT(setCpIdx(int)));
	connect( &parameters, SIGNAL(cpChangedIdx(int)), this, SLOT(setCpIdx(int)) );

	QLabel* branchswLabel = new QLabel("SWITCH");
	branchswLabel->setToolTip("Switches to another branch at the bifurcation point.");
	branchsw = new QComboBox();
	for( int i = 0; i < parameters.branchswSize(); ++i )
	{
		branchsw->addItem( parameters.branchswString( i ).c_str() );
	}
	systemGrid->addWidget( branchswLabel, 3, 4, Qt::AlignHCenter | Qt::AlignBottom);
	systemGrid->addWidget( branchsw, 4, 4 );
	// this is the bidirectional connection
	connect( branchsw, SIGNAL(currentIndexChanged(int)), &parameters, SLOT(setBranchSWIdx(int)));
	connect( &parameters, SIGNAL(branchswChangedIdx(int)), this, SLOT(setBranchSWIdx(int)) );

	eqnsLabel = new QLabel("NEQNS");
	eqnsLabel->setToolTip( QString("NPARX: Number of additional parameters to be used in the continuation.\n"
	                               "NEQNS: Number of equations and variables that define the system.") );
	eqns = new QSpinBox();
	eqns->setRange( 0, 0xffff );
	systemGrid->addWidget( pttypeLabel, 3, 1, 1, 2, Qt::AlignHCenter | Qt::AlignBottom );
	systemGrid->addWidget( cpLabel, 3, 3, Qt::AlignHCenter | Qt::AlignBottom );
	systemGrid->addWidget( pttype, 4, 1, 1, 2 );
	systemGrid->addWidget( cp, 4, 3 );
	systemGrid->addWidget( eqnsLabel, 5, 0, Qt::AlignLeft | Qt::AlignBottom );
	systemGrid->addWidget( eqns, 6, 0, Qt::AlignTop );

	table = new EqnVarTableView( &parameters );
	systemGrid->addWidget( table, 5, 1, 2, 4, Qt::AlignVCenter );
	// this has to reconfigure the table
	connect( eqns, SIGNAL(valueChanged(int)), &parameters, SLOT(setNEqns(int)) );
	connect( &parameters, SIGNAL(neqnsChanged(int)), this, SLOT(setNEqns(int)) );
	
	// setting NINT, NDEG, NMUL, STAB, NMAT
	QLabel* nintLabel = new QLabel("NINT");
	QLabel* ndegLabel = new QLabel("NDEG");
	QLabel* nmulLabel = new QLabel("NMUL");
	QLabel* stabLabel = new QLabel("STAB");
	QLabel* nmatLabel = new QLabel("NMAT");
	nintLabel->setToolTip( QString("Number of collocation intervals.") );
	ndegLabel->setToolTip( QString("The degree of the piecewise polynomial.") );
	nmulLabel->setToolTip( QString("Number of Floquet multipliers to be computed.") );
	stabLabel->setToolTip( QString("Checked when stability information (the Floquet multipliers) is computed.") );
	nmatLabel->setToolTip( QString("An integer N (the smallest) such that tau<sub>max</sub>/T < N") );
	nint = new QSpinBox(); nint->setRange( 0, 0xffff );
	ndeg = new QSpinBox(); ndeg->setRange( 0, 7 );
	nmul = new QSpinBox(); nmul->setRange( 0, 0xffff );
	stab = new QCheckBox();
	nmat = new QSpinBox(); nmat->setRange( 0, 0xffff );
	
	numericsGrid->addWidget( nintLabel, 0, 0, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( ndegLabel, 0, 1, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( nmulLabel, 0, 2, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( stabLabel, 0, 4, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( nmatLabel, 0, 3, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( nint, 1, 0 );
	numericsGrid->addWidget( ndeg, 1, 1 );
	numericsGrid->addWidget( nmul, 1, 2 );
	numericsGrid->addWidget( stab, 1, 4, Qt::AlignHCenter );
	numericsGrid->addWidget( nmat, 1, 3 );
	connect( nint, SIGNAL(valueChanged(int)), &parameters, SLOT(setNInt(int)) );
	connect( ndeg, SIGNAL(valueChanged(int)), &parameters, SLOT(setNDeg(int)) );
	connect( nmul, SIGNAL(valueChanged(int)), &parameters, SLOT(setNMul(int)) );
	connect( stab, SIGNAL(stateChanged(int)), &parameters, SLOT(setStab(int)) );
	connect( nmat, SIGNAL(valueChanged(int)), &parameters, SLOT(setNMat(int)) );
	connect( &parameters, SIGNAL(nintChanged(int)), this, SLOT(setNInt(int)) );
	connect( &parameters, SIGNAL(ndegChanged(int)), this, SLOT(setNDeg(int)) );
	connect( &parameters, SIGNAL(nmulChanged(int)), this, SLOT(setNMul(int)) );
	connect( &parameters, SIGNAL(stabChanged(bool)), this, SLOT(setStab(bool)) );
	connect( &parameters, SIGNAL(nmatChanged(int)), this, SLOT(setNMat(int)) );

	QLabel* nint1Label = new QLabel("NINT1");
	QLabel* nint2Label = new QLabel("NINT2");
	QLabel* ndeg1Label = new QLabel("NDEG1");
	QLabel* ndeg2Label = new QLabel("NDEG2");
	nint1 = new QSpinBox(); nint1->setRange( 0, 0xffff );
	nint2 = new QSpinBox(); nint1->setRange( 0, 0xffff );
	ndeg1 = new QSpinBox(); ndeg1->setRange( 0, 0xffff );
	ndeg2 = new QSpinBox(); ndeg2->setRange( 0, 0xffff );
	torusGrid->addWidget( nint1Label, 0, 0, Qt::AlignHCenter | Qt::AlignBottom );
	torusGrid->addWidget( nint2Label, 0, 1, Qt::AlignHCenter | Qt::AlignBottom );
	torusGrid->addWidget( ndeg1Label, 0, 2, Qt::AlignHCenter | Qt::AlignBottom );
	torusGrid->addWidget( ndeg2Label, 0, 3, Qt::AlignHCenter | Qt::AlignBottom );
	torusGrid->addWidget( nint1, 1, 0 );
	torusGrid->addWidget( nint2, 1, 1 );
	torusGrid->addWidget( ndeg1, 1, 2 );
	torusGrid->addWidget( ndeg2, 1, 3 );
	connect( nint1, SIGNAL(valueChanged(int)), &parameters, SLOT(setNInt1(int)) );
	connect( nint2, SIGNAL(valueChanged(int)), &parameters, SLOT(setNInt2(int)) );
	connect( ndeg1, SIGNAL(valueChanged(int)), &parameters, SLOT(setNDeg1(int)) );
	connect( ndeg2, SIGNAL(valueChanged(int)), &parameters, SLOT(setNDeg2(int)) );
	connect( &parameters, SIGNAL(nint1Changed(int)), this, SLOT(setNInt1(int)) );
	connect( &parameters, SIGNAL(nint2Changed(int)), this, SLOT(setNInt2(int)) );
	connect( &parameters, SIGNAL(ndeg1Changed(int)), this, SLOT(setNDeg1(int)) );
	connect( &parameters, SIGNAL(ndeg2Changed(int)), this, SLOT(setNDeg2(int)) );

	QDoubleValidator* dbValid = new QDoubleValidator( -DBL_MAX, DBL_MAX, 16, this );

	QLabel* stepsLabel = new QLabel("STEPS");
	QLabel* dsLabel = new QLabel("DS");
	QLabel* dsMinLabel = new QLabel("DSMIN");
	QLabel* dsMaxLabel = new QLabel("DSMAX");
	QLabel* dsStartLabel = new QLabel("DSSTART");
	stepsLabel->setToolTip( QString("Number of continuation steps to be taken.") );
	dsLabel->setToolTip( QString("The default arc-length step size.") );
	dsMinLabel->setToolTip( QString("The minimal arc-length step size.") );
	dsMaxLabel->setToolTip( QString("The maximal arc-length step size.") );
	dsStartLabel->setToolTip( QString("The amplitude of the initial guess along the critical eigenvector at branch switching.") );
	steps = new QSpinBox(); steps->setRange( 0, 0xffff );
	ds = new QLineEdit(); ds->setValidator( dbValid );
	dsMin = new QLineEdit(); dsMin->setValidator( dbValid );
	dsMax = new QLineEdit(); dsMax->setValidator( dbValid );
	dsStart = new QLineEdit(); dsStart->setValidator( dbValid );
	numericsGrid->addWidget( stepsLabel, 2, 0, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( dsLabel, 2, 1, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( dsMinLabel, 2, 2, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( dsMaxLabel, 2, 3, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( dsStartLabel, 2, 4, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( steps, 3, 0 );
	numericsGrid->addWidget( ds, 3, 1 );
	numericsGrid->addWidget( dsMin, 3, 2 );
	numericsGrid->addWidget( dsMax, 3, 3 );
	numericsGrid->addWidget( dsStart, 3, 4 );
	connect( steps, SIGNAL(valueChanged(int)), &parameters, SLOT(setSteps(int)) );
	connect( ds, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setDs(const QString &)) );
	connect( dsMin, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setDsMin(const QString &)) );
	connect( dsMax, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setDsMax(const QString &)) );
	connect( dsStart, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setDsStart(const QString &)) );
	connect( ds, SIGNAL(editingFinished()), &parameters, SLOT(editedDs()) );
	connect( dsMin, SIGNAL(editingFinished()), &parameters, SLOT(editedDsMin()) );
	connect( dsMax, SIGNAL(editingFinished()), &parameters, SLOT(editedDsMax()) );
	connect( dsStart, SIGNAL(editingFinished()), &parameters, SLOT(editedDsStart()) );
	connect( &parameters, SIGNAL(stepsChanged(int)), this, SLOT(setSteps(int)) );
	connect( &parameters, SIGNAL(dsChanged(const QString&)), ds, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(dsMinChanged(const QString&)), dsMin, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(dsMaxChanged(const QString&)), dsMax, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(dsStartChanged(const QString&)), dsStart, SLOT(setText(const QString&)) );

	QLabel* epsCLabel = new QLabel("EPSC");
	QLabel* epsRLabel = new QLabel("EPSR");
	QLabel* epsSLabel = new QLabel("EPSS");
	QLabel* cpMinLabel = new QLabel("CPMIN");
	QLabel* cpMaxLabel = new QLabel("CPMAX");
	epsCLabel->setToolTip( QString("Convergence tolerance in continuation.") );
	epsRLabel->setToolTip( QString("Convergence tolerance when converging to a solution.") );
	epsSLabel->setToolTip( QString("Convergence tolerance when finding a special (bifurcation) point.") );
	cpMinLabel->setToolTip( QString("The minimal value of the continuation parameter.") );
	cpMaxLabel->setToolTip( QString("The maximal value of the continuation parameter.") );
	epsC = new QLineEdit(); epsC->setValidator( dbValid );
	epsR = new QLineEdit(); epsR->setValidator( dbValid );
	epsS = new QLineEdit(); epsS->setValidator( dbValid );
	cpMin = new QLineEdit(); cpMin->setValidator( dbValid );
	cpMax = new QLineEdit(); cpMax->setValidator( dbValid );
	numericsGrid->addWidget( epsCLabel, 4, 0, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( epsRLabel, 4, 1, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( epsSLabel, 4, 2, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( cpMinLabel, 4, 3, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( cpMaxLabel, 4, 4, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( epsC, 5, 0 );
	numericsGrid->addWidget( epsR, 5, 1 );
	numericsGrid->addWidget( epsS, 5, 2 );
	numericsGrid->addWidget( cpMin, 5, 3 );
	numericsGrid->addWidget( cpMax, 5, 4 );
	connect( epsC, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setEpsC(const QString &)) );
	connect( epsR, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setEpsR(const QString &)) );
	connect( epsS, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setEpsS(const QString &)) );
	connect( cpMin, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setCpMin(const QString &)) );
	connect( cpMax, SIGNAL(textChanged(const QString&)), &parameters, SLOT(setCpMax(const QString &)) );
	connect( epsC, SIGNAL(editingFinished()), &parameters, SLOT(editedEpsC()) );
	connect( epsR, SIGNAL(editingFinished()), &parameters, SLOT(editedEpsR()) );
	connect( epsS, SIGNAL(editingFinished()), &parameters, SLOT(editedEpsS()) );
	connect( cpMin, SIGNAL(editingFinished()), &parameters, SLOT(editedCpMin()) );
	connect( cpMax, SIGNAL(editingFinished()), &parameters, SLOT(editedCpMax()) );
	connect( &parameters, SIGNAL(epsCChanged(const QString&)), epsC, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(epsRChanged(const QString&)), epsR, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(epsSChanged(const QString&)), epsS, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(cpMinChanged(const QString&)), cpMin, SLOT(setText(const QString&)) );
	connect( &parameters, SIGNAL(cpMaxChanged(const QString&)), cpMax, SLOT(setText(const QString&)) );

	QLabel* nitCLabel = new QLabel("NITC");
	QLabel* nitRLabel = new QLabel("NITR");
	QLabel* nitSLabel = new QLabel("NITS");
	QLabel* nsymLabel = new QLabel("NSYM");
	nitCLabel->setToolTip( QString("Maximal number of iteration steps in continuation.") );
	nitRLabel->setToolTip( QString("Maximal number of iteration when converging to a solution.") );
	nitSLabel->setToolTip( QString("Maximal number of iteration steps when finding a special (bifurcation) point.") );
	nsymLabel->setToolTip( QString("Number of rotational symmetric complex dimensions.") );

	nitC = new QSpinBox(); nitC->setRange( 0, 0xffff );
	nitR = new QSpinBox(); nitR->setRange( 0, 0xffff );
	nitS = new QSpinBox(); nitS->setRange( 0, 0xffff );
	nsym = new QSpinBox(); nsym->setRange( 0, 0xffff );
	numericsGrid->addWidget( nitCLabel, 6, 0, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( nitRLabel, 6, 1, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( nitSLabel, 6, 2, Qt::AlignHCenter | Qt::AlignBottom );
	symmetryGrid->addWidget( nsymLabel, 0, 0, Qt::AlignHCenter | Qt::AlignBottom );
	numericsGrid->addWidget( nitC, 7, 0 );
	numericsGrid->addWidget( nitR, 7, 1 );
	numericsGrid->addWidget( nitS, 7, 2 );
	symmetryGrid->addWidget( nsym, 2, 0 );

	sym = new SYMTableView( &parameters );
	
	symmetryGrid->addWidget( sym, 3, 0, 2, 5 );
	connect( nsym, SIGNAL(valueChanged(int)), &parameters, SLOT(setNSym(int)) );
	connect( nitC, SIGNAL(valueChanged(int)), &parameters, SLOT(setNItC(int)) );
	connect( nitR, SIGNAL(valueChanged(int)), &parameters, SLOT(setNItR(int)) );
	connect( nitS, SIGNAL(valueChanged(int)), &parameters, SLOT(setNItS(int)) );
	connect( &parameters, SIGNAL(nsymChanged(int)), this, SLOT(setNSym(int)) );
	connect( &parameters, SIGNAL(nitCChanged(int)), this, SLOT(setNItC(int)) );
	connect( &parameters, SIGNAL(nitRChanged(int)), this, SLOT(setNItR(int)) );
	connect( &parameters, SIGNAL(nitSChanged(int)), this, SLOT(setNItS(int)) );
	
	// connecting exceptions
	connect( &compThread, SIGNAL(exceptionOccured(const pddeException&)), this, SLOT(externalException(const pddeException&)) );
	connect( &parameters, SIGNAL(exceptionOccured(const pddeException&)), this, SLOT(externalException(const pddeException&)) );
	// text output
	connect( &compThread, SIGNAL(printToScreen(const std::string&)), this, SLOT(terminalTextAppend(const std::string&)) );

	createActions();
	createMenus();
	createToolBars();
	createStatusBar();

	readSettings();

	setCurrentFile("");
}

void MainWindow::run()
{
	if( !compThread.isRunning() )
	{
		terminalText.clear();
		compThread.setConstants(parameters);
		compThread.setStopFlag(false);
		compThread.start();
	}
}

void MainWindow::stop()
{
	compThread.setStopFlag(true);
}

void MainWindow::setSysName()
{
	QString fileName = QFileDialog::getOpenFileName(this);
	if (!fileName.isEmpty())
	{
		sysname->setText( fileName );
	}
}

void MainWindow::setInputFile()
{
	QString fileName = QFileDialog::getOpenFileName(this, "Open input file", QString(), "v4 MAT files (*.mat);;All files (*)");
	if (!fileName.isEmpty())
	{
		inputFile->setText( QDir::current().relativeFilePath(fileName) );
	}
}

void MainWindow::setOutputFile()
{
	QString fileName = QFileDialog::getSaveFileName(this, "Save output file as", QString(), "v4 MAT files (*.mat);;All files (*)");
	if (!fileName.isEmpty())
	{
		outputFile->setText( QDir::current().relativeFilePath(fileName) );
	}
}

void MainWindow::setPointType( )
{
	if( parameters.getPointType() == SolUser )
	{
		eqnsLabel->setText(QString("NEQNS"));
		branchsw->setEnabled(true);
	} else 
	{
		eqnsLabel->setText(QString("NPARX"));
		branchsw->setEnabled(false);
	}
}

void MainWindow::inputPlotDestroyed()
{
// 	std::cout<<"plot destroyed\n";
	delete inputPlotWindow;
	delete inputData;
	inputPlotWindow = 0;
	inputData = 0;
}

void MainWindow::inputPlot()
{
	if( inputPlotWindow == 0 )
	{
		try{ inputData = new mat4Data( inputFile->text().toStdString() ); }
		catch( pddeException ex )
		{
			delete inputData;
			QMessageBox::critical( this, "MainWindow::inputPlot()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
			return;
		}
		if( inputData->isTorus() ) { delete inputData; return; }
		inputPlotWindow = new plotWindow( inputData );
		QDialog *inputPlotDialog = new QDialog( this );
		QVBoxLayout *inputPlotLayout = new QVBoxLayout();
		inputPlotLayout->addWidget( inputPlotWindow );
		inputPlotLayout->setMargin(0);
		inputPlotDialog->setLayout( inputPlotLayout );
		inputPlotDialog->setWindowTitle("Plot - input data");
		inputPlotDialog->setWindowFlags( Qt::Window );
		inputPlotDialog->show();
		inputPlotDialog->raise();
		connect( inputPlotDialog, SIGNAL(finished(int)), this, SLOT(inputPlotDestroyed()) );
	}
}

void MainWindow::outputPlotDestroyed()
{
// 	std::cout<<"plot destroyed\n";
	delete outputPlotWindow;
	delete outputData;
	outputPlotWindow = 0;
	outputData = 0;
}

void MainWindow::outputPlot()
{
	if( outputPlotWindow == 0 )
	{
		try
		{
			outputData = new mat4Data( outputFile->text().toStdString() );
		}
		catch( pddeException ex )
		{
			delete outputData;
			QMessageBox::critical( this, "MainWindow::outputPlot()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
			return;
		}
		if( outputData->isTorus() ) { delete outputData; return; }
		outputPlotWindow = new plotWindow( outputData );
		QDialog *outputPlotDialog = new QDialog( this );
		QVBoxLayout *outputPlotLayout = new QVBoxLayout();
		outputPlotLayout->addWidget( outputPlotWindow );
		outputPlotLayout->setMargin(0);
		outputPlotDialog->setLayout( outputPlotLayout );
		outputPlotDialog->setWindowTitle("Plot - output data");
		outputPlotDialog->setWindowFlags( Qt::Window );
		outputPlotDialog->show();
		outputPlotDialog->raise();
		connect( outputPlotDialog, SIGNAL(finished(int)), this, SLOT(outputPlotDestroyed()) );
	}
}

void MainWindow::terminalViewDestroyed()
{
// 	std::cout<<"TERM DESTROYED\n";
// 	delete terminalDialog;
	terminalDialog = 0;
}

void MainWindow::terminalView()
{
	if( terminalDialog == 0 )
	{
		terminalDialog = new screenDialog( this );
		terminalDialog->setWindowTitle("terminal view");
		terminalDialog->setText( terminalText );
		terminalDialog->setWindowFlags( Qt::Window );
		terminalDialog->show();
		terminalDialog->raise();
		connect( &compThread, SIGNAL(printToScreen(const std::string&)), terminalDialog, SLOT(append(const std::string&)) );
		connect( terminalDialog, SIGNAL(finished(int)), this, SLOT(terminalViewDestroyed()) );
	}
}

void MainWindow::compileSystem()
{
	QString fileName = QFileDialog::getOpenFileName(this, "Open system definition", QString(), "C++ source (*.cpp)");
	if ( !fileName.isEmpty() && compilerProcess == 0 )
	{
		QString newfile( fileName );
		newfile.replace( QString(".cpp"), QString(".so") );
		QFileInfo compiler( CMAKE_CXX_COMPILER );
		compilerProcess = new QProcess();
		compilerProcess->setProcessChannelMode( QProcess::MergedChannels );
		QString cmdLine( QString( "%1 %2 %3 %4 -I%5 %6 -o %7" )
			.arg(compiler.fileName())
			.arg(CMAKE_CXX_FLAGS)
			.arg(CMAKE_SHARED_LIBRARY_C_FLAGS)
			.arg(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS)
			.arg("\""+QDir::cleanPath(executableDir+"/../include/")+"\"")
			.arg("\""+fileName+"\"").arg("\""+newfile+"\"") );
		std::cout<<"The command line:\n\t"<<cmdLine.toStdString()<<"\n";
		compilerProcess->start( cmdLine );
		bool finished = compilerProcess->waitForFinished();
		QErrorMessage *errorPopup = new QErrorMessage( this );
		if ( finished )
		{
			if( compilerProcess->exitCode() == 0 )
			{
				if( !compilerProcess->readAll().isEmpty() )
				{
					errorPopup->setWindowTitle( "Successful compilation" );
					errorPopup->showMessage( QString( compilerProcess->readAll().constData() ) );
					errorPopup->show();
					errorPopup->raise();
				}
				sysname->setText(newfile);
			} else
			{
				errorPopup->setWindowTitle( "Compilation failed" );
				errorPopup->showMessage( QString( compilerProcess->readAll().constData() ) );
				errorPopup->show();
				errorPopup->raise();
			}
		} else
		{
			errorPopup->setWindowTitle( "Compilation timed out" );
			errorPopup->showMessage( compilerProcess->errorString() );
			errorPopup->show();
			errorPopup->raise();
		}
		delete compilerProcess;
		compilerProcess = 0;
	}
}

void MainWindow::closeEvent(QCloseEvent *event)
{
	writeSettings();
	event->accept();
}

void MainWindow::newFile()
{
	setCurrentFile("");
}

void MainWindow::open()
{
	QString fileName = QFileDialog::getOpenFileName(this, "Open constants file", QString(), "Constants files (*.xml);;All files (*)");
	if (!fileName.isEmpty())
	{
		QDir::setCurrent( QFileInfo(fileName).absolutePath() );
		loadFile(fileName);
	}
}

bool MainWindow::save()
{
	if (curFile.isEmpty()) {
		return saveAs();
	} else {
		return saveFile(curFile);
	}
}

bool MainWindow::saveAs()
{
	QString fileName = QFileDialog::getSaveFileName(this);
	if (fileName.isEmpty())
		return false;

	return saveFile(fileName);
}

void MainWindow::about()
{
	QMessageBox::about(this, tr("About PDDE-CONT"),
		QString(tr(
			"<h3>%1: A continuation software for delay-differential equations</h3>"
			"<p>Version %2 (%3), Copyright (c) 2002-2006 Robert Szalai</p>"
			"<p>Available at <a href=http://www.mm.bme.hu/~szalai/pdde>"
			"http://www.mm.bme.hu/~szalai/pdde</a></p>"
			"<p>This program is free software; you can redistribute it and/or "
			"modify it under the terms of the "
			"<a href=http://www.gnu.org/copyleft/gpl.html>GNU General Public License</a> "
			"as published by the Free Software Foundation; version 2 "
			"of the License.</p>"
			"<p>This program is distributed in the hope that it will be useful, "
			"but WITHOUT ANY WARRANTY; without even the implied warranty of "
			"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the "
			"GNU General Public License for more details.</p>"
			"<p>You should have received a copy of the GNU General Public License "
			"along with this program; if not, write to the Free Software "
			"Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.</p>"
		)).arg(PACKAGE_NAME).arg(PACKAGE_VERSION).arg(PACKAGE_REVISION) );
}

void MainWindow::createActions()
{
	runAct = new QAction(/*QIcon(":/images/new.png"),*/ tr("&Run"), this);
	runAct->setShortcut(tr("Ctrl+R"));
	runAct->setStatusTip(tr("Start the computation"));
	connect(runAct, SIGNAL(triggered()), this, SLOT(run()));

	stopAct = new QAction(/*QIcon(":/images/new.png"),*/ tr("&Stop"), this);
	stopAct->setShortcut(tr("Ctrl+T"));
	stopAct->setStatusTip(tr("Stop the computation"));
	connect(stopAct, SIGNAL(triggered()), this, SLOT(stop()));

	terminalAct = new QAction(/*QIcon(":/images/new.png"),*/ tr("&Text"), this);
	terminalAct->setShortcut(tr("Ctrl+V"));
	terminalAct->setStatusTip(tr("View textual output"));
	connect(terminalAct, SIGNAL(triggered()), this, SLOT(terminalView()));

	openAct = new QAction(/*QIcon(":/images/open.png"),*/ tr("&Open..."), this);
	openAct->setShortcut(tr("Ctrl+O"));
	openAct->setStatusTip(tr("Open an existing file"));
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

	saveAct = new QAction(/*QIcon(":/images/save.png"),*/ tr("&Save"), this);
	saveAct->setShortcut(tr("Ctrl+S"));
	saveAct->setStatusTip(tr("Save the document to disk"));
	connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

	saveAsAct = new QAction(tr("Save &As..."), this);
	saveAsAct->setStatusTip(tr("Save the document under a new name"));
	connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));

	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

	aboutAct = new QAction(tr("&About"), this);
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
	QSettings settings("Trolltech", "Application Example");
	QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
	QSize size = settings.value("size", QSize(400, 400)).toSize();
	resize(size);
	move(pos);
}

void MainWindow::writeSettings()
{
	QSettings settings("Trolltech", "Application Example");
	settings.setValue("pos", pos());
	settings.setValue("size", size());
}


void MainWindow::loadFile(const QString &fileName)
{
	try{ parameters.loadXmlFile( fileName.toStdString() ); }
	catch( pddeException ex ){ externalException( ex ); }
	setCurrentFile(fileName);
	statusBar()->showMessage(tr("File loaded"), 2000);
}

bool MainWindow::saveFile(const QString &fileName)
{
	parameters.saveXmlFile( fileName.toStdString() );
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

	setWindowTitle(tr("%1[*] - %2").arg(shownName).arg(tr("PDDE-CONT")));
}

QString MainWindow::strippedName(const QString &fullFileName)
{
	return QFileInfo(fullFileName).fileName();
}
