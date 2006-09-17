// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "plotwindow.h"
#include <vector>
#include <cmath>

#include <QGraphicsView>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QMessageBox>
#include <QLayout>
#include <QToolBar>
#include <QAction>
#include <QToolButton>
#include <QFileDialog>
#include <QDockWidget>
#include <QListWidget>

plotWindow::plotWindow( const QString& fname, QWidget *parent ) :
	QMainWindow(parent), data(0)
{
	// open data file
	openFile( fname );
	if( !data ) return;
	if( data->isTorus() ) return;
	//
	QGraphicsView *plot = new QGraphicsView;
	plot->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
	plot->setScene( &plotdata );
// 	QToolBar *toolbar = addToolBar("default");
	QWidget *centralWidget = new QWidget;
	QVBoxLayout *centralLayout = new QVBoxLayout;
	QHBoxLayout *topContainerLayout = new QHBoxLayout;
	QGridLayout *topLayout = new QGridLayout;
	centralWidget->setLayout( centralLayout );
	const int margin = centralLayout->margin();
	centralLayout->setMargin(0);
	centralLayout->addLayout( topContainerLayout );
	centralLayout->addWidget( plot );
	topContainerLayout->setMargin(margin);
	topContainerLayout->addLayout( topLayout );
	topContainerLayout->addStretch();
	this->setCentralWidget( centralWidget );
	
	QLabel *matfileLabel = new QLabel( "MAT file" );
	QLabel *xvarLabel = new QLabel( "X Coordinate" );
	QLabel *yvarLabel = new QLabel( "Y Coordinate" );
	QLabel *ptlabelLabel = new QLabel( "LABEL" );
	QLabel *dimLabel = new QLabel( "DIM" );

	matfile = new QLineEdit;
	matfile->setReadOnly(true);
	matfile->setText( fname );
	QAction *matfileAct = new QAction(QIcon(":/res/images/cr22-action-fileopen.png"), tr("&Open"), this);
	QToolButton *matfileButton = new QToolButton();
	matfileButton->setDefaultAction( matfileAct );
	connect( matfileAct, SIGNAL(triggered()), this, SLOT(open()));
	
	xvar = new QComboBox;
	yvar = new QComboBox;
	ptlabel = new QSpinBox;
	dim = new QSpinBox;
	
	plotSelect = new QListWidget;
	plotSelect->setSelectionMode(QAbstractItemView::SingleSelection);
	QDockWidget *plotDock = new QDockWidget;
	QVBoxLayout *dockLayout = new QVBoxLayout;
	QWidget     *dockWidget = new QWidget;
	QAction     *removePlotAct = new QAction(QIcon(":/res/images/cr22-action-eraser.png"), "Remove Selected", this);
	QToolButton *removePlot = new QToolButton;
	dockLayout->setMargin(0);
	removePlot->setDefaultAction( removePlotAct );
	connect( removePlotAct, SIGNAL(triggered()), this, SLOT(removePlot()) );
	plotDock->setWidget( dockWidget );
	dockWidget->setLayout( dockLayout );
	dockLayout->addWidget( removePlot );
	dockLayout->addWidget( plotSelect );
	addDockWidget( Qt::LeftDockWidgetArea, plotDock );

	xvarMap.push_back("None");
	xvarMap.push_back("Label");
	xvarMap.push_back("Mesh");
	xvarMap.push_back("RealMultiplier");
	for( unsigned int i = XParameter0; i < XParameter0 + data->getNPar(); ++i ) {
		xvarMap.push_back(QString("P ") + QString::number( i - XParameter0 ));
	}
	for( unsigned int i = 0; i < xvarMap.size(); ++i ) xvar->insertItem( i, xvarMap.at(i) );
	
	yvarMap.push_back("None");
	yvarMap.push_back("L2Norm");
	yvarMap.push_back("Amplitude");
	yvarMap.push_back("ImagMultiplier");
	yvarMap.push_back("AbsMultiplier");
	yvarMap.push_back("Profile");
	for( unsigned int i = YParameter0; i < YParameter0 + data->getNPar(); ++i ) {
		yvarMap.push_back(QString("P ") + QString::number( i - YParameter0 ));
	}
	for( unsigned int i = 0; i < yvarMap.size(); ++i ) yvar->insertItem( i, yvarMap.at(i) );

	ptlabel->setRange( 0, data->getNCols()-1 );
	dim->setRange( 0, data->getNDim()-1 );
	
	QAction *addPlotAct = new QAction(QIcon(":/res/images/cr22-action-add.png"), tr("Add Plot"), this);
	QToolButton *addPlotButton = new QToolButton();
	addPlotButton->setDefaultAction( addPlotAct );
	connect( addPlotAct, SIGNAL(triggered()), this, SLOT(addPlot()) );
	
	QAction *clearAllPlotAct = new QAction(QIcon(":/res/images/cr22-action-remove.png"), tr("Clear All"), this);
	QToolButton *clearAllPlotButton = new QToolButton();
	clearAllPlotButton->setDefaultAction( clearAllPlotAct );
	connect( clearAllPlotAct, SIGNAL(triggered()), this, SLOT(clearPlot()) );
	
	topLayout->addWidget( matfileLabel, 0, 0 );
	topLayout->addWidget( xvarLabel, 0, 2 );
	topLayout->addWidget( yvarLabel, 0, 3 );
	topLayout->addWidget( ptlabelLabel, 0, 4 );
	topLayout->addWidget( dimLabel, 0, 5 );
	topLayout->addWidget( matfile, 1, 0 );
	topLayout->addWidget( matfileButton, 1, 1 );
	topLayout->addWidget( xvar, 1, 2 );
	topLayout->addWidget( yvar, 1, 3 );
	topLayout->addWidget( ptlabel, 1, 4 );
	topLayout->addWidget( dim, 1, 5 );
	topLayout->addWidget( addPlotButton, 1, 6 );
	topLayout->addWidget( clearAllPlotButton, 1, 7 );
	
	plot->setMinimumSize(plot->mapFromScene(plotdata.sceneRect()).boundingRect().size()*1.1+
	   QSize(2*plot->frameWidth(),2*plot->frameWidth()) );
	resize(minimumSizeHint());
}

plotWindow::~plotWindow()
{
	delete data;
}

void plotWindow::openFile( const QString& fileName )
{
	const mat4Data *t_data = 0;
	try{ t_data = new mat4Data( fileName.toStdString() ); }
	catch( pddeException ex )
	{
		QMessageBox::critical( this, "PlotWindow::openFile()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
		return;
	}
	if( t_data )
	{
		delete data;
		data = t_data;
	}
	QFileInfo fi( fileName );
	shortFileName = fi.fileName();
}

void plotWindow::open( )
{
	QString fileName = QFileDialog::getOpenFileName(this, "Open data", QString(), "v4 MAT files (*.mat);;All files (*)");
	if( !fileName.isEmpty() )
	{
		openFile( fileName );
		matfile->setText( QDir::current().relativeFilePath(fileName) );
	}
}

void plotWindow::addPlot()
{
	if( data )
	{
		if( !data->isTorus() )
		{
			try{
				plotdata.addPlot( data, (PlotXVariable)xvar->currentIndex(), (PlotYVariable)yvar->currentIndex(), ptlabel->value(), dim->value(), "r" );
			}
			catch( pddeException ex ){
				QMessageBox::critical( this, "plotWindow::addPlot()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
			}
			plotSelect->addItem(QString( "%1 : (%2, %3, L%4, D%5)")
				.arg(shortFileName).arg(xvarMap.at(xvar->currentIndex()))
				.arg(yvarMap.at(yvar->currentIndex()))
				.arg(ptlabel->value()).arg(dim->value()));
		}
	}
}

void plotWindow::removePlot()
{
	if(plotSelect->currentRow() != -1)
	{
		plotdata.clear(plotSelect->currentRow()+1);
		delete plotSelect->takeItem(plotSelect->currentRow());
	}
}

void plotWindow::clearPlot()
{
	plotdata.clearAll();
	plotSelect->clear();
}
