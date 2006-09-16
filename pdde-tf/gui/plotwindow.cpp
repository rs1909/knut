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

plotWindow::plotWindow( const QString& fname, QWidget *parent ) : QMainWindow(parent), data(0)
{
	QGraphicsView *plot = new QGraphicsView;
	plot->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
	// open data file
	try{ data = new mat4Data( fname.toStdString() ); }
	catch( pddeException ex )
	{
		delete data; data = 0;
		QMessageBox::critical( this, "PlotWindow::PlotWindow()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
		return;
	}
	if( data->isTorus() ) return;
	//
	plot->setScene( &plotdata );
	QToolBar *toolbar = addToolBar("default");
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
	QAction *matfileAct = new QAction( tr("&Open"), this );
	QToolButton *matfileButton = new QToolButton();
	matfileButton->setDefaultAction( matfileAct );
	xvar = new QComboBox;
	yvar = new QComboBox;
	ptlabel = new QSpinBox;
	dim = new QSpinBox;
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
	connect( matfileAct, SIGNAL(triggered()), this, SLOT(open()));

	std::vector<QString> xvarMap( XParameter0 + data->getNPar() );
	xvarMap[XNone] = "None";
	xvarMap[XLabel] = "Label";
	xvarMap[XMesh] = "Mesh";
	xvarMap[XRealMultiplier] = "RealMultiplier";
	for( unsigned int i = XParameter0; i < xvarMap.size(); ++i ) xvarMap[i] = QString("P ") + QString::number( i - XParameter0 );
	for( unsigned int i = 0; i < xvarMap.size(); ++i ) xvar->insertItem( i, xvarMap.at(i) );
	
	std::vector<QString> yvarMap( XParameter0 + data->getNPar() );
	yvarMap[YNone] = "None";
	yvarMap[YL2Norm] = "L2Norm";
	yvarMap[YAmplitude] = "Amplitude";
	yvarMap[YImagMultiplier] = "ImagMultiplier";
	yvarMap[YAbsMultiplier] = "AbsMultiplier";
	yvarMap[YProfile] = "Profile";
	for( unsigned int i = YParameter0; i < xvarMap.size(); ++i ) yvarMap[i] = QString("P ") + QString::number( i - YParameter0 );
	for( unsigned int i = 0; i < yvarMap.size(); ++i ) yvar->insertItem( i, yvarMap.at(i) );

	ptlabel->setRange( 0, data->getNCols()-1 );
	dim->setRange( 0, data->getNDim()-1 );
	
	QAction *addnewplot = toolbar->addAction("Add Plot");
	QAction *clearallplot = toolbar->addAction("Clear All");
	connect( addnewplot, SIGNAL(triggered()), this, SLOT(addPlot()) );
	connect( clearallplot, SIGNAL(triggered()), this, SLOT(clearPlot()) );

	plot->setMinimumSize(plot->mapFromScene(plotdata.sceneRect()).boundingRect().size()+
	   QSize(2*plot->frameWidth(),2*plot->frameWidth()) );
}

plotWindow::~plotWindow()
{
	delete data;
}

void plotWindow::open()
{
	QString fileName = QFileDialog::getOpenFileName(this, "Open data", QString(), "v4 MAT files (*.mat);;All files (*)");
	if( !fileName.isEmpty() )
	{
		const mat4Data *t_data;
		try{ t_data = new mat4Data( fileName.toStdString() ); }
		catch( pddeException ex )
		{
			delete t_data; t_data = 0;
			QMessageBox::critical( this, "PlotWindow::open()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
			return;
		}
		if( t_data )
		{
			delete data; data = t_data;
			matfile->setText( QDir::current().relativeFilePath(fileName) );
		}
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
		}
	}
}

void plotWindow::clearPlot()
{
	plotdata.clear();
	plotdata.PlotPaint();
}
