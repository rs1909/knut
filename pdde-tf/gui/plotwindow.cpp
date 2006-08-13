#include "plotwindow.h"
#include <vector>
#include <cmath>

plotWindow::plotWindow( const mat4Data* d, QWidget *parent ) : QMainWindow(parent), data(d)
{
	QGraphicsView *plot = new QGraphicsView;
// 	plot->setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
// 	plot->setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
// 	plot->setFrameRect
// 	PlotData *plotdata = new PlotData;

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
	
	QLabel *xvarLabel = new QLabel( "X Coordinate" );
	QLabel *yvarLabel = new QLabel( "Y Coordinate" );
	QLabel *ptlabelLabel = new QLabel( "LABEL" );
	QLabel *dimLabel = new QLabel( "DIM" );

	xvar = new QComboBox;
	yvar = new QComboBox;
	ptlabel = new QSpinBox;
	dim = new QSpinBox;
	topLayout->addWidget( xvarLabel, 0, 0 );
	topLayout->addWidget( yvarLabel, 0, 1 );
	topLayout->addWidget( ptlabelLabel, 0, 2 );
	topLayout->addWidget( dimLabel, 0, 3 );
	topLayout->addWidget( xvar, 1, 0 );
	topLayout->addWidget( yvar, 1, 1 );
	topLayout->addWidget( ptlabel, 1, 2 );
	topLayout->addWidget( dim, 1, 3 );

	std::vector<QString> xvarMap( XParameter0 + data->getNPar() );
	xvarMap[XNone] = "None";
	xvarMap[XLabel] = "Label";
	xvarMap[XMesh] = "Mesh";
	xvarMap[XRealMultiplier] = "RealMultiplier";
	for( int i = XParameter0; i < xvarMap.size(); ++i ) xvarMap[i] = QString("P ") + QString::number( i - XParameter0 );
	for( int i = 0; i < xvarMap.size(); ++i ) xvar->insertItem( i, xvarMap.at(i) );
	
	std::vector<QString> yvarMap( XParameter0 + data->getNPar() );
	yvarMap[YNone] = "None";
	yvarMap[YL2Norm] = "L2Norm";
	yvarMap[YAmplitude] = "Amplitude";
	yvarMap[YImagMultiplier] = "ImagMultiplier";
	yvarMap[YAbsMultiplier] = "AbsMultiplier";
	yvarMap[YProfile] = "Profile";
	for( int i = YParameter0; i < xvarMap.size(); ++i ) yvarMap[i] = QString("P ") + QString::number( i - YParameter0 );
	for( int i = 0; i < yvarMap.size(); ++i ) yvar->insertItem( i, yvarMap.at(i) );

	ptlabel->setRange( 0, data->getNCols()-1 );
	dim->setRange( 0, data->getNDim()-1 );
	
	QAction *addnewplot = toolbar->addAction("Add Plot");
	QAction *clearallplot = toolbar->addAction("Clear All");
	connect( addnewplot, SIGNAL(triggered()), this, SLOT(addPlot()) );
	connect( clearallplot, SIGNAL(triggered()), this, SLOT(clearPlot()) );

// 	plotdata.addPlot( *data, XMesh, YProfile, 0, 0, "b" );
// 	plotdata.clear();

// 	Vector v1(1000), v2(1000);
// 	for( int i=0; i < 1000; i++ ){ v1(i) = i*2*M_PI/1000.0; v2(i) = sin(v1(i)); }
// 	plotdata->addPlot( v1, v2, "r--" );
// 	for( int i=0; i < 1000; i++ ){ v1(i) = i*2*M_PI/1000.0; v2(i) = 1.1*sin(v1(i)+1); }
// 	plotdata->addPlot( v1, v2, "g--" );
// 	for( int i=0; i < 1000; i++ ){ v1(i) = i*2*M_PI/1000.0; v2(i) = 1.2*sin(v1(i)+2); }
// 	plotdata->addPlot( v1, v2, "c:" );
// 	plotdata->addPlot( v2, "y-." );
	plot->setMinimumSize(plot->mapFromScene(plotdata.sceneRect()).boundingRect().size()+
	   QSize(2*plot->frameWidth(),2*plot->frameWidth()) );
}

void plotWindow::addPlot()
{
	if( !data->isTorus() )
	{
		try{
			plotdata.addPlot( *data, (PlotXVariable)xvar->currentIndex(), (PlotYVariable)yvar->currentIndex(), ptlabel->value(), dim->value(), "b" );
		}
		catch( pddeException ex ){
			QMessageBox::critical( this, "plotWindow::addPlot()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
		}
	}
}

void plotWindow::clearPlot()
{
	plotdata.clear();
	plotdata.PlotPaint();
}
