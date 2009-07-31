// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "plotwindow.h"
#include "mainwindow.h"
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
#include <QColorDialog>
#include <QPrintDialog>
#include <QPrinter>
#include <QSplitter>
#include <QSvgGenerator>

plotWindow::plotWindow(const QString& fname, QWidget *parent) :
    QMainWindow(parent), data(0)
{
  plot = new QGraphicsView(this);
  plot->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform | QPainter::TextAntialiasing);
  plot->setScene(&plotdata);
  QWidget *centralWidget = new QWidget;
  QVBoxLayout *centralLayout = new QVBoxLayout;
  QHBoxLayout *topContainerLayout = new QHBoxLayout;
  QGridLayout *topLayout = new QGridLayout;
  centralWidget->setLayout(centralLayout);
  const int margin = centralLayout->margin();
  centralLayout->setMargin(0);
  centralLayout->addLayout(topContainerLayout);
  centralLayout->addWidget(plot);
  topContainerLayout->setMargin(margin);
  topContainerLayout->addLayout(topLayout);
  topContainerLayout->addStretch();

  QLabel *matfileLabel = new QLabel("MAT file");
  QLabel *xvarLabel = new QLabel("X Coordinate");
  QLabel *yvarLabel = new QLabel("Y Coordinate");
  QLabel *ptlabelLabel = new QLabel("LABEL");
  QLabel *dimLabel = new QLabel("DIM");

  matfile = new QLineEdit;
  matfile->setReadOnly(true);
  QAction *matfileAct = new QAction(QIcon(":/res/images/cr16-action-fileopen.png"), tr("&Open"), this);
  QToolButton *matfileButton = new QToolButton();
  matfileButton->setDefaultAction(matfileAct);
  connect(matfileAct, SIGNAL(triggered()), this, SLOT(open()));

  xvar = new QComboBox;
  yvar = new QComboBox;
  ptlabel = new QSpinBox;
  dim = new QSpinBox;

  plotSelect = new QListWidget;
  plotSelect->setSelectionMode(QAbstractItemView::SingleSelection);
  QVBoxLayout *dockLayout = new QVBoxLayout;
  QHBoxLayout *buttonsLayout = new QHBoxLayout;
  QWidget     *dockWidget = new QWidget;
  QAction     *removePlotAct = new QAction(QIcon(":/res/images/cr16-action-eraser.png"), "Remove Selected", this);
  QToolButton *removePlotButton = new QToolButton;
  removePlotButton->setDefaultAction(removePlotAct);
  connect(removePlotAct, SIGNAL(triggered()), this, SLOT(removePlot()));

  QAction     *colorizePlotAct = new QAction(QIcon(":/res/images/cr16-action-colorize.png"), "Change Color", this);
  QToolButton *colorizePlotButton = new QToolButton;
  colorizePlotButton->setDefaultAction(colorizePlotAct);
  connect(colorizePlotAct, SIGNAL(triggered()), this, SLOT(colorizePlot()));

  dockWidget->setLayout(dockLayout);
  dockLayout->setMargin(0);
  dockLayout->addSpacing(margin);
  dockLayout->addLayout(buttonsLayout);
  dockLayout->addWidget(plotSelect);
  buttonsLayout->addWidget(removePlotButton);
  buttonsLayout->addWidget(colorizePlotButton);

  QAction *addPlotAct = new QAction(QIcon(":/res/images/cr16-action-add.png"), tr("Add Plot"), this);
  QToolButton *addPlotButton = new QToolButton();
  addPlotButton->setDefaultAction(addPlotAct);
  connect(addPlotAct, SIGNAL(triggered()), this, SLOT(addPlot()));

  QAction *clearAllPlotAct = new QAction(QIcon(":/res/images/cr16-action-remove.png"), tr("Clear All"), this);
  QToolButton *clearAllPlotButton = new QToolButton();
  clearAllPlotButton->setDefaultAction(clearAllPlotAct);
  connect(clearAllPlotAct, SIGNAL(triggered()), this, SLOT(clearPlot()));

  QAction *printAct = new QAction(QIcon(":/res/images/cr16-action-fileprint.png"), tr("Print"), this);
  QToolButton *printButton = new QToolButton();
  printButton->setDefaultAction(printAct);
  connect(printAct, SIGNAL(triggered()), this, SLOT(print()));

  QAction *exportSvgAct = new QAction(QIcon(":/res/images/cr16-action-svg-export.png"), tr("Export to SVG"), this);
  QToolButton *exportSvgButton = new QToolButton();
  exportSvgButton->setDefaultAction( exportSvgAct );
  connect( exportSvgAct, SIGNAL(triggered()), this, SLOT(exportSvg()) );

  QSplitter *splitter = new QSplitter;
  splitter->addWidget(dockWidget);
  splitter->addWidget(centralWidget);
  this->setCentralWidget(splitter);

  topLayout->addWidget(matfileLabel, 0, 0);
  topLayout->addWidget(xvarLabel, 0, 2);
  topLayout->addWidget(yvarLabel, 0, 3);
  topLayout->addWidget(ptlabelLabel, 0, 4);
  topLayout->addWidget(dimLabel, 0, 5);
  topLayout->addWidget(matfile, 1, 0);
  topLayout->addWidget(matfileButton, 1, 1);
  topLayout->addWidget(xvar, 1, 2);
  topLayout->addWidget(yvar, 1, 3);
  topLayout->addWidget(ptlabel, 1, 4);
  topLayout->addWidget(dim, 1, 5);
  topLayout->addWidget(addPlotButton, 1, 6);
  topLayout->addWidget(clearAllPlotButton, 1, 7);
  topLayout->addWidget(printButton, 1, 8);
  topLayout->addWidget(exportSvgButton, 1, 9 );

  plot->setMinimumSize(plot->mapFromScene(plotdata.sceneRect()).boundingRect().size()*1.1 +
                       QSize(2*plot->frameWidth(), 2*plot->frameWidth()));
                       
  // open data file
  openFile(fname);
  if (!data) return;
  if (data->isTorus()) return;
  matfile->setText(QDir::current().relativeFilePath(fname));
}

plotWindow::~plotWindow()
{
  delete data;
}

void plotWindow::openFile(const QString& fileName)
{
  const mat4Data *t_data = 0;
  try
  {
    t_data = new mat4Data(fileName.toStdString());
  }
  catch (knutException ex)
  {
    MainWindow::showException(this, ex);
    return;
  }
  if (t_data)
  {
    delete data;
    data = t_data;
   
    data->lock(); 
    ptlabel->setRange(0, data->getNCols() - 1);
    dim->setRange(0, data->getNDim() - 1);

    const int xidx = xvar->currentIndex();
    std::vector<std::string> parNames;
    data->getParNames(parNames);
    xvarMap.clear();
    xvar->clear();
    xvarMap.push_back("None");
    xvarMap.push_back("Label");
    xvarMap.push_back("Mesh");
    xvarMap.push_back("RealMultiplier");
    for (int i = XParameter0; i < XParameter0 + data->getNPar(); ++i)
    {
      if (i - XParameter0 < (int)parNames.size()) xvarMap.push_back(QString(parNames[i - XParameter0].c_str()));
    }
    for (unsigned int i = 0; i < xvarMap.size(); ++i)
      xvar->insertItem(static_cast<int>(i), xvarMap.at(i));
    if ((xidx != -1) && (xidx < xvar->count())) xvar->setCurrentIndex(xidx);

    const int yidx = yvar->currentIndex();  
    yvarMap.clear();
    yvar->clear();
    yvarMap.push_back("None");
    yvarMap.push_back("L2Norm");
    yvarMap.push_back("Amplitude");
    yvarMap.push_back("ImagMultiplier");
    yvarMap.push_back("AbsMultiplier");
    yvarMap.push_back("Profile");
    for (int i = YParameter0; i < YParameter0 + data->getNPar(); ++i)
    {
      if (i - YParameter0 < (int)parNames.size()) yvarMap.push_back(QString(parNames[i - YParameter0].c_str()));
    }
    for (unsigned int i = 0; i < yvarMap.size(); ++i)
      yvar->insertItem(static_cast<int>(i), yvarMap.at(i));
    if ((yidx != -1) && (yidx < yvar->count())) yvar->setCurrentIndex(yidx);
    data->unlock();
  }
  QFileInfo fi(fileName);
  shortFileName = fi.fileName();
}

void plotWindow::open()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open data", QString(), "v4 MAT files (*.mat);;All files (*)");
  if (!fileName.isEmpty())
  {
    openFile(fileName);
    matfile->setText(QDir::current().relativeFilePath(fileName));
  }
}

void plotWindow::addPlot()
{
  if (data)
  {
    data->lock();
    // make sure that the data is consistent
    const_cast<mat4Data*>(data)->initHeaders("unnamed file");
    if (!data->isTorus())
    {
      bool added = false;
      try
      {
        added = plotdata.addPlot(data, (PlotXVariable)xvar->currentIndex(),
                                 (PlotYVariable)yvar->currentIndex(), ptlabel->value(), dim->value());
      }
      catch (knutException ex)
      {
        MainWindow::showException(this, ex);
      }
      if (added)
      {
        plotSelect->addItem(QString("%1 : (%2, %3, L%4, D%5)")
          .arg(shortFileName).arg(xvarMap.at(static_cast<unsigned int>(xvar->currentIndex())))
          .arg(yvarMap.at(static_cast<unsigned int>(yvar->currentIndex())))
          .arg(ptlabel->value()).arg(dim->value()));
      }
    }
    data->unlock();
  }
}

void plotWindow::removePlot()
{
  if (plotSelect->currentRow() != -1)
  {
    plotdata.clear(plotSelect->currentRow() + 1);
    delete plotSelect->takeItem(plotSelect->currentRow());
  }
}

void plotWindow::colorizePlot()
{
  if (plotSelect->currentRow() != -1)
  {
    QColor newcolor = QColorDialog::getColor(plotdata.getColor(plotSelect->currentRow() + 1));
    if (newcolor.isValid()) plotdata.setColor(plotSelect->currentRow() + 1, newcolor);
  }
}

void plotWindow::clearPlot()
{
  plotdata.clearAll();
  plotSelect->clear();
}

void plotWindow::print()
{
  QPrinter printer;
  printer.setOrientation(QPrinter::Landscape);
  if (QPrintDialog(&printer).exec() == QDialog::Accepted)
  {
    QRectF sourceRect = plotdata.itemsBoundingRect();
    QPainter painter(&printer);
    painter.setRenderHint(QPainter::Antialiasing);
    plotdata.render(&painter, QRectF(printer.pageRect()), sourceRect);
  }
}

void plotWindow::exportSvg()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Export to SVG", "plot.svg", "SVG files (*.svg);;All files (*)");
  if( !fileName.isEmpty() )
  {
    QRectF sourceRect = plotdata.itemsBoundingRect();
    QRectF targetRect = QRectF(0.0, 0.0, sourceRect.size().width(), sourceRect.size().height());
    QSvgGenerator printer;
    printer.setFileName(fileName);
    printer.setSize(sourceRect.size().toSize());
    printer.setViewBox(targetRect);
    QPainter painter(&printer);
    plotdata.render(&painter, targetRect, sourceRect);
    painter.end();
  }
}
