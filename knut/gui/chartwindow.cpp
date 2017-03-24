// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "chartwindow.h"
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

#include <QtCharts/QChartView>
#include <QtGui/QMouseEvent>

KnutChartView::KnutChartView(QChart *chart, QWidget *parent) :
    QChartView(chart, parent),
    m_isTouching(false)
{
    setRubberBand(QChartView::RectangleRubberBand);
}

bool KnutChartView::viewportEvent(QEvent *event)
{
    if (event->type() == QEvent::TouchBegin) {
        // By default touch events are converted to mouse events. So
        // after this event we will get a mouse event also but we want
        // to handle touch events as gestures only. So we need this safeguard
        // to block mouse events that are actually generated from touch.
        m_isTouching = true;

        // Turn off animations when handling gestures they
        // will only slow us down.
        chart()->setAnimationOptions(QChart::NoAnimation);
    }
    return QChartView::viewportEvent(event);
}

void KnutChartView::mousePressEvent(QMouseEvent *event)
{
    if (m_isTouching)
        return;
    QChartView::mousePressEvent(event);
}

void KnutChartView::mouseMoveEvent(QMouseEvent *event)
{
    if (m_isTouching)
        return;
    QChartView::mouseMoveEvent(event);
}

void KnutChartView::mouseReleaseEvent(QMouseEvent *event)
{
    if (m_isTouching)
        m_isTouching = false;

    // Because we disabled animations when touch event was detected
    // we must put them back on.
    chart()->setAnimationOptions(QChart::SeriesAnimations);

    QChartView::mouseReleaseEvent(event);
}

void KnutChartView::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Plus:
        chart()->zoomIn();
        break;
    case Qt::Key_Minus:
        chart()->zoomOut();
        break;
    case Qt::Key_Left:
        chart()->scroll(-10, 0);
        break;
    case Qt::Key_Right:
        chart()->scroll(10, 0);
        break;
    case Qt::Key_Up:
        chart()->scroll(0, 10);
        break;
    case Qt::Key_Down:
        chart()->scroll(0, -10);
        break;
    default:
        QGraphicsView::keyPressEvent(event);
        break;
    }
}

plotWindow::plotWindow(QWidget *parent) :
    QSplitter(parent) //, data(0)
{
  setupPlotWindow();
}

plotWindow::plotWindow(const QString& filename, QWidget *parent) :
    QSplitter(parent), dataFileInfo(filename)
{
  setupPlotWindow();
  matfile->setText(dataFileInfo.fileName());
}

plotWindow::~plotWindow()
{
//  if (data) emit closeFile(data);
}

void plotWindow::init(Var cp)
{
  emit openFile(dataFileInfo.absoluteFilePath());
  bool found = false;
  for (unsigned int i=0; i < xvarMap.size(); ++i)
  {
    if ((xvarMap[i].value-XParameter0) == (cp-VarPAR0))
    {
      xvar->setCurrentIndex(i);
      found = true;
    }
  }
  if (found)
  {
    yvar->setCurrentIndex(2); // Amplitude
    emit requestPlot (dataFileInfo.absoluteFilePath());
  }
}

void plotWindow::setupPlotWindow()
{
  plot = new KnutChartView(plotdata.getChart(), this);
  plot->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform | QPainter::TextAntialiasing);
//   plot->setChart(plotdata.getChart());
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
  

  // There is no easy way to change axes in QCharts
//   QComboBox* xyLog = new QComboBox;
//   xyLog->insertItem(0, "Linear");
//   xyLog->insertItem(1, "SemiLogX");
//   xyLog->insertItem(2, "SemiLogY");
//   xyLog->insertItem(3, "LogLog");
//   xyLog->setCurrentIndex(0);
//   connect( xyLog, SIGNAL(currentIndexChanged(int)), &plotdata, SLOT(setAxes(int)) );

  // Qcharts is automatically resizable
//   QSpinBox* plotsize = new QSpinBox;
//   plotsize->setMinimum(100);
//   plotsize->setMaximum(1280);
//   connect(plotsize, SIGNAL(valueChanged(int)), &plotdata, SLOT(setXSize(int)));
  
  this->addWidget(dockWidget);
  this->addWidget(centralWidget);

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
  topLayout->addWidget(addPlotButton, 2, 1);
  topLayout->addWidget(clearAllPlotButton, 2, 2);
  topLayout->addWidget(printButton, 2, 3);
  topLayout->addWidget(exportSvgButton, 2, 4 );
//   topLayout->addWidget(plotsize, 2, 5 );
//   topLayout->addWidget(xyLog, 2, 0 );

//   plot->setMinimumSize(plot->mapFromScene(plotdata.sceneRect()).boundingRect().size()*1.1 +
//                        QSize(2*plot->frameWidth(), 2*plot->frameWidth()));
}

void plotWindow::initPlot(const KNDataFile* mat)
{
  if (mat == nullptr) { std::cout << "Can't set data\n"; return; }
//  data = mat;
  mat->lockRead();
  QFileInfo fi(QString::fromStdString(mat->getFileName()));
  if (fi.isSymLink()) fi = QFileInfo(fi.symLinkTarget());
  dataFileInfo = fi;
  shortFileName = fi.fileName();

  matfile->setText(QDir::current().relativeFilePath(mat->getFileName().c_str()));
  ptlabel->setRange(0, mat->getNCols() - 1);
  dim->setRange(0, mat->getNDim() - 1);

  const int xidx = xvar->currentIndex();
  std::vector<std::string> parNames;
  mat->getParNames(parNames);
  xvarMap.clear();
  xvarMap.push_back(XTranslate(XNone, "None"));
  xvarMap.push_back(XTranslate(XLabel, "Label"));
  xvarMap.push_back(XTranslate(XMesh, "Mesh"));
  xvarMap.push_back(XTranslate(XRealMultiplier, "RealMultiplier"));
  for (size_t i = 0; i < mat->getNPar(); ++i)
  {
    if (i < parNames.size())
    {
      xvarMap.push_back(XTranslate((PlotXVariable)(XParameter0 + i), QString::fromStdString(parNames[i])));
    }
  }
  xvar->clear();
  for (size_t i = 0; i < xvarMap.size(); ++i)
  {
    xvar->insertItem(static_cast<int>(i), xvarMap.at(i).name);
  }
  if ((xidx != -1) && (xidx < xvar->count())) xvar->setCurrentIndex(xidx);

  const int yidx = yvar->currentIndex();  
  yvarMap.clear();
  yvarMap.push_back(YTranslate(YNone, "None"));
  yvarMap.push_back(YTranslate(YL2Norm, "L2Norm"));
  yvarMap.push_back(YTranslate(YAmplitude, "Amplitude"));
  yvarMap.push_back(YTranslate(YPoincare, "Map"));
  yvarMap.push_back(YTranslate(YImagMultiplier, "ImagMultiplier"));
  yvarMap.push_back(YTranslate(YAbsMultiplier, "AbsMultiplier"));
  yvarMap.push_back(YTranslate(YProfile, "Profile"));
  for (size_t i = 0; i < mat->getNPar(); ++i)
  {
    if (i < parNames.size())
    {
      yvarMap.push_back(YTranslate((PlotYVariable)(YParameter0 + i), QString::fromStdString(parNames[i])));
    }
  }
  yvar->clear();
  for (size_t i = 0; i < yvarMap.size(); ++i)
  {
    yvar->insertItem(static_cast<int>(i), yvarMap.at(i).name);
  }
  if ((yidx != -1) && (yidx < yvar->count())) yvar->setCurrentIndex(yidx);

  mat->unlock();
}

void plotWindow::open()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open data", QString(), "v4 MAT files (*.mat);;All files (*)");
  if (!fileName.isEmpty())
  {
    emit openFile(fileName);
  }
}

void plotWindow::addPlot(const KNDataFile* data)
{
  // Need to query data
  if (data)
  {
//     std::cout << "plotWindow::addPlot 1\n";
    data->lockRead();
    QFileInfo fi(QString::fromStdString(data->getFileName()));
    data->unlock();
    if (fi.isSymLink()) fi = QFileInfo(fi.symLinkTarget());
    if (dataFileInfo != fi) return; // exit if not the same file as before
//     std::cout << "plotWindow::addPlot 2\n";

    data->lockRead();
    // make sure that the data is consistent
    const_cast<KNDataFile*>(data)->initHeaders();
    if (!data->isTorus())
    {
      bool added = false;
      try
      {
        added = plotdata.addPlot(data, xvarMap.at(xvar->currentIndex()).value,
                                 yvarMap.at(yvar->currentIndex()).value, ptlabel->value(), dim->value());
      }
      catch (KNException ex)
      {
        MainWindow::showException(this, ex);
      }
      if (added)
      {
      	QListWidgetItem* item = new QListWidgetItem(QString("%1 : (%2, %3, L%4, D%5)")
          .arg(shortFileName).arg(xvarMap.at(static_cast<unsigned int>(xvar->currentIndex())).name)
          .arg(yvarMap.at(static_cast<unsigned int>(yvar->currentIndex())).name)
          .arg(ptlabel->value()).arg(dim->value()),
          plotSelect);
        item->setForeground(QBrush(plotdata.getColor(plotSelect->count())));
        plotSelect->addItem(item);
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
    QRectF sourceRect = plot -> sceneRect();
    QPainter painter(&printer);
    painter.setRenderHint(QPainter::Antialiasing);
    plot -> render(&painter, QRectF(printer.pageRect()), sourceRect.toRect());
  }
}

void plotWindow::exportSvg()
{
  QString fileName = QFileDialog::getSaveFileName(plot, "Export to SVG", "plot.svg", "SVG files (*.svg);;All files (*)");
  if( !fileName.isEmpty() )
  {
    QRectF sourceRect = plot -> sceneRect();
    QRectF targetRect = QRectF(0.0, 0.0, sourceRect.size().width(), sourceRect.size().height());
    QSvgGenerator printer;
    printer.setFileName(fileName);
    printer.setSize(sourceRect.size().toSize());
    printer.setViewBox(targetRect);
    QPainter painter(&printer);
    plot -> render(&painter, targetRect, sourceRect.toRect());
    painter.end();
  }
}

void plotWindow::updatePlot(const KNDataFile* data)
{
//   std::cout << "plotWindow::updatePlot\n";
  data->lockRead();
  plotdata.updatePlot(data);
  data->unlock();
  emit updated();
}
