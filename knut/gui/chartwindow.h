// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include "datachart.h"
#include "pointtype.h"

#include <QSplitter>
#include <QFileInfo>
#include <QCloseEvent>
#include <QtCharts/QChartView>

class QLineEdit;
class QComboBox;
class QSpinBox;
class QListWidget;

using namespace QtCharts;


#include <QtCharts/QChartView>
#include <QtWidgets/QRubberBand>

QT_CHARTS_USE_NAMESPACE

class KnutChartView : public QChartView
{
public:
    KnutChartView(QChart *chart, QWidget *parent = nullptr);

protected:
    bool viewportEvent(QEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

private:
    bool m_isTouching;
};

template <class T> class Translate
{
  public:
  Translate(T v, const QString& n) : value(v), name(n) { }
  Translate(const Translate& cp) : value(cp.value), name(cp.name) {}
  Translate & operator=(Translate& cp) { value = cp.value; name = cp.name; return *this; }
  T value;
  QString name;
};

typedef Translate<PlotXVariable> XTranslate;
typedef Translate<PlotYVariable> YTranslate;

class plotWindow : public QSplitter
{
    Q_OBJECT
  public:
    plotWindow(const QString& filename, QWidget *parent = nullptr);
    plotWindow(QWidget *parent = nullptr);
    ~plotWindow() override;
    void init(Var cp);
  protected:
    void closeEvent(QCloseEvent *event) override
    {
      emit windowClosed();
      event->accept();
    }
  private:
    void setupPlotWindow();
    // variables
    DataChart plotdata;
    KnutChartView *plot;
    // gui elements
    QLineEdit *matfile;
    QComboBox *xvar;
    QComboBox *yvar;
    std::vector<XTranslate> xvarMap;
    std::vector<YTranslate> yvarMap;
    QSpinBox  *ptlabel;
    QSpinBox  *dim;
    QListWidget *plotSelect;
    std::list<QString> plotList;
    QString   shortFileName;
    QFileInfo dataFileInfo;
  private slots:
    // this is connected to the addplpot button
    void addPlot() { emit requestPlot (dataFileInfo.absoluteFilePath()); }
    void clearPlot();
    void open();
    void removePlot();
    void colorizePlot();
    void print();
    void exportSvg();
  public slots:
    // this is when it is called for the first time with the same data file
    void initPlot(const KNDataFile* dataFile);
    // add a plot with the same data file
    void addPlot(const KNDataFile* mat);
    // is called when the computing thread made a step
    void updatePlot(const KNDataFile* dataFile);
  signals:
    void requestPlot(const QString& fileName);
    void windowClosed();
    // gets emitted when a new file is opened
    // this makes the system open the file
    void openFile(const QString& filename);
    // called when plotting is finished
    void updated();
};

#endif
