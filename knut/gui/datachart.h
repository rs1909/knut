// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2017 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include <QObject>
#include <QString>
#include <QList>
#include <QFileInfo>
#include <map>

#ifndef DATACHART_H
#define DATACHART_H

class KNDataFile;
#include <QColor>
#include <QtCharts/QChart>
#include <QtCharts/QXYSeries>
#include <QtCharts/QValueAxis>

using namespace QtCharts;

enum PlotXVariable
{
  XNone,
  // only one row of data is used
  XMesh,           // -> YProfile
  XRealMultiplier, // -> YImagMultiplier
  XSeparator,
  // From here all rows are necessary
  // can be combined with any of the X > XSeparator
  XLabel,
  XParameter0
};

enum PlotYVariable
{
  YNone,
  // only one row of data is used
  YImagMultiplier, // -> XRealMultiplier
  YProfile,        // - > XMesh
  YAbsMultiplier,
  YSeparator,
  // From here all rows are necessary
  // can be combined with any of the Y > YSeparator
  YL2Norm,
  YAmplitude,
  YPoincare,
  YParameter0
};

class DataChart : public QObject
{
    Q_OBJECT

  public:
    DataChart(QObject *parent = 0);
    ~DataChart();
    bool addPlot(const KNDataFile* mat, 
      PlotXVariable x, PlotYVariable y, size_t pt, size_t dim);
    void clearAll();
    void clear(size_t n);
    size_t nplots();
    QColor getColor(size_t n);
    void setColor(size_t n, QColor& Color);
    void updatePlot(const KNDataFile* mat);
    QChart* getChart() { return chart; }
    
  private:
    class OnePlot
    {
      public:
        OnePlot ( QList<QXYSeries*>* graph_,
                    QString filename_,
                    PlotXVariable xvar_,
                    PlotYVariable yvar_, size_t pt_, size_t dim_ ) :
            graph(graph_), file(filename_), xvar(xvar_), yvar(yvar_), pt(pt_), dim(dim_) { }
        QList<QXYSeries*>* graph;
        QFileInfo file;
        PlotXVariable xvar;
        PlotYVariable yvar;
        size_t pt;
        size_t dim;
    };

    void addToChart(OnePlot& lst);

    QChart* chart;
    QValueAxis *hAxis;
    QValueAxis *vAxis;
    QList< OnePlot > Graph;
    std::map<PlotXVariable,QString> XCoordMap;
    std::map<PlotYVariable,QString> YCoordMap;

};

#endif // DATACHART_H
