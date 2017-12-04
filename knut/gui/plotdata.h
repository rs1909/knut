// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef PLOTDATA_H
#define PLOTDATA_H

#include <QPen>
#include <QPolygon>
#include <QPainterPath>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QPainter>
#include <QGraphicsRectItem>
#include <QLinkedList>
#include <QFileInfo>
class QGraphicsSceneMouseEvent;

#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <variant>

#include "matrix.h"
#include "mat4data.h"

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

enum PlotType
{
  PlotLineType,
  PlotCircleType,
  PlotPolygonType
};

enum PlotDataType
{
  PlotBasicData,
  PlotStability,
  PlotAuxiliary
};

enum PlotMarkerStyle
{
  PlotMarkerCircle,
  PlotMarkerSquare,
  PlotMarkerTriangle,
  PlotMarkerCross
};

class PlotPolyLine
{
  public:
    PlotPolyLine() { }
    // need to avoid at all costs
//    PlotPolyLine(const PlotPolyLine& pl)
    PlotPolyLine(const QPen& pen_) : pen(pen_), item(nullptr) { }
    
    QPen               pen;
    QPainterPath       path;
    std::unique_ptr<QGraphicsPathItem> item;
};

class PlotLine : public std::list<PlotPolyLine>
{
  public:
    PlotLine() { }
    PlotLine(const QPen& p) : pen(p) { }
    void clear()
    {
      std::list<PlotPolyLine>::clear();
    }
    QPen pen;
};

class PlotCircle
{
  public:
    PlotCircle(const QPen& pen_, const QRectF& point_, bool scale_ = false) 
      : pen(pen_), point(point_), scale(scale_) { }
    QPen                 pen;
    QRectF               point;
    QRectF               scaledPoint;
    std::list<QPointF> pos;
    std::list<std::unique_ptr<QGraphicsEllipseItem> > item;
    bool scale;
};

// use std::variant instead
class PlotPolygon
{
  public:
    PlotPolygon(const QPen& pen_, const QPolygonF& point_) : pen(pen_), point(point_)
    { }
    QPen                   pen;
    QPolygonF              point;
    std::list<QPointF>   pos;
    std::list<std::unique_ptr<QGraphicsPolygonItem> > item;
};

// union PlotItemUnion
// {
//   std::shared_ptr<PlotLine>    line;
//   std::shared_ptr<PlotCircle>  circle;
//   std::shared_ptr<PlotPolygon> polygon;
//   PlotItemUnion() : line(nullptr) {}
//   ~PlotItemUnion() {}
// };

class PlotItem
{
  public:
    PlotItem(const KNDataFile* mat, PlotDataType type,
             PlotXVariable xv, PlotYVariable yv, size_t pt, size_t dim)
     : dataType(type),
     sourcePath(QFileInfo(QString::fromStdString(mat->getFileName())).canonicalFilePath()),
     varX(xv), varY(yv), point(pt), dimension(dim) {}
    bool isFrom(const KNDataFile* mat) { return sourcePath == QFileInfo(QString::fromStdString(mat->getFileName())).canonicalFilePath(); }
    bool isUnitCircle() { return (varX == XRealMultiplier) && (varY == YImagMultiplier); } 
    std::variant<PlotLine, PlotCircle, PlotPolygon> data;
    PlotType      type;
    PlotDataType  dataType;
    KNVector      x;
    KNVector      y;
    size_t        number;
    bool          principal;
    QString       sourcePath;
    const PlotXVariable varX;
    const PlotYVariable varY;
    const size_t  point;
    const size_t  dimension;
};

struct ViewBox
{
  qreal xmax {1e6}, xmin {-1e-6}, ymax {1e6}, ymin {1e-6};
  size_t xticks {1}, yticks {1};
};

class PlotData : public QGraphicsScene
{
    Q_OBJECT

  public:
    PlotData(QObject *parent = nullptr);
    ~PlotData() override;
    bool addPlot(const KNDataFile* mat, 
      PlotXVariable x, PlotYVariable y, size_t pt, size_t dim);
    void clearAll();
    void clear(size_t n);
    size_t nplots();
    QColor getColor(size_t n);
    void setColor(size_t n, QColor& Color);
    void updatePlot(const KNDataFile* mat);
    
  public slots:
    void setAxes(int b)
    {
      bool xl = (b & 1) != 0;
      bool yl = (b & 2) != 0;
      if ((xlog != xl)||(ylog != yl))
      {
       xlog = xl;
       ylog = yl;
       dataToGraphics(Graph.begin(), Graph.end());
      }
    }
    void setXSize(int size);
       
  protected:
    void mousePressEvent(QGraphicsSceneMouseEvent * event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent * event) override;
    void mouseMoveEvent(QGraphicsSceneMouseEvent * event) override;
    void keyPressEvent(QKeyEvent * event) override;

  private:
    void addPlotLine(std::list<std::unique_ptr<PlotItem>>::iterator it, const QPen& pen, bool p, bool s = true);
    void addPlotPoint(std::list<std::unique_ptr<PlotItem>>::iterator it, const QPen& pen, 
                      PlotMarkerStyle type, bool p, qreal radius = 3.0, bool scale = false);
    void dataToGraphics(std::list<std::unique_ptr<PlotItem>>::const_iterator begin,
                        std::list<std::unique_ptr<PlotItem>>::const_iterator end);
    void getScale(qreal& transx, qreal& transy, qreal& scale);
    QPointF intersect(QPointF p1, QPointF p2);
    inline bool contains(double x, double y);
    bool crossbox(QPointF p1, QPointF p2, QPointF& i1, QPointF& i2);
    void rescaleData(std::list<std::unique_ptr<PlotItem>>::const_iterator begin,
                     std::list<std::unique_ptr<PlotItem>>::const_iterator end);
    QGraphicsRectItem* makeBox();
    void PlotPaint(std::list<std::unique_ptr<PlotItem>>::const_iterator begin,
                   std::list<std::unique_ptr<PlotItem>>::const_iterator end, bool zoom);
//     void replot();
    void clearAxes();
    void labelColor();
    inline double xcoord(double x) { if (xlog) { if (x>0) return log10(fabs(x)); else return 0; } else return x; }
    inline double ycoord(double y) { if (ylog) { if (y>0) return log10(fabs(y)); else return 0; } else return y; }

    // geometry
    std::list<ViewBox> ZoomHistory;
    std::list<ViewBox>::iterator currZoom;
    bool xlog;
    bool ylog;
    static const ViewBox defaultBox;
    // relative geometry
    qreal AspectRatio;
    qreal plotXSize, plotYSize;
    const int   FontSize;
    // mouse
    QPointF           mouseBegin;
    QPointF           mouseEnd;
    QPointF           mouseMove;
    QGraphicsRectItem selection;
    // data
    std::list<std::unique_ptr<PlotItem>> Graph;
    // unitcircle
    std::unique_ptr<QGraphicsEllipseItem> unitCircleItem;
    int unitCircleCount;

    // the plot box
    std::unique_ptr<QGraphicsRectItem>              Box;
    QPainterPath                     clipBox;
    std::vector<std::unique_ptr<QGraphicsTextItem> >  HText;
    std::vector<std::unique_ptr<QGraphicsTextItem> >  VText;
    std::vector<std::unique_ptr<QGraphicsLineItem> >  TopTicks;
    std::vector<std::unique_ptr<QGraphicsLineItem> >  BottomTicks;
    std::vector<std::unique_ptr<QGraphicsLineItem> >  LeftTicks;
    std::vector<std::unique_ptr<QGraphicsLineItem> >  RightTicks;
    
    // for the axes
    std::vector<QString> XCoordText;
    std::vector<QString> YCoordText;
    std::vector<std::unique_ptr<QGraphicsTextItem> > XCoordTextItems;
    std::vector<std::unique_ptr<QGraphicsTextItem> > YCoordTextItems;
    std::map<PlotXVariable,QString> XCoordMap;
    std::map<PlotYVariable,QString> YCoordMap;
};

#endif
