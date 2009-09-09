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
#include <QGraphicsRectItem>
class QGraphicsSceneMouseEvent;

#include <vector>
#include <list>

#include "matrix.h"
#include "mat4data.h"

enum PlotXVariable
{
  XNone,
  XLabel,
  XMesh,
  XRealMultiplier,
  XParameter0
};

enum PlotYVariable
{
  YNone,
  YL2Norm,
  YAmplitude,
  YImagMultiplier,
  YAbsMultiplier,
  YProfile,
  YParameter0
};

enum PlotType
{
  PlotLineType,
  PlotCircleType,
  PlotPolygonType
};

enum PlotMarkerStyle
{
  PlotMarkerCircle,
  PlotMarkerSquare,
  PlotMarkerTriangle,
  PlotMarkerCross
};

class PlotLine
{
  public:
    PlotLine(const QPen& p) : pen(p)
    { }
    QPen           pen;
    QPainterPath   line;
    QGraphicsPathItem *item;
};

class PlotCircle
{
  public:
    PlotCircle(const QPen& pen_, const QRectF& point_) : pen(pen_), point(point_)
    { }
    QPen                 pen;
    QRectF               point;
    std::vector<QPointF> pos;
    std::vector<QGraphicsEllipseItem*> item;
};

class PlotPolygon
{
  public:
    PlotPolygon(const QPen& pen_, const QPolygonF& point_) : pen(pen_), point(point_)
    { }
    QPen                   pen;
    QPolygonF              point;
    std::vector<QPointF>   pos;
    std::vector<QGraphicsPolygonItem*> item;
};

union PlotItemUnion
{
  PlotLine*    line;
  PlotCircle*  circle;
  PlotPolygon* polygon;
};

class PlotItem
{
  public:
    PlotItemUnion data;
    PlotType      type;
    Vector        x;
    Vector        y;
    unsigned int  number;
    bool          principal;
};

struct ViewBox
{
  qreal xmax, xmin, ymax, ymin;
  unsigned int xticks, yticks;
};

class PlotData : public QGraphicsScene
{
    Q_OBJECT

  public:
    PlotData(QObject *parent = 0);
    ~PlotData();
    void makeBox();
    void PlotPaint();
    bool addPlot(const mat4Data* data, 
      PlotXVariable x, PlotYVariable y, int pt, int dim);
    void clearAll();
    void clear(unsigned int n);
    int  nplots();
    QColor getColor(unsigned int n);
    void setColor(unsigned int n, QColor& Color);

  protected:
    void mousePressEvent(QGraphicsSceneMouseEvent * event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent * event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent * event);
    void keyPressEvent(QKeyEvent * event);

  private:
    void addPlotLine(std::list<PlotItem>::iterator& it, const QPen& pen, bool p);
    void addPlotPoint(std::list<PlotItem>::iterator& it, const QPen& pen, PlotMarkerStyle type, bool p);
    void dataToGraphics();
    void getScale(qreal& transx, qreal& transy, qreal& scale);
//    void plotStyle(QPen& pen, const char* style);
    QPointF intersect(QPointF p1, QPointF p2);
    inline bool contains(double x, double y);
    bool crossbox(QPointF p1, QPointF p2, QPointF& i1, QPointF& i2);
    void rescaleData();
    void replot();
    void clearAxes();
    void labelColor();

    // geometry
    std::list<ViewBox> ZoomHistory;
    std::list<ViewBox>::iterator currZoom;
    // relative geometry
    const qreal AspectRatio;
    const qreal plotXSize, plotYSize;
    const int   FontSize;
    // mouse
    QPointF           mouseBegin;
    QPointF           mouseEnd;
    QPointF           mouseMove;
    QGraphicsRectItem selection;
    // data
    std::list<PlotItem> Graph;

    // the plot box
    QGraphicsRectItem               *Box;
    std::vector<QGraphicsTextItem*>  HText;
    std::vector<QGraphicsTextItem*>  VText;
    std::vector<QGraphicsLineItem*>  TopTicks;
    std::vector<QGraphicsLineItem*>  BottomTicks;
    std::vector<QGraphicsLineItem*>  LeftTicks;
    std::vector<QGraphicsLineItem*>  RightTicks;
    
    // for the axes
    std::vector<QString> XCoordText;
    std::vector<QString> YCoordText;
    std::vector<QGraphicsTextItem*> XCoordTextItems;
    std::vector<QGraphicsTextItem*> YCoordTextItems;
    std::map<PlotXVariable,QString> XCoordMap;
    std::map<PlotYVariable,QString> YCoordMap;
};

#endif