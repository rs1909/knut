// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
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

class PlotLine
{
 public:
	PlotLine( const QPen& p ) : pen(p) { }
	QPen           pen;
	QPainterPath   line;
	QGraphicsPathItem *item;
};

class PlotCircle
{
 public:
	PlotCircle( const QPen& pen_, const QRectF& point_ ) : pen(pen_), point(point_) { }
	QPen                 pen;
	QRectF               point;
	std::vector<QPointF> pos;
	std::vector<QGraphicsEllipseItem*> item;
};

class PlotPolygon
{
 public:
	PlotPolygon( const QPen& pen_, const QPolygonF& point_ ) : pen(pen_), point(point_) { }
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
	int           number;
};

struct ViewBox
{
	qreal xmax, xmin, ymax, ymin;
	int xticks, yticks;
};

class PlotData : public QGraphicsScene
{
	Q_OBJECT

public:
	PlotData( QObject *parent = 0 );
	~PlotData( );
	void makeBox( );
	void PlotPaint( );
	void addPlot( const mat4Data* data, PlotXVariable x, PlotYVariable y, int pt, int d, const char* style );
	void clearAll();
	void clear( int n );
	int  nplots();

protected:
	void mousePressEvent ( QGraphicsSceneMouseEvent * event );
	void mouseReleaseEvent ( QGraphicsSceneMouseEvent * event );
	void mouseMoveEvent ( QGraphicsSceneMouseEvent * event );
	void keyPressEvent ( QKeyEvent * event );
	bool event( QEvent* event );

private:
	void addPlotLine( std::list<PlotItem>::iterator& it, const char* style );
	void addPlotPoint( std::list<PlotItem>::iterator& it, const char* style, int type );
	void dataToGraphics( );
	void getScale( qreal& transx, qreal& transy, qreal& scale );
	void plotStyle( QPen& pen, const char* style );
	QPointF intersect( QPointF p1, QPointF p2 );
	inline bool contains( double x, double y );
	bool crossbox( QPointF p1, QPointF p2, QPointF& i1, QPointF& i2 );
	void rescaleData();
	void replot();

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
};

#endif
