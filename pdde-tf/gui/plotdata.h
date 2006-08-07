
#ifndef PLOTDATA_H
#define PLOTDATA_H

#include <QObject>
#include <QWidget>
#include <QPen>
#include <QPolygon>
#include <QList>
#include <QPainterPath>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QFocusEvent>
#include <QGraphicsSceneMouseEvent>
#define PDDESYS_H
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

struct PlotLine
{
	PlotLine( const QPen& p ) : pen(p) { }
	QPen         pen;
	QPainterPath line;
};

struct PlotCircle
{
	PlotCircle( const QPen& pen_, const QRectF& point_ ) : pen(pen_), point(point_) { }
	QPen             pen;
	QRectF           point;
	QVector<QPointF> pos;
};

struct PlotPolygon
{
	PlotPolygon( const QPen& pen_, const QPolygonF& point_ ) : pen(pen_), point(point_) { }
	QPen               pen;
	QPolygonF          point;
	QVector<QPointF>   pos;
};

union PlotItemUnion
{
	PlotLine*    line;
	PlotCircle*  circle;
	PlotPolygon* polygon;
};

struct PlotItem
{
	PlotItemUnion data;
	PlotType      type;
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
	void addPlot( const Vector& v );
	void addPlot( const Vector& v, const char* style );
	void addPlot( const Vector& v1, const Vector& v2 );
	void addPlot( const Vector& v1, const Vector& v2, const char* style );
	void addPlot( const mat4Data& data, PlotXVariable x, PlotYVariable y, int pt, int d, const char* style );
	void clear();

protected:
// 	void paintEvent( QPaintEvent *event );
// 	void resizeEvent ( QResizeEvent * event );
// 	void mousePressEvent ( QGraphicsSceneMouseEvent * event );
// 	void mouseReleaseEvent ( QGraphicsSceneMouseEvent * event );
// 	void mouseMoveEvent ( QGraphicsSceneMouseEvent * event );
// 	void keyPressEvent ( QKeyEvent * event );
// 	void focusInEvent ( QFocusEvent * focusEvent ) { std::cout<<"FocusIn\n"; std::cout.flush(); }
// 	void focusOutEvent ( QFocusEvent * focusEvent ) { std::cout<<"FocusOut\n"; std::cout.flush(); }
	bool event( QEvent* event );

private:
	void addPlotLine( const char* style );
	void addPlotCircle( const char* style );
	void addPlotSquare( const char* style );
	void addPlotUpTriangle( const char* style );
	void dataToGraphics( int n = 1 );
	void getScale( qreal& transx, qreal& transy, qreal& scale );
	void plotStyle( QPen& pen, const char* style );
	QPointF intersect( QPointF p1, QPointF p2 );
	bool contains( QPointF p );
	bool crossbox( QPointF p1, QPointF p2, QPointF& i1, QPointF& i2 );
	void rescaleData();
	void replot();

	QGraphicsScene* PlotScene;
	// geometry
	QList<ViewBox> ZoomHistory;
	QList<ViewBox>::iterator currZoom;
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
	QList<Vector>   DataX;
	QList<Vector>   DataY;
	QList<PlotType> DataType;
	
	QList<PlotItem> Graph;
	
	QVector<QGraphicsTextItem*>    HText;
	QVector<QGraphicsTextItem*>    VText;
	QVector<QGraphicsLineItem*>    TopTicks;
	QVector<QGraphicsLineItem*>    BottomTicks;
	QVector<QGraphicsLineItem*>    LeftTicks;
	QVector<QGraphicsLineItem*>    RightTicks;
	QVector<QGraphicsPathItem*>    PathItems;
	QVector<QGraphicsEllipseItem*> CircleItems;
	QVector<QGraphicsPolygonItem*> PolygonItems;
};

#endif
