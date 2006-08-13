#include <QtGui>

#include "plotdata.h"
#include <cmath>
#include <cfloat>

PlotData::PlotData(QObject *parent) :
	QGraphicsScene(parent),
	AspectRatio( 4.0/3.0 ),
	plotXSize  ( 540 ), plotYSize ( plotXSize / AspectRatio ),
	FontSize ( 12 )
{
	const ViewBox cvb = { 1.0, -1.0, 1.0, -1.0, 1, 1 };
	ZoomHistory.append( cvb );
	currZoom = ZoomHistory.begin();
	this->addItem( &selection );
	selection.setVisible(false);
	makeBox();
}

PlotData::~PlotData( )
{
	clear();
}

void PlotData::clear()
{
	QList<PlotItem>::iterator i;
	QList<Vector>::iterator j1;
	QList<Vector>::iterator j2;
	for( i = Graph.begin(), j1 = DataX.begin(), j2=DataY.begin(); i != Graph.end(); i++, j1++, j2++ )
	{
		if( (*i).type == PlotLineType )
		{
			delete (*i).data.line;
		}
		else if( (*i).type == PlotCircleType )
		{
			delete (*i).data.circle;
		}
		else if( (*i).type == PlotPolygonType )
		{
			delete (*i).data.polygon;
		}
		else
		{
			std::cout<<"Something wrong\n";
		}
	}
	DataX.clear();
	DataY.clear();
	Graph.clear();
	const ViewBox cvb = *currZoom;
	ZoomHistory.clear();
	ZoomHistory.append( cvb );
	currZoom = ZoomHistory.begin();
}

/// this computes the minimum and maximum value of an axis
static inline void adjustAxis( qreal& min, qreal& max, int& numTicks )
{
	const int MinTicks = 4;
	qreal grossStep = (max - min) / MinTicks;
	qreal step = pow(10,floor(log10(grossStep)));

	if( 5 * step < grossStep ) step *= 5;
	else if( 2 * step < grossStep ) step *= 2;

	numTicks = (int)(ceil(max / step) - floor(min / step));
	min = floor(min / step) * step;
	max = ceil(max / step) * step;
}

void PlotData::plotStyle( QPen& pen, const char* style )
{
	while( *style != '\0' )
	{
		switch( style[0] )
		{
			case 'b': pen.setColor("blue"); break;
			case 'g': pen.setColor("green"); break;
			case 'r': pen.setColor("red"); break;
			case 'c': pen.setColor("cyan"); break;
			case 'm': pen.setColor("magenta"); break;
			case 'y': pen.setColor("yellow"); break;
			case 'k': pen.setColor("black"); break;
			case '-': 
				switch( style[1] )
				{
					case '-': pen.setStyle(Qt::DashLine); style++; break;
					case '.': pen.setStyle(Qt::DashDotLine); style++; break;
					default:  pen.setStyle(Qt::SolidLine); break;
				}
				break;
			case ':': pen.setStyle(Qt::DotLine); break;
			default: std::cout<<"style Error\n"; break;
		}
		style++;
	}
}

void PlotData::addPlot( const Vector& v1 )
{
	const Vector v0;
	addPlot( v0, v1, "b" );
}

void PlotData::addPlot( const Vector& v1, const char* style )
{
	const Vector v0;
	addPlot( v0, v1, style );
}

void PlotData::addPlot( const Vector& v1, const Vector& v2 )
{
	addPlot( v1, v2, "b" );
}

void PlotData::addPlot( const Vector& v1, const Vector& v2, const char* style )
{
	if( ( v1.Size() != 0 )&&( v1.Size() != v2.Size() ) ) { std::cout<<"error\n"; exit(-1); }
	
	DataX.push_back( Vector() );
	DataY.push_back( Vector() );
	if( v1.Size() != 0 ) DataX.last().Init(v1);
	else
	{
		DataX.last().Init( v2.Size() );
		for( int r=0; r < v2.Size(); r++ ) DataX.last()(r) = r;
	}
	DataY.last().Init(v2);
	addPlotLine( style );
	dataToGraphics();
}

void PlotData::addPlotLine( const char* style )
{
	PlotItem t_pl = { {0}, PlotLineType };
	Graph.push_back( t_pl );
	Graph.last().data.line = new PlotLine( QPen( QColor( "blue" ) ) );
	Graph.last().data.line->pen.setWidthF( 1 );
	plotStyle( Graph.last().data.line->pen, style );
}

void PlotData::addPlotCircle( const char* style )
{
	PlotItem t_pl = { {0}, PlotCircleType };
	Graph.push_back( t_pl );
	Graph.last().data.circle = new PlotCircle( QPen( QColor( "blue" ) ), QRectF( -3.0, -3.0, 6.0 ,6.0 ) );
	Graph.last().data.circle->pen.setWidthF( 1 );
	plotStyle( Graph.last().data.circle->pen, style );
}

void PlotData::addPlotSquare( const char* style )
{
	PlotItem t_pl = { {0}, PlotPolygonType };
	Graph.push_back( t_pl );
	QPolygonF pl(4);
	pl[0] = QPointF( -3.0, -3.0 );
	pl[1] = QPointF( -3.0, 3.0 );
	pl[2] = QPointF( 3.0, 3.0 );
	pl[3] = QPointF( 3.0, -3.0 );
	Graph.last().data.polygon = new PlotPolygon( QPen( QColor( "blue" ) ), pl );
	Graph.last().data.polygon->pen.setWidthF( 1 );
	plotStyle( Graph.last().data.polygon->pen, style );
}

void PlotData::addPlotUpTriangle( const char* style )
{
	PlotItem t_pl = { {0}, PlotPolygonType };
	Graph.push_back( t_pl );
	QPolygonF pl(3);
	pl[0] = QPointF( -3.0, -3.0 );
	pl[1] = QPointF( 0.0, 3.0 );
	pl[2] = QPointF( 3.0, -3.0 );
	Graph.last().data.polygon = new PlotPolygon( QPen( QColor( "blue" ) ), pl );
	Graph.last().data.polygon->pen.setWidthF( 1 );
	plotStyle( Graph.last().data.polygon->pen, style );
}

// this only called by addPlot...
void PlotData::dataToGraphics( )
{
	if( ZoomHistory.begin() == currZoom )
	{
		ViewBox cvb = { -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, 1, 1 };
		QList<Vector>::const_iterator j1;
		QList<Vector>::const_iterator j2;
		for( j1 = DataX.begin(), j2 = DataY.begin(); j1 != DataX.end(); j1++, j2++ )
		{
			for( int k = 0; k < (*j1).Size(); k++ )
			{
				if( (*j1)(k) > cvb.xmax ) cvb.xmax = (*j1)(k);
				if( (*j1)(k) < cvb.xmin ) cvb.xmin = (*j1)(k);
				if( (*j2)(k) > cvb.ymax ) cvb.ymax = (*j2)(k);
				if( (*j2)(k) < cvb.ymin ) cvb.ymin = (*j2)(k);
			}
		}
		adjustAxis( cvb.xmin, cvb.xmax, cvb.xticks );
		adjustAxis( cvb.ymin, cvb.ymax, cvb.yticks );
		*currZoom = cvb;
	}
	rescaleData();
	PlotPaint( );
}

void PlotData::addPlot( const mat4Data& data, PlotXVariable x, PlotYVariable y, int pt, int dim, const char* style )
{
	int xadded = 0;
	int yadded = 0;
	// add stability
	QVector<int> bifidx;
	QVector<int> biftype;
	if( x == XLabel || x >= XParameter0 )
	{
		int k, k_p = 0;
		do{
			k = data.getNextBifurcation( k_p, 1 );
			if( k != -1 )
			{
				bifidx.push_back(k);
				biftype.push_back( data.getBifurcationType( k, 1 ) );
				k_p = k;
			}
		}while( k != -1 );
	}
	// mindig az Y utan kovetkezik csak addPlot...()
	if( x >= XParameter0 )
	{
		DataX.push_back( Vector() );
		DataX.last().Init( data.getNPoints() );
		for( int i = 0; i < data.getNPoints(); i++ )
		{
			DataX.last()(i) = data.getPar( i, x-XParameter0 );
		}
		++xadded;
		for( int i = 0; i < bifidx.size(); i++ )
		{
			DataX.push_back( Vector() );
			DataX.last().Init( 1 );
			DataX.last()(0) = data.getPar( bifidx[i], x-XParameter0 );
			++xadded;
		}
	}
	if( x == XLabel && y != YAbsMultiplier && y != YProfile )
	{
		DataX.push_back( Vector() );
		DataX.last().Init( data.getNPoints() );
		for( int i = 0; i < data.getNPoints(); i++ )
		{
			DataX.last()(i) = i;
		}
		++xadded;
		for( int i = 0; i < bifidx.size(); i++ )
		{
			DataX.push_back( Vector() );
			DataX.last().Init( 1 );
			DataX.last()(0) = bifidx[i];
			++xadded;
		}
	}
	if( x == XMesh && y == YProfile )
	{
		DataX.push_back( Vector() );
		DataY.push_back( Vector() );
		DataX.last().Init( data.getMeshLength() );
		DataY.last().Init( data.getMeshLength() );
		for( int i = 0; i < data.getMeshLength(); i++ )
		{
			DataX.last()(i) = data.getMesh( pt, i );
			DataY.last()(i) = data.getProfile( pt, dim, i );
		}
		++xadded; ++yadded;
		addPlotLine( style );
		dataToGraphics();
	}
	if( y >= YParameter0 && (x != XMesh && x != XRealMultiplier) )
	{
		DataY.push_back( Vector() );
		DataY.last().Init( data.getNPoints() );
		for( int i = 0; i < data.getNPoints(); i++ )
		{
			DataY.last()(i) = data.getPar( i, y-YParameter0 );
		}
		++yadded;
		addPlotLine( style );
		for( int i = 0; i < bifidx.size(); i++ )
		{
			DataY.push_back( Vector() );
			DataY.last().Init( 1 );
			DataY.last()(0) = data.getPar( bifidx[i], y-YParameter0 );
			++yadded;
			addPlotCircle( style );
		}
		dataToGraphics();
	}
	if( y == YAmplitude && (x != XMesh && x != XRealMultiplier) )
	{
		DataY.push_back( Vector() );
		DataY.last().Init( data.getNPoints() );
		for( int i = 0; i < data.getNPoints(); i++ )
		{
			double min = DBL_MAX;
			double max = -DBL_MAX;
			for( int j = 0; j < data.getMeshLength(); ++j )
			{
				for( int k = 0; k < data.getNDim(); ++k )
				{
					if( min > data.getProfile( i, k, j ) ) min = data.getProfile( i, k, j );
					if( max < data.getProfile( i, k, j ) ) max = data.getProfile( i, k, j );
				}
			}
			DataY.last()(i) = max - min;
		}
		++yadded;
		addPlotLine( style );
		for( int i = 0; i < bifidx.size(); i++ )
		{
			DataY.push_back( Vector() );
			DataY.last().Init( 1 );
			double min = DBL_MAX;
			double max = -DBL_MAX;
			for( int j = 0; j < data.getMeshLength(); ++j )
			{
				for( int k = 0; k < data.getNDim(); ++k )
				{
					if( min > data.getProfile( bifidx[i], k, j ) ) min = data.getProfile( bifidx[i], k, j );
					if( max < data.getProfile( bifidx[i], k, j ) ) max = data.getProfile( bifidx[i], k, j );
				}
			}
			DataY.last()(0) = max - min;
			++yadded;
			addPlotCircle( style );
		}
		dataToGraphics();
	}
	if( x == XRealMultiplier && y == YImagMultiplier )
	{
		// plot the unit circel!!!
		const int segments = 100;
		DataX.push_back( Vector() );
		DataY.push_back( Vector() );
		DataX.last().Init( segments );
		DataY.last().Init( segments );
		for( int i = 0; i < segments; ++i )
		{
			DataX.last()(i) = sin( i*2.0*M_PI/(segments-1) );
			DataY.last()(i) = cos( i*2.0*M_PI/(segments-1) );
		}
		addPlotLine( "k" );
		// plotting the multipliers
		DataX.push_back( Vector() );
		DataY.push_back( Vector() );
		DataX.last().Init( data.getNMul() );
		DataY.last().Init( data.getNMul() );
		for( int i = 0; i < data.getNMul(); i++ )
		{
			DataX.last()(i) = data.getMulRe( pt, i );
			DataY.last()(i) = data.getMulIm( pt, i );
		}
		++xadded; ++yadded;
		addPlotCircle( style );
		dataToGraphics();
	}
	if( x == XLabel && y == YAbsMultiplier )
	{
		for( int r = 0; r < data.getNMul(); r++ )
		{
			DataX.push_back( Vector() );
			DataY.push_back( Vector() );
			DataX.last().Init( data.getNPoints() );
			DataY.last().Init( data.getNPoints() );
			for( int i = 0; i < data.getNPoints(); i++ )
			{
				DataX.last()(i) = i;
				DataY.last()(i) = sqrt(data.getMulRe(i,r)*data.getMulRe(i,r)+data.getMulIm(i,r)*data.getMulIm(i,r));
			}
			addPlotCircle( style );
		}
		dataToGraphics();
	}
	if( xadded != yadded )
	{
		for( int i = 0; i < xadded; ++i ) DataX.erase( DataX.end()-1 );
		for( int i = 0; i < yadded; ++i ) DataY.erase( DataY.end()-1 );
	}
}

///
/// The inlines
///

static inline qreal itr( qreal xx, qreal x1, qreal x2, qreal y1, qreal y2 )
{
	return ( (x2-xx)*y1 + (xx-x1)*y2 ) / ( x2-x1 );
}

static inline int intpos( qreal x, qreal min, qreal max )
{
	if( x < min ) return -1;
	else if( x < max ) return 0;
	else return 1;
}

static inline bool crossing( int x1, int y1, int x2, int y2 )
{
	if( (x1*x2 == 1) || (y1*y2 == 1) ) return false;
	else return true;
}

///
/// End inlines
///

// /-----------------|
// |  1  |  4  |  7  |
// |-----+-----+-----|
// |  2  |  5  |  8  |
// |-----+-----+-----|
// |  3  |  6  |  9  |
// \-----------------|

/// finds out how the viewbox is intersected if
/// p1 is inside the box,
/// p2 is outside the box
QPointF PlotData::intersect( QPointF p1, QPointF p2 )
{
	// p1 is in the box, p2 is out of the box
	const ViewBox cvb = { plotXSize, 0.0, plotYSize, 0.0, 0, 0 };
	qreal ycr, xcr;
	if( p2.x() > cvb.xmax )
	{
		ycr = itr( cvb.xmax, p1.x(), p2.x(), p1.y(), p2.y() );
		if( ycr > cvb.ymax ) { ycr = cvb.ymax; xcr = itr( cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x() ); }
		else if( ycr > cvb.ymin ) xcr = cvb.xmax;
		else { ycr = cvb.ymin; xcr = itr( cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x() ); }
	}
	else if( p2.x() > cvb.xmin )
	{
		if( p2.y() > cvb.ymax ) { xcr = itr( cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x() ); ycr = cvb.ymax; }
		else if( p2.y() > cvb.ymin ) { xcr = p2.x(), ycr = p2.y(); }
		else { xcr = itr( cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x() ); ycr = cvb.ymin; }
	}else
	{
		ycr = itr( cvb.xmin, p1.x(), p2.x(), p1.y(), p2.y() );
		if( ycr > cvb.ymax ) { ycr = cvb.ymax; xcr = itr( cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x() ); }
		else if( ycr > cvb.ymin ) xcr = cvb.xmin;
		else { ycr = cvb.ymin; xcr = itr( cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x() ); }
	}
	return QPointF( xcr, ycr );
}

inline bool PlotData::contains( double x, double y )
{
	return (x >= 0.0) && (y >= 0.0) && (x <= plotXSize) && (y <= plotYSize);
}

bool PlotData::crossbox( QPointF p1, QPointF p2, QPointF& i1, QPointF& i2 )
{
	const ViewBox cvb = { plotXSize, 0.0, plotYSize, 0.0, 0, 0 };
	
	qreal ycrmin = itr( cvb.xmin, p1.x(), p2.x(), p1.y(), p2.y() );
	qreal ycrmax = itr( cvb.xmax, p1.x(), p2.x(), p1.y(), p2.y() );
	qreal xcrmin = itr( cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x() );
	qreal xcrmax = itr( cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x() );
	int id = 0;
	qreal x[4], y[4]; // security reasons
	if( (ycrmin > cvb.ymin) && (ycrmin < cvb.ymax) )
	{
		x[id] = cvb.xmin; y[id] = ycrmin; ++id;
	}
	if( (ycrmax > cvb.ymin) && (ycrmax < cvb.ymax) )
	{
		x[id] = cvb.xmax; y[id] = ycrmax; ++id;
	}
	if( (xcrmin > cvb.xmin) && (xcrmin < cvb.xmax) )
	{
		y[id] = cvb.ymin; x[id] = xcrmin; ++id;
	}
	if( (xcrmax > cvb.xmin) && (xcrmax < cvb.xmax) )
	{
		y[id] = cvb.ymax; x[id] = xcrmax; ++id;
	}
	if( id == 0 ) return false;
	else if( id == 2 )
	{
		i1 = QPointF( x[0], y[0] );
		i2 = QPointF( x[1], y[1] );
		return true;
	}else{ std::cout<<"GEBASZ van "<<id<<"\n"; return false; }
}

void PlotData::rescaleData()
{
	ViewBox cvb = *currZoom;

	const double xscale = plotXSize/(cvb.xmax-cvb.xmin);
	const double yscale = plotYSize/(cvb.ymax-cvb.ymin);
	
	// rescaling all the data
	QList<PlotItem>::iterator i;
	QList<Vector>::iterator j1;
	QList<Vector>::iterator j2;
	for( i = Graph.begin(), j1 = DataX.begin(), j2=DataY.begin(); i != Graph.end(); i++, j1++, j2++ )
	{
		if( (*i).type == PlotLineType )
		{
			(*i).data.line->line = QPainterPath();
			int x = 0, y = 0, prx = 0, pry = 0;
			bool pr = true;
			if( (*j1).Size() != (*j2).Size() )
			{
				std::cout<<"DataX DataY Sizes differ\n";
				return;
			}
			for( int k = 0; k < (*j1).Size(); k++ )
			{
				x = intpos( (*j1)(k), cvb.xmin, cvb.xmax );
				y = intpos( (*j2)(k), cvb.ymin, cvb.ymax );
				if( crossing( prx, pry, x, y ) )
				{
					if( (x == 0) && (y == 0) )
					{
						if( !pr )
						{
							// the previous was not in the viewport
							(*i).data.line->line.moveTo( intersect( QPointF( xscale*((*j1)(k)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k)) ), 
																				QPointF( xscale*((*j1)(k-1)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k-1)) ) ) );
						}
						if( k == 0 )
						{
							(*i).data.line->line.moveTo( xscale*((*j1)(k)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k)) );
						}else
						{
							(*i).data.line->line.lineTo( xscale*((*j1)(k)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k)) );
						}
						pr = true;
					}else
					{
						if( k != 0 )
						{
							if( pr )
							{
								(*i).data.line->line.lineTo( intersect( QPointF( xscale*((*j1)(k-1)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k-1)) ), 
																	QPointF( xscale*((*j1)(k)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k)) ) ) );
							}else
							{
								QPointF pt1, pt2;
								if( crossbox( QPointF( xscale*((*j1)(k-1)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k-1)) ),
												QPointF( xscale*((*j1)(k)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k)) ),
												pt1, pt2 ) )
								{
									(*i).data.line->line.moveTo( pt1 );
									(*i).data.line->line.lineTo( pt2 );
								}
							}
						}
						pr = false;
					}
				}
				prx = x; pry = y;
			}
		}
		if( (*i).type == PlotCircleType )
		{
			(*i).data.circle->pos.clear();
			for( int k = 0; k < (*j1).Size(); k++ )
			{
				const QPointF pt = QPointF( xscale*((*j1)(k)-cvb.xmin), yscale*(cvb.ymax-(*j2)(k)) );
				if( contains( pt.x(), pt.y() ) ) (*i).data.circle->pos.push_back( pt );
			}
		}
		if( (*i).type == PlotPolygonType )
		{
			(*i).data.polygon->pos.clear();
			for( int k = 0; k < (*j1).Size(); k++ )
			{
				const QPointF pt = QPointF(  );
				if( contains( pt.x(), pt.y() ) ) (*i).data.polygon->pos.push_back( pt );
			}
		}
	}
}

bool PlotData::event( QEvent* ev )
{
	QGraphicsSceneMouseEvent *event;
	if( (event=dynamic_cast<QGraphicsSceneMouseEvent*>(ev)) != 0)
	{
		if( event->type() == QEvent::GraphicsSceneMouseMove )
		{
			if( event->buttons() == Qt::RightButton )
			{
				mouseMove = QPointF( event->scenePos() );
				if( sceneRect().contains( mouseMove ) )
				{
					selection.setRect( QRectF(mouseBegin.x(), mouseBegin.y(), mouseMove.x() - mouseBegin.x() , mouseMove.y() - mouseBegin.y() ).normalized() );
					selection.setVisible(true);
					selection.update();
				}
				event->accept();
				return true;
			}
		}
		if( event->type() == QEvent::GraphicsSceneMousePress )
		{
			if( event->button() == Qt::RightButton )
			{
				mouseBegin = QPointF( event->scenePos() );
				mouseMove = mouseBegin;
				event->accept();
				return true;
			}
		}
		if( event->type() == QEvent::GraphicsSceneMouseRelease )
		{
			if( event->button() == Qt::RightButton )
			{
				mouseEnd = QPointF( event->scenePos() );
				if ( mouseEnd == mouseBegin ) return true;
				
				ViewBox cvb = *currZoom, newvb;
				const double xscale = plotXSize/(cvb.xmax-cvb.xmin);
				const double yscale = plotYSize/(cvb.ymax-cvb.ymin);
		
				newvb.xmin = qMin( mouseBegin.x(), mouseEnd.x() )/xscale + cvb.xmin;
				newvb.xmax = qMax( mouseBegin.x(), mouseEnd.x() )/xscale + cvb.xmin;
				newvb.ymin = cvb.ymax - qMax( mouseBegin.y(), mouseEnd.y() )/yscale;
				newvb.ymax = cvb.ymax - qMin( mouseBegin.y(), mouseEnd.y() )/yscale;
				adjustAxis( newvb.xmin, newvb.xmax, newvb.xticks );
				adjustAxis( newvb.ymin, newvb.ymax, newvb.yticks );
				// remove the next zoom levels
				if( currZoom != ZoomHistory.end() ) ZoomHistory.erase(currZoom+1,ZoomHistory.end());
				// add the new level
				ZoomHistory.append( newvb );
				currZoom = ZoomHistory.end()-1;
				dataToGraphics();
				selection.setVisible(false);
				update( selection.boundingRect().normalized() );
				event->accept();
				return true;
			}
		}
	}
	QKeyEvent *key;
	if( (key=dynamic_cast<QKeyEvent*>(ev)) != 0 )
	{
		if( key->type() == QEvent::KeyPress )
		{
			if( key->key() == Qt::Key_P )
			{
				if( ZoomHistory.constBegin() != currZoom )
				{
					--currZoom;
					dataToGraphics();
				}
				key->accept();
				return true;
			}
			if( key->key() == Qt::Key_N )
			{
				if( ZoomHistory.constEnd()-1 != currZoom )
				{
					++currZoom;
					dataToGraphics();
				}
				key->accept();
				return true;
			}
		}
	}
	std::cout.flush();
	return false;
}

void PlotData::makeBox( )
{
	ViewBox cvb = *currZoom;
	// drawing the box
	this->addRect( QRectF( 0.0, 0.0, plotXSize, plotYSize ), QPen( QBrush(Qt::SolidPattern), 2.0 ) );
	
	// drawing the ticks and tickmarks
	for( int i=0; i < HText.size(); i++ )
	{
		this->removeItem( BottomTicks[i] );
		this->removeItem( TopTicks[i] );
		this->removeItem( HText[i] );
		delete BottomTicks[i];
		delete TopTicks[i];
		delete HText[i];
	}
	BottomTicks.clear();
	TopTicks.clear();
	HText.clear();
	for( int i=0; i < cvb.xticks+1; i++ )
	{
		BottomTicks.push_back( new QGraphicsLineItem );
		BottomTicks[i]->setLine( plotXSize * i / cvb.xticks, 0.0, plotXSize * i/cvb.xticks, 5.0 );
		BottomTicks[i]->setPen( QPen( QBrush(Qt::SolidPattern), 1.0 ) );
		TopTicks.push_back( new QGraphicsLineItem );
		TopTicks[i]->setLine( plotXSize * i / cvb.xticks, plotYSize, plotXSize * i/cvb.xticks, plotYSize - 5.0 );
		TopTicks[i]->setPen( QPen( QBrush(Qt::SolidPattern), 1.0 ) );
		HText.push_back( new QGraphicsTextItem );
		HText[i]->setPlainText( QString::number( cvb.xmin + (cvb.xmax-cvb.xmin)*i/cvb.xticks ) );
		HText[i]->setFont( QFont("Arial", FontSize) );
		QRectF b = HText[i]->boundingRect().normalized();
		HText[i]->setPos( plotXSize * i/cvb.xticks - b.width()/2.0, plotYSize /*- b.height()*/ );
		this->addItem( TopTicks[i] );
		this->addItem( BottomTicks[i] );
		this->addItem( HText[i] );
	}
	for( int i=0; i < VText.size(); i++ )
	{
		this->removeItem( LeftTicks[i] );
		this->removeItem( RightTicks[i] );
		this->removeItem( VText[i] );
		delete LeftTicks[i];
		delete RightTicks[i];
		delete VText[i];
	}
	LeftTicks.clear();
	RightTicks.clear();
	VText.clear();
	for( int i=0; i < cvb.yticks+1; i++ )
	{
		LeftTicks.push_back( new QGraphicsLineItem );
		LeftTicks[i]->setLine( 0.0, plotYSize * i/cvb.yticks, 5.0, plotYSize * i/cvb.yticks );
		LeftTicks[i]->setPen( QPen( QBrush(Qt::SolidPattern), 1.0 ) );
		RightTicks.push_back( new QGraphicsLineItem );
		RightTicks[i]->setLine( plotXSize, plotYSize * i/cvb.yticks, plotXSize - 5.0, plotYSize * i/cvb.yticks );
		RightTicks[i]->setPen( QPen( QBrush(Qt::SolidPattern), 1.0 ) );
		VText.push_back( new QGraphicsTextItem );
		VText[i]->setPlainText( QString::number( cvb.ymin + (cvb.ymax-cvb.ymin)*i/cvb.yticks ) );
		VText[i]->setFont( QFont("Arial", FontSize) );
		QRectF b = VText[i]->boundingRect().normalized();
		VText[i]->setPos( -b.width(), plotYSize - plotYSize*i/cvb.yticks - b.height()/2.0 );
		this->addItem( LeftTicks[i] );
		this->addItem( RightTicks[i] );
		this->addItem( VText[i] );
	}
}

void PlotData::PlotPaint( )
{
	// drawing the box
	makeBox( );
	// drawing the lines
	for( int i = 0; i < PathItems.size(); ++i )
	{
		this->removeItem( PathItems[i] );
		delete PathItems[i];
	}
	PathItems.clear();
	for( int i = 0; i < CircleItems.size(); ++i )
	{
		this->removeItem( CircleItems[i] );
		delete CircleItems[i];
	}
	CircleItems.clear();
	for( int i = 0; i < PolygonItems.size(); ++i )
	{
		this->removeItem( PolygonItems[i] );
		delete PolygonItems[i];
	}
	PolygonItems.clear();
	for( QList<PlotItem>::const_iterator i = Graph.constBegin(); i != Graph.constEnd(); i++ )
	{
		if( (*i).type == PlotLineType )
		{
			QGraphicsPathItem *pt = this->addPath( (*i).data.line->line, (*i).data.line->pen );
			PathItems.push_back( pt );
		}
		if( (*i).type == PlotCircleType )
		{
// 			std::cout<<"circles "<<(*i).data.circle->pos.size()<<"\n";
			for( int j = 0; j < (*i).data.circle->pos.size(); ++j )
			{
				QGraphicsEllipseItem *pt = this->addEllipse( (*i).data.circle->point, (*i).data.circle->pen );
				pt->setPos( (*i).data.circle->pos[j] );
				CircleItems.push_back( pt );
			}
		}
		if( (*i).type == PlotPolygonType )
		{
// 			std::cout<<"polygons "<<(*i).data.polygon->pos.size()<<"\n";
			for( int j = 0; j < (*i).data.polygon->pos.size(); ++j )
			{
				QGraphicsPolygonItem *pt = this->addPolygon( (*i).data.polygon->point, (*i).data.polygon->pen );
				pt->setPos( (*i).data.polygon->pos[j] );
				PolygonItems.push_back( pt );
			}
		}
	}
	this->clearFocus();
}
