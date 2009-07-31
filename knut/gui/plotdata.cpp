// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "plotdata.h"
#include "ncolloc.h"

#include <cmath>
#include <cfloat>

#include <QPen>
#include <QPolygon>
#include <QPainterPath>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>

PlotData::PlotData(QObject *parent) :
    QGraphicsScene(parent),
    AspectRatio(4.0 / 3.0),
    plotXSize(540), plotYSize(plotXSize / AspectRatio),
    FontSize(12), Box(0)
{
  const ViewBox cvb =
    {
      1.0, -1.0, 1.0, -1.0, 1, 1
    };
  ZoomHistory.push_back(cvb);
  currZoom = ZoomHistory.begin();
  addItem(&selection);
  selection.setVisible(false);
  makeBox();
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XNone,"None"));
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XLabel,"Label"));
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XMesh,"t/Period"));
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XRealMultiplier,"Re(Floquet Multiplier)"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YNone,"None"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YL2Norm,"L2Norm"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YAmplitude,"max(x(.))-min(x(.))"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YImagMultiplier,"Im(Floquet Multiplier)"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YAbsMultiplier,"Abs(Floquet Multiplier)"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YProfile,"X_%1"));
}

PlotData::~PlotData()
{
  clearAll();
}

int PlotData::nplots()
{
  unsigned int number = 0, it = 0;
  for (std::list<PlotItem>::iterator i = Graph.begin(); i != Graph.end(); ++i)
  {
    if (i->number != number)
    {
      number = i->number;
      ++it;
    }
  }
  return it;
}

void PlotData::clear(unsigned int n)
{
  unsigned int number = 0, it = 0;
  std::list<PlotItem>::iterator i = Graph.begin();
  while (i != Graph.end())
  {
    if (i->number != number)
    {
      number = i->number;
      ++it;
    }
    if (it == n)
    {
      if ((*i).type == PlotLineType)
      {
        delete(*i).data.line->item;
        delete(*i).data.line;
      }
      else if ((*i).type == PlotCircleType)
      {
        for (unsigned int j = 0; j < (*i).data.circle->item.size(); ++j) delete(*i).data.circle->item[j];
        delete(*i).data.circle;
      }
      else if ((*i).type == PlotPolygonType)
      {
        for (unsigned int j = 0; j < (*i).data.polygon->item.size(); ++j) delete(*i).data.polygon->item[j];
        delete(*i).data.polygon;
      }
      else
      {
        std::cout << "Something wrong\n";
      }
      i = Graph.erase(i);
    }
    else
    {
      ++i;
    }
  }
  if (n-1 < XCoordTextItems.size()) { delete XCoordTextItems[n-1]; XCoordTextItems.erase(XCoordTextItems.begin()+n-1); }
  else std::cout<<"GB1\n";
  if (n-1 < YCoordTextItems.size()) { delete YCoordTextItems[n-1]; YCoordTextItems.erase(YCoordTextItems.begin()+n-1); }
  else std::cout<<"GB2\n";
  if (n-1 < XCoordText.size()) XCoordText.erase(XCoordText.begin()+n-1);
  else std::cout<<"GB3\n";
  if (n-1 < YCoordText.size()) YCoordText.erase(YCoordText.begin()+n-1);
  else std::cout<<"GB4\n";
  makeBox();
  labelColor();
  //   update();
}

QColor PlotData::getColor(unsigned int n)
{
  unsigned int number = 0, it = 0;
  for (std::list<PlotItem>::const_iterator i = Graph.begin(); i != Graph.end(); ++i)
  {
    if (i->number != number)
    {
      number = i->number;
      ++it;
    }
    if (it == n)
    {
      if ((*i).type == PlotLineType)
      {
        return (*i).data.line->item->pen().color();
      }
      else if ((*i).type == PlotCircleType)
      {
        return (*i).data.circle->item[0]->pen().color();
      }
      else if ((*i).type == PlotPolygonType)
      {
        return (*i).data.polygon->item[0]->pen().color();
      }
      else
      {
        std::cout << "Something wrong\n";
      }
    }
  }
  return QColor(Qt::black);
}

void PlotData::setColor(unsigned int n, QColor& color)
{
  unsigned int number = 0, it = 0;
  for (std::list<PlotItem>::iterator i = Graph.begin(); i != Graph.end(); ++i)
  {
    if (i->number != number)
    {
      number = i->number;
      ++it;
    }
    if (it == n)
    {
      if ((*i).type == PlotLineType)
      {
        QPen pen = (*i).data.line->item->pen();
        pen.setColor(color);
        (*i).data.line->pen = pen;
        (*i).data.line->item->setPen(pen);
      }
      else if ((*i).type == PlotCircleType)
      {
        for (unsigned int j = 0; j < (*i).data.circle->item.size(); ++j)
        {
          QPen pen = (*i).data.circle->item[j]->pen();
          pen.setColor(color);
          (*i).data.circle->pen = pen;
          (*i).data.circle->item[j]->setPen(pen);
        }
      }
      else if ((*i).type == PlotPolygonType)
      {
        for (unsigned int j = 0; j < (*i).data.polygon->item.size(); ++j)
        {
          QPen pen = (*i).data.polygon->item[j]->pen();
          pen.setColor(color);
          (*i).data.polygon->pen = pen;
          (*i).data.polygon->item[j]->setPen(pen);
        }
      }
      else
      {
        std::cout << "Something wrong\n";
      }
    }
  }
  labelColor();
//   update();
}

void PlotData::clearAxes()
{
  for (unsigned int i = 0; i < XCoordTextItems.size(); ++i)
  {
    removeItem(XCoordTextItems[i]);
    delete XCoordTextItems[i];
  }
  XCoordTextItems.clear();
  for (unsigned int i = 0; i < YCoordTextItems.size(); ++i)
  {
    removeItem(YCoordTextItems[i]);
    delete YCoordTextItems[i];
  }
  YCoordTextItems.clear();
}

void PlotData::clearAll()
{
  std::list<PlotItem>::iterator i;
  for (i = Graph.begin(); i != Graph.end(); i++)
  {
    if ((*i).type == PlotLineType)
    {
      delete(*i).data.line->item;
      delete(*i).data.line;
    }
    else if ((*i).type == PlotCircleType)
    {
      for (unsigned int j = 0; j < (*i).data.circle->item.size(); ++j) delete(*i).data.circle->item[j];
      delete(*i).data.circle;
    }
    else if ((*i).type == PlotPolygonType)
    {
      for (unsigned int j = 0; j < (*i).data.polygon->item.size(); ++j) delete(*i).data.polygon->item[j];
      delete(*i).data.polygon;
    }
    else
    {
      std::cout << "Something wrong\n";
    }
  }
  Graph.clear();
  clearAxes();
  XCoordText.clear();
  YCoordText.clear();
  const ViewBox cvb = *currZoom;
  ZoomHistory.clear();
  ZoomHistory.push_back(cvb);
  currZoom = ZoomHistory.begin();
//   update();
}

/// this computes the minimum and maximum value of an axis
static inline void adjustAxis(qreal& min, qreal& max, unsigned int& numTicks)
{
  const unsigned int MinTicks = 4;
  qreal grossStep = (max - min) / MinTicks;
  qreal step = pow(10, floor(log10(grossStep)));

  if (5 * step < grossStep) step *= 5;
  else if (2 * step < grossStep) step *= 2;

  numTicks = (unsigned int)(ceil(max / step) - floor(min / step));
  min = floor(min / step) * step;
  max = ceil(max / step) * step;
}

void PlotData::addPlotLine(std::list<PlotItem>::iterator& it, const QPen& pen, bool principal)
{
  it->type = PlotLineType;
  it->data.line = new PlotLine(pen);
  it->data.line->item = 0;
  it->data.line->pen.setWidthF(1);
  it->principal = principal;
}

void PlotData::addPlotPoint(std::list<PlotItem>::iterator& it, const QPen& pen, PlotMarkerStyle type, bool principal)
{
  switch (type)
  {
    case PlotMarkerCircle:  // CIRCLE
      {
        it->type = PlotCircleType;
        it->data.circle = new PlotCircle(pen, QRectF(-3.0, -3.0, 6.0 , 6.0));
        it->data.circle->pen.setWidthF(1);
      }
      break;
    case PlotMarkerSquare:  // SQUARE
      {
        it->type = PlotPolygonType;
        QPolygonF pl(4);
        pl[0] = QPointF(-3.0, -3.0);
        pl[1] = QPointF(-3.0, 3.0);
        pl[2] = QPointF(3.0, 3.0);
        pl[3] = QPointF(3.0, -3.0);
        it->data.polygon = new PlotPolygon(pen, pl);
        it->data.polygon->pen.setWidthF(1);
      }
      break;
    case PlotMarkerTriangle: // TRIANGLE
      {
        it->type = PlotPolygonType;
        QPolygonF pl(3);
        pl[0] = QPointF(-3.0, -3.0);
        pl[1] = QPointF(0.0, 3.0);
        pl[2] = QPointF(3.0, -3.0);
        it->data.polygon = new PlotPolygon(pen, pl);
        it->data.polygon->pen.setWidthF(1);
      }
      break;
    case PlotMarkerCross: // CROSS
      {
        it->type = PlotPolygonType;
        QPolygonF pl(6);
        pl[0] = QPointF(-3.0, 3.0);
        pl[1] = QPointF(3.0, -3.0);
        pl[2] = QPointF(0.0, 0.0);
        pl[3] = QPointF(-3.0, -3.0);
        pl[4] = QPointF(3.0, 3.0);
        pl[5] = QPointF(0.0, 0.0);
        it->data.polygon = new PlotPolygon(pen, pl);
        it->data.polygon->pen.setWidthF(1);
      }
      break;
    default:
      std::cout << "Another Serious BUG\n";
      abort();
      break;
  }
  it->principal = principal;
}

// this only called by addPlot...
void PlotData::dataToGraphics()
{
  if (ZoomHistory.begin() == currZoom)
  {
    ViewBox cvb = { -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, 1, 1 };
    std::list<PlotItem>::const_iterator it;
    for (it = Graph.begin(); it != Graph.end(); it++)
    {
      for (int k = 0; k < it->x.Size(); k++)
      {
        if (it->x(k) > cvb.xmax) cvb.xmax = it->x(k);
        if (it->x(k) < cvb.xmin) cvb.xmin = it->x(k);
        if (it->y(k) > cvb.ymax) cvb.ymax = it->y(k);
        if (it->y(k) < cvb.ymin) cvb.ymin = it->y(k);
      }
    }
    adjustAxis(cvb.xmin, cvb.xmax, cvb.xticks);
    adjustAxis(cvb.ymin, cvb.ymax, cvb.yticks);
    *currZoom = cvb;
  }
  rescaleData();
  PlotPaint();
  labelColor();
}

bool PlotData::addPlot(const mat4Data* data, PlotXVariable x, PlotYVariable y, 
  int pt, int dim)
{
  int xadded = 0;
  int yadded = 0;
  std::list<PlotItem>::iterator start = --Graph.end();
  int startnumber = 0;
  if (Graph.begin() != Graph.end()) startnumber = Graph.rbegin()->number;
  QColor plcolor, stabcolor(Qt::black);
  plcolor.setHsv((240 + startnumber*40)%360, 255, 255);
  // mindig az Y utan kovetkezik csak addPlot...()
  if (x >= XParameter0 && y != YAbsMultiplier && y != YProfile)
  {
    Graph.push_back(PlotItem());
    Graph.rbegin()->x.Init(data->getNPoints());
    Graph.rbegin()->y.Init(data->getNPoints());
    for (int i = 0; i < data->getNPoints(); i++)
    {
      if (x != XParameter0) 
        Graph.rbegin()->x(i) = data->getPar(i, x - XParameter0);
      else
        Graph.rbegin()->x(i) = data->getPar(i, x - XParameter0)/data->getPar(i,data->getNPar()-ParEnd+ParPeriod);
    }
    ++xadded;
  }
  if (x == XLabel && y != YAbsMultiplier && y != YProfile)
  {
    Graph.push_back(PlotItem());
    Graph.rbegin()->x.Init(data->getNPoints());
    Graph.rbegin()->y.Init(data->getNPoints());
    for (int i = 0; i < data->getNPoints(); i++)
    {
      Graph.rbegin()->x(i) = i;
    }
    ++xadded;
  }
  if (x == XMesh && y == YProfile)
  {
    const int ndeg = data->getNDeg();
    const int nint = data->getNInt();
    Graph.push_back(PlotItem());
    Graph.rbegin()->x.Init(ndeg*nint + 1);
    Graph.rbegin()->y.Init(ndeg*nint + 1);
    for (int i = 0; i < nint; i++)
    {
      for (int j = 0; j < ndeg; j++)
      {
        Graph.rbegin()->x(j + ndeg*i) = data->getMesh(pt, i) + data->getElem(pt, j) * (data->getMesh(pt, i + 1) - data->getMesh(pt, i));
        Graph.rbegin()->y(j + ndeg*i) = data->getProfile(pt, dim, j + ndeg * i);
      }
    }
    Graph.rbegin()->x(ndeg*nint) = data->getMesh(pt, nint);
    Graph.rbegin()->y(ndeg*nint) = data->getProfile(pt, dim, ndeg * nint);
    ++xadded;
    ++yadded;
    addPlotLine(--Graph.end(), QPen(plcolor), true);
  }
  if (y >= YParameter0 && (x != XMesh && x != XRealMultiplier))
  {
    for (int i = 0; i < data->getNPoints(); i++)
    {
      Graph.rbegin()->y(i) = data->getPar(i, y - YParameter0);
    }
    ++yadded;
    addPlotLine(--Graph.end(), QPen(plcolor), true);
  }
  if (y == YAmplitude && (x != XMesh && x != XRealMultiplier))
  {
    for (int i = 0; i < data->getNPoints(); i++)
    {
      const int ndeg = data->getNDeg();
      const int nint = data->getNInt();
      double nrm = 0.0;
      for (int p = 0; p < data->getNDim(); ++p)
      {
        double min = DBL_MAX;
        double max = -DBL_MAX;
        for (int j = 0; j < nint; ++j)
        {
          for (int k = 0; k < ndeg; ++k)
          {
            if (min > data->getProfile(i, p, k + ndeg*j)) min = data->getProfile(i, p, k + ndeg * j);
            if (max < data->getProfile(i, p, k + ndeg*j)) max = data->getProfile(i, p, k + ndeg * j);
          }
        }
        nrm += (max - min) * (max - min);
      }
      Graph.rbegin()->y(i) = sqrt(nrm);
    }
    ++yadded;
    addPlotLine(--Graph.end(), QPen(plcolor), true);
  }
  if (y == YL2Norm && (x != XMesh && x != XRealMultiplier))
  {
    Vector elem(false);
    const_cast<mat4Data*>(data)->getElemRef(0, elem);
    Matrix metric(elem.Size(), elem.Size());
    NColloc::getMetric(metric, elem);
    Vector prof(false), msh(false);
    for (int i = 0; i < data->getNPoints(); i++)
    {
      const_cast<mat4Data*>(data)->getMeshRef(i, msh);
      const_cast<mat4Data*>(data)->getProfileRef(i, prof);
      Graph.rbegin()->y(i) = NColloc::integrate(prof, prof, metric, msh, data->getNDim());
    }
    ++yadded;
    addPlotLine(--Graph.end(), QPen(plcolor), true);
  }
  if (x == XRealMultiplier && y == YImagMultiplier)
  {
    // plot the unit circel!!!
    const int segments = 100;
    Graph.push_back(PlotItem());
    Graph.rbegin()->x.Init(segments);
    Graph.rbegin()->y.Init(segments);
    for (int i = 0; i < segments; ++i)
    {
      Graph.rbegin()->x(i) = sin(i * 2.0 * M_PI / (segments - 1));
      Graph.rbegin()->y(i) = cos(i * 2.0 * M_PI / (segments - 1));
    }
    addPlotLine(--Graph.end(), QPen(stabcolor), false);
    // plotting the multipliers
    Graph.push_back(PlotItem());
    Graph.rbegin()->x.Init(data->getNMul());
    Graph.rbegin()->y.Init(data->getNMul());
    for (int i = 0; i < data->getNMul(); i++)
    {
      Graph.rbegin()->x(i) = data->getMulRe(pt, i);
      Graph.rbegin()->y(i) = data->getMulIm(pt, i);
    }
    ++xadded;
    ++yadded;
    addPlotPoint(--Graph.end(), QPen(plcolor), PlotMarkerCross, true);
  }
  if ((x == XLabel || x >= XParameter0) && y == YAbsMultiplier)
  {
    for (int r = 0; r < data->getNMul(); r++)
    {
      Graph.push_back(PlotItem());
      Graph.rbegin()->x.Init(data->getNPoints());
      Graph.rbegin()->y.Init(data->getNPoints());
      for (int i = 0; i < data->getNPoints(); i++)
      {
        if (x >= XParameter0) Graph.rbegin()->x(i) = data->getPar(i, x - XParameter0);
        else Graph.rbegin()->x(i) = i;
        Graph.rbegin()->y(i) = sqrt(data->getMulRe(i, r) * data->getMulRe(i, r) + data->getMulIm(i, r) * data->getMulIm(i, r));
      }
      ++xadded;
      ++yadded;
//       addPlotPoint(--Graph.end(), QPen(plcolor), 0);
      addPlotLine(--Graph.end(), QPen(plcolor), true);
    }
  }
  // add stability
  std::vector<int>     bifidx;
  std::vector<BifType> biftype;
  if ((x == XLabel || x >= XParameter0) && xadded != 0 && xadded == yadded)
  {
    int k, k_p = 0;
    do
    {
      k = data->getNextBifurcation(k_p);
      if (k != -1)
      {
        bifidx.push_back(k);
        biftype.push_back(data->getBifurcationType(k));
        k_p = k;
      }
    }
    while (k != -1);
    std::list<PlotItem>::const_iterator it = --Graph.end();
    if (it->x.Size() == data->getNPoints())
    {
      for (unsigned int i = 0; i < bifidx.size(); i++)
      {
        Graph.push_back(PlotItem());
        Graph.rbegin()->x.Init(1);
        Graph.rbegin()->y.Init(1);
        Graph.rbegin()->x(0) = (it->x(bifidx[i] - 1) + it->x(bifidx[i])) / 2.0;
        if (y == YAbsMultiplier) Graph.rbegin()->y(0) = 1.0;
        else Graph.rbegin()->y(0) = (it->y(bifidx[i] - 1) + it->y(bifidx[i])) / 2.0;
        ++xadded;
        ++yadded;
        switch (biftype[i])
        {
          case BifLP:
            addPlotPoint(--Graph.end(), QPen(stabcolor), PlotMarkerSquare, false);
            break;
          case BifPD:
            addPlotPoint(--Graph.end(), QPen(stabcolor), PlotMarkerTriangle, false);
            break;
          case BifNS:
            addPlotPoint(--Graph.end(), QPen(stabcolor), PlotMarkerCircle, false);
            break;
          default:
            addPlotPoint(--Graph.end(), QPen(stabcolor), PlotMarkerCross, false);
            break;
        } 
      }
    }
  }
  for (std::list<PlotItem>::iterator it = ++start; it != Graph.end(); ++it)
  {
    it->number = startnumber + 1;
  }
  if (xadded != yadded)
  {
    std::cout << "bad number of X and Y coordinates\n";
    return false;
  }
  else
  {
    if (xadded != 0)
    {
      std::vector<std::string> parNames;
      data->getParNames(parNames);
      if (x >= XParameter0) XCoordText.push_back(parNames.at(x-XParameter0).c_str());
      else XCoordText.push_back(XCoordMap[x]);
      if (y >= YParameter0) YCoordText.push_back(parNames.at(y-YParameter0).c_str());
      else if (y == YProfile) YCoordText.push_back(QString("X_%1(t)").arg(dim));
      else YCoordText.push_back(YCoordMap[y]);
      
      dataToGraphics();
      labelColor();
      return true;
    }
    else
    {
      return false;
    }
  }
}

///
/// The inlines
///

static inline qreal itr(qreal xx, qreal x1, qreal x2, qreal y1, qreal y2)
{
  return ((x2 -xx)*y1 + (xx - x1)*y2) / (x2 - x1);
}

static inline int intpos(qreal x, qreal min, qreal max)
{
  if (x < min) return -1;
  else if (x < max) return 0;
  else return 1;
}

static inline bool crossing(int x1, int y1, int x2, int y2)
{
  if ((x1*x2 == 1) || (y1*y2 == 1)) return false;
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
QPointF PlotData::intersect(QPointF p1, QPointF p2)
{
  // p1 is in the box, p2 is out of the box
  const ViewBox cvb =
    {
      plotXSize, 0.0, plotYSize, 0.0, 0, 0
    };
  qreal ycr, xcr;
  if (p2.x() > cvb.xmax)
  {
    ycr = itr(cvb.xmax, p1.x(), p2.x(), p1.y(), p2.y());
    if (ycr > cvb.ymax)
    {
      ycr = cvb.ymax;
      xcr = itr(cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x());
    }
    else if (ycr > cvb.ymin) xcr = cvb.xmax;
    else
    {
      ycr = cvb.ymin;
      xcr = itr(cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x());
    }
  }
  else if (p2.x() > cvb.xmin)
  {
    if (p2.y() > cvb.ymax)
    {
      xcr = itr(cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x());
      ycr = cvb.ymax;
    }
    else if (p2.y() > cvb.ymin)
    {
      xcr = p2.x(), ycr = p2.y();
    }
    else
    {
      xcr = itr(cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x());
      ycr = cvb.ymin;
    }
  }
  else
  {
    ycr = itr(cvb.xmin, p1.x(), p2.x(), p1.y(), p2.y());
    if (ycr > cvb.ymax)
    {
      ycr = cvb.ymax;
      xcr = itr(cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x());
    }
    else if (ycr > cvb.ymin) xcr = cvb.xmin;
    else
    {
      ycr = cvb.ymin;
      xcr = itr(cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x());
    }
  }
  return QPointF(xcr, ycr);
}

inline bool PlotData::contains(double x, double y)
{
  return (x >= 0.0) && (y >= 0.0) && (x <= plotXSize) && (y <= plotYSize);
}

bool PlotData::crossbox(QPointF p1, QPointF p2, QPointF& i1, QPointF& i2)
{
  const ViewBox cvb =
    {
      plotXSize, 0.0, plotYSize, 0.0, 0, 0
    };

  qreal ycrmin = itr(cvb.xmin, p1.x(), p2.x(), p1.y(), p2.y());
  qreal ycrmax = itr(cvb.xmax, p1.x(), p2.x(), p1.y(), p2.y());
  qreal xcrmin = itr(cvb.ymin, p1.y(), p2.y(), p1.x(), p2.x());
  qreal xcrmax = itr(cvb.ymax, p1.y(), p2.y(), p1.x(), p2.x());
  int id = 0;
  qreal x[4], y[4]; // security reasons
  if ((ycrmin > cvb.ymin) && (ycrmin < cvb.ymax))
  {
    x[id] = cvb.xmin;
    y[id] = ycrmin;
    ++id;
  }
  if ((ycrmax > cvb.ymin) && (ycrmax < cvb.ymax))
  {
    x[id] = cvb.xmax;
    y[id] = ycrmax;
    ++id;
  }
  if ((xcrmin > cvb.xmin) && (xcrmin < cvb.xmax))
  {
    y[id] = cvb.ymin;
    x[id] = xcrmin;
    ++id;
  }
  if ((xcrmax > cvb.xmin) && (xcrmax < cvb.xmax))
  {
    y[id] = cvb.ymax;
    x[id] = xcrmax;
    ++id;
  }
  if (id == 0) return false;
  else if (id == 2)
  {
    i1 = QPointF(x[0], y[0]);
    i2 = QPointF(x[1], y[1]);
    return true;
  }
  else
  {
    std::cout << "GEBASZ van " << id << "\n";
    return false;
  }
}

void PlotData::rescaleData()
{
  ViewBox cvb = *currZoom;

  const double xscale = plotXSize / (cvb.xmax - cvb.xmin);
  const double yscale = plotYSize / (cvb.ymax - cvb.ymin);

  // rescaling all the data
  std::list<PlotItem>::iterator i;
  for (i = Graph.begin(); i != Graph.end(); i++)
  {
    if ((*i).type == PlotLineType)
    {
      delete(*i).data.line->item;
      (*i).data.line->item = 0;
      (*i).data.line->line = QPainterPath();
      int x = 0, y = 0, prx = 0, pry = 0;
      bool pr = true;
      if (i->x.Size() != i->y.Size())
      {
        std::cout << "DataX DataY Sizes differ\n";
        return;
      }
      for (int k = 0; k < i->x.Size(); k++)
      {
        x = intpos(i->x(k), cvb.xmin, cvb.xmax);
        y = intpos(i->y(k), cvb.ymin, cvb.ymax);
        if (crossing(prx, pry, x, y))
        {
          if ((x == 0) && (y == 0))
          {
            if (!pr)
            {
              // the previous was not in the viewport
              (*i).data.line->line.moveTo(intersect(QPointF(xscale*(i->x(k) - cvb.xmin), yscale*(cvb.ymax - i->y(k))),
                                                    QPointF(xscale*(i->x(k - 1) - cvb.xmin), yscale*(cvb.ymax - i->y(k - 1)))));
            }
            if (k == 0)
            {
              (*i).data.line->line.moveTo(xscale*(i->x(k) - cvb.xmin), yscale*(cvb.ymax - i->y(k)));
            }
            else
            {
              (*i).data.line->line.lineTo(xscale*(i->x(k) - cvb.xmin), yscale*(cvb.ymax - i->y(k)));
            }
            pr = true;
          }
          else
          {
            if (k != 0)
            {
              if (pr)
              {
                (*i).data.line->line.lineTo(intersect(QPointF(xscale*(i->x(k - 1) - cvb.xmin), yscale*(cvb.ymax - i->y(k - 1))),
                                                      QPointF(xscale*(i->x(k) - cvb.xmin), yscale*(cvb.ymax - i->y(k)))));
              }
              else
              {
                QPointF pt1, pt2;
                if (crossbox(QPointF(xscale*(i->x(k - 1) - cvb.xmin), yscale*(cvb.ymax - i->y(k - 1))),
                             QPointF(xscale*(i->x(k) - cvb.xmin), yscale*(cvb.ymax - i->y(k))),
                             pt1, pt2))
                {
                  (*i).data.line->line.moveTo(pt1);
                  (*i).data.line->line.lineTo(pt2);
                }
              }
            }
            pr = false;
          }
        }
        prx = x;
        pry = y;
      }
    }
    if ((*i).type == PlotCircleType)
    {
      (*i).data.circle->pos.clear();
      for (int k = 0; k < i->x.Size(); k++)
      {
        const QPointF pt = QPointF(xscale * (i->x(k) - cvb.xmin), yscale * (cvb.ymax - i->y(k)));
        if (contains(pt.x(), pt.y()))(*i).data.circle->pos.push_back(pt);
      }
      for (unsigned int p = (*i).data.circle->pos.size(); p < (*i).data.circle->item.size(); ++p) delete(*i).data.circle->item[p];
      (*i).data.circle->item.resize((*i).data.circle->pos.size(), 0);
    }
    if ((*i).type == PlotPolygonType)
    {
      (*i).data.polygon->pos.clear();
      for (int k = 0; k < i->x.Size(); k++)
      {
        const QPointF pt = QPointF(xscale * (i->x(k) - cvb.xmin), yscale * (cvb.ymax - i->y(k)));
        if (contains(pt.x(), pt.y()))(*i).data.polygon->pos.push_back(pt);
      }
      for (unsigned int p = (*i).data.polygon->pos.size(); p < (*i).data.polygon->item.size(); ++p) delete(*i).data.polygon->item[p];
      (*i).data.polygon->item.resize((*i).data.polygon->pos.size(), 0);
    }
  }
}

void PlotData::mousePressEvent(QGraphicsSceneMouseEvent * event)
{
  if (event->button() == Qt::LeftButton)
  {
    mouseBegin = QPointF(event->scenePos());
    mouseMove = mouseBegin;
    event->accept();
  }
}

void PlotData::mouseReleaseEvent(QGraphicsSceneMouseEvent * event)
{
  if (event->button() == Qt::LeftButton)
  {
    mouseEnd = QPointF(event->scenePos());
    if (mouseEnd == mouseBegin)
    {
      event->accept();
      return;
    }

    ViewBox cvb = *currZoom, newvb;
    const double xscale = plotXSize / (cvb.xmax - cvb.xmin);
    const double yscale = plotYSize / (cvb.ymax - cvb.ymin);

    newvb.xmin = qMin(mouseBegin.x(), mouseEnd.x()) / xscale + cvb.xmin;
    newvb.xmax = qMax(mouseBegin.x(), mouseEnd.x()) / xscale + cvb.xmin;
    newvb.ymin = cvb.ymax - qMax(mouseBegin.y(), mouseEnd.y()) / yscale;
    newvb.ymax = cvb.ymax - qMin(mouseBegin.y(), mouseEnd.y()) / yscale;
    adjustAxis(newvb.xmin, newvb.xmax, newvb.xticks);
    adjustAxis(newvb.ymin, newvb.ymax, newvb.yticks);
    // remove the next zoom levels
    ZoomHistory.erase(++currZoom, ZoomHistory.end());
    --currZoom;
    // add the new level
    ZoomHistory.push_back(newvb);
    currZoom = --ZoomHistory.end();
    dataToGraphics();
    selection.setVisible(false);
    update(selection.boundingRect().normalized());
    event->accept();
  }
}

void PlotData::mouseMoveEvent(QGraphicsSceneMouseEvent * event)
{
  if (event->buttons() == Qt::LeftButton)
  {
    mouseMove = QPointF(event->scenePos());
    if (sceneRect().contains(mouseMove))
    {
      selection.setRect(QRectF(mouseBegin.x(), mouseBegin.y(), mouseMove.x() - mouseBegin.x() , mouseMove.y() - mouseBegin.y()).normalized());
      selection.setVisible(true);
//       selection.update();
    }
    event->accept();
  }
}

void PlotData::keyPressEvent(QKeyEvent * key)
{
  if (key->key() == Qt::Key_P)
  {
    if (ZoomHistory.begin() != currZoom)
    {
      --currZoom;
      dataToGraphics();
    }
    key->accept();
  }
  if (key->key() == Qt::Key_N)
  {
    if (--ZoomHistory.end() != currZoom)
    {
      ++currZoom;
      dataToGraphics();
    }
    key->accept();
  }
}

void PlotData::makeBox()
{
  ViewBox cvb = *currZoom;
  // drawing the box
  if (Box == 0) Box = addRect(QRectF(0.0, 0.0, plotXSize, plotYSize), QPen(QBrush(Qt::SolidPattern), 2.0));

  // drawing the ticks and tickmarks
  for (unsigned int i = 0; i < HText.size(); i++)
  {
    removeItem(BottomTicks[i]);
    removeItem(TopTicks[i]);
    removeItem(HText[i]);
    delete BottomTicks[i];
    delete TopTicks[i];
    delete HText[i];
  }
  BottomTicks.clear();
  TopTicks.clear();
  HText.clear();
  qreal highest = 0;
  for (unsigned int i = 0; i < cvb.xticks + 1; i++)
  {
    BottomTicks.push_back(new QGraphicsLineItem);
    BottomTicks[i]->setLine(plotXSize * i / cvb.xticks, 0.0, plotXSize * i / cvb.xticks, 5.0);
    BottomTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    TopTicks.push_back(new QGraphicsLineItem);
    TopTicks[i]->setLine(plotXSize * i / cvb.xticks, plotYSize, plotXSize * i / cvb.xticks, plotYSize - 5.0);
    TopTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    HText.push_back(new QGraphicsTextItem);
    HText[i]->setPlainText(QString::number(cvb.xmin + (cvb.xmax - cvb.xmin)*i / cvb.xticks));
    HText[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = HText[i]->boundingRect().normalized();
    HText[i]->setPos(plotXSize * i / cvb.xticks - b.width() / 2.0, plotYSize /*- b.height()*/);
    addItem(TopTicks[i]);
    addItem(BottomTicks[i]);
    addItem(HText[i]);
    if (highest < b.height()) highest = b.height();
  }
  for (unsigned int i = 0; i < VText.size(); i++)
  {
    removeItem(LeftTicks[i]);
    removeItem(RightTicks[i]);
    removeItem(VText[i]);
    delete LeftTicks[i];
    delete RightTicks[i];
    delete VText[i];
  }
  LeftTicks.clear();
  RightTicks.clear();
  VText.clear();
  qreal widest = 0;
  for (unsigned int i = 0; i < cvb.yticks + 1; i++)
  {
    LeftTicks.push_back(new QGraphicsLineItem);
    LeftTicks[i]->setLine(0.0, plotYSize * i / cvb.yticks, 5.0, plotYSize * i / cvb.yticks);
    LeftTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    RightTicks.push_back(new QGraphicsLineItem);
    RightTicks[i]->setLine(plotXSize, plotYSize * i / cvb.yticks, plotXSize - 5.0, plotYSize * i / cvb.yticks);
    RightTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    VText.push_back(new QGraphicsTextItem);
    VText[i]->setPlainText(QString::number(cvb.ymin + (cvb.ymax - cvb.ymin)*i / cvb.yticks));
    VText[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = VText[i]->boundingRect().normalized();
    VText[i]->setPos(-b.width(), plotYSize - plotYSize*i / cvb.yticks - b.height() / 2.0);
    addItem(LeftTicks[i]);
    addItem(RightTicks[i]);
    addItem(VText[i]);
    if (widest < b.width()) widest = b.width();
  }
  // remove previous texts
  clearAxes();
  // Add labels
  qreal sumwidth = 0;
  qreal sumheight = 0;
  unsigned int startmove = 0;
  for (unsigned int i = 0; i < YCoordText.size(); ++i)
  {
    YCoordTextItems.push_back(new QGraphicsTextItem);
    YCoordTextItems[i]->setPlainText(YCoordText[i]);
    YCoordTextItems[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = YCoordTextItems[i]->boundingRect().normalized();
    YCoordTextItems[i]->setPos(-widest - b.height() - sumheight, plotYSize/2.0 + sumwidth + b.width());
    YCoordTextItems[i]->setTransform(QTransform().rotate(-90));
    sumwidth += b.width();
    if (sumwidth > plotYSize)
    {
      for (unsigned int i = startmove; i < YCoordTextItems.size(); ++i)
      {
        YCoordTextItems[i]->moveBy(0.0, -sumwidth/2.0);
        addItem(YCoordTextItems[i]);
      }
      sumheight += highest;
      startmove = YCoordTextItems.size();
      sumwidth = 0;
    }
  }
  for (unsigned int i = startmove; i < YCoordTextItems.size(); ++i)
  {
    YCoordTextItems[i]->moveBy(0.0, -sumwidth/2.0);
    addItem(YCoordTextItems[i]);
  }
  sumwidth = 0;
  sumheight = 0;
  startmove = 0;
  for (unsigned int i = 0; i < XCoordText.size(); ++i)
  {
    XCoordTextItems.push_back(new QGraphicsTextItem);
    XCoordTextItems[i]->setPlainText(XCoordText[i]);
    XCoordTextItems[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = XCoordTextItems[i]->boundingRect().normalized();
    XCoordTextItems[i]->setPos(plotXSize / 2.0 + sumwidth, plotYSize + highest + sumheight);
    sumwidth += b.width();
    if (sumwidth > plotXSize)
    {
      for (unsigned int i = startmove; i < XCoordTextItems.size(); ++i)
      {
        XCoordTextItems[i]->moveBy(-sumwidth/2.0, 0.0);
        addItem(XCoordTextItems[i]);
      }
      sumheight += highest;
      startmove = XCoordTextItems.size();
      sumwidth = 0;
    }
  }
  for (unsigned int i = startmove; i < XCoordText.size(); ++i)
  {
    XCoordTextItems[i]->moveBy(-sumwidth/2.0, 0.0);
    addItem(XCoordTextItems[i]);
  }
}

void PlotData::labelColor()
{
  unsigned int prenumber = 0, it = 0;
  for (std::list<PlotItem>::iterator i = Graph.begin(); i != Graph.end(); i++)
  {
    if ((*i).principal && prenumber != (*i).number)
    {
      QColor textColor;
      if ((*i).type == PlotLineType)
      {      
        textColor = (*i).data.line->pen.color();
      }
      if ((*i).type == PlotCircleType)
      {
        textColor = (*i).data.circle->pen.color();
      }
      if ((*i).type == PlotPolygonType)
      {
        textColor = (*i).data.polygon->pen.color();
      }
      unsigned int num = it++;
      if (num < XCoordTextItems.size() && num < YCoordTextItems.size() && 
          num < XCoordText.size() && num < YCoordText.size()) // safety measures
      {
        XCoordTextItems[num]->setDefaultTextColor(textColor);
        YCoordTextItems[num]->setDefaultTextColor(textColor);
//        std::cout << "Color @" << num << " " << textColor.red() << " " << textColor.green() << " " << textColor.blue() << "\n";
      }
//      else std::cout<<"GB5\n";
    }
    prenumber = (*i).number;
  }
}

void PlotData::PlotPaint()
{
  // drawing the box
  makeBox();
  // drawing the lines
  for (std::list<PlotItem>::iterator i = Graph.begin(); i != Graph.end(); i++)
  {
    QColor textColor;
    if ((*i).type == PlotLineType)
    {
      delete(*i).data.line->item;
      (*i).data.line->item = addPath((*i).data.line->line, (*i).data.line->pen);
    }
    if ((*i).type == PlotCircleType)
    {
      for (unsigned int j = 0; j < (*i).data.circle->pos.size(); ++j)
      {
        delete(*i).data.circle->item[j];
        (*i).data.circle->item[j] = addEllipse((*i).data.circle->point, (*i).data.circle->pen);
        (*i).data.circle->item[j]->setPos((*i).data.circle->pos[j]);
      }
    }
    if ((*i).type == PlotPolygonType)
    {
      for (unsigned int j = 0; j < (*i).data.polygon->pos.size(); ++j)
      {
        delete(*i).data.polygon->item[j];
        (*i).data.polygon->item[j] = addPolygon((*i).data.polygon->point, (*i).data.polygon->pen);
        (*i).data.polygon->item[j]->setPos((*i).data.polygon->pos[j]);
      }
    }
  }
//   update();
}
