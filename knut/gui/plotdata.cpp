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
#include <memory>

#include <QPen>
#include <QPolygon>
#include <QPainterPath>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>

const ViewBox PlotData::defaultBox;

PlotData::PlotData(QObject *parent) :
    QGraphicsScene(parent),
    xlog(false), ylog(false),
    AspectRatio(4.0 / 3.0),
    plotXSize(480), plotYSize(plotXSize / AspectRatio),
    FontSize(12), unitCircleCount(0)
{
  const ViewBox cvb = defaultBox;
  setSceneRect(-0.1*plotXSize, -0.1*plotYSize, 1.2*plotXSize, 1.2*plotYSize);
  ZoomHistory.push_back(cvb);
  currZoom = ZoomHistory.begin();
  addItem(&selection);
//  std::cout << "addItem " << __LINE__ << "\n";
  selection.setVisible(false);
  clipBox.addRect(QRectF(0.0, 0.0, plotXSize, plotYSize));
  Box.reset(makeBox());
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XNone,"None"));
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XLabel,"Label"));
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XMesh,"t/Period"));
  XCoordMap.insert(std::pair<PlotXVariable,QString>(XRealMultiplier,"Re(Floquet Multiplier)"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YNone,"None"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YL2Norm,"L2Norm"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YAmplitude,"max(x(.))-min(x(.))"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YPoincare,"Map"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YImagMultiplier,"Im(Floquet Multiplier)"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YAbsMultiplier,"Abs(Floquet Multiplier)"));
  YCoordMap.insert(std::pair<PlotYVariable,QString>(YProfile,"X_%1"));
}

PlotData::~PlotData()
{
  clearAll();
}

void PlotData::setXSize(int size)
{
  plotXSize = size;
  plotYSize = plotXSize / AspectRatio;
  dataToGraphics(Graph.begin(), Graph.end());
//   std::cout << "GRAPH RESCALE " << size << "\n";
}

size_t PlotData::nplots()
{
  size_t number = 0, it = 0;
  for (const auto& ii : Graph)
  {
    if (ii->number != number)
    {
      number = ii->number;
      ++it;
    }
  }
  return it;
}

void PlotData::clear(size_t n)
{
  size_t number = 0, it = 0;
  auto i = Graph.begin();
  while (i != Graph.end())
  {
    if ((*i)->number != number)
    {
      number = (*i)->number;
      ++it;
    }
    if (it == n)
    {
      i = Graph.erase(i);
    }
    else
    {
      ++i;
    }
  }
  if (n-1 < XCoordTextItems.size()) { XCoordTextItems[n-1].reset(nullptr); XCoordTextItems.erase(XCoordTextItems.begin()+n-1); }
  else std::cout<<"GB1\n";
  if (n-1 < YCoordTextItems.size()) { YCoordTextItems[n-1].reset(nullptr); YCoordTextItems.erase(YCoordTextItems.begin()+n-1); }
  else std::cout<<"GB2\n";
  if (n-1 < XCoordText.size()) XCoordText.erase(XCoordText.begin()+n-1);
  else std::cout<<"GB3\n";
  if (n-1 < YCoordText.size()) YCoordText.erase(YCoordText.begin()+n-1);
  else std::cout<<"GB4\n";
//   Box = makeBox();
  labelColor();
}

QColor PlotData::getColor(size_t n)
{
  size_t number = 0, it = 0;
  for (const auto& ii : Graph)
  {
    if (ii->number != number)
    {
      number = ii->number;
      ++it;
    }
    if (it == n)
    {
      if (ii->type == PlotLineType)
      {
        return std::get<PlotLine>(ii->data).pen.color();
      }
      else if (ii->type == PlotCircleType)
      {
        return std::get<PlotCircle>(ii->data).item.front()->pen().color();
      }
      else if (ii->type == PlotPolygonType)
      {
        return std::get<PlotPolygon>(ii->data).item.front()->pen().color();
      }
      else
      {
        std::cout << "Something wrong\n";
      }
    }
  }
  return QColor(Qt::black);
}

void PlotData::setColor(size_t n, QColor& color)
{
  size_t number = 0, it = 0;
  for (auto& ii : Graph)
  {
    if (ii->number != number)
    {
      number = ii->number;
      ++it;
    }
    if (it == n)
    {
      if (ii->type == PlotLineType)
      {
        std::get<PlotLine>(ii->data).pen.setColor(color);
        for (PlotLine::iterator it = std::get<PlotLine>(ii->data).begin(); it != std::get<PlotLine>(ii->data).end(); ++it)
        {
          it->pen.setColor(color);
          it->item->setPen(it->pen);
        }
      }
      else if (ii->type == PlotCircleType)
      {
        for (auto& it : std::get<PlotCircle>(ii->data).item)
        {
          QPen pen = it->pen();
          pen.setColor(color);
          std::get<PlotCircle>(ii->data).pen = pen;
          it->setPen(pen);
        }
      }
      else if (ii->type == PlotPolygonType)
      {
        for (auto& it : std::get<PlotPolygon>(ii->data).item)
        {
          QPen pen = it->pen();
          pen.setColor(color);
          std::get<PlotPolygon>(ii->data).pen = pen;
          it->setPen(pen);
        }
      }
      else
      {
        std::cout << "Something wrong\n";
      }
    }
  }
  labelColor();
}

void PlotData::clearAxes()
{
  for (size_t i = 0; i < XCoordTextItems.size(); ++i)
  {
    removeItem(XCoordTextItems[i].get());
    XCoordTextItems[i].reset(nullptr);
  }
  XCoordTextItems.clear();
  for (size_t i = 0; i < YCoordTextItems.size(); ++i)
  {
    removeItem(YCoordTextItems[i].get());
    YCoordTextItems[i].reset(nullptr);
  }
  YCoordTextItems.clear();
}

void PlotData::clearAll()
{
//   for (auto& ii : Graph)
//   {
//     if (ii->type == PlotLineType)
//     {
//       for (PlotLine::iterator it = std::get<PlotLine>(ii->data).begin(); it != std::get<PlotLine>(ii->data).end(); ++it)
//       {
//         it->item.reset(nullptr);
//       }
// //       ii->data.line.reset();
// //       (*i).data.line = nullptr;
//     }
//     else if (ii->type == PlotCircleType)
//     {
//       for (auto& it : std::get<PlotCircle>(ii->data).item)
//       {
//         it.reset(nullptr);
//       }
//     }
//     else if (ii->type == PlotPolygonType)
//     {
//       for (auto& it : std::get<PlotPolygon>(ii->data).item)
//       {
//         it.reset(nullptr);
//       }
//     }
//     else
//     {
//       std::cout << "Something wrong\n";
//     }
//   }
  Graph.clear();
  clearAxes();
  XCoordText.clear();
  YCoordText.clear();
  const ViewBox cvb = defaultBox;
  ZoomHistory.clear();
  ZoomHistory.push_back(cvb);
  currZoom = ZoomHistory.begin();
  unitCircleCount = 0;
  if (unitCircleItem != 0)
  {
    removeItem(unitCircleItem.get());
    unitCircleItem.reset(nullptr);
//     unitCircleItem = 0;
  }
}

/// this computes the minimum and maximum value of an axis
static inline void adjustAxis(qreal& min, qreal& max, size_t& numTicks)
{
  // keep invalid axis
  if (numTicks == 1) return;
  if (min == max) return;
  const size_t MinTicks = 4;
  qreal grossStep = (max - min) / MinTicks;
  qreal step = pow(10, floor(log10(grossStep)));

  if (5 * step < grossStep) step *= 5;
  else if (2 * step < grossStep) step *= 2;

  numTicks = static_cast<size_t>(fabs(ceil(max / step) - floor(min / step)));
  min = floor(min / step) * step;
  max = ceil(max / step) * step;
}

void PlotData::addPlotLine(std::list<std::unique_ptr<PlotItem>>::iterator it, const QPen& pen, bool principal, bool stab)
{
  (*it)->type = PlotLineType;
  (*it)->data = PlotLine(pen);
  std::get<PlotLine>((*it)->data).pen.setWidthF(1);
  if (stab) std::get<PlotLine>((*it)->data).pen.setStyle(Qt::SolidLine);
  else std::get<PlotLine>((*it)->data).pen.setStyle(Qt::DashLine);
  (*it)->principal = principal;
}

void PlotData::addPlotPoint(std::list<std::unique_ptr<PlotItem>>::iterator it, const QPen& pen,
                            PlotMarkerStyle type, bool principal, qreal radius, bool scale)
{
  switch (type)
  {
    case PlotMarkerCircle:  // CIRCLE
      {
        (*it)->type = PlotCircleType;
        (*it)->data = PlotCircle(pen, QRectF(-radius, -radius, 2*radius, 2*radius), scale);
        std::get<PlotCircle>((*it)->data).pen.setWidthF(1);
      }
      break;
    case PlotMarkerSquare:  // SQUARE
      {
        (*it)->type = PlotPolygonType;
        QPolygonF pl(4);
        pl[0] = QPointF(-radius, -radius);
        pl[1] = QPointF(-radius, radius);
        pl[2] = QPointF(radius, radius);
        pl[3] = QPointF(radius, -radius);
        (*it)->data = PlotPolygon(pen, pl);
        std::get<PlotPolygon>((*it)->data).pen.setWidthF(1);
      }
      break;
    case PlotMarkerTriangle: // TRIANGLE
      {
        (*it)->type = PlotPolygonType;
        QPolygonF pl(3);
        pl[0] = QPointF(-radius, radius/sqrt(3));
        pl[1] = QPointF(0.0, -radius*2/sqrt(3));
        pl[2] = QPointF(radius, radius/sqrt(3));
        (*it)->data = PlotPolygon(pen, pl);
        std::get<PlotPolygon>((*it)->data).pen.setWidthF(1);
      }
      break;
    case PlotMarkerCross: // CROSS
      {
        (*it)->type = PlotPolygonType;
        QPolygonF pl(6);
        pl[0] = QPointF(-radius, radius);
        pl[1] = QPointF(radius, -radius);
        pl[2] = QPointF(0.0, 0.0);
        pl[3] = QPointF(-radius, -radius);
        pl[4] = QPointF(radius, radius);
        pl[5] = QPointF(0.0, 0.0);
        (*it)->data = PlotPolygon(pen, pl);
        std::get<PlotPolygon>((*it)->data).pen.setWidthF(1);
      }
      break;
    default:
      std::cout << "Another Serious BUG\n";
      abort();
      break;
  }
  (*it)->principal = principal;
}

// this only called by addPlot...
// begin and end refers to what was added to the graphics
void PlotData::dataToGraphics(std::list<std::unique_ptr<PlotItem>>::const_iterator begin,
                              std::list<std::unique_ptr<PlotItem>>::const_iterator end)
{
  bool toZoom = false;
  // if current is the top of zoom levels
  if (ZoomHistory.begin() == currZoom)
  {
    ViewBox cvb = *currZoom;
    bool touched = false;
    std::list<std::unique_ptr<PlotItem>>::const_iterator it;
//    std::cout << "ENTER ymin " << cvb.ymin << ", ymax" << cvb.ymax << " ticks " << cvb.yticks << "\n";
    for (it = begin; it != end; it++)
    {
      for (size_t k = 0; k < (*it)->x.size(); k++)
      {
        if (xcoord((*it)->x(k)) > cvb.xmax) { cvb.xmax = xcoord((*it)->x(k)); toZoom = true; }
        if (xcoord((*it)->x(k)) < cvb.xmin) { cvb.xmin = xcoord((*it)->x(k)); toZoom = true; }
        if (ycoord((*it)->y(k)) > cvb.ymax) { cvb.ymax = ycoord((*it)->y(k)); toZoom = true; }
        if (ycoord((*it)->y(k)) < cvb.ymin) { cvb.ymin = ycoord((*it)->y(k)); toZoom = true; }
//        std::cout << "y=" << ycoord(it->y(k)) << "\n";
        touched = true;
      }
    }
//    std::cout << "MID ymin " << cvb.ymin << ", ymax" << cvb.ymax << " ticks " << cvb.yticks << "\n";
    if (touched)
    {
      if (cvb.xmax == cvb.xmin ) { cvb.xmin *= 0.95; cvb.xmax *= 1.05; }
      if (cvb.ymax == cvb.ymin ) { cvb.ymin *= 0.95; cvb.ymax *= 1.05; }
      // make them valid axes, in case it was invalid
      cvb.xticks = std::max(cvb.xticks, (size_t)2u);
      cvb.yticks = std::max(cvb.yticks, (size_t)2u);
      adjustAxis(cvb.xmin, cvb.xmax, cvb.xticks);
      adjustAxis(cvb.ymin, cvb.ymax, cvb.yticks);
    }
    *currZoom = cvb;
//    std::cout << "EXIT ymin " << cvb.ymin << ", ymax" << cvb.ymax << " ticks " << cvb.yticks << "\n";
  }
  if (toZoom)
  {
    rescaleData(Graph.begin(), Graph.end());
    PlotPaint(Graph.begin(), Graph.end(), true);
  } else
  {
    rescaleData(begin, end);
    PlotPaint(begin, end, false);
  }
  labelColor();
}

// Need to simplify this to less cases
// Draw a line until stability changes
// 1. find out if it is stable start with the solid line
// 2. find out where is the next zero crossing point and draw line until
// 3. repeat antuil the end
// 4. put bifurcation points
// Cases
// 1. x > XSeparator, y > YSeparator
//    - We draw stability in this case with dashed lines
// 2. everything else doesn't change.
bool PlotData::addPlot(const KNDataFile* data, PlotXVariable x, PlotYVariable y,
  size_t pt, size_t dim)
{
  size_t xadded = 0;
  size_t yadded = 0;
  std::list<std::unique_ptr<PlotItem>>::iterator start = --Graph.end();
  size_t startnumber = 0;
  if (!Graph.empty()) startnumber = (*Graph.rbegin())->number;
  QColor plcolor, stabcolor(Qt::black);
  plcolor.setHsv((240 + startnumber*40)%360, 255, 255);
  // get stability data
  bool stab_ini = data->getUnstableMultipliers(0) == 0;
  std::vector<size_t> stabidx; // stability change + 2: the beginning and the end are also included
  std::vector<size_t> bifidx;  // bifurcation points
  std::vector<BifType> biftype; // the type of bifurcation points
  {
    size_t k_p = 0, k;
    bool stab = false;
    stabidx.push_back(0);
    do
    {
      k = data->getNextBifurcation(k_p, &stab);
      if (k < data->getNPoints())
      {
        bifidx.push_back(k);
        biftype.push_back(data->getBifurcationType(k));
        if (stab) stabidx.push_back(k);
        k_p = k;
      }
    } while (k < data->getNPoints());
    if (data->getNPoints() > 0) stabidx.push_back(data->getNPoints()-1);
    else stabidx.push_back(0);
  }
  // Putting in the data
  if (x > XSeparator && y > YSeparator)
  {
    double pxend = 0.0, pyend = 0.0;
    for (size_t b = 1; b < stabidx.size(); ++b)
    {
      std::cout << "npoints " << data->getNPoints() << " b " << stabidx[b] << " b-1 " << stabidx[b-1] << "\n";
      Graph.push_back(std::make_unique<PlotItem>(data, PlotBasicData, x, y, pt, dim));
      int bskip = 0, eskip = 0;
      if (b == 1) bskip = 0; else bskip = 1;
      if (b == stabidx.size()-1) eskip = 0; else eskip = 1;
      std::cout << "Allocating size " << stabidx[b] - stabidx[b-1] + bskip + eskip << "\n";
      (*Graph.rbegin())->x.init(stabidx[b] - stabidx[b-1] + bskip + eskip);
      (*Graph.rbegin())->y.init(stabidx[b] - stabidx[b-1] + bskip + eskip);

      // X Coordinate
      if (x == XLabel)
      {
        for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j) (*Graph.rbegin())->x(j) = i;
      } else
      if (x >= XParameter0)
      {
        for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
        {
          if (x != XParameter0)
            (*Graph.rbegin())->x(j) = data->getPar(i, x - XParameter0);
          else
            (*Graph.rbegin())->x(j) =
              data->getPar(i, x - XParameter0)/data->getPar(i, data->getNPar()+VarPeriod-VarEnd);
        }
      } else std::cout << "Bad X coord\n";
      ++xadded;
      // Y Coordinate
      if (y == YL2Norm)
      {
        KNVector elem(false);
        const_cast<KNDataFile*>(data)->getElemRef(0, elem);
        KNMatrix metric(elem.size(), elem.size());
        KNDdeBvpCollocation::getMetric(metric, elem);
        KNVector prof(false), msh(false);
        for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
        {
          const_cast<KNDataFile*>(data)->getMeshRef(i, msh);
          const_cast<KNDataFile*>(data)->getProfileRef(i, prof);
          (*Graph.rbegin())->y(j) = sqrt(KNDdeBvpCollocation::integrate(prof, prof, metric, msh, data->getNDim()));
        }
      } else
      if (y == YAmplitude)
      {
        const size_t ndeg = data->getNDeg();
        const size_t nint = data->getNInt();
//        std::cout << "NINT " << nint << " NDEG " << ndeg << "\n";
        for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
        {
          double nrm = 0.0;
          size_t p = dim;
//          for (int p = 0; p < data->getNDim(); ++p)
//          {
            double min = data->getProfile(i, p, 0);
            double max = data->getProfile(i, p, 0);
            for (size_t l = 0; l < nint; ++l)
            {
              for (size_t k = 0; k < ndeg; ++k)
              {
                if (min > data->getProfile(i, p, k + ndeg*l)) min = data->getProfile(i, p, k + ndeg*l);
                if (max < data->getProfile(i, p, k + ndeg*l)) max = data->getProfile(i, p, k + ndeg*l);
              }
            }
            nrm += (max - min) * (max - min);
//          }
//          std::cout << "min " << min << " max " << max << " nrm " << nrm << " j=" << j <<"\n";
          if (nint == 0) (*Graph.rbegin())->y(j) = data->getProfile(i, p, 0);
          else (*Graph.rbegin())->y(j) = sqrt(nrm);
        }
      } else
      if (y == YPoincare)
      {
	// displays the 0-th point, which is also the last point of the periodic orbit
	// only works for forced systems
	// TODO: take into account the period, so that each cross section is counted
        for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
        {
          (*Graph.rbegin())->y(j) = data->getProfile(i, dim, 0);
        }
      } else
      if (y >= YParameter0)
      {
        for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
        {
          if (y != YParameter0)
            (*Graph.rbegin())->y(j) = data->getPar(i, y - YParameter0);
          else
            (*Graph.rbegin())->y(j) =
              data->getPar(i, y - YParameter0)/data->getPar(i,data->getNPar()+VarPeriod-VarEnd);
        }
      } else std::cout << "Bad Y coord\n";
      ++yadded;
      // set the beginning and the end
      const size_t end = (*Graph.rbegin())->x.size()-1;
      if (bskip != 0)
      {
        (*Graph.rbegin())->x(0) = pxend;
        (*Graph.rbegin())->y(0) = pyend;
      }
      if (eskip != 0)
      {
        // The end is incorrectly set.
        pxend = ((*Graph.rbegin())->x(end) + (*Graph.rbegin())->x(end-1))/2;
        pyend = ((*Graph.rbegin())->y(end) + (*Graph.rbegin())->y(end-1))/2;
        (*Graph.rbegin())->x(end) = pxend;
        (*Graph.rbegin())->y(end) = pyend;
      }
      addPlotLine(--Graph.end(), QPen(plcolor), true, stab_ini);
      stab_ini = !stab_ini;
    }
  }
  // Extra bits
  // The profile
  if (x == XMesh && y == YProfile)
  {
    const size_t ndeg = data->getNDeg();
    const size_t nint = data->getNInt();
    Graph.push_back(std::make_unique<PlotItem>(data, PlotBasicData, x, y, pt, dim));
    (*Graph.rbegin())->x.init(ndeg*nint + 1);
    (*Graph.rbegin())->y.init(ndeg*nint + 1);
    for (size_t i = 0; i < nint; i++)
    {
      for (size_t j = 0; j < ndeg; j++)
      {
        (*Graph.rbegin())->x(j + ndeg*i) = data->getMesh(pt, i) + data->getElem(pt, j) * (data->getMesh(pt, i + 1) - data->getMesh(pt, i));
        (*Graph.rbegin())->y(j + ndeg*i) = data->getProfile(pt, dim, j + ndeg * i);
      }
    }
    (*Graph.rbegin())->x(ndeg*nint) = data->getMesh(pt, nint);
    (*Graph.rbegin())->y(ndeg*nint) = data->getProfile(pt, dim, ndeg * nint);
    ++xadded;
    ++yadded;
    addPlotLine(--Graph.end(), QPen(plcolor), true);
  }
  // Unit circle
  if (x == XRealMultiplier && y == YImagMultiplier)
  {
    // plotting the multipliers
    Graph.push_back(std::make_unique<PlotItem>(data, PlotBasicData, x, y, pt, dim));
    (*Graph.rbegin())->x.init(data->getNMul());
    (*Graph.rbegin())->y.init(data->getNMul());
    for (size_t i = 0; i < data->getNMul(); i++)
    {
      (*Graph.rbegin())->x(i) = data->getMulRe(pt, i);
      (*Graph.rbegin())->y(i) = data->getMulIm(pt, i);
    }
    ++xadded;
    ++yadded;
    addPlotPoint(--Graph.end(), QPen(plcolor), PlotMarkerCross, true);
    ++unitCircleCount;
  }
  // multipliers
  if ((x == XLabel || x >= XParameter0) && y == YAbsMultiplier)
  {
    for (size_t r = 0; r < data->getNMul(); r++)
    {
      Graph.push_back(std::make_unique<PlotItem>(data, PlotBasicData, x, y, pt, dim));
      (*Graph.rbegin())->x.init(data->getNPoints());
      (*Graph.rbegin())->y.init(data->getNPoints());
      for (size_t i = 0; i < data->getNPoints(); i++)
      {
        if (x >= XParameter0) (*Graph.rbegin())->x(i) = data->getPar(i, x - XParameter0);
        else (*Graph.rbegin())->x(i) = i;
        (*Graph.rbegin())->y(i) = sqrt(data->getMulRe(i, r) * data->getMulRe(i, r) + data->getMulIm(i, r) * data->getMulIm(i, r));
      }
      ++xadded;
      ++yadded;
      addPlotLine(--Graph.end(), QPen(plcolor), true);
    }
  }
  // add stability
  if ((x > XSeparator && y > YSeparator) || y == YAbsMultiplier)
  {
    auto itc = Graph.end();
    for (size_t i = 0; i < xadded; ++i) --itc;
    for (size_t i = 0; i < bifidx.size(); ++i)
    {
      Graph.push_back(std::make_unique<PlotItem>(data, PlotStability, x, y, pt, dim));
      (*Graph.rbegin())->x.init(1);
      (*Graph.rbegin())->y.init(1);

      if (y != YAbsMultiplier)
      {
        // find out in which stability region is our point
        size_t k = 1, k2;
        while ((bifidx[i] > stabidx[k]) && k < stabidx.size()) ++k;
        const bool stbif = (bifidx[i] == stabidx[k]);
        if (k > 1) k2 = bifidx[i] - stabidx[k-1] + 1;
        else k2 = bifidx[i] - stabidx[k-1];

        auto it = itc;
        for (size_t p = 0; p < k-1; ++p) ++it;

        if (stbif)
        {
          (*Graph.rbegin())->x(0) = (*it)->x((*it)->x.size()-1);
          (*Graph.rbegin())->y(0) = (*it)->y((*it)->y.size()-1);
        } else
        {
          (*Graph.rbegin())->x(0) = ((*it)->x(k2 - 1) + xcoord((*it)->x(k2))) / 2.0;
          (*Graph.rbegin())->y(0) = ((*it)->y(k2 - 1) + (*it)->y(k2)) / 2.0;
        }
      } else
      {
        (*Graph.rbegin())->x(0) = ((*itc)->x(bifidx[i] - 1) + (*itc)->x(bifidx[i])) / 2.0;
        (*Graph.rbegin())->y(0) = 1.0;
      }
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
  for (auto it = ++start; it != Graph.end(); ++it)
  {
    (*it)->number = startnumber + 1;
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

      auto it = Graph.end();
      for (size_t i = 0; i < xadded; ++i){ --it; if ((*it)->x.size() != (*it)->y.size()) std::cout << "bad number of X and Y coordinates\n"; }
      dataToGraphics(it, Graph.end());
      labelColor();
      return true;
    }
    else
    {
      return false;
    }
  }
}

struct rcStruc
{
  rcStruc(size_t n, PlotXVariable x, PlotYVariable y, size_t p, size_t d) :
    num(n), pt(p), dim(d), X(x), Y(y) {}
  const size_t num;
  const size_t pt;
  const size_t dim;
  const PlotXVariable X;
  const PlotYVariable Y;
};

void PlotData::updatePlot(const KNDataFile* mat)
{
  size_t number = 0;
  bool dirty = false;

  std::list<rcStruc> lst;

  size_t ct_num = 0;
  size_t ct_it = 0;
  auto first = Graph.end();
  auto last = Graph.end();
  for (auto i = Graph.begin(); i != Graph.end(); i++)
  {
    if (ct_num != (*i)->number) { ct_num = (*i)->number; ++ct_it; }
    // only update if all the data is necessary
    if ((*i)->isFrom(mat) && (*i)->varX > XSeparator && (*i)->varY > YSeparator)
    {
      // only select one of the, but recreate all of them
      if (number != (*i)->number)
      {
        number = (*i)->number;
        lst.push_back(rcStruc(ct_it, (*i)->varX, (*i)->varY, (*i)->point, (*i)->dimension));
      } else dirty = true;
      if (first == Graph.end()) first = i;
      last = i;
    } else dirty = true;
  }

  if (dirty)
  {
    size_t it = 0;
    for (std::list<rcStruc>::iterator i = lst.begin(); i != lst.end(); i++)
    {
      clear(i->num - it);
      ++it;
    }
  } else clearAll();
  for (std::list<rcStruc>::iterator i = lst.begin(); i != lst.end(); i++)
  {
    addPlot(mat, i->X, i->Y, i->pt, i->dim);
  }
}

// The intersection point is defined by
// A + alpha*(B-A) == C + beta*(D-C)
// if 0 < alpha <= 1 AND 0 < beta <= 1 -> intersect in the middle
static inline bool intersect(qreal& alpha, qreal& beta,
                             QPointF& itPoint,
                             const QPointF& A, const QPointF& B,
                             const QPointF& C, const QPointF& D)
{
  const qreal Vx = B.x() - A.x();
  const qreal Vy = B.y() - A.y();
  const qreal Wx = D.x() - C.x();
  const qreal Wy = D.y() - C.y();
  const qreal VpW = -Vy*Wx + Vx*Wy;
  const qreal AmCx = A.x() - C.x();
  const qreal AmCy = A.y() - C.y();
  if (VpW != 0.0)
  {
    alpha = (-Wy*AmCx + Wx*AmCy)/VpW;
    beta = (-Vy*AmCx + Vx*AmCy)/VpW;
    itPoint = QPointF(A.x() + alpha*Vx, A.y() + alpha*Vy);
    return true;
  }
  return false;
}

enum horizPosition
{
  Left = 1,
  HMiddle = 2,
  Right = 3
};

enum vertPosition
{
  Bottom = 1,
  VMiddle = 2,
  Top = 3
};

static inline horizPosition horizPoint(const QPointF& BottomRight, const QPointF& TopLeft, const QPointF& A)
{
  if (A.x() < TopLeft.x()) return Left;
  else if (A.x() <= BottomRight.x()) return HMiddle;
  else return Right;
}

static inline vertPosition vertPoint(const QPointF& BottomRight, const QPointF& TopLeft, const QPointF& A)
{
  if (A.y() < BottomRight.y()) return Bottom;
  else if (A.y() <= TopLeft.y()) return VMiddle;
  else return Top;
}

// A is the orevoius point, B is the next point
// returns the point to be plotted 'toPlot'
static inline void pointOutside(PlotLine& ppath,
                                const QPointF& BottomLeft, const QPointF& BottomRight,
                                const QPointF& TopLeft,    const QPointF& TopRight,
                                const QPointF& A,          const QPointF& B)
{
  // Handle the easy cases quickly
  const horizPosition Ahoriz = horizPoint (BottomRight, TopLeft, A);
  const horizPosition Bhoriz = horizPoint (BottomRight, TopLeft, B);
  const vertPosition Avert = vertPoint (BottomRight, TopLeft, A);
  const vertPosition Bvert = vertPoint (BottomRight, TopLeft, B);

  if (Ahoriz == HMiddle && Avert == VMiddle && Bhoriz == HMiddle && Bvert == VMiddle)
  {
    ppath.back().path.lineTo(B);
    return;
  }
  else if (Ahoriz == Bhoriz && Ahoriz != HMiddle) return;
  else if (Avert == Bvert && Avert != VMiddle) return;
  else if (Ahoriz == Bhoriz && Avert == Bvert) return;

  // when it can intersect
  qreal alpha[4];
  qreal beta[4];
  QPointF itPoint[4];
  bool it[4]{false,false,false,false};
  it[0] = intersect (alpha[0], beta[0], itPoint[0], BottomLeft,  BottomRight, A, B);
  it[1] = intersect (alpha[1], beta[1], itPoint[1], BottomRight, TopRight, A, B);
  it[2] = intersect (alpha[2], beta[2], itPoint[2], TopRight,    TopLeft, A, B);
  it[3] = intersect (alpha[3], beta[3], itPoint[3], TopLeft,     BottomLeft, A, B);
  int oneSideIt = 0;
  int twoSideIt = 0;
  int oneSidePt[4]{0,0,0,0};
  int twoSidePt[4]{0,0,0,0};
//   int otherSidePt[4];
  for (int i = 0; i < 4; ++i)
  {
    if (it[i])
    {
      if (alpha[i] > 0 && alpha[i] <= 1 && beta[i] > 0 && beta[i] <= 1) { twoSidePt[twoSideIt] = i; ++twoSideIt; }
      else if (alpha[i] > 0 && alpha[i] <= 1) { oneSidePt[oneSideIt] = i; ++oneSideIt; }
    }
  }
  if (oneSideIt == 2)
  {
    if (beta[oneSidePt[0]]*beta[oneSidePt[1]] < 0)
    {
      ppath.back().path.lineTo(B);
//       std::cout << "Should have been handled already";
    }
  } else
  if (twoSideIt == 1)
  {
    if (beta[oneSidePt[0]] > 1)
    {
      ppath.push_back(PlotPolyLine(ppath.pen));
      ppath.back().path.moveTo(itPoint[twoSidePt[0]]);
      ppath.back().path.lineTo(B);
//       std::cout << "-in-";
    } else
    {
      ppath.back().path.lineTo(itPoint[twoSidePt[0]]);
//       std::cout << "-out-";
    }
  } else
  if (twoSideIt == 2)
  {
    ppath.push_back(PlotPolyLine(ppath.pen));
    ppath.back().path.moveTo(itPoint[twoSidePt[0]]);
    ppath.back().path.lineTo(itPoint[twoSidePt[1]]);
//     std::cout << "-cross-";
  }
//   std::cout.flush();
}

void PlotData::rescaleData(std::list<std::unique_ptr<PlotItem>>::const_iterator begin,
                           std::list<std::unique_ptr<PlotItem>>::const_iterator end)
{
  ViewBox cvb = *currZoom;

  if (cvb.xmax == cvb.xmin) return;
  if (cvb.ymax == cvb.ymin) return;

  const double xscale = plotXSize / (cvb.xmax - cvb.xmin);
  const double yscale = plotYSize / (cvb.ymax - cvb.ymin);

  const QPointF BottomLeft(0,0), BottomRight(plotXSize,0),
                TopLeft(0,plotYSize), TopRight(plotXSize,plotYSize);
  // rescaling all the data
  for (auto i = begin; i != end; i++)
  {
    if ((*i)->x.size() == 0) continue;
    if ((*i)->type == PlotLineType)
    {
      auto& ppath = std::get<PlotLine>((*i)->data);
      QPointF prevPoint(xscale*(xcoord((*i)->x(0)) - cvb.xmin), yscale*(cvb.ymax - ycoord((*i)->y(0))));

      // initializing ppath
      ppath.clear();
      ppath.push_back(PlotPolyLine(ppath.pen));
      ppath.back().path.moveTo(prevPoint);

      for (size_t k = 1; k < (*i)->x.size(); k++)
      {
        QPointF currentPoint(xscale*(xcoord((*i)->x(k)) - cvb.xmin), yscale*(cvb.ymax - ycoord((*i)->y(k))));
        pointOutside(ppath, BottomLeft, BottomRight, TopLeft, TopRight, prevPoint, currentPoint);
        prevPoint = currentPoint;
      }
    }
    if ((*i)->type == PlotCircleType)
    {
//       std::cout << "Scaling PlotCircle\n";
      std::get<PlotCircle>((*i)->data).pos.clear();
      for (size_t k = 0; k < (*i)->x.size(); k++)
      {
        const QPointF pt = QPointF(xscale * (xcoord((*i)->x(k)) - cvb.xmin), yscale * (cvb.ymax - ycoord((*i)->y(k))));
        QRectF& rect = std::get<PlotCircle>((*i)->data).point;
        QRectF& scaledRect = std::get<PlotCircle>((*i)->data).scaledPoint;
        if (std::get<PlotCircle>((*i)->data).scale)
        {
          scaledRect.setLeft(xscale*rect.left());
          scaledRect.setRight(xscale*rect.right());
          scaledRect.setBottom(yscale*rect.bottom());
          scaledRect.setTop(yscale*rect.top());
          if (Box->rect().intersects(scaledRect)) std::get<PlotCircle>((*i)->data).pos.push_back(pt);
        } else
        {
          scaledRect = rect;
          if (Box->rect().contains(pt.x(), pt.y())) std::get<PlotCircle>((*i)->data).pos.push_back(pt);
        }
      }
//      for (auto& it : std::get<PlotCircle>((*i)->data).item) it.reset(nullptr);
      std::get<PlotCircle>((*i)->data).item.resize(std::get<PlotCircle>((*i)->data).pos.size());
    }
    if ((*i)->type == PlotPolygonType)
    {
      std::get<PlotPolygon>((*i)->data).pos.clear();
      for (size_t k = 0; k < (*i)->x.size(); k++)
      {
        const QPointF pt = QPointF(xscale * (xcoord((*i)->x(k)) - cvb.xmin), yscale * (cvb.ymax - ycoord((*i)->y(k))));
        if (Box->rect().contains(pt.x(), pt.y())) std::get<PlotPolygon>((*i)->data).pos.push_back(pt);
      }
//      for (auto& it : std::get<PlotPolygon>((*i)->data).item) it.reset(nullptr);
      std::get<PlotPolygon>((*i)->data).item.resize(std::get<PlotPolygon>((*i)->data).pos.size());
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
    if ((cvb.xmax == cvb.xmin) || (cvb.ymax == cvb.ymin))
    {
      event->accept();
      return;
    }
    if ((cvb.xticks > 1)&&(cvb.yticks > 1))
    {
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
      dataToGraphics(Graph.begin(), Graph.end());
    }
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
      dataToGraphics(Graph.begin(), Graph.end());
    }
    key->accept();
  }
  if (key->key() == Qt::Key_N)
  {
    if (--ZoomHistory.end() != currZoom)
    {
      ++currZoom;
      dataToGraphics(Graph.begin(), Graph.end());
    }
    key->accept();
  }
}

QGraphicsRectItem* PlotData::makeBox()
{
  ViewBox cvb = *currZoom;
  // drawing the box
  QGraphicsRectItem * newBox = addRect(QRectF(0.0, 0.0, plotXSize, plotYSize), QPen(QBrush(Qt::SolidPattern), 2.0));
  newBox->setFlags(newBox->flags() | QGraphicsItem::ItemClipsChildrenToShape);
  // drawing the ticks and tickmarks
  for (size_t i = 0; i < HText.size(); i++)
  {
    removeItem(BottomTicks[i].get());
    removeItem(TopTicks[i].get());
    removeItem(HText[i].get());
//     BottomTicks[i].reset(nullptr);
//     TopTicks[i].reset(nullptr);
//     HText[i].reset(nullptr);
  }
  BottomTicks.resize(cvb.xticks + 1);
  TopTicks.resize(cvb.xticks + 1);
  HText.resize(cvb.xticks + 1);
  qreal highest = 0;
//  std::cout << "xmin " << cvb.xmin << ", xmax" << cvb.xmax << " ticks " << cvb.xticks << "\n";
  for (size_t i = 0; i < cvb.xticks + 1; i++)
  {
    BottomTicks[i] = std::make_unique<QGraphicsLineItem>();
    BottomTicks[i]->setLine(plotXSize * i / cvb.xticks, 0.0, plotXSize * i / cvb.xticks, 5.0);
    BottomTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    TopTicks[i] = std::make_unique<QGraphicsLineItem>();
    TopTicks[i]->setLine(plotXSize * i / cvb.xticks, plotYSize, plotXSize * i / cvb.xticks, plotYSize - 5.0);
    TopTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    HText[i] = std::make_unique<QGraphicsTextItem>();
    HText[i]->setPlainText(QString::number(cvb.xmin + (cvb.xmax - cvb.xmin)*i / cvb.xticks));
    HText[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = HText[i]->boundingRect().normalized();
    HText[i]->setPos(plotXSize * i / cvb.xticks - b.width() / 2.0, plotYSize /*- b.height()*/);
    addItem(TopTicks[i].get());
    addItem(BottomTicks[i].get());
    addItem(HText[i].get());
//     std::cout << "addItem " << i << " " << __LINE__ << "" << TopTicks[i].get() << " " << BottomTicks[i].get() << " " << HText[i].get() << "\n";
    if (highest < b.height()) highest = b.height();
  }
  for (size_t i = 0; i < VText.size(); i++)
  {
    removeItem(LeftTicks[i].get());
    removeItem(RightTicks[i].get());
    removeItem(VText[i].get());
//     LeftTicks[i].reset(nullptr);
//     RightTicks[i].reset(nullptr);
//     VText[i].reset(nullptr);
  }
  LeftTicks.resize(cvb.yticks + 1);
  RightTicks.resize(cvb.yticks + 1);
  VText.resize(cvb.yticks + 1);
  qreal widest = 0;
//  std::cout << "ymin " << cvb.ymin << ", ymax" << cvb.ymax << " ticks " << cvb.yticks << "\n";
  for (size_t i = 0; i < cvb.yticks + 1; i++)
  {
    LeftTicks[i] = std::make_unique<QGraphicsLineItem>();
    LeftTicks[i]->setLine(0.0, plotYSize * i / cvb.yticks, 5.0, plotYSize * i / cvb.yticks);
    LeftTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    RightTicks[i] = std::make_unique<QGraphicsLineItem>();
    RightTicks[i]->setLine(plotXSize, plotYSize * i / cvb.yticks, plotXSize - 5.0, plotYSize * i / cvb.yticks);
    RightTicks[i]->setPen(QPen(QBrush(Qt::SolidPattern), 1.0));
    VText[i] = std::make_unique<QGraphicsTextItem>();
    VText[i]->setPlainText(QString::number(cvb.ymin + (cvb.ymax - cvb.ymin)*i / cvb.yticks));
    VText[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = VText[i]->boundingRect().normalized();
    VText[i]->setPos(-b.width(), plotYSize - plotYSize*i / cvb.yticks - b.height() / 2.0);
    addItem(LeftTicks[i].get());
    addItem(RightTicks[i].get());
    addItem(VText[i].get());
//     std::cout << "addItem " << __LINE__ << "\n";
    if (widest < b.width()) widest = b.width();
  }
  // remove previous texts
  clearAxes();
  // Add labels
  qreal sumwidth = 0;
  qreal sumheight = 0;
  size_t startmove = 0;
  for (size_t i = 0; i < YCoordText.size(); ++i)
  {
    YCoordTextItems.push_back(std::make_unique<QGraphicsTextItem>());
    YCoordTextItems[i]->setPlainText(YCoordText[i]);
    YCoordTextItems[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = YCoordTextItems[i]->boundingRect().normalized();
    YCoordTextItems[i]->setPos(-widest - b.height() - sumheight, plotYSize/2.0 + sumwidth + b.width());
    YCoordTextItems[i]->setTransform(QTransform().rotate(-90));
    sumwidth += b.width();
    if (sumwidth > plotYSize)
    {
      for (size_t i = startmove; i < YCoordTextItems.size(); ++i)
      {
        YCoordTextItems[i]->moveBy(0.0, -sumwidth/2.0);
        addItem(YCoordTextItems[i].get());
//         std::cout << "addItem " << __LINE__ << "\n";
      }
      sumheight += highest;
      startmove = YCoordTextItems.size();
      sumwidth = 0;
    }
  }
  for (size_t i = startmove; i < YCoordTextItems.size(); ++i)
  {
    YCoordTextItems[i]->moveBy(0.0, -sumwidth/2.0);
    addItem(YCoordTextItems[i].get());
//     std::cout << "addItem " << __LINE__ << "\n";
  }
  sumwidth = 0;
  sumheight = 0;
  startmove = 0;
  for (size_t i = 0; i < XCoordText.size(); ++i)
  {
    XCoordTextItems.push_back(std::make_unique<QGraphicsTextItem>());
    XCoordTextItems[i]->setPlainText(XCoordText[i]);
    XCoordTextItems[i]->setFont(QFont("Helvetica", FontSize));
    QRectF b = XCoordTextItems[i]->boundingRect().normalized();
    XCoordTextItems[i]->setPos(plotXSize / 2.0 + sumwidth, plotYSize + highest + sumheight);
    sumwidth += b.width();
    if (sumwidth > plotXSize)
    {
      for (size_t i = startmove; i < XCoordTextItems.size(); ++i)
      {
        XCoordTextItems[i]->moveBy(-sumwidth/2.0, 0.0);
        addItem(XCoordTextItems[i].get());
//         std::cout << "addItem " << __LINE__ << "\n";
      }
      sumheight += highest;
      startmove = XCoordTextItems.size();
      sumwidth = 0;
    }
  }
  for (size_t i = startmove; i < XCoordText.size(); ++i)
  {
    XCoordTextItems[i]->moveBy(-sumwidth/2.0, 0.0);
    addItem(XCoordTextItems[i].get());
//     std::cout << "addItem " << __LINE__ << "\n";
  }

  // add unit circle if necessary
  if (unitCircleItem != 0)
  {
    removeItem(unitCircleItem.get());
    unitCircleItem.reset(nullptr);
    unitCircleItem = 0;
  }
  if (unitCircleCount > 0)
  {
    if ((cvb.xmax == cvb.xmin) || (cvb.ymax == cvb.ymin)) return newBox;
    const double xscale = plotXSize / (cvb.xmax - cvb.xmin);
    const double yscale = plotYSize / (cvb.ymax - cvb.ymin);
    const QPointF pos = QPointF(xscale * (xcoord(0.0) - cvb.xmin), yscale * (cvb.ymax - ycoord(0.0)));
    QRectF scaledRect(-xscale, -yscale, 2*xscale, 2*yscale);
    unitCircleItem = std::make_unique<QGraphicsEllipseItem>(scaledRect, newBox);
    unitCircleItem->setPos(pos);
    unitCircleItem->setPen(QPen(QBrush(Qt::SolidPattern), 1));
  }
  return newBox;
}

void PlotData::labelColor()
{
  size_t prenumber = 0, it = 0;
  for (auto i = Graph.begin(); i != Graph.end(); i++)
  {
    if ((*i)->principal && prenumber != (*i)->number)
    {
      QColor textColor;
      if ((*i)->type == PlotLineType)
      {
        textColor = std::get<PlotLine>((*i)->data).pen.color();
      }
      if ((*i)->type == PlotCircleType)
      {
        textColor = std::get<PlotCircle>((*i)->data).pen.color();
      }
      if ((*i)->type == PlotPolygonType)
      {
        textColor = std::get<PlotPolygon>((*i)->data).pen.color();
      }
      size_t num = it++;
      if (num < XCoordTextItems.size() && num < YCoordTextItems.size() &&
          num < XCoordText.size() && num < YCoordText.size()) // safety measures
      {
        XCoordTextItems[num]->setDefaultTextColor(textColor);
        YCoordTextItems[num]->setDefaultTextColor(textColor);
      }
    }
    prenumber = (*i)->number;
  }
}

void PlotData::PlotPaint(std::list<std::unique_ptr<PlotItem>>::const_iterator begin,
                         std::list<std::unique_ptr<PlotItem>>::const_iterator end, bool zoom)
{
  // drawing the box
  QGraphicsRectItem* newBox = makeBox();
  // drawing the lines
  for (auto i = begin; i != end; i++)
  {
    QColor textColor;
    if ((*i)->type == PlotLineType)
    {
      for (auto& it : std::get<PlotLine>((*i)->data))
      {
//         delete it->item;
//         it.item = std::make_unique<QGraphicsPathItem>(it.path, newBox);
        // this is managed by unique_ptr, hence no need for parent
        if(it.item) 
        {
//          std::cout << "remove a line " << it.item.get() << " scene " << it.item->scene() << "\n";
          removeItem(it.item.get());
        }
        it.item = std::make_unique<QGraphicsPathItem>(it.path);
        it.item->setPen(it.pen);
//        std::cout << "<adding a line " << it.item.get() << " scene " << it.item.get()->scene() << "\n";
        addItem(it.item.get());
//        std::cout << "<adding a line " << it.item.get() << " scene " << it.item->scene() << ">\n";
      }
    }
    if ((*i)->type == PlotCircleType)
   {
//      std::cout << "adding a circle ";
      // need to iterate through both pos and item at the same time
      {
        auto itp = std::get<PlotCircle>((*i)->data).pos.begin();
        auto iti = std::get<PlotCircle>((*i)->data).item.begin();
        for ( ; itp != std::get<PlotCircle>((*i)->data).pos.end(); itp++, iti++)
        {
          // this is managed by unique_ptr, hence no need for parent
          auto item = std::make_unique<QGraphicsEllipseItem>(std::get<PlotCircle>((*i)->data).scaledPoint);
          item->setPen(std::get<PlotCircle>((*i)->data).pen);
          item->setPos(*itp);
          if (*iti) removeItem(iti->get());
//          std::cout << "adding a circle " << item.get() << "\n";
          addItem(item.get());
          *iti = std::move(item);
        }
      }
    }
    if ((*i)->type == PlotPolygonType)
    {
//      std::cout << "adding a polygon ";
      // need to iterate through both pos and item at the same time
      {
        auto itp = std::get<PlotPolygon>((*i)->data).pos.begin();
        auto iti = std::get<PlotPolygon>((*i)->data).item.begin();
        for ( ; itp != std::get<PlotPolygon>((*i)->data).pos.end(); itp++, iti++)
        {
          // this is managed by unique_ptr, hence no need for parent
          auto item = std::make_unique<QGraphicsPolygonItem>(std::get<PlotPolygon>((*i)->data).point);
          item->setPen(std::get<PlotPolygon>((*i)->data).pen);
          item->setPos(*itp);
          if (*iti) removeItem(iti->get());
//          std::cout << "adding a polygon " << item.get() << "\n";
          addItem(item.get());
          *iti = std::move(item);
        }
      }
    }
  }
  Box.reset(newBox);
  labelColor();
}
