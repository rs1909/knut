// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include <QtCharts/QChart>
#include <QtCharts/QXYSeries>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <vector>
#include <list>

#include "datachart.h"
#include "pointtype.h"
#include "mat4data.h"
#include "ncolloc.h"

using namespace QtCharts;

DataChart::DataChart(QObject *parent) : QObject(parent)
{
  chart = new QChart ();
  chart -> setTheme (QChart::ChartThemeHighContrast);
  chart -> setAnimationOptions (QChart::NoAnimation);
  hAxis = new QValueAxis ();
  hAxis -> setTitleVisible(false);
  vAxis = new QValueAxis ();
  vAxis -> setTitleVisible(false);
  chart -> addAxis (vAxis, Qt::AlignLeft);
  chart -> addAxis (hAxis, Qt::AlignBottom);
    
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

DataChart::~DataChart()
{
  clearAll();
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
DataChart::OnePlot* DataChart::addPartPlot(const KNDataFile* data, DataChart::OnePlot* opl, PlotXVariable x, PlotYVariable y,
  size_t pt, size_t dim)
{
  // sanity check
  if ( data->getNDim() <= dim )
  {
    std::cout << "DataChart::addPlot: wrong requested dimension. NDIM=" << data->getNDim() << " d=" << dim << "\n";
    return nullptr;
  }
  if ( data->getNPoints() <= pt )
  {
    std::cout << "DataChart::addPlot: wrong requested point.\n";
    return nullptr;
  }
  size_t xadded = 0;
  size_t yadded = 0;

  // if there is a new plot set colour by number
  QColor plcolor; // stabcolor(Qt::black);
  const size_t startnumber = Graph.size();
  plcolor.setHsv((240 + startnumber*40)%360, 255, 255);

  // if there is a plot here already, redifine colour
  if (opl != nullptr)
  {
    if (opl->graph != nullptr)
    {
      if (opl->graph->size() != 0) plcolor = opl->graph->at(0)->color();
    }
  }

  // get stability data
  bool stab_ini = data->getUnstableMultipliers(0) == 0;
  std::vector<size_t> stabidx;
  std::vector<size_t> bifidx;
  std::vector<BifType> biftype;
  {
    size_t k_p = 0, k;
    bool stab = false;
    stabidx.push_back(0);
    do
    {
      k = data->getNextBifurcation(k_p, &stab);
      if (k < data->getNCols())
      {
        bifidx.push_back(k);
        biftype.push_back(data->getBifurcationType(k));
        if (stab) stabidx.push_back(k);
        k_p = k;
      }
    } while (k < data->getNCols());
    stabidx.push_back(data->getNPoints()-1);
  }
  // Putting in the data
  QList<QXYSeries * > * theGraph = nullptr;
  // This is where we have to decide about each case
  if (x > XSeparator && y > YSeparator)
  {
    size_t bstart = 1;
    QList<QXYSeries * > *newGraph = nullptr;
    if (opl != nullptr)
    {
      QList<QXYSeries * > *oplGraph = opl->graph;
      int oplPointNum = 0;
      // check the number of points
      for (int k = 0; k < oplGraph->size(); ++k)
      {
        oplPointNum += oplGraph->at(k)->count();
      }
      oplPointNum -= 2*(oplGraph->size()) - 2;
      // We need to update the Graphs
      newGraph = oplGraph;
//       if ((oplGraph->size()) + 1 < stabidx.size()) std::cout << "Creating some graph : ";
//       else std::cout << "Reusing graph : ";
      // add whatever is mising
      for (size_t b = (oplGraph->size()) + 1; b < stabidx.size(); ++b) newGraph->push_back(new QLineSeries());
      if (oplPointNum <= data->getNPoints())
      {
//         std::cout << "APPEND\n";
        bstart = oplGraph->size();
      } else
      {
//         std::cout << "CLEAR\n";
        for (int k = 0; k < newGraph->size(); ++k) (*newGraph)[k]->clear();
        bstart = 0;
      }
    }
    if (newGraph == nullptr)
    {
//       std::cout << "Creating graph\n";
      // We need a new graph;
      newGraph = new QList<QXYSeries*>;
      for (size_t b = 1; b < stabidx.size(); ++b)
      {
        QLineSeries* line = new QLineSeries(); 
        if (b == 1) line -> setName(QString("%1 : d%2").arg(data->getFileName().c_str()).arg(dim));
        newGraph->push_back (line);
      }
    }
    
    // now add in the points
    QPointF pEnd (0.0, 0.0); // end of previous line
    for (size_t b = 1; b < stabidx.size(); ++b)
    {
      if (b >= bstart)
      {
  //      std::cout << "npoints " << data->getNPoints() << " b " << stabidx[b] << " b-1 " << stabidx[b-1] << "\n";
        int bskip = 0, eskip = 0;
        if (b == 1) bskip = 0; else bskip = 1;
        if (b == stabidx.size()-1) eskip = 0; else eskip = 1;
        // Creating a new line
  //       QLineSeries* lineN = new QLineSeries();
  //       lineN -> setName(QString("Name=%1 : dim=%2").arg(data->getFileName().c_str()).arg(dim));
        // creating storage for the line
        QVector<QPointF> newLine (stabidx[b] - stabidx[b-1] + bskip + eskip);
        // X Coordinate
        if (x == XLabel)
        {
          for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j) newLine[j].setX(i);
        } else
        if (x >= XParameter0)
        {
          for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
          {
            if (x != XParameter0)
              newLine[j].setX ( data->getPar(i, x - XParameter0) );
            else
              newLine[j].setX(
                data->getPar(i, x - XParameter0)/data->getPar(i, data->getNPar()+VarPeriod-VarEnd) );
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
            newLine[j].setY( sqrt(KNDdeBvpCollocation::integrate(prof, prof, metric, msh, data->getNDim())) );
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
            if (nint == 0) newLine[j].setY( data->getProfile(i, p, 0) );
            else newLine[j].setY( sqrt(nrm) );
          }
        } else
        if (y == YPoincare)
        {
      // displays the 0-th point, which is also the last point of the periodic orbit
      // only works for forced systems
      // TODO: take into account the period, so that each cross section is counted
          for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
          {
            newLine[j].setY( data->getProfile(i, dim, 0) );
          }
        } else
        if (y >= YParameter0)
        {
          for (size_t i = stabidx[b-1], j = bskip; i < stabidx[b]+eskip; ++i, ++j)
          {
            if (y != YParameter0)
              newLine[j].setY( data->getPar(i, y - YParameter0) );
            else
              newLine[j].setY(
                data->getPar(i, y - YParameter0)/data->getPar(i,data->getNPar()+VarPeriod-VarEnd) );
          }
        } else std::cout << "Bad Y coord\n";
        ++yadded;
        // set the beginning and the end
        const size_t end = newLine.size()-1;
        if (bskip != 0)
        {
          // need to attach the start to the previous line segment
          if (bstart > 1)
          {
            // if there are lines from previous call
            const int ct = newGraph->at(b-2)->count();
            if (ct > 0) newLine[0] = newGraph->at(b-2)->at(ct-1);
          } else
          {
            // if there are lines from this loop 
            newLine[0] = pEnd;
          }
        }
        if (eskip != 0)
        {
          // The end is incorrectly set.
          pEnd = (newLine[end] + newLine[end-1])/2;
          newLine[end] = pEnd;
        }
        newGraph->at(b-1)->clear();
        for (int k = newGraph->at(b-1)->count(); k < newLine.size(); ++k)
        {
          newGraph->at(b-1)->append(newLine[k]);
        }
      }
      QPen plpen(QBrush(), 1, stab_ini ? Qt::SolidLine : Qt::DashLine);
      newGraph->at(b-1) -> setPen (plpen);
      newGraph->at(b-1) -> setColor (plcolor);

      stab_ini = !stab_ini;
    }
    // newGraph ends
    theGraph = newGraph;
  }
  // Extra bits
  // TODO convert these to update as well TODO
  // The profile
  if (x == XMesh && y == YProfile)
  {
    QList<QXYSeries * > *newGraph = new QList<QXYSeries*>;
    const size_t ndeg = data->getNDeg();
    const size_t nint = data->getNInt();
    // Creating a new plot
    QLineSeries* newLine = new QLineSeries();
    newLine -> setName(QString("%1 : d%2").arg(data->getFileName().c_str()).arg(dim));

    for (size_t i = 0; i < nint; i++)
    {
      for (size_t j = 0; j < ndeg; j++)
      {
        newLine -> append (
            data->getMesh(pt, i) + data->getElem(pt, j) * (data->getMesh(pt, i + 1) - data -> getMesh(pt, i)),
            data->getProfile(pt, dim, j + ndeg * i) );
      }
    }
    newLine -> append (data->getMesh(pt, nint), data->getProfile(pt, dim, ndeg * nint));
    ++xadded;
    ++yadded;
    newLine -> setColor (plcolor);
    QPen plpen(QBrush(), 1, stab_ini ? Qt::SolidLine : Qt::DashLine);
    newLine -> setPen (plpen);
    // just a single graph, no stability
    newGraph -> push_back (newLine);
    theGraph = newGraph;
  }
  // Unit circle
  // TODO convert these to update as well TODO
  if (x == XRealMultiplier && y == YImagMultiplier)
  {
    QList<QXYSeries * > *newGraph = new QList<QXYSeries*>;
    QScatterSeries* newLine = new QScatterSeries();
    newLine -> setName(QString("%1 : d%2").arg(data->getFileName().c_str()).arg(dim));
    // plotting the multipliers
    for (size_t i = 0; i < data->getNMul(); i++)
    {
      newLine -> append (data->getMulRe(pt, i), data->getMulIm(pt, i) );
    }
    ++xadded;
    ++yadded;
    newLine -> setColor (plcolor);
    // just a single graph no stability
    newGraph -> push_back (newLine);
    theGraph = newGraph;
  }
  // multipliers
  if ((x == XLabel || x >= XParameter0) && y == YAbsMultiplier)
  {
    QList<QXYSeries * > *newGraph = nullptr;
    int stpoint = 0;
    bool grNew = true;
    if (opl != nullptr)
    {
      if (opl->graph->size() == data->getNMul())
      {
        // not to allocate, only clear here!
        grNew = false;
        if (opl->graph->at(0)->count() <= data->getNPoints())
        {
          stpoint = opl->graph->at(0)->count();
          newGraph = opl->graph;
//           std::cout << "use old np=" << stpoint << " nnew= " << data->getNPoints() << "\n";
        } else
        {
          newGraph = opl->graph;
          for (size_t r = 0; r < opl->graph->size(); r++) opl->graph->at(r)->clear();
        }
      }
    }
    // allocate new only if reset
    if (grNew)
    {
//       std::cout << "Allocating new lines\n";
      newGraph = new QList<QXYSeries*>;
      for (size_t r = 0; r < data->getNMul(); r++)
      {
        newGraph -> push_back (new QLineSeries());
      }
//       std::cout << "Creating graph\n";
    }
//     else std::cout << "Reusing graph\n";
    // fill up with data    
    for (size_t r = 0; r < data->getNMul(); r++)
    {
      QXYSeries* newLine = newGraph->at(r);
      if (r == 0) newLine -> setName(QString("%1 : d%2").arg(data->getFileName().c_str()).arg(dim));
      for (size_t i = stpoint; i < data->getNPoints(); i++)
      {
        qreal xpos;
        if (x >= XParameter0) xpos = data->getPar(i, x - XParameter0);
        else xpos = i;
        
        newLine -> append (xpos, sqrt(data->getMulRe(i, r) * data->getMulRe(i, r) + data->getMulIm(i, r) * data->getMulIm(i, r)));
      }
      ++xadded;
      ++yadded;
      newLine -> setColor (plcolor);
    }
    // NMul number of graps
    theGraph = newGraph;
  }
  // add stability
  // This assumes that theGraph is set up
  QScatterSeries* stabLine = nullptr;
  if ((x > XSeparator && y > YSeparator) || ((x == XLabel || x >= XParameter0) && y == YAbsMultiplier))
  {
    if (opl != nullptr)
    {
      // if there is a new graph we need to create new stabLine
      if ((opl->bifpoints != nullptr) && (opl->graph == theGraph))
      {
        if (bifidx.size() >= opl->bifpoints->count())
        {
          stabLine = opl->bifpoints;
//           std::cout << "Reusing stabline\n";
        } else
        {
//           std::cout << "SHOULD Clear and Reuse stabline\n";
          stabLine = opl->bifpoints;
//           for (int k = 0; k < stabLine->count(); ++k)
//           {
//             std::cout << "P(" << stabLine->at(k).x() << "," << stabLine->at(k).y() << ")\n";
//           }
          stabLine->clear(); // Bug in QtCharts
        }
      }
    }
    if (stabLine == nullptr)
    {
//       std::cout << "Creating stabline\n";
      stabLine = new QScatterSeries(this); 
    }
//     std::cout << "IN: stabLine->count() " << stabLine->count() << "bifidx.size() " << bifidx.size() << "\n";
    for (size_t i = stabLine->count(); i < bifidx.size(); ++i)
    {
      if (y != YAbsMultiplier)
      {
        // find out in which stability region is our point
        size_t k = 1, k2;
        // stabidx[k] is where the point is
        // k-1 is the index of the line
        while ((bifidx[i] > stabidx[k]) && k < stabidx.size()) ++k;
        // true if bifurcation at the end of line
        const bool stbif = (bifidx[i] == stabidx[k]);
        // the index within part of the series
        if (k > 1) k2 = bifidx[i] - stabidx[k-1] + 1;
        else k2 = bifidx[i] - stabidx[k-1];

        if (stbif)
        {
          int end = theGraph->at(k-1)->count() - 1;
          stabLine -> append (theGraph->at(k-1)->at(end));
//           std::cout << "Add at end of line" << theGraph->at(k-1)->at(end).x() << "," << theGraph->at(k-1)->at(end).y() << "\n";
        } else
        {
          stabLine -> append ((theGraph->at(k-1)->at(k2) + theGraph->at(k-1)->at(k2-1))/2.0);
//           std::cout << "Add at middle of line\n";
        }
      } else
      {
        stabLine -> append ((theGraph->at(0)->at(bifidx[i]).x() + theGraph->at(0)->at(bifidx[i]-1).x())/2.0, 1.0);
      }
      ++xadded;
      ++yadded;
      switch (biftype[i])
      {
        case BifLP:
          stabLine -> setMarkerShape (QScatterSeries::MarkerShapeCircle);
          break;
        case BifPD:
          stabLine -> setMarkerShape (QScatterSeries::MarkerShapeCircle);
          break;
        case BifNS:
          stabLine -> setMarkerShape (QScatterSeries::MarkerShapeRectangle);
          break;
        default:
          stabLine -> setMarkerShape (QScatterSeries::MarkerShapeRectangle);
          break;
      }
      stabLine -> setColor (plcolor);
    }
//     std::cout << "OUT: stabLine->count() " << stabLine->count() << "bifidx.size() " << bifidx.size() << "\n";
  }
  if (xadded != yadded)
  {
    std::cout << "bad number of X and Y coordinates\n";
    return nullptr;
  }
  // TODO this is definitely leaking memory
  OnePlot *oplret = new OnePlot(theGraph, stabLine, data->getFileName().c_str(), x, y, pt, dim); 
  return oplret;

//   else
//   {
//     if (xadded != 0)
//     {
//       std::vector<std::string> parNames;
//       data->getParNames(parNames);
//       if (x >= XParameter0) XCoordText.push_back(parNames.at(x-XParameter0).c_str());
//       else XCoordText.push_back(XCoordMap[x]);
//       if (y >= YParameter0) YCoordText.push_back(parNames.at(y-YParameter0).c_str());
//       else if (y == YProfile) YCoordText.push_back(QString("X_%1(t)").arg(dim));
//       else YCoordText.push_back(YCoordMap[y]);
// 
//       std::list<PlotItem>::const_iterator it = Graph.end();
//       for (size_t i = 0; i < xadded; ++i){ --it; if (it->x.size() != it->y.size()) std::cout << "bad number of X and Y coordinates\n"; }
//       dataToGraphics(it, Graph.end());
//       labelColor();
//       return true;
//     }
//     else
//     {
//       return false;
//     }
//   }
}

bool DataChart::addPlot(const KNDataFile* data, PlotXVariable x, PlotYVariable y,
  size_t pt, size_t dim)
{
  OnePlot *opl = addPartPlot(data, nullptr, x, y, pt, dim);
  if (opl != nullptr)
  {
    Graph.push_back (opl);
    addToChart (opl);
    return true;
  } else
  {
    return false;
  }
}

void DataChart::addToChart(OnePlot* plot)
{
//   std::cout << "Adding " << plot->file.fileName().toStdString() << "\n";
  foreach( QXYSeries* mem, *(plot->graph) )
  {
//     std::cout << "addSeries (mem)\n";
    chart->addSeries (mem); 
//     std::cout << "attachAxis (hAxis)\n";
    mem -> attachAxis (hAxis);
//     std::cout << "attachAxis (vAxis)\n";
    mem -> attachAxis (vAxis);
//     std::cout << "show ()\n";
    mem -> show();
  }
  if (plot->bifpoints != nullptr)
  {
//     std::cout << "addSeries (mem)\n";
    chart->addSeries (plot->bifpoints);
//     std::cout << "attachAxis (hAxis)\n";
    plot->bifpoints -> attachAxis (hAxis);
//     std::cout << "attachAxis (vAxis)\n";
    plot->bifpoints -> attachAxis (vAxis);
//     std::cout << "show ()\n";
    plot->bifpoints -> show();
  }
}

void DataChart::clearAll()
{
  foreach( OnePlot* lst, Graph )
  {
    auto gr = lst->graph;
    for( int k = 0; k < gr->size(); ++k )
    {
      (*gr)[k] -> detachAxis(vAxis);
      (*gr)[k] -> detachAxis(hAxis);
      chart -> removeSeries ((*gr)[k]);
      delete (*gr)[k];
      (*gr)[k] = nullptr;
    }
    chart -> removeSeries (lst -> bifpoints);
    delete lst -> bifpoints;
    lst -> bifpoints = nullptr;
    delete lst->graph;
    delete lst;
  }
  Graph.clear();
}

void DataChart::clear(size_t n)
{
//   std::cout << "Clear n = " << n << " size = " << Graph.size() << "\n";
  if (n < Graph.size() )
  {
    auto gr = Graph[n]->graph;
    for( int k = 0; k < gr->size(); ++k )
    {
      (*gr)[k] -> detachAxis(vAxis);
      (*gr)[k] -> detachAxis(hAxis);
      chart -> removeSeries ((*gr)[k]);
      delete (*gr)[k];
      (*gr)[k] = nullptr;
    }
    chart -> removeSeries (Graph.at(n)->bifpoints);
    delete Graph[n]->bifpoints;
    Graph[n]->bifpoints = nullptr;
    delete Graph[n] -> graph;
    Graph[n] -> graph = nullptr;
    delete Graph[n];
    Graph[n] = nullptr;
    Graph.removeAt(n);
  }
}

size_t DataChart::nplots()
{
 return Graph.size();
}

QColor DataChart::getColor(size_t n)
{
  if (n < Graph.size() )
    return Graph.at(n)->graph->front()->color();
  else
    return QColor();
}

void DataChart::setColor(size_t n, QColor& Color)
{
  foreach( QXYSeries* mem, *(Graph.at(n)->graph) )
  {
    mem -> setColor (Color);
  }
  Graph.at(n) -> bifpoints -> setColor (Color);
}

void DataChart::updatePlot(const KNDataFile* mat)
{
  QFileInfo finf(mat->getFileName().c_str());
  // list of items to remove
  QList<int> toRemove;
  int gsize = Graph.size();
  for (int k = 0; k < gsize; ++k)
  {
    if (Graph.at(k)->file == finf)
    {
      // updates the plot
      // this always creates a new OnePlot, but not a graph
      int szgraph = Graph.at(k)->graph->size();
//       std::cout << "addPartPlot\n";
      auto gr = addPartPlot(mat, Graph[k], Graph[k]->xvar, Graph[k]->yvar, Graph[k]->pt, Graph[k]->dim);
      // only removes after the new is added
      if (gr->graph != nullptr)
      {
        if (gr->graph != Graph[k]->graph)
        {
            Graph.push_back (gr);
            addToChart (gr);
            clear (k);
        } else
        {
            for (int k = szgraph; k < gr->graph->size(); ++k)
            {
            chart->addSeries (gr->graph->at(k)); 
            gr->graph->at(k) -> attachAxis (hAxis);
            gr->graph->at(k) -> attachAxis (vAxis);
            gr->graph->at(k) -> show();          
            }
            // the stability must be the same, so no updating is necessary. (It migh be empty)
        }
      }
      // No need to delete these
    }
  }

  // now we can remove the items
  foreach( int id, toRemove )
  {
    Graph.removeAt (id);
  }
}
