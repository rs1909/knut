// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "plotdata.h"

#include <QMainWindow>
#include <QSharedPointer>
class QLineEdit;
class QComboBox;
class QSpinBox;
class QListWidget;

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

class plotWindow : public QMainWindow
{
    Q_OBJECT
  public:
    plotWindow(const QSharedPointer<const mat4Data>& mat, QWidget *parent = 0);
    plotWindow(QWidget *parent = 0);
    ~plotWindow();
    bool isDataSet()
    {
      return data != 0;
    }
  private:
    void setupPlotWindow();
    // variables
    QSharedPointer<const mat4Data> data;
    PlotData plotdata;
    QGraphicsView *plot;
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
  private slots:
    void addPlot();
    void clearPlot();
    void open();
    void removePlot();
    void colorizePlot();
    void print();
    void exportSvg();
  public slots:
    // sets the data file
    void setData(const QSharedPointer<const mat4Data>& mat);
    // is called when the computing thread made a step
    void updatePlot(const QSharedPointer<const mat4Data>& data);
  signals:
    // gets emitted when a new file is opened
    // this makes the system open the file and the will called
    // the setData slot to add the plot
    void openFile(const QString& filename);
};
