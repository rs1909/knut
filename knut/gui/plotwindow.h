// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "plotdata.h"
#include "pointtype.h"

#include <QSplitter>
#include <QFileInfo>
#include <QCloseEvent>

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

class plotWindow : public QSplitter
{
    Q_OBJECT
  public:
    plotWindow(const QString& filename, QWidget *parent = 0);
    plotWindow(QWidget *parent = 0);
    ~plotWindow();
    void init(Var cp);
    bool isDataSet()
    {
      return data != 0;
    }
    bool isThisData(const KNDataFile* dataFile)
    {
      QFileInfo f1(QString::fromStdString(data->getFileName()));
      QFileInfo f2(QString::fromStdString(dataFile->getFileName()));
//      std::cout << f1.canonicalFilePath().toStdString() << " v.s.\n"
//                << f2.canonicalFilePath().toStdString() << "\n" << (int)(f1 == f2) << "\n";
      return (f1.canonicalFilePath() == f2.canonicalFilePath());
    }
  protected:
  	void closeEvent(QCloseEvent *event)
 	{
 	  emit windowClosed();
      event->accept();
    }
  private:
    void setupPlotWindow();
    // variables
    const KNDataFile* data;
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
    QFileInfo dataFileInfo;
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
    void setData(const KNDataFile* dataFile);
    // is called when the computing thread made a step
    void updatePlot(const KNDataFile* dataFile);
  signals:
    void windowClosed();
    // gets emitted when a new file is opened
    // this makes the system open the file and the will called
    // the setData slot to add the plot
    void openFile(const QString& filename);
    // if the file is no longer in use, the signal is emitted
    void closeFile(const KNDataFile* dataFile);
    // called when plotting is finished
    void updated();
};
