#include "plotdata.h"
#include <QtGui>

class plotWindow : public QMainWindow
{
	Q_OBJECT
 public:
	plotWindow( const mat4Data* d, QWidget *parent = 0 );
	~plotWindow() {}
 private:
	const mat4Data *data;
	PlotData plotdata;
	// gui elements
	QComboBox *xvar;
	QComboBox *yvar;
	QSpinBox  *ptlabel;
	QSpinBox  *dim;
 private slots:
	void addPlot();
	void clearPlot();
};
