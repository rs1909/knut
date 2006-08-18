// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "plotdata.h"

#include <QMainWindow>
class QComboBox;
class QSpinBox;

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
