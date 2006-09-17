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
class QLineEdit;
class QComboBox;
class QSpinBox;
class QListWidget;

class plotWindow : public QMainWindow
{
	Q_OBJECT
 public:
	plotWindow( const QString& fname, QWidget *parent = 0 );
	~plotWindow();
	bool isPlotting() { return data != 0; }
 private:
	void openFile( const QString& fname );
	// variables
	const mat4Data *data;
	PlotData plotdata;
	// gui elements
	QLineEdit *matfile;
	QComboBox *xvar;
	QComboBox *yvar;
	std::vector<QString> xvarMap;
	std::vector<QString> yvarMap;
	QSpinBox  *ptlabel;
	QSpinBox  *dim;
	QListWidget *plotSelect;
	std::list<QString> plotList;
	QString   shortFileName;
 private slots:
	void addPlot();
	void clearPlot();
	void open();
	void removePlot( );
};
