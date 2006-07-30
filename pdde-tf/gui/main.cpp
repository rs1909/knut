#include <QApplication>

#include "mainwindow.h"
#include "point.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	MainWindow mainWin;
	mmappedPointData( "dummy.data", 112, 2, 4, 100, 4, 12 );
	mainWin.show();
	return app.exec();
}

