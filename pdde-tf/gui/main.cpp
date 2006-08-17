#include <QApplication>
#include <string>

#include "mainwindow.h"
#include "point.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	MainWindow mainWin( app.applicationDirPath() );
	qRegisterMetaType<pddeException>("pddeException");
	qRegisterMetaType<std::string>("std::string");
	mainWin.show();
	return app.exec();
}

