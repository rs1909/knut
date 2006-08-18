// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "constqtgui.h"
#include <string>
#include <QThread>

class MThread : public QThread
{
	Q_OBJECT
 public:
	MThread( const NConstantsQtGui& constants, QObject* parent = 0 ) : QThread(parent), stopFlag(false)
	   { params = new NConstantsQtGui( constants, this ); }
	~MThread()
	   { delete params; }
	void setConstants( const NConstantsQtGui& constants )
	   { delete params; params = new NConstantsQtGui( constants, this ); }
	void run();
	void setStopFlag( bool flag ) { stopFlag = flag; }
 public slots:

 signals:
	void exceptionOccured( const pddeException& ex );
	void printToScreen( const std::string& str );
 private:
	
	NConstantsQtGui* params;
	bool             stopFlag;

};
