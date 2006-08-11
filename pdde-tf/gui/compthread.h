#include <QThread>
#include "constants.h"

class MThread : public QThread
{
	Q_OBJECT
 public:
	MThread( const NConstants& constants, QObject* parent = 0 ) : QThread(parent), stopFlag(false)
	   { params = new NConstants( constants ); }
	~MThread()
	   { delete params; }
	void setConstants( const NConstants& constants )
	   { delete params; params = new NConstants( constants ); }
	void run();
	void setStopFlag( bool flag ) { stopFlag = flag; }
 public slots:

 signals:
	void exceptionOccured( const pddeException& ex );
	void printToScreen( const std::string& str );
 private:
	
	NConstants* params;
	bool        stopFlag;

};
