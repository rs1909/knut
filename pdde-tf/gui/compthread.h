#include <QThread>
#include "constants.h"

class MThread : public QThread
{
	Q_OBJECT
 public:
	MThread( const NConstants& constants, QObject* parent = 0 ) : QThread(parent)
	   { params = new NConstants( constants ); }
	~MThread()
	   { delete params; }
	void setConstants( const NConstants& constants )
	   { delete params; params = new NConstants( constants ); }
	void run();
 public slots:

 private:
	
	NConstants* params;
	
};
