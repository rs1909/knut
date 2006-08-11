#include <QDialog>
#include <QTextEdit>
#include <string>
#include <QVBoxLayout>
#include <iostream>

class screenDialog : public QDialog
{
	Q_OBJECT
 public:
	screenDialog( QWidget *parent = 0 ) : QDialog(parent)
	{
		display = new QTextEdit;
		QVBoxLayout *layout = new QVBoxLayout;
		layout->addWidget( display );
		layout->setMargin(0);
		this->setLayout( layout );
	}
 public slots:
	void setTitle( const std::string& str ) { display->setDocumentTitle( QString( str.c_str() )); }
	void append( const std::string& str ) { display->append( QString( str.c_str() )); }
	void setText( const std::string& str ) { display->setPlainText( QString( str.c_str() )); }
 private:
	QTextEdit *display;
};
