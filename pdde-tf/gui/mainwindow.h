#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QCheckBox>

#include "system.h"
#include "paramview.h"
#include "compthread.h"
#include "screendialog.h"

class QAction;
class QMenu;
class QTextEdit;
class QVBoxLayout;
class QLineEdit;
class QSpinBox;
class QGridLayout;
class QComboBox;
class QPushButton;
class QLabel;
class mat4Data;
class plotWindow;

class EqnVarTableView;


class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();

public slots:

	void setOutputFileText( const std::string& st )
	{
		outputFile->blockSignals(true);
		int cpos = outputFile->cursorPosition();
		outputFile->setText(st.c_str());
		outputFile->setCursorPosition(cpos);
		outputFile->blockSignals(false);
	}
	void setInputFileText( const std::string& st )
	{
		inputFile->blockSignals(true);
		int cpos = inputFile->cursorPosition();
		inputFile->setText(st.c_str());
		inputFile->setCursorPosition(cpos);
		inputFile->blockSignals(false);
	}
	void setSysNameText( const std::string& st )
	{
		sysname->blockSignals(true);
		int cpos = sysname->cursorPosition();
		sysname->setText(st.c_str());
		sysname->setCursorPosition(cpos);
		sysname->blockSignals(false);
	}
	void setLabel( int i ) { label->blockSignals(true); label->setValue(i); label->blockSignals(false); }
	void setPointTypeIdx( int i ) { pttype->blockSignals(true); pttype->setCurrentIndex(i); pttype->blockSignals(false); }
	void setCpIdx( int i ) { cp->blockSignals(true); cp->setCurrentIndex(i); cp->blockSignals(false); }
	void setBranchSWIdx( int i ) { branchsw->blockSignals(true); branchsw->setCurrentIndex(i); branchsw->blockSignals(false); }
	void setNEqns( int i ) { eqns->blockSignals(true); eqns->setValue(i); eqns->blockSignals(false); }
	void setNInt( int i ) { nint->blockSignals(true); nint->setValue(i); nint->blockSignals(false); }
	void setNDeg( int i ) { ndeg->blockSignals(true); ndeg->setValue(i); ndeg->blockSignals(false); }
	void setNMul( int i ) { nmul->blockSignals(true); nmul->setValue(i); nmul->blockSignals(false); }
	void setStab( bool b )
	{
		stab->blockSignals(true); 
		if(b) stab->setCheckState(Qt::Checked); else stab->setCheckState(Qt::Unchecked);
		stab->blockSignals(false);
	}
	void setNMat( int i ) { nmat->blockSignals(true); nmat->setValue(i); nmat->blockSignals(false); }
	void setNInt1( int i ) { nint1->blockSignals(true); nint1->setValue(i); nint1->blockSignals(false); }
	void setNInt2( int i ) { nint2->blockSignals(true); nint2->setValue(i); nint2->blockSignals(false); }
	void setNDeg1( int i ) { ndeg1->blockSignals(true); ndeg1->setValue(i); ndeg1->blockSignals(false); }
	void setNDeg2( int i ) { ndeg2->blockSignals(true); ndeg2->setValue(i); ndeg2->blockSignals(false); }
	void setSteps( int i ) { steps->blockSignals(true); steps->setValue(i); steps->blockSignals(false); }
	void setNItC( int i ) { nitC->blockSignals(true); nitC->setValue(i); nitC->blockSignals(false); }
	void setNItR( int i ) { nitR->blockSignals(true); nitR->setValue(i); nitR->blockSignals(false); }
	void setNItS( int i ) { nitS->blockSignals(true); nitS->setValue(i); nitS->blockSignals(false); }
	void setNSym( int i ) { nsym->blockSignals(true); nsym->setValue(i); nsym->blockSignals(false); }

	void externalException( const pddeException& ex )
	{
		QMessageBox::critical( this, "MainWindow::externalException()", QString( "%1:%2 %3" ).arg(ex.file.c_str()).arg(ex.line).arg(ex.message.message.c_str()), QMessageBox::Ok, 0, 0 );
	}
	
protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void run();
	void stop();
	void newFile();
	void open();
	bool save();
	bool saveAs();
	void about();
	
	void setSysName();
	void setInputFile();
	void setOutputFile();
	void setPointType();
	void setupCp()
	{
		cp->blockSignals(true);
		cp->clear();
		cp->blockSignals(false);
		for( int i = 0; i < parameters.cpSize(); ++i ) cp->addItem( parameters.cpString( i ).c_str() );
	}
	void inputPlot();
	void inputPlotDestroyed();
	void outputPlot();
	void outputPlotDestroyed();
	void terminalView();
	void terminalViewDestroyed();
	void terminalTextAppend(const std::string& str) { terminalText.append(QString(str.c_str()).append("\n")); }
	void compileSystem();

private:
	inline bool inputAssert( std::istream& is );

	void createActions();
	void createMenus();
	void createToolBars();
	void createStatusBar();
	void readSettings();
	void writeSettings();
	bool maybeSave();
	void loadFile(const QString &fileName);
	bool saveFile(const QString &fileName);
	void setCurrentFile(const QString &fileName);
	QString strippedName(const QString &fullFileName);

	// all the parameters
	NConstantsQtGui parameters;
	MThread      compThread;
	QWidget     *paramsWidget;
	QGridLayout *paramsGrid;
	QLabel      *eqnsLabel;
	
	//the textual output
	QString      terminalText;
	
	// these contain the constants
	QLineEdit *inputFile;
	mat4Data  *inputData;
	plotWindow *inputPlotWindow;
	QLineEdit *outputFile;
	mat4Data  *outputData;
	plotWindow *outputPlotWindow;
	screenDialog* terminalDialog;
	QLineEdit *sysname;
	QSpinBox  *label;
	QComboBox *pttype;
	QComboBox *cp;
	QSpinBox  *eqns;
	QComboBox *branchsw;

	QSpinBox  *nint;
	QSpinBox  *ndeg;
	QSpinBox  *nmul;
	QCheckBox *stab;
	QSpinBox  *nmat;
	QSpinBox  *nint1;
	QSpinBox  *nint2;
	QSpinBox  *ndeg1;
	QSpinBox  *ndeg2;
	QSpinBox  *steps;
	QLineEdit *cpMin;
	QLineEdit *cpMax;
	QLineEdit *ds;
	QLineEdit *dsMin;
	QLineEdit *dsMax;
	QLineEdit *dsStart;
	QLineEdit *epsC;
	QLineEdit *epsR;
	QLineEdit *epsS;
	QSpinBox  *nitC;
	QSpinBox  *nitR;
	QSpinBox  *nitS;
	SYMTableView *sym;
	QSpinBox  *nsym;

	EqnVarTableView* table;

	QString  curFile;
	QProcess *compilerProcess;

	QMenu    *fileMenu;
	QMenu    *helpMenu;
	QToolBar *fileToolBar;
	QAction *runAct;
	QAction *stopAct;
	QAction *terminalAct;
	QAction *openAct;
	QAction *saveAct;
	QAction *saveAsAct;
	QAction *exitAct;
	QAction *aboutAct;
	QAction *aboutQtAct;
};

#endif
