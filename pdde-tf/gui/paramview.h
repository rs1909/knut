#ifndef PARAMVIEW_H
#define PARAMVIEW_H

#include <QAbstractTableModel>
#include <QItemDelegate>
#include <QVector>
#include <QtGui>
#include <iostream>

#include "pointtype.h"
#include "constqtgui.h"

class QTableView;

//---------------------------------------------------------------
// this is the model for the table
//---------------------------------------------------------------

class ParamsModel : public QAbstractTableModel
{
	Q_OBJECT

 public:
	ParamsModel( NConstantsQtGui* params, QObject *parent_ = 0 )
	 : QAbstractTableModel(parent_), parameters(params) { }

	int rowCount(const QModelIndex &parent = QModelIndex()) const 
	 { if( parameters->getPointType() == SolUser ) return 2; else return 1; }
	int columnCount(const QModelIndex &parent = QModelIndex()) const { return parameters->getNEqns(); };

	QVariant data(const QModelIndex &index, int role) const;
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

	Qt::ItemFlags flags(const QModelIndex &index) const;
	bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
	
 public slots:
	void dataUpdated() { reset(); }
 private:
	NConstantsQtGui* parameters;
};


class BoxDelegate : public QItemDelegate
{
	Q_OBJECT

 public:
	BoxDelegate( NConstantsQtGui* params, QObject *parent = 0 ) 
	 : QItemDelegate(parent), parameters(params) { }

	QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;

	void setEditorData(QWidget *editor, const QModelIndex &index) const;
	void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;

	void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;

 private:

	NConstantsQtGui* parameters;
};

class EqnVarTableView : public QTableView
{
	Q_OBJECT

 public:
	EqnVarTableView( NConstantsQtGui* params, QWidget* parent_ = 0 ) : QTableView( parent_ )
	{
		parameters = params;
		model = new ParamsModel( params );
		delegate = new BoxDelegate( params );
		setItemDelegate(static_cast<QAbstractItemDelegate*>(delegate));
		setModel( model );
		resetSize( );
		connect( parameters, SIGNAL(sysnameChanged(const std::string&)), model, SLOT(dataUpdated()) );
		connect( parameters, SIGNAL(neqnsChanged(int)), model, SLOT(dataUpdated()) );
		connect( parameters, SIGNAL(pointTypeChangedIdx(int)), model, SLOT(dataUpdated()) );
		connect( parameters, SIGNAL(pointTypeChangedIdx(int)), this, SLOT(resetSize()) );
		connect( parameters, SIGNAL(neqnsChanged(int)), this, SLOT(resetSize()) );
	}

 private slots:

	void resetSize( )
	{
		resizeRowToContents(0);
		if( parameters->getPointType() == SolUser ) resizeRowToContents(1);
		horizontalHeader()->setDefaultAlignment( Qt::AlignLeft );
		setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
		setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
		if( parameters->getPointType() == SolUser )
		{
			setMaximumSize( columnWidth(0)+2*frameWidth()+2000,
			                rowHeight(0)+rowHeight(1)+2*frameWidth() +
			                horizontalHeader()->frameRect().height() );
			setMinimumSize( columnWidth(0)+2*frameWidth(),
			                rowHeight(0)+rowHeight(1)+2*frameWidth() +
			                horizontalHeader()->frameRect().height() );
		}else
		{
			setMaximumSize( columnWidth(0)+2*frameWidth()+2000,
			                rowHeight(0)+2*frameWidth() +
			                horizontalHeader()->frameRect().height() );
			setMinimumSize( columnWidth(0)+2*frameWidth(),
			                rowHeight(0)+2*frameWidth() +
			                horizontalHeader()->frameRect().height() );
		}
	}

 private:
	
	NConstantsQtGui*  parameters;
	ParamsModel* model;
	BoxDelegate* delegate;

};

//-----------------------------------------------------
// SYM model
//-----------------------------------------------------


class SYMModel : public QAbstractTableModel
{
	Q_OBJECT

 public:
	SYMModel( NConstantsQtGui* params, QObject *parent_ = 0 )
	 : QAbstractTableModel(parent_), parameters(params) { }

	int rowCount(const QModelIndex &parent = QModelIndex()) const { return 2; }
	int columnCount(const QModelIndex &parent = QModelIndex()) const { return parameters->getNSym(); };

	QVariant data(const QModelIndex &index, int role) const;
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

	Qt::ItemFlags flags(const QModelIndex &index) const;
	bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
	
 public slots:
	
	void dataUpdated(  )
	{
		reset();
	}
	
 private:
	NConstantsQtGui* parameters;
};

class SYMDelegate : public QItemDelegate
{
	Q_OBJECT

 public:
	SYMDelegate( NConstantsQtGui* params, QObject *parent_ = 0 ) 
	 : QItemDelegate(parent_), parameters(params) { }

	QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;

	void setEditorData(QWidget *editor, const QModelIndex &index) const;
	void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;

	void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;

 private:

	NConstantsQtGui* parameters;

};


class SYMTableView : public QTableView
{
	Q_OBJECT

 public:
	SYMTableView( NConstantsQtGui* params, QWidget* parent_ = 0 ) : QTableView( parent_ ), parameters(params)
	{
		model = new SYMModel( parameters );
		delegate = new SYMDelegate( parameters );
		setItemDelegate(static_cast<QAbstractItemDelegate*>(delegate));
		setModel( model );
		resetSize();
		connect( parameters, SIGNAL(nsymChanged(int)), model, SLOT(dataUpdated()) );
		connect( parameters, SIGNAL(nsymChanged(int)), this, SLOT(resetSize()) );
	}

 public slots:
	
	void resetSize()
	{
		resizeRowToContents(0);
		resizeRowToContents(1);
		horizontalHeader()->setDefaultAlignment( Qt::AlignLeft );
		setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
		setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
		setMaximumSize( columnWidth(0)+2*frameWidth()+2000,
								rowHeight(0)+rowHeight(1)+2*frameWidth() +
								horizontalHeader()->frameRect().height() );
		setMinimumSize( columnWidth(0)+2*frameWidth(),
								rowHeight(0)+rowHeight(1)+2*frameWidth() +
								horizontalHeader()->frameRect().height() );
	}

 private:

	SYMModel*    model;
	SYMDelegate* delegate;
	NConstantsQtGui*  parameters;
};


#endif
