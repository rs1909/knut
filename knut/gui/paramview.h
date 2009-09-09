// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef PARAMVIEW_H
#define PARAMVIEW_H

#include "pointtype.h"
#include "constqtgui.h"

#include <QAbstractTableModel>
#include <QTableView>
#include <QItemDelegate>
#include <QVector>

class ParamsModel : public QAbstractTableModel
{
    Q_OBJECT

  public:
    ParamsModel(NConstantsQtGui* params, QObject *parent_ = 0)
        : QAbstractTableModel(parent_), parameters(params)
    { }

    int rowCount(const QModelIndex &/*parent = QModelIndex()*/) const
    {
      if (parameters->getPointType() == SolUser) return 2;
      else return 1;
    }
    int columnCount(const QModelIndex &/*parent = QModelIndex()*/) const
    {
      if (parameters->getPointType() == SolUser) return parameters->getEqnsNumSize();
      else return parameters->getParxNumSize();
    };

    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

    Qt::ItemFlags flags(const QModelIndex &index) const;
    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);

  public slots:
    void dataUpdated()
    {
      reset();
    }
  private:
    NConstantsQtGui* parameters;
};


class BoxDelegate : public QItemDelegate
{
    Q_OBJECT

  public:
    BoxDelegate(NConstantsQtGui* params, QObject *parent = 0)
        : QItemDelegate(parent), parameters(params)
    { }

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
    EqnVarTableView(NConstantsQtGui* params, QWidget* parent_ = 0);
    ~EqnVarTableView();

  private slots:

    void setConstant(const char* name)
    {
      bool call = false;
      if (!strcmp(name,"pointType")) call = true;
      else if (!strcmp(name,"parxNum")) call = true;
      else if (!strcmp(name,"eqnsNum")) call = true;
      else if (!strcmp(name,"varsNum")) call = true;
      else if (!strcmp(name,"translationMaps")) call = true;
      if (call)
      {
        model->dataUpdated();
        this->resetSize();
      }
    }
    
    void resetSize();

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
    SYMModel(NConstantsQtGui* params, QObject *parent_ = 0)
        : QAbstractTableModel(parent_), parameters(params)
    { }

    int rowCount(const QModelIndex &/*parent = QModelIndex()*/) const
    {
      return 2;
    }
    int columnCount(const QModelIndex &/*parent = QModelIndex()*/) const
    {
      return parameters->getSymReSize();
    };

    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

    Qt::ItemFlags flags(const QModelIndex &index) const;
    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);

  public slots:

    void dataUpdated()
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
    SYMDelegate(NConstantsQtGui* params, QObject *parent_ = 0)
        : QItemDelegate(parent_), parameters(params)
    { }

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
    SYMTableView(NConstantsQtGui* params, QWidget* parent_ = 0);
    ~SYMTableView();

  private slots:

    void setConstant(const char* name)
    {
      bool call = false;
      if (!strcmp(name,"symIm")) call = true;
      if (call)
      {
        model->dataUpdated();
        this->resetSize();
      }
    }

    void resetSize();

  private:

    SYMModel*    model;
    SYMDelegate* delegate;
    NConstantsQtGui*  parameters;
};


#endif