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
    ParamsModel(NConstantsQtGui* params, QObject *parent_ = nullptr)
        : QAbstractTableModel(parent_), parameters(params)
    { }

    int rowCount(const QModelIndex &parent = QModelIndex()) const override
    {
      if (parameters->getPointType() == PtType::SolUser) return 2;
      else return 1;
    }
    int columnCount(const QModelIndex &parent = QModelIndex()) const override
    {
      if (parameters->getPointType() == PtType::SolUser) return static_cast<int>(parameters->getEqnsSize());
      else return static_cast<int>(parameters->getParxSize());
    };

    QVariant data(const QModelIndex &index, int role) const override;
    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

    Qt::ItemFlags flags(const QModelIndex &index) const override
    {
      if (!index.isValid()) return Qt::NoItemFlags;
      else return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
    }

  public slots:
    void dataUpdated()
    {
      beginResetModel();
      endResetModel();
    }
  private:
    NConstantsQtGui* parameters;
};


class BoxDelegate : public QItemDelegate
{
    Q_OBJECT

  public:
    BoxDelegate(NConstantsQtGui* params, QObject *parent = nullptr)
        : QItemDelegate(parent), parameters(params)
    { }

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;

  private:

    NConstantsQtGui* parameters;
};

class EqnVarTableView : public QTableView
{
    Q_OBJECT

  public:
    EqnVarTableView(NConstantsQtGui* params, QWidget* parent_ = nullptr);
    ~EqnVarTableView() override;

  public slots:
    void dataUpdated(const std::vector<Var>& val)
    {
      model->dataUpdated();
      this->resetSize();
      emit sizeChanged(val.size());
    }
    void dataUpdated(const std::vector<Eqn>& val)
    {
      model->dataUpdated();
      this->resetSize();
      emit sizeChanged(val.size());
    }
    void dataUpdated(int typeidx)
    {
      model->dataUpdated();
      this->resetSize();
    }
    // set size of the widget
    void resetSize();
    QSize minimumSizeHint() const override;
    
  signals:
    void sizeChanged(int);
    
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
    SYMModel(NConstantsQtGui* params, QObject *parent_ = nullptr)
        : QAbstractTableModel(parent_), parameters(params)
    { }

    int rowCount(const QModelIndex &parent = QModelIndex()) const override
    {
      return 2;
    }
    int columnCount(const QModelIndex &parent = QModelIndex()) const override
    {
      return static_cast<int>(parameters->getSymReSize());
    };

    QVariant data(const QModelIndex &index, int role) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

    Qt::ItemFlags flags(const QModelIndex &index) const override;
    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;

  public slots:

    void dataUpdated()
    {
      beginResetModel();
      endResetModel();
    }

  private:
    NConstantsQtGui* parameters;
};

class SYMDelegate : public QItemDelegate
{
    Q_OBJECT

  public:
    SYMDelegate(NConstantsQtGui* params, QObject *parent_ = nullptr)
        : QItemDelegate(parent_), parameters(params)
    { }

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    void setEditorData(QWidget *editor, const QModelIndex &index) const override;
    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;

    void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

  private:

    NConstantsQtGui* parameters;

};


class SYMTableView : public QTableView
{
    Q_OBJECT

  public:
    SYMTableView(NConstantsQtGui* params, QWidget* parent_ = nullptr);
    ~SYMTableView() override;

  public slots:

    void dataUpdated(const std::vector<size_t>& val)
    {
      model->dataUpdated();
      this->resetSize();
      emit sizeChanged(val.size());
    }
    void resetSize();
  signals:
    void sizeChanged(int);
  private:

    SYMModel*    model;
    SYMDelegate* delegate;
    NConstantsQtGui*  parameters;
};


#endif
