// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "paramview.h"
#include <QComboBox>
#include <QSpinBox>
#include <QHeaderView>

QVariant ParamsModel::data(const QModelIndex &index, int role) const
{
  if (!index.isValid()) return QVariant();

  if (index.column() >= parameters->getNEqns()) return QVariant();
  if (index.row() >= 2) return QVariant();

  if (role == Qt::DisplayRole)
  {
    if (index.row() == 0)
    {
      if (parameters->getPointType() == SolUser)
        return QVariant(parameters->findEqnsString(parameters->getEqns(index.column())).c_str());
      else return QVariant(parameters->findParXString(parameters->getParX(index.column())).c_str());
    }
    if (index.row() == 1)
    {
      if (parameters->getPointType() == SolUser)
        return QVariant(parameters->findVarsString(parameters->getVars(index.column())).c_str());
    }
  }
  return QVariant();
}

Qt::ItemFlags ParamsModel::flags(const QModelIndex &index) const
{
  if (!index.isValid())
    return Qt::ItemIsEnabled;

  return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}

bool ParamsModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
  if (index.isValid() && role == Qt::EditRole)
  {
    if (index.row() == 0)
    {
      if (parameters->getPointType() == SolUser) parameters->setEqnsIdx(index.column(), value.toInt());
      else parameters->setParXIdx(index.column(), value.toInt());
    }
    if ((index.row() == 1) && (parameters->getPointType() == SolUser)) parameters->setVarsIdx(index.column(), value.toInt());
    emit dataChanged(index, index);
    return true;
  }
  return false;
}

QVariant ParamsModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (role == Qt::DisplayRole)
  {
    if (parameters->getPointType() == SolUser)
    {
      if (section == 0 && orientation == Qt::Vertical) return QVariant("Eqn");
      if (section == 1 && orientation == Qt::Vertical) return QVariant("Var");
      else return QVariant(section);
    }
    else
    {
      if (section == 0 && orientation == Qt::Vertical) return QVariant("PARX");
      else return QVariant(section);
    }
  }
  return QVariant();
}

//---------------------------
// BoxDelegate
//-----------------------------

QWidget *BoxDelegate::createEditor(QWidget *parent,
                                   const QStyleOptionViewItem &/* option */,
                                   const QModelIndex& index) const
{
  QComboBox *editor = new QComboBox(parent);

  if (parameters->getPointType() == SolUser)
  {
    if (index.row() == 0)
    {
      for (int i = 0; i < parameters->eqnsSize(); ++i) editor->addItem(parameters->eqnsString(i).c_str());
    }
    if (index.row() == 1)
    {
      for (int i = 0; i < parameters->varsSize(); ++i) editor->addItem(parameters->varsString(i).c_str());
    }
  }
  else
  {
    if (index.row() == 0)
    {
      for (int i = 0; i < parameters->parxSize(); ++i) editor->addItem(parameters->parxString(i).c_str());
    }
  }

  editor->installEventFilter(const_cast<BoxDelegate*>(this));

  return editor;
}

void BoxDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
  int value = index.model()->data(index, Qt::DisplayRole).toInt();

  QComboBox *cbox = static_cast<QComboBox*>(editor);
  cbox->setCurrentIndex(value);
}

void BoxDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
  QComboBox *cbox = static_cast<QComboBox*>(editor);

  if (parameters->getPointType() == SolUser)
  {
    if (index.row() == 0)  model->setData(index, cbox->currentIndex());
    if (index.row() == 1)  model->setData(index, cbox->currentIndex());
  }
  else
  {
    if (index.row() == 0)  model->setData(index, cbox->currentIndex());
  }
}

void BoxDelegate::updateEditorGeometry(QWidget *editor,
                                       const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
  editor->setGeometry(option.rect);
}


QVariant SYMModel::data(const QModelIndex &index, int role) const
{
  if (!index.isValid()) return QVariant();

  if (index.column() >= parameters->getNSym()) return QVariant();
  if (index.row() >= 2) return QVariant();

  if (role == Qt::DisplayRole)
  {
    if (index.row() == 0)
    {
      return QVariant(parameters->getSymRe(index.column()));
    }
    if (index.row() == 1)
    {
      return QVariant(parameters->getSymIm(index.column()));
    }
  }
  return QVariant();
}

Qt::ItemFlags SYMModel::flags(const QModelIndex &index) const
{
  if (!index.isValid())
    return Qt::ItemIsEnabled;

  return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}

bool SYMModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
  if (index.isValid() && role == Qt::EditRole)
  {
    if (index.row() == 0) parameters->setSymRe(index.column(), value.toInt());
    if (index.row() == 1) parameters->setSymIm(index.column(), value.toInt());
    emit dataChanged(index, index);
    return true;
  }
  return false;
}

QVariant SYMModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (role == Qt::DisplayRole)
  {
    if (section == 0 && orientation == Qt::Vertical) return QVariant("Re");
    if (section == 1 && orientation == Qt::Vertical) return QVariant("Im");
    else return QVariant(section);
  }
  return QVariant();
}

//---------------------------
// SYMDelegate
//-----------------------------

QWidget *SYMDelegate::createEditor(QWidget *parent,
                                   const QStyleOptionViewItem & /*option*/,
                                   const QModelIndex& /*index*/) const
{
  QSpinBox *editor = new QSpinBox(parent);
  editor -> setRange(0, parameters->getNDim() - 1);
  editor -> installEventFilter(const_cast<SYMDelegate*>(this));

  return editor;
}

void SYMDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
  int value = index.model()->data(index, Qt::DisplayRole).toInt();

  QSpinBox *sbox = static_cast<QSpinBox*>(editor);
  sbox->setValue(value);
}

void SYMDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
  QSpinBox *sbox = static_cast<QSpinBox*>(editor);

  if (index.row() < 2)
  {
    int value = sbox->value();
    model->setData(index, value);
  }
}

void SYMDelegate::updateEditorGeometry(QWidget *editor,
                                       const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
  editor->setGeometry(option.rect);
}

EqnVarTableView::EqnVarTableView(NConstantsQtGui* params, QWidget* parent_) : QTableView(parent_)
{
  parameters = params;
  model = new ParamsModel(params);
  delegate = new BoxDelegate(params);
  setItemDelegate(static_cast<QAbstractItemDelegate*>(delegate));
  setModel(model);
  resetSize();
  connect(parameters, SIGNAL(sysnameChanged(const std::string&)), model, SLOT(dataUpdated()));
  connect(parameters, SIGNAL(neqnsChanged(int)), model, SLOT(dataUpdated()));
  connect(parameters, SIGNAL(pointTypeChangedIdx(int)), model, SLOT(dataUpdated()));
  connect(parameters, SIGNAL(pointTypeChangedIdx(int)), this, SLOT(resetSize()));
  connect(parameters, SIGNAL(neqnsChanged(int)), this, SLOT(resetSize()));
}

void EqnVarTableView::resetSize()
{
  resizeRowToContents(0);
  if (parameters->getPointType() == SolUser) resizeRowToContents(1);
  horizontalHeader()->setDefaultAlignment(Qt::AlignLeft);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  if (parameters->getPointType() == SolUser)
  {
    setMaximumSize(columnWidth(0) + 2*frameWidth() + 2000,
                   rowHeight(0) + rowHeight(1) + 2*frameWidth() +
                   horizontalHeader()->frameRect().height());
    setMinimumSize(columnWidth(0) + 2*frameWidth(),
                   rowHeight(0) + rowHeight(1) + 2*frameWidth() +
                   horizontalHeader()->frameRect().height());
  }
  else
  {
    setMaximumSize(columnWidth(0) + 2*frameWidth() + 2000,
                   rowHeight(0) + 2*frameWidth() +
                   horizontalHeader()->frameRect().height());
    setMinimumSize(columnWidth(0) + 2*frameWidth(),
                   rowHeight(0) + 2*frameWidth() +
                   horizontalHeader()->frameRect().height());
  }
}

SYMTableView::SYMTableView(NConstantsQtGui* params, QWidget* parent_) : QTableView(parent_), parameters(params)
{
  model = new SYMModel(parameters);
  delegate = new SYMDelegate(parameters);
  setItemDelegate(static_cast<QAbstractItemDelegate*>(delegate));
  setModel(model);
  resetSize();
  connect(parameters, SIGNAL(nsymChanged(int)), model, SLOT(dataUpdated()));
  connect(parameters, SIGNAL(nsymChanged(int)), this, SLOT(resetSize()));
}

void SYMTableView::resetSize()
{
  resizeRowToContents(0);
  resizeRowToContents(1);
  horizontalHeader()->setDefaultAlignment(Qt::AlignLeft);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setMaximumSize(columnWidth(0) + 2*frameWidth() + 2000,
                 rowHeight(0) + rowHeight(1) + 2*frameWidth() +
                 horizontalHeader()->frameRect().height());
  setMinimumSize(columnWidth(0) + 2*frameWidth(),
                 rowHeight(0) + rowHeight(1) + 2*frameWidth() +
                 horizontalHeader()->frameRect().height());
}
