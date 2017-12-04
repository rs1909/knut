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

  if (index.row() >= 2) return QVariant();

  if (role == Qt::DisplayRole)
  {
    if (index.row() == 0)
    {
      if (parameters->getPointType() == SolUser)
      {
        // Find the equation NAME based on its index in the list
        return QVariant(parameters->EqnTable.CIndexToTypeName(parameters->EqnTable.TypeToCIndex(parameters->getEqns(index.column()))).c_str());
      }
      else
      {
        // Find the parameter NAME based on its index in the list
        return QVariant(parameters->VarTable.CIndexToTypeName(parameters->VarTable.TypeToCIndex(parameters->getParx(index.column()))).c_str());
      }
    }
    if (index.row() == 1)
    {
      if (parameters->getPointType() == SolUser)
      {
        // Find the variable NAME based on its index in the list
        return QVariant(parameters->VarTable.CIndexToTypeName(parameters->VarTable.TypeToCIndex(parameters->getVars(index.column()))).c_str());
      }
    }
  }
  return QVariant();
}

bool ParamsModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
  if (index.isValid() && role == Qt::EditRole)
  {
    if (index.row() == 0)
    {
      if (parameters->getPointType() == SolUser)
      {
        // set the equation at position index.column() by its index
        parameters->setEqns(index.column(), parameters->EqnTable.CIndexToType(value.toUInt()));
      } else
      {
        // set the extra parameter at position index.column() by its index
        parameters->setParx(index.column(), parameters->VarTable.CIndexToType(value.toUInt()));
      }
    }
    if ((index.row() == 1) && (parameters->getPointType() == SolUser))
    {
      // set the variable at position index.column() by its index
      parameters->setVars(index.column(), parameters->VarTable.CIndexToType(value.toUInt()));
    }
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
      for (size_t i = 0; i < parameters->EqnTable.size(); ++i) editor->addItem(parameters->EqnTable.CIndexToTypeName(i).c_str());
      editor->setCurrentIndex(parameters->EqnTable.TypeToCIndex(parameters->getEqns(index.column())));
    }
    if (index.row() == 1)
    {
      for (size_t i = 0; i < parameters->VarTable.size(); ++i) editor->addItem(parameters->VarTable.CIndexToTypeName(i).c_str());
      editor->setCurrentIndex(parameters->VarTable.TypeToCIndex(parameters->getVars(index.column())));
    }
  }
  else
  {
    if (index.row() == 0)
    {
      for (size_t i = 0; i < parameters->VarTable.size(); ++i) editor->addItem(parameters->VarTable.CIndexToTypeName(i).c_str());
      editor->setCurrentIndex(parameters->VarTable.TypeToCIndex(parameters->getParx(index.column())));
    }
  }

  editor->installEventFilter(const_cast<BoxDelegate*>(this));

  return editor;
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

QVariant SYMModel::data(const QModelIndex &index, int role) const
{
  if (!index.isValid()) return QVariant();

  if (index.column() >= static_cast<int>(parameters->getSymReSize())) return QVariant();
  if (index.row() >= 2) return QVariant();

  if (role == Qt::DisplayRole)
  {
    if (index.row() == 0)
    {
      return QVariant(static_cast<int>(parameters->getSymRe(index.column())));
    }
    if (index.row() == 1)
    {
      return QVariant(static_cast<int>(parameters->getSymIm(index.column())));
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
}

EqnVarTableView::~EqnVarTableView()
{
  delete delegate;
  delete model;
}

void EqnVarTableView::resetSize()
{
  resizeRowsToContents();
//   resizeColumnsToContents(); // Do not resize, it gets too small
  horizontalHeader()->setDefaultAlignment(Qt::AlignLeft);
  horizontalHeader()->setSectionResizeMode(QHeaderView::Fixed);
  verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  
  QSize minSize = QTableView::minimumSizeHint();
  setMinimumSize(minSize);
  setMaximumSize(minSize);
}

QSize EqnVarTableView::minimumSizeHint() const
{
  int tableWidth = 2*frameWidth() + verticalHeader()->frameRect().width() + 
    2*verticalHeader()->frameWidth();
  int tableHeight = 2*frameWidth() + horizontalHeader()->frameRect().height() +
    2*horizontalHeader()->frameWidth();
  for (int i=0; i < model->columnCount(); ++i) tableWidth += columnWidth(i);
  for (int i=0; i < model->rowCount(); ++i) tableHeight += rowHeight(i);
  
  return QSize(tableWidth, tableHeight);
}

SYMTableView::SYMTableView(NConstantsQtGui* params, QWidget* parent_) : QTableView(parent_), parameters(params)
{
  model = new SYMModel(parameters);
  delegate = new SYMDelegate(parameters);
  setItemDelegate(static_cast<QAbstractItemDelegate*>(delegate));
  setModel(model);
  resetSize();
}

SYMTableView::~SYMTableView()
{
  delete delegate;
  delete model;
}

void SYMTableView::resetSize()
{
  resizeRowsToContents();
//   resizeColumnsToContents(); // Do not resize, it gets too small
  horizontalHeader()->setDefaultAlignment(Qt::AlignLeft);
  horizontalHeader()->setSectionResizeMode(QHeaderView::Fixed);
  verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  
  int tableWidth = 2*frameWidth() + verticalHeader()->frameRect().width() + 
    2*verticalHeader()->frameWidth();
  int tableHeight = 2*frameWidth() + horizontalHeader()->frameRect().height() +
    2*horizontalHeader()->frameWidth();
  for (int i=0; i < model->columnCount(); ++i) tableWidth += columnWidth(i);
  for (int i=0; i < model->rowCount(); ++i) tableHeight += rowHeight(i);
  
  setMinimumSize(tableWidth, tableHeight);
  setMaximumSize(tableWidth, tableHeight);
}
