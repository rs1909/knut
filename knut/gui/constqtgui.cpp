// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2011 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "constqtgui.h"
#include <QSpinBox>
#include <QComboBox>

void NConstantsQtGui::constantChanged(const char* name)
{
  for (unsigned int k=0; k<connectionList.size(); k++)
  {
    if (!strcmp(connectionList [k].name, name))
    {
      std::string type = cnames.findType(name);
      const void * value = cnames.findValue(name);
      QObject* object = connectionList [k].object;
      const char * slotName = connectionList [k].slotName;
      object->blockSignals(true);
      if (type == "int")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, *static_cast<const int *>(value)));
      } else if (type == "size_t")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, static_cast<const int>(*static_cast<const size_t *>(value))));
      } else if (type == "char")
      {
        std::cout << name << "=" << *static_cast<const char *>(value) << "\n";
      } else if (type == "bool")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(bool, *static_cast<const bool *>(value)));
      } else if (type == "double")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(const QString&, QString::number(*static_cast<const double *>(value))));
      } else if (type == "std::string")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(const QString&, QString(static_cast<const std::string *>(value)->c_str())));
      } else if (type == "PtType")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, PtTypeTable.TypeToCIndex(*static_cast<const PtType *>(value))));
      } else if (type == "Eqn")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, EqnTable.TypeToCIndex(*static_cast<const Eqn *>(value))));
      } else if (type == "Var")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, VarTable.TypeToCIndex(*static_cast<const Var *>(value))));
      } else if (type == "BranchSW")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, BranchSWTable.TypeToCIndex(*static_cast<const BranchSW *>(value))));
      } else if (type == "BifType")
      {
        QMetaObject::invokeMethod (object, slotName, Q_ARG(int, BifTypeTable.TypeToCIndex(*static_cast<const BifType *>(value))));
      }else if (type == "vector<Var>")
      {
        if (!strcmp(QSpinBox::staticMetaObject.className(), object->metaObject()->className()))
        {
          QMetaObject::invokeMethod (object, slotName, Q_ARG(int, static_cast<const std::vector<Var>*>(value)->size()) );
        } else
        QMetaObject::invokeMethod (object, slotName, Q_ARG(const std::vector<Var>&, *static_cast<const std::vector<Var>*>(value)));
      } else if (type == "vector<Eqn>")
      {
        if (!strcmp(QSpinBox::staticMetaObject.className(), object->metaObject()->className()))
        {
          QMetaObject::invokeMethod (object, slotName, Q_ARG(int, static_cast<const std::vector<Eqn>*>(value)->size()) );
        } else
          QMetaObject::invokeMethod (object, slotName, Q_ARG(const std::vector<Eqn>&, *static_cast<const std::vector<Eqn>*>(value)));
      } else if (type == "vector<int>")
      {
        if (!strcmp(QSpinBox::staticMetaObject.className(), object->metaObject()->className()))
        {
          QMetaObject::invokeMethod (object, slotName, Q_ARG(int, static_cast<const std::vector<int>*>(value)->size()) );
        } else
          QMetaObject::invokeMethod (object, slotName, Q_ARG(const std::vector<int>&, *static_cast<const std::vector<int>*>(value)));
      } else if (type == "vector<size_t>")
      {
        if (!strcmp(QSpinBox::staticMetaObject.className(), object->metaObject()->className()))
        {
          QMetaObject::invokeMethod (object, slotName, Q_ARG(int, static_cast<const int>(static_cast<const std::vector<size_t>*>(value)->size())) );
        } else
          QMetaObject::invokeMethod (object, slotName, Q_ARG(const std::vector<size_t>&, *static_cast<const std::vector<size_t>*>(value)));
      } else
      {
        std::cout << type << " is not a registered type.\n";
      }
      object->blockSignals(false);
    } 
  }
}

void NConstantsQtGui::registerCallback(const char * type, const char * name, QObject* object, const char * signal, const char * slot)
{
  const tuple_t tp = {type, name, object, signal, slot};
  connectionList.push_back(tp);
  if (signal == 0) return;
  if (!strncmp(type,"int",12))
  {
    connect(object, signal, this, SLOT(slotToInt(int)));
  } else if (!strncmp(type,"size_t",12))
  {
    connect(object, signal, this, SLOT(slotToSizeT(int)));
  } else if (!strncmp(type,"bool",12))
  {
    connect(object, signal, this, SLOT(slotToBool(bool)));
  } else if (!strncmp(type,"double",12))
  {
    connect(object, signal, this, SLOT(slotToDouble(const QString &)));
  } else if (!strncmp(type,"PtType",12))
  {
    connect(object, signal, this, SLOT(slotToPtType(int)));
  } else if (!strncmp(type,"Eqn",12))
  {
    connect(object, signal, this, SLOT(slotToEqn(int)));
  } else if (!strncmp(type,"Var",12))
  {
    connect(object, signal, this, SLOT(slotToVar(int)));
  } else if (!strncmp(type,"BranchSW",12))
  {
    connect(object, signal, this, SLOT(slotToBranchSW(int)));
  } else if (!strncmp(type,"BifType",12))
  {
    connect(object, signal, this, SLOT(slotToBifType(int)));
  } else if (!strncmp(type,"std::string",12))
  {
    connect(object, signal, this, SLOT(slotToString(const QString &)));
  }
}

/// Note that there is a deficiency; 
/// if there are two signals connecting from the same object, they cannot be told apart.
template<typename TP> NConstantsQtGui::callbackPtr<TP> NConstantsQtGui::findConnection(QObject *object)
{
  int found = 0;
  const char *name = 0;
  for (unsigned int k=0; k<connectionList.size(); k++)
  {
    if (connectionList [k].object == object)
    {
      name = connectionList [k].name; 
      found++;
    } 
  }
  if (found == 1) return callbackPtr<TP>(cnames.findFun(name));
  return 0;
}

void NConstantsQtGui::slotToInt(int val)
{
  callbackPtr<int> ptr = findConnection<int>(sender());
  if (ptr.value) ptr.value(*this, val);
}

void NConstantsQtGui::slotToSizeT(int val)
{
  callbackPtr<size_t> ptr = findConnection<size_t>(sender());
  if (ptr.value) ptr.value(*this, static_cast<size_t>(val));
}

void NConstantsQtGui::slotToBool(bool val)
{
  callbackPtr<bool> ptr = findConnection<bool>(sender());
  if (ptr.value) ptr.value(*this, val);
}

void NConstantsQtGui::slotToDouble(const QString& val)
{
  callbackPtr<double> ptr = findConnection<double>(sender());
  if (ptr.value) ptr.value(*this, val.toDouble());
}

void NConstantsQtGui::slotToPtType(int val)
{
  callbackPtr<PtType> ptr = findConnection<PtType>(sender());
  if (ptr.value) ptr.value(*this, PtTypeTable.CIndexToType(val));  
}

void NConstantsQtGui::slotToEqn(int val)
{
  callbackPtr<Eqn> ptr = findConnection<Eqn>(sender());
  if (ptr.value) ptr.value(*this, EqnTable.CIndexToType(val)); 
}

void NConstantsQtGui::slotToBranchSW(int val)
{
  callbackPtr<BranchSW> ptr = findConnection<BranchSW>(sender());
  if (ptr.value) ptr.value(*this, BranchSWTable.CIndexToType(val)); 
}

void NConstantsQtGui::slotToBifType(int val)
{
  callbackPtr<BifType> ptr = findConnection<BifType>(sender());
  if (ptr.value) ptr.value(*this, BifTypeTable.CIndexToType(val)); 
}

void NConstantsQtGui::slotToVar(int val)
{
  callbackPtr<Var> ptr = findConnection<Var>(sender());
  if (ptr.value) ptr.value(*this, VarTable.CIndexToType(val));
//   if ( VarToCIndex(CIndexToVar(val)) != val ) std::cout << "wrong conversion\n";
}

void NConstantsQtGui::slotToString(const QString& val)
{
  callbackPtr<std::string> ptr = findConnection<std::string>(sender());
  if (ptr.value) ptr.value(*this, val.toStdString());
}
