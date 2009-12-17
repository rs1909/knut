// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "compthread.h"

MThread::MThread(const NConstants& constants, QObject* parent) 
  : QThread(parent), BaseComp(constants)
{
}

MThread::~MThread()
{
}

void MThread::run()
{
  BaseComp::run();
}

void MThread::print(std::ostringstream& str)
{
  emit printToScreen(str.str()); str.str("");
}

void MThread::raiseException(const knutException& ex)
{
  emit exceptionOccured(ex);
}

void MThread::setData(mat4Data* data)
{
  output = QSharedPointer<mat4Data>(data);
}

mat4Data& MThread::data()
{
  return *output;
}

const QSharedPointer<mat4Data>& MThread::dataPointer()
{
  return output;
}

void MThread::deleteData()
{
  output.clear();
}

void MThread::dataUpdated()
{
  emit dataChanged( qSharedPointerConstCast<mat4Data>(output) );
}

// SLOT:
void MThread::consolePrint(const std::string& str)
{
  std::cout << str;
}
