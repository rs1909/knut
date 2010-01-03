// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "compthread.h"

MThread::MThread(const KNConstants& constants, QObject* parent) 
  : QThread(parent), KNAbstractContinuation(constants)
{
}

MThread::~MThread()
{
}

void MThread::run()
{
  KNAbstractContinuation::run();
}

void MThread::print(std::ostringstream& str)
{
  emit printToScreen(str.str()); str.str("");
}

void MThread::raiseException(const KNException& ex)
{
  emit exceptionOccured(ex);
}

void MThread::setData(KNDataFile* data)
{
  output = QSharedPointer<KNDataFile>(data);
}

KNDataFile& MThread::data()
{
  return *output;
}

const QSharedPointer<KNDataFile>& MThread::dataPointer()
{
  return output;
}

void MThread::deleteData()
{
  output.clear();
}

void MThread::dataUpdated()
{
  emit dataChanged( qSharedPointerConstCast<KNDataFile>(output) );
}

// SLOT:
void MThread::consolePrint(const std::string& str)
{
  std::cout << str;
}
