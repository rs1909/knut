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
  : QThread(parent), KNAbstractContinuation(constants), output(0), changeQueued(false)
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
  output = data;
//  std::cout << "MThread: new data\n";
  emit dataSet(output);
}

KNDataFile& MThread::data()
{
  return *output;
}

const KNDataFile* MThread::dataPointer()
{
  return output;
}

// deleteData -> dataDeleteReq -> <MainWindow> -> dataDeleteAck
void MThread::deleteData()
{
  if (output) emit dataDeleteReq();
}

void MThread::dataUpdated()
{
  changeLock.lock();
  volatile bool queued = changeQueued;
  changeLock.unlock();
  if (!queued)
  {
    changeLock.lock();
    changeQueued = true;
    changeLock.unlock();
//    std::cout << "+req " << queued << "\n";
    emit dataChanged( const_cast<const KNDataFile*>(output) );
  } else
  {
//    std::cout << "0req " << queued << "\n";
  }
}

// SLOTS:

void MThread::dataDeleteAck()
{
//  std::cout << "MThread: delete data\n";
  delete output;
  output = 0;
}

void MThread::dataChangedAck()
{
//  std::cout << "-ack\n";
  changeLock.lock();
  changeQueued = false;
  changeLock.unlock();
}

void MThread::consolePrint(const std::string& str)
{
  std::cout << str;
}
