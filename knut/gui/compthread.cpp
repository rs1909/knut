// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "compthread.h"
#include <QCoreApplication>

MThread::MThread(const KNConstants& constants, QObject* parent) 
  : QObject(parent), KNAbstractContinuation(constants), output(0), changeQueued(false)
{
}

MThread::~MThread()
{
}

void MThread::process()
{
  KNAbstractContinuation::run();
  emit finished();
}

void MThread::printStream()
{
  emit printToScreen(screenout.str()); screenout.str("");
}

void MThread::raiseException(const KNException& ex)
{
  emit exceptionOccured(ex);
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
//   std::cout << "MThread::deleteData in\n";
  if (output) emit dataDeleteReq();
  // wait for the slot to be activated, output turns zero
  do {
    QCoreApplication::processEvents(QEventLoop::WaitForMoreEvents);
  } while (output != 0);
//   std::cout << "MThread::deleteData out\n";
}

void MThread::dataDeleteAck()
{
  output = 0;
}

void MThread::dataUpdated()
{
  QCoreApplication::processEvents ();
  if (!changeQueued)
  {
    changeQueued = true;
//     std::cout << "+req " << changeQueued << "\n";
    emit dataChanged( const_cast<const KNDataFile*>(output) );
  } else
  {
//     std::cout << "0req " << changeQueued << "\n";
  }
}

// SLOTS:

// this is called to notify that the result is plotted already
void MThread::dataChangedAck()
{
//   std::cout << "-ack\n";
  changeQueued = false;
}

void MThread::consolePrint(const std::string& str)
{
  std::cout << str;
}

void MThread::dataCreated(KNDataFile* dataFile)
{
//   std::cout << "MThread::dataCreated\n";
  output = dataFile;
}

void MThread::createDataLC (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms)
{
//   std::cout << "MThread::createDataLC in\n";
  output = 0;
  emit createDataRequestLC (fileName, ndim, npar, prms);
  // wait for the slot to be activated, output turns nonzero
  do {
    QCoreApplication::processEvents(QEventLoop::WaitForMoreEvents);
//     std::cout << "MThread::createDataLC loop " << output << "\n";
  } while (output == 0);
//   std::cout << "MThread::createDataLC out\n";
}

void MThread::createDataTR (const std::string& fileName, size_t ndim, size_t npar, KNConstants* prms)
{
  output = 0;
  emit createDataRequestTR (fileName, ndim, npar, prms);
  do {
    QCoreApplication::processEvents(QEventLoop::WaitForMoreEvents);
  } while (output == 0);
}
