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

MThread::MThread(QObject* parent)
  : QObject(parent), output(0), changeQueued(false)
{
}

MThread::~MThread()
{
}

void MThread::setConstants (const KNConstants& prms)
{
  params = new KNConstants(prms);
  try {
    sys = new KNSystem (prms.getSysName ());
  }
  catch (KNException& ex)
  {
    delete params;
    params = nullptr;
    sys = nullptr;
    emit exceptionOccured (ex);
  }
}

const QEventLoop::ProcessEventsFlags MThread::waitFlag = QEventLoop::WaitForMoreEvents | QEventLoop::ExcludeUserInputEvents | QEventLoop::ExcludeSocketNotifiers;

void MThread::process()
{
  if ((sys != nullptr) && (params != nullptr))
  {
    KNDataFile *inputData = nullptr;
    if (params->getLabel() != 0)
    {
      inputData = new KNDataFile (params->getInputFile());
    }
    try {
      KNAbstractContinuation::run(sys, params, inputData);
    }
    catch (KNException& ex)
    {
      emit exceptionOccured (ex);
      delete inputData;
    }
    delete inputData;
  }
  emit finished();
}

void MThread::printStream()
{
//  std::cout << "MThread::printStream in\n";
  emit printToScreen(screenout.str()); screenout.str("");
//  std::cout << "MThread::printStream out\n";
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
//  std::cout << "MThread::deleteData in\n";
  if (output) emit dataDeleteReq();
  // wait for the slot to be activated, output turns zero
  int it = 0;
  do {
    if (it == 0) QCoreApplication::processEvents(waitFlag);
    else QCoreApplication::processEvents(waitFlag, 100);
    ++it;
//    std::cout << "MThread::deleteData loop " << output << " " << it << "\n";
  } while (output != 0);
//  std::cout << "MThread::deleteData out\n";
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
//    std::cout << "+req " << changeQueued << "\n";
    emit dataChanged( const_cast<const KNDataFile*>(output) );
  } else
  {
//     std::cout << "0req " << changeQueued << "\n";
  }
}

// SLOTS:

void MThread::stopReq(bool flag)
{
  KNAbstractContinuation::setStopFlag(flag);
}

// this is called to notify that the result is plotted already
void MThread::dataChangedAck()
{
//  std::cout << "-ack\n";
  changeQueued = false;
}

void MThread::consolePrint(const std::string& str)
{
  std::cout << str;
}

void MThread::dataCreated(KNDataFile* dataFile)
{
//  std::cout << "MThread::dataCreated\n";
  output = dataFile;
}

void MThread::createData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms)
{
//  std::cout << "MThread::createDataLC in\n";
  output = 0;
  emit createDataRequest (fileName, t, ndim, npar, prms);
  // wait for the slot to be activated, output turns nonzero
  int it = 0;
  do {
    if (it == 0) QCoreApplication::processEvents(waitFlag);
    else QCoreApplication::processEvents(waitFlag, 100);
    ++it;
//    std::cout << "MThread::createDataLC loop " << output << " " << it << "\n";
  } while (output == 0);
//  std::cout << "MThread::createDataLC out\n";
}
