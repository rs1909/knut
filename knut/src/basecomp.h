#ifndef BASECOMP_H
#define BASECOMP_H

#include "constants.h"
#include "mat4data.h"
#include <iostream>

class KNAbstractContinuation
{
  public:
    KNAbstractContinuation(const KNConstants& constants) : stopFlag(false)
    {
      params = new KNConstants(constants);
    }
    virtual ~KNAbstractContinuation()
    {
      delete params;
    }
    void setConstants(const KNConstants& constants)
    {
      delete params;
      params = new KNConstants(constants);
    }
    void run(const char* branchFile);
    void run() { run(0); }
    void setStopFlag(bool flag)
    {
      stopFlag = flag;
    }
    virtual std::ostream& outStream() { return screenout; };
    virtual void printStream() = 0;
    virtual void clearLastLine() = 0;
    virtual void storeCursor() = 0;
    virtual void raiseException(const KNException& ex) = 0;
    virtual void setData(KNDataFile* data) = 0;
    virtual KNDataFile& data() = 0;
    virtual void deleteData() = 0;
    virtual void dataUpdated() = 0;
  protected:
    KNConstants* params;
    bool        stopFlag;
    // this is used within the thread
    std::ostringstream screenout;
};

class KNCliContinuation : public KNAbstractContinuation
{
  public:
    KNCliContinuation(const KNConstants& constants) : KNAbstractContinuation(constants), output(0), charsPrinted(0) { }
    ~KNCliContinuation() { }
    void printStream() { std::cout << screenout.str(); charsPrinted += screenout.str().size(); screenout.str(""); std::cout.flush();  }
    virtual void storeCursor() { charsPrinted = 0; }
    virtual void clearLastLine()
    {
      for (size_t k=0; k<charsPrinted; ++k) std::cout << "\b";
    }
    static void printException(const KNException& ex)
    {
      std::cerr << ex.getMessage().str() << " This has occurred in file '" << ex.getFile() << "' at line " << ex.getLine() << ".\n";
    }
    void raiseException(const KNException& ex)
    {
      printException(ex);
      exit(-1);
    }
    void setData(KNDataFile* data) { output = data; }
    KNDataFile& data() { return *output; }
    void deleteData() { delete output; output = 0; }
    void dataUpdated() { }
  private:
    KNDataFile* output;
    size_t charsPrinted;
};

#endif
