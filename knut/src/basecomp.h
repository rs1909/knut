#ifndef BASECOMP_H
#define BASECOMP_H

#include "constants.h"
#include "mat4data.h"

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
    virtual void print(std::ostringstream& str) = 0;
    virtual void raiseException(const KNException& ex) = 0;
    virtual void setData(KNDataFile* data) = 0;
    virtual KNDataFile& data() = 0;
    virtual void deleteData() = 0;
    virtual void dataUpdated() = 0;
  protected:
    KNConstants* params;
    bool        stopFlag;
};

class KNCliContinuation : public KNAbstractContinuation
{
  public:
    KNCliContinuation(const KNConstants& constants) : KNAbstractContinuation(constants), output(0) { }
    ~KNCliContinuation() { }
    void print(std::ostringstream& str) { std::cout<<str.str(); str.str(""); }
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
};

#endif
