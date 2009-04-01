#ifndef BASECOMP_H
#define BASECOMP_H

#include "constants.h"

class BaseComp
{
  public:
    BaseComp(const NConstants& constants) : stopFlag(false)
    {
      params = new NConstants(constants);
    }
    virtual ~BaseComp()
    {
      delete params;
    }
    void setConstants(const NConstants& constants)
    {
      delete params;
      params = new NConstants(constants);
    }
    void run(const char* branchFile);
    void run() { run(0); }
    void setStopFlag(bool flag)
    {
      stopFlag = flag;
    }
    virtual void print(std::ostringstream& str) = 0;
    virtual void raiseException(const knutException& ex) = 0;
  protected:
    NConstants* params;
    bool        stopFlag;
};

class CLIComp : public BaseComp
{
  public:
    CLIComp(const NConstants& constants) : BaseComp(constants) { }
    ~CLIComp() { }
    void print(std::ostringstream& str) { std::cout<<str.str(); str.str(""); }
    static void printException(const knutException& ex)
    {
      std::cerr << ex.getMessage().str() << " This has occured in file '" << ex.getFile() << "' at line " << ex.getLine() << ".\n";
    }
    void raiseException(const knutException& ex)
    {
      printException(ex);
      exit(-1);
    }
};

#endif
