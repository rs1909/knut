#ifndef BASECOMP_H
#define BASECOMP_H

#include "constants.h"
#include "mat4data.h"
#include <iostream>

enum class DataType : int { ST, LC, TR };

class KNAbstractContinuation
{
  public:
    KNAbstractContinuation ();
//     KNAbstractContinuation(const KNConstants& constants);
    virtual ~KNAbstractContinuation();
    
    // sets the constants file pointer
//     void setConstants(const KNConstants& constants);
    // run the continuation 
    void run(KNExprSystem* sys, KNConstants* constants, 
             const KNAbstractData* inputData);
//     void run(const char* branchFile);
//     void run() { run(0); }
    // sets the flag to stop computation
    void setStopFlag(bool flag);
    bool getStopFlag () { return stopFlag;}

    // these are all called from the thread, hence should not be called from outside.
    virtual std::ostream& outStream() { return screenout; };
    // screenout should be printed to the output
    virtual void printStream() = 0;
    // clear the last line on the screen
    virtual void clearLastLine() = 0;
    // store the current position on the screen so that on can delete anything after that
    virtual void storeCursor() = 0;
    // notify the user that an exception has occured
    virtual void raiseException(const KNException& ex) = 0;
    // transforms error message into a stream
    static void printException(std::ostream& err, const KNException& ex)
    {
      err << ex.getMessage().str() << " This has occurred in file '" << ex.getFile() << "' at line " << ex.getLine() << ".\n";
    }
    // create a data file for steady states
    virtual void createData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms) = 0;
    // the thread requests the data file
    virtual KNAbstractData& data() = 0;
    // the thread wants to close the data file
    virtual void deleteData() = 0;
    // the thread notifies when data should be closed
    virtual void dataUpdated() = 0;
    // create a data file, to be used to implement the virtual function
    static KNDataFile* createDataStatic (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms);
  protected:
//     KNConstants* params;
    bool        stopFlag;
    // this is used within the thread
    std::ostringstream screenout;
};

inline KNAbstractContinuation::KNAbstractContinuation () : stopFlag(false) {}

// inline KNAbstractContinuation::KNAbstractContinuation(const KNConstants& constants) : stopFlag(false)
// {
//   params = new KNConstants(constants);
// }

inline KNAbstractContinuation::~KNAbstractContinuation()
{
//   delete params;
}

// void inline KNAbstractContinuation::setConstants(const KNConstants& constants)
// {
//   delete params;
//   params = new KNConstants(constants);
// }

void inline KNAbstractContinuation::setStopFlag(bool flag)
{
  stopFlag = flag;
}

class KNCliContinuation : public KNAbstractContinuation
{
  public:
//     KNCliContinuation(const KNConstants& constants) : KNAbstractContinuation(constants), output(0), charsPrinted(0) { }
    KNCliContinuation () : output(0), charsPrinted(0) { }
    ~KNCliContinuation() { }
    void printStream() { std::cout << screenout.str(); charsPrinted += screenout.str().size(); screenout.str(""); std::cout.flush();  }
    virtual void storeCursor() { charsPrinted = 0; }
    virtual void clearLastLine()
    {
      for (size_t k=0; k<charsPrinted; ++k) std::cout << "\b";
    }
    static void printException(const KNException& ex)
    {
      KNAbstractContinuation::printException(std::cerr, ex);
    }
    void raiseException(const KNException& ex)
    {
      printException(ex);
      exit(-1);
    }
    void createData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms)
    {
      output = KNAbstractContinuation::createDataStatic (fileName, t, ndim, npar, prms);
    }
    KNAbstractData& data() { return *output; }
    void deleteData() { delete output; output = 0; }
    void dataUpdated() { }
  private:
    KNDataFile* output;
    size_t charsPrinted;
};

#endif
