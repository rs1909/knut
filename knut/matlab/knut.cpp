// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005, 2012 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "basecomp.h" // includes exprsyste.h
#include "knerror.h"
#include <unistd.h>
#include <string>
#define __GXX_EXPERIMENTAL_CXX0X__ 1
#include <functional>
#include <algorithm>
// #include "constants.h" included in basecomp
// #include <iostream> included in basecomp
extern "C"
{
#include <signal.h>
#include <stdio.h>
// just to make the compiler accept the header file
#define __STDC_UTF_16__
#include "mex.h"
}

#ifdef MATLAB_MEX_FILE
extern "C" bool utIsInterruptPending();
#else
#include <octave/oct.h>
#include <octave/unwind-prot.h>
#endif

#include "mxdata.h"

static void ReplaceStringInPlace(std::string& subject, const std::string& search, const std::string& replace) 
{
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != std::string::npos)
  {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
}

class KNMatlabConstantsReader : public KNAbstractConstantsReader
{
  public:
    KNMatlabConstantsReader();
    KNMatlabConstantsReader(const mxArray *options) : options_in(options) {}
    virtual bool getSystem(std::string& strout);
    virtual bool getTextField(const char* field, std::string& strout);
    virtual bool getIndexField(const char* field, size_t& out);
    virtual bool getDoubleField(const char* field, double& out);
    virtual bool getStringListField(const char* field, std::vector<std::string>& out);
    virtual bool getIndexListField(const char* field, std::vector<size_t>& out);
  private:
    const mxArray *options_in;
};

class KNMatlabConstantsWriter : public KNAbstractConstantsWriter
{
  public:
    KNMatlabConstantsWriter() {}
    virtual void setSystem(const std::string& strout);
    virtual void setTextField(const char* fieldname, const std::string& str);
    virtual void setDoubleField(const char* fieldname, double val);
    virtual void setIndexField(const char* fieldname, size_t val);
    virtual void setStringListField(const char* field, const std::vector<std::string>& vec);
    virtual void setIndexListField(const char* field, const std::vector<size_t>& vec);
    mxArray *write();
  private:
    std::vector<std::string> fields;
    std::vector<mxArray *> cells;
};

class KNMatlabConstants : public KNConstants
{
  public:
    void fromMatlab (const mxArray *options);
    mxArray * toMatlab();
};

class KNMatlabContinuation : public KNAbstractContinuation
{
  public:
    KNMatlabContinuation () : output(nullptr), charsPrinted(0) { }
    ~KNMatlabContinuation () { delete output; }
    void printStream () 
    {
      mexPrintf (screenout.str().c_str());
      charsPrinted += screenout.str().size(); screenout.str("");
      mexEvalString("drawnow");
    #ifdef MATLAB_MEX_FILE
      if (utIsInterruptPending() && !getStopFlag()) /* check for a Ctrl-C event */
      {
        setStopFlag (true);
        mexPrintf ("Ctrl-C Detected. END\n\n");
      }
    #else
      OCTAVE_QUIT;
    #endif
    }
    virtual void storeCursor () { charsPrinted = 0; }
    virtual void clearLastLine ()
    {
      for (size_t k=0; k<charsPrinted; ++k) mexPrintf ("\b");
    }
    void raiseException (const KNException& ex)
    {
      std::ostringstream err;
      KNAbstractContinuation::printException (err, ex);
      mexErrMsgIdAndTxt ("MATLAB:Knut", err.str().c_str());
      std::cout <<"MATLAB:Knut" << err.str().c_str() << "\n";
    }
    void createData (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms)
    {
      delete output;
      output = KNMatlabData::createDataStatic (t, ndim, npar, prms);
    }
    KNMatlabData& data () { return *output; }
    void deleteData () { /*delete output; output = 0;*/ }
    void dataUpdated () 
    {
      // this is the place to put a callback to MATLAB
//       mexPrintf ("C");
    }
    bool validData () { return output != nullptr; }
  private:
    KNMatlabData* output;
    size_t charsPrinted;
};

static bool getText(const mxArray *elem, std::string& strout)
{
  if (elem)
  {
    const size_t len = mxGetM(elem)*mxGetN(elem)+1;
    char* str = new char[len];
    int ecd = mxGetString(elem, str, len);
    if (ecd == 0) { strout = str; return true; }
  }
  return false;
}

// API:
// 'get'
// 'set' +1 parameter
// 'run'
// 'cfile' +1 parameter
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{  
  static KNMatlabConstants* params = nullptr;
  static KNMatlabContinuation* comp = nullptr;
//   static std::string* vfString = nullptr;
  try {
//     P_ERROR_X1(nlhs == 0, "does not provide output.");
    P_ERROR_X1(nrhs > 0, "at least a command.");
    
    std::string cmd;
    getText(prhs[0], cmd);
    if (cmd.compare("set") == 0)
    {
      P_ERROR_X1(nrhs == 2, "constants structure is required");
      if (!params) params = new KNMatlabConstants;
      else mexPrintf ("Not a new instance.\n");
      params->fromMatlab (prhs[1]);
    } else if (cmd.compare("get") == 0)
    {
      if (params)
      {
        P_ERROR_X1(nlhs == 1, "output struct is required");
        mxDestroyArray(plhs[0]);
        plhs[0] = params->toMatlab();
      }
    } else if (cmd.compare("run") == 0)
    {
      P_ERROR_X1 (params != 0, "No constants are specified");
      P_ERROR_X1(nlhs == 1, "output struct is required");
      if (!comp) comp = new KNMatlabContinuation;
      else mexPrintf ("Not the first computation.\n");
      KNExprSystem* sys = nullptr;
      if (params->getSystem ().empty())
      {
        sys = new KNSystem (params->getSysName (), false);
      } else
      {
        sys = new KNExprSystem (params->getSystem (), false);
      }
      std::string tmp;
      sys->toString (tmp);
      // making matlab string
      tmp.insert (tmp.begin(), '\'');
      tmp.insert (tmp.begin(), '[');
      tmp.push_back ('\'');
      tmp.push_back (']');
      ReplaceStringInPlace (tmp, "\n", "'.""..\n'");
      mexPrintf ("%s\n", tmp.c_str());
      // E
      KNMatlabData* inputData = 0;
      if (nrhs == 2) inputData = new KNMatlabData (prhs[1]);
      comp -> run (sys, params, inputData);
      if (comp -> validData ()) plhs[0] = comp -> data() . getData();
      delete sys;
    } else if (cmd.compare("cfile") == 0)
    {
      P_ERROR_X1(nrhs == 2, "filename is missing");
      std::string fname;
      getText(prhs[1], fname);
      if (!params) params = new KNMatlabConstants;
      else mexPrintf ("Not a new instance.\n");
      params->loadXmlFileV5(fname);
    } else if (cmd.compare("steps") == 0)
    {
      if (comp)
      {
        if (comp -> validData ())
        {
          plhs[0] = mxCreateDoubleScalar (comp -> data() . getNCols());
        }
      }
    }
  }
  catch (KNException ex)
  {
    std::ostringstream err;
    KNAbstractContinuation::printException (err, ex);
    mexErrMsgIdAndTxt ("MATLAB:Knut", err.str().c_str());
    delete params;
  }
}

bool KNMatlabConstantsReader::getSystem(std::string& strout)
{
  return getTextField ("system", strout);
}

bool KNMatlabConstantsReader::getTextField(const char* field, std::string& strout)
{
  mxArray *elem = mxGetField(options_in, 0, field);
  if (elem)
  {
    if (!getText(elem, strout)) strout.clear();
    return true;
  }
  else return false;
}

bool KNMatlabConstantsReader::getIndexField(const char* field, size_t& out)
{
  mxArray *elem = mxGetField(options_in, 0, field);
  if (elem)
  {
    out = static_cast<size_t>(mxGetScalar(elem));
//     mexPrintf ("%s is being set to %d\n", field, out);
    return true;
  }
  else return false;
}

bool KNMatlabConstantsReader::getDoubleField(const char* field, double& out)
{
  mxArray *elem = mxGetField(options_in, 0, field);
  if (elem)
  {
    out = mxGetScalar(elem);
//     mexPrintf ("%s is being set to %lf\n", field, out);
    return true;
  }
  else return false;
}

bool KNMatlabConstantsReader::getStringListField(const char* field, std::vector<std::string>& out)
{
  mxArray *elem = mxGetField(options_in, 0, field);
  if (elem)
  {
    if (mxIsCell(elem))
    {
      size_t len_cell = mxGetNumberOfElements(elem);
      out.resize(len_cell);
      for (size_t k = 0; k < len_cell; k++)
      {
        mxArray* celem = mxGetCell(elem, k);
        if (celem)
        {
          getText(celem, out[k]);
        }
      }
      return true;
    }
  }
  return false;
}

bool KNMatlabConstantsReader::getIndexListField(const char* field, std::vector<size_t>& out)
{
  mxArray *elem = mxGetField(options_in, 0, field);
  if (elem)
  {
    if (mxIsCell(elem))
    {
      size_t len_cell = mxGetNumberOfElements(elem);
      out.resize(len_cell);
      for (size_t k = 0; k < len_cell; k++)
      {
        mxArray* celem = mxGetCell(elem, k);
        if (celem)
        {
          out[k] = static_cast<size_t>(mxGetScalar(celem));
        }
      }
      return true;
    }
  }
  return false;
}

void KNMatlabConstantsWriter::setSystem(const std::string& str)
{
  setTextField("system", str);
}

void KNMatlabConstantsWriter::setTextField(const char* fieldname, const std::string& str)
{
  if (!str.empty())
  {
    mxArray *fieldval = mxCreateString(str.c_str());
    fields.push_back(fieldname);
    cells.push_back(fieldval);
  }
}

void KNMatlabConstantsWriter::setDoubleField(const char* fieldname, double val)
{
  mxArray *fieldval = mxCreateDoubleScalar (val);
  fields.push_back(fieldname);
  cells.push_back(fieldval);
}

void KNMatlabConstantsWriter::setIndexField(const char* fieldname, size_t val)
{
  mxArray *fieldval = mxCreateDoubleScalar (static_cast<double>(val));
  fields.push_back(fieldname);
  cells.push_back(fieldval);
}

void KNMatlabConstantsWriter::setStringListField(const char* fieldname, const std::vector<std::string>& vec)
{
  mxArray* cell = mxCreateCellMatrix(1, vec.size());
  for (size_t k = 0; k < vec.size(); k++)
  {
    mxArray *el = mxCreateString(vec[k].c_str());
    mxSetCell(cell, k, el);
  }
  fields.push_back(fieldname);
  cells.push_back(cell);
}

void KNMatlabConstantsWriter::setIndexListField(const char* fieldname, const std::vector<size_t>& vec)
{
  mxArray* cell = mxCreateCellMatrix(1, vec.size());
  for (size_t k = 0; k < vec.size(); k++)
  {
    mxArray *el = mxCreateDoubleScalar(vec[k]);
    mxSetCell(cell, k, el);
  }
  fields.push_back(fieldname);
  cells.push_back(cell);
}

mxArray * KNMatlabConstantsWriter::write()
{
  const char* fnames[fields.size()];
  for (size_t k = 0; k < fields.size(); k++) fnames[k] = fields[k].c_str();
  mxArray *opts = mxCreateStructMatrix(1, 1, fields.size(), fnames);
  for (size_t k = 0; k < fields.size(); k++) mxSetField(opts, 0, fields[k].c_str(), cells[k]);
  return opts;
}

void KNMatlabConstants::fromMatlab (const mxArray *options)
{
  KNMatlabConstantsReader reader(options);
  read(reader);
}

mxArray * KNMatlabConstants::toMatlab()
{
  KNMatlabConstantsWriter writer;
  write(writer);
  return writer.write();
}
