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
// #include "constants.h" included in basecomp
// #include <iostream> included in basecomp
extern "C"
{
#include <signal.h>
#include <stdio.h>
// just to make the compiler accept the header file
#define __STDC_UTF_16__
#include "mex.h"

bool utIsInterruptPending();
  
}
// #include "vf.h"

class KNMatlabConstantsReader : public KNAbstractConstantsReader
{
  public:
    KNMatlabConstantsReader();
    KNMatlabConstantsReader(const mxArray *options) : options_in(options) {}
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

// class KNMatlabVectorField : public VectorField
// {
// public:
//   mxArray * toMatlab();
//   void fromMatlab(const mxArray *);
// };

class KNMatlabContinuation : public KNCliContinuation
{
  public:
    KNMatlabContinuation () : output(0), charsPrinted(0) { }
    ~KNMatlabContinuation () { }
    void printStream () 
    {
      mexPrintf (screenout.str().c_str());
      charsPrinted += screenout.str().size(); screenout.str("");
      mexEvalString("drawnow");
      if (utIsInterruptPending()) /* check for a Ctrl-C event */
      {
	setStopFlag(true);
	mexPrintf("Ctrl-C Detected. END\n\n");
      }
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
    KNDataFile& data () { return *output; }
    void deleteData () { delete output; output = 0; }
    void dataUpdated () { }
  private:
    KNDataFile* output;
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

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{  
  static KNMatlabConstants* params = nullptr;
  static std::string* vfString = nullptr;
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
//       delete params;
    } else if (cmd.compare("get") == 0)
    {
      if (params)
      {
// 	mexPrintf ("Getting someting.\n");
	P_ERROR_X1(nlhs == 1, "output struct is required");
// 	mexPrintf ("ptr b %d.\n", plhs[0]);
	mxDestroyArray(plhs[0]);
	plhs[0] = params->toMatlab();
// 	mexPrintf ("ptr e %d.\n", plhs[0]);
      }
    } else if (cmd.compare("run") == 0)
    {
      P_ERROR_X1 (params != 0, "No constants are specified");
      KNCliContinuation comp;
      KNSystem* sys = new KNSystem (params->getSysName ());
//       comp.run(sys, params);
      mexPrintf ("%s\n", params->getSysName ().c_str());
      std::string tmp;
      sys->toString (tmp);
      mexPrintf ("%s\n", tmp.c_str());
      comp.run (sys, params);
      mexPrintf ("END\n");
      delete sys;
    } else if (cmd.compare("cfile") == 0)
    {
      P_ERROR_X1(nrhs == 2, "filename is missing");
      std::string fname;
      getText(prhs[1], fname);
      if (!params) params = new KNMatlabConstants;
      else mexPrintf ("Not a new instance.\n");
      params->loadXmlFileV5(fname);
    } else if (cmd.compare("system") == 0)
    {
      if (params)
      {
// 	P_ERROR_X1(nlhs == 1, "output struct is required");
//         KNMatlabVectorField vf;
// 	vf.ReadXML(params->getSysName());
// 	int pserr = vf.ProcessSymbols();
// 	plhs[0] = vf.toMatlab();
      }
    } else if (cmd.compare("load") == 0)
    {
      P_ERROR_X1(nrhs == 2, "output struct is required");
      std::string vfexpr;
      getText(prhs[1], vfexpr);
      KNExprSystem vf (vfexpr, false);
      std::string tmp;
      vf.toString (tmp);
      mexPrintf ("%s\n", tmp.c_str());
    }
  }
  catch (KNException ex)
  {
    std::ostringstream err;
    KNAbstractContinuation::printException (err, ex);
    mexErrMsgIdAndTxt ("MATLAB:Knut", err.str().c_str());
//     std::cout <<"MATLAB:Knut" << err.str().c_str() << "\n";
    delete params;
//     exit(-1);
  }
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
// 	  mexPrintf ("%s[%d] is being set to %s\n", field, k, out[k].c_str());
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
// 	  mexPrintf ("%s[%d] is being set to %lf\n", field, k, out[k]);
	}
      }
      return true;
    }
  }
  return false;
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
  for (int k = 0; k < vec.size(); k++)
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
  for (int k = 0; k < vec.size(); k++)
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

class fieldWriter
{
public:
  void addField(const char *field, mxArray *src)
  {
    if (src)
    {
      fieldname.push_back(field);
      source.push_back(src);
    }
  }
  mxArray *write()
  {
    if (fieldname.size() > 0)
    {
      const char * fields[fieldname.size()];
      for (size_t k = 0; k < fieldname.size(); k++) fields[k] = fieldname[k];
      mxArray *res = mxCreateStructMatrix(1, 1, fieldname.size(), fields);
      for (size_t k = 0; k < fieldname.size(); k++) mxSetField(res, 0, fieldname[k], source[k]);
      return res;
    }
    return 0;
  }
private:
  std::vector<const char *> fieldname;
  std::vector<mxArray *> source;
};

static void makeCellField(fieldWriter& writer, const char *field, size_t sz, std::function<const std::string& (size_t)> func)
{
  mxArray* cell = mxCreateCellMatrix(1, sz);
  bool nonempty = false;
  for (size_t k = 0; k < sz; k++) { mxSetCell(cell, k, mxCreateString(func(k).c_str())); nonempty |= !func(k).empty(); }
  if (nonempty) writer.addField(field, cell);
}

// mxArray * KNMatlabVectorField::toMatlab()
// {
//   fieldWriter writer;
//     
//   fieldWriter name;
//   makeCellField(name, "name", 1, [this] (size_t k) { return this->Name(); } );
//   makeCellField(name, "description", 1, [this] (size_t k) { return this->Description(); } );
//   makeCellField(name, "independentvariable", 1, [this] (size_t k) { return this->IndependentVariable; } );
//   writer.addField( "name", name.write() );
//   
//   fieldWriter constant_writer;
//   makeCellField(constant_writer, "name", Constants.size(), [this] (size_t k) { return this->Constants[k]->Name(); } );
//   makeCellField(constant_writer, "value", Constants.size(), [this] (size_t k) { return this->Constants[k]->Value(); } );
//   makeCellField(constant_writer, "description", Constants.size(), [this] (size_t k) { return this->Constants[k]->Description(); } );
//   makeCellField(constant_writer, "latex", Constants.size(), [this] (size_t k) { return this->Constants[k]->Latex(); } );
//   writer.addField( "constant", constant_writer.write() );
//     
//   fieldWriter parameter;
//   makeCellField(parameter, "name", Parameters.size(), [this] (size_t k) { return this->Parameters[k]->Name(); } );
//   makeCellField(parameter, "defaultvalue", Parameters.size(), [this] (size_t k) { return this->Parameters[k]->DefaultValue(); } );
//   makeCellField(parameter, "description", Parameters.size(), [this] (size_t k) { return this->Parameters[k]->Description(); } );
//   makeCellField(parameter, "latex", Parameters.size(), [this] (size_t k) { return this->Parameters[k]->Latex(); } );
//   writer.addField( "parameter", parameter.write() );
//   
//   fieldWriter expression;
//   makeCellField(expression, "name", Expressions.size(), [this] (size_t k) { return this->Expressions[k]->Name(); } );
//   makeCellField(expression, "formula", Expressions.size(), [this] (size_t k) { return this->Expressions[k]->Formula(); } );
//   makeCellField(expression, "description", Expressions.size(), [this] (size_t k) { return this->Expressions[k]->Description(); } );
//   makeCellField(expression, "latex", Expressions.size(), [this] (size_t k) { return this->Expressions[k]->Latex(); } );
//   writer.addField( "expression", expression.write() );
// 
//   fieldWriter statevariable;
//   makeCellField(statevariable, "name", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->Name(); } );
//   makeCellField(statevariable, "formula", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->Formula(); } );
//   makeCellField(statevariable, "description", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->Description(); } );
//   makeCellField(statevariable, "latex", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->Latex(); } );
//   makeCellField(statevariable, "periodfrom", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->PeriodicFrom(); } );
//   makeCellField(statevariable, "periodto", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->PeriodicTo(); } );
//   makeCellField(statevariable, "defaultinitialcondition", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->DefaultInitialCondition(); } );
//   makeCellField(statevariable, "defaulthistory", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->DefaultHistory(); } );
//   makeCellField(statevariable, "mass", StateVariables.size(), [this] (size_t k) { return this->StateVariables[k]->Mass(); } );
//   writer.addField( "statevariable", statevariable.write() );
//   
//   fieldWriter function;
//   makeCellField(function, "name", Functions.size(), [this] (size_t k) { return this->Functions[k]->Name(); } );
//   makeCellField(function, "formula", Functions.size(), [this] (size_t k) { return this->Functions[k]->Formula(); } );
//   makeCellField(function, "description", Functions.size(), [this] (size_t k) { return this->Functions[k]->Description(); } );
//   writer.addField( "function", function.write() );
//   
//   return writer.write();
// }

template<class T> static void resize_vec(std::vector<T*>& vec, size_t sz )
{
  const size_t szp = vec.size();
  for (int k = sz; k < vec.size(); k++) delete vec[k];
  if (szp != sz) vec.resize(sz);
  for (int k = szp; k < sz; k++) vec[k] = new T("");
}

static size_t getCellField(mxArray* source, const char *field, 
			   std::function<void (size_t k)> rsfunc, 
			   std::function<void (size_t k, const char *)> func)
{
  if (source)
  {
    mxArray* cell = mxGetField(source, 0, field);
    if ((cell != 0) && (mxIsCell(cell)))
    {
      const size_t len = mxGetM(cell)*mxGetN(cell);
      if (rsfunc) rsfunc(len);
      for (size_t k = 0; k < len; k++)
      {
	mxArray* elem = mxGetCell(cell, k);
	if (elem)
	{
	  const size_t lens = mxGetM(elem)*mxGetN(elem);
	  char * str = new char[lens+2];
	  int c = mxGetString(elem, str, lens+1);
	  if (!c) func(k, str);
	  delete[] str;
	}
      }
      return len;
    }
  }
  return 0;
}

// void KNMatlabVectorField::fromMatlab(const mxArray *vf)
// {
//   mxArray *elem = mxGetField(vf, 0, "name");
//   if (elem) 
//   {
//     getCellField(elem, "name", 0, [this] (size_t k, const char * str) { this->Name(str); } );
//     getCellField(elem, "description", 0, [this] (size_t k, const char * str) { this->Description(str); } );
//     getCellField(elem, "independentvariable", 0, [this] (size_t k, const char * str) { this->IndependentVariable = str; } );
//   }
//   
//   elem = mxGetField(vf, 0, "constant");
//   if (elem)
//   {
//     getCellField(elem, "name", [this](size_t k){resize_vec<Constant>(this->Constants,k);}, [this] (size_t k, const char * str) { this->Constants[k]->Name(str); } );
//     getCellField(elem, "value", [this](size_t k){resize_vec<Constant>(this->Constants,k);}, [this] (size_t k, const char * str) { this->Constants[k]->Value(str); } );
//     getCellField(elem, "description", [this](size_t k){resize_vec<Constant>(this->Constants,k);}, [this] (size_t k, const char * str) { this->Constants[k]->Description(str); } );
//     getCellField(elem, "latex", [this](size_t k){resize_vec<Constant>(this->Constants,k);}, [this] (size_t k, const char * str) { this->Constants[k]->Latex(str); } );
//   }
// 
//   elem = mxGetField(vf, 0, "parameter");
//   if (elem)
//   {
//     getCellField(elem, "name", [this](size_t k){resize_vec<Parameter>(this->Parameters,k);}, [this] (size_t k, const char * str) { this->Parameters[k]->Name(str); } );
//     getCellField(elem, "defaultvalue", [this](size_t k){resize_vec<Parameter>(this->Parameters,k);}, [this] (size_t k, const char * str) { this->Parameters[k]->DefaultValue(str); });
//     getCellField(elem, "description", [this](size_t k){resize_vec<Parameter>(this->Parameters,k);}, [this] (size_t k, const char * str) { this->Parameters[k]->Description(str); });
//     getCellField(elem, "latex", [this](size_t k){resize_vec<Parameter>(this->Parameters,k);}, [this] (size_t k, const char * str) { this->Parameters[k]->Latex(str); });
//   }
//   
//   elem = mxGetField(vf, 0, "expression");
//   if (elem)
//   {
//     getCellField(elem, "name", [this](size_t k){resize_vec<Expression>(this->Expressions,k);}, [this] (size_t k, const char * str) { this->Expressions[k]->Name(str); } );
//     getCellField(elem, "formula", [this](size_t k){resize_vec<Expression>(this->Expressions,k);}, [this] (size_t k, const char * str) { this->Expressions[k]->Formula(str); } );
//     getCellField(elem, "description", [this](size_t k){resize_vec<Expression>(this->Expressions,k);}, [this] (size_t k, const char * str) { this->Expressions[k]->Description(str); } );
//     getCellField(elem, "latex", [this](size_t k){resize_vec<Expression>(this->Expressions,k);}, [this] (size_t k, const char * str) { this->Expressions[k]->Latex(str); } );
//   }
// 
//   elem = mxGetField(vf, 0, "statevariable");
//   if (elem)
//   {
//     getCellField(elem, "name", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->Name(str); } );
//     getCellField(elem, "formula", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->Formula(str); } );
//     getCellField(elem, "description", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->Description(str); } );
//     getCellField(elem, "latex", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->Latex(str); } );
//     getCellField(elem, "periodfrom", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->PeriodicFrom(str); } );
//     getCellField(elem, "periodto", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->PeriodicTo(str); } );
//     getCellField(elem, "defaultinitialcondition", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->DefaultInitialCondition(str); } );
//     getCellField(elem, "defaulthistory", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->DefaultHistory(str); } );
//     getCellField(elem, "mass", [this](size_t k){resize_vec<StateVariable>(this->StateVariables,k);}, [this] (size_t k, const char * str) { this->StateVariables[k]->Mass(str); } );
//   }
//   
//   elem = mxGetField(vf, 0, "function");
//   if (elem)
//   {
//     getCellField(elem, "name", [this](size_t k){resize_vec<Function>(this->Functions,k);}, [this] (size_t k, const char * str) { this->Functions[k]->Name(str); } );
//     getCellField(elem, "formula", [this](size_t k){resize_vec<Function>(this->Functions,k);}, [this] (size_t k, const char * str) { this->Functions[k]->Formula(str); } );
//     getCellField(elem, "description", [this](size_t k){resize_vec<Function>(this->Functions,k);}, [this] (size_t k, const char * str) { this->Functions[k]->Description(str); } );
//     getCellField(elem, "latex", [this](size_t k){resize_vec<Function>(this->Functions,k);}, [this] (size_t k, const char * str) { this->Functions[k]->Latex(str); } );
//   }
// }
