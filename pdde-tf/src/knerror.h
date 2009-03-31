// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef ERROR_H
#define ERROR_H

#ifdef _MSC_VER
# define snprintf sprintf_s
#endif

#include <string>

class knutMessage
{
  public:
    knutMessage(const char *str) : message(str)
    { }
    knutMessage(const std::string& str) : message(str)
    { }
    knutMessage(const knutMessage& m) : message(m.message)
    { }
    knutMessage& operator<<(const std::string& str)
    {
      message += std::string(" ") + str;
      return *this;
    }
    knutMessage& operator<<(int i)
    {
      char buf[12+1];
      buf[12] = '\0';
      snprintf(buf, 12, " %d", i);
      message += std::string(buf);
      return *this;
    }
    knutMessage& operator<<(double i)
    {
      char buf[12+1];
      buf[12] = '\0';
      snprintf(buf, 12, " %lf", i);
      message += std::string(buf);
      return *this;
    }
    std::string message;
};

class knutException
{
  public:
    knutException() : file(""), line(0), message("")
    { }
    knutException(const std::string& f, int l, const knutMessage& m) :
        file(f), line(l), message(m)
    { }
    std::string file;
    int         line;
    knutMessage message;
};

#define P_MESSAGE(msg) throw(knutException( __FILE__, __LINE__, knutMessage(msg) ))
#define P_MESSAGE1(msg,arg1) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1 ))
#define P_MESSAGE2(msg,arg1,arg2) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1<<arg2 ))
#define P_MESSAGE3(msg,arg1,arg2,arg3) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1<<arg2<<arg3 ))

#ifdef DEBUG
#  define P_ASSERT(cond) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\ninternal error: " )); }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\ninternal error: "<<"\n\t"<<arg1 )); }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\ninternal error: "<<"\n\t"<<arg1<<arg2 )); }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\ninternal error: "<<"\n\t"<<arg1<<arg2<<arg3 )); }while(0)
#else
#  define P_ASSERT(cond) do{ }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2, arg3) do{ }while(0)
#endif

#define P_ERROR(cond) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\nerror: " )); }while(0)
#define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\nerror: "<<"\n\t"<<arg1 )); }while(0)
#define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\nerror: "<<arg1<<arg2 )); }while(0)
#define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\nerror: "<<arg1<<arg2<<arg3 )); }while(0)
#define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("(")<<#cond<<")"<<"\nerror: "<<arg1<<arg2<<arg3<<arg4 )); }while(0)
#endif
