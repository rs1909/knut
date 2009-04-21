// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef KNERROR_H
#define KNERROR_H

#ifdef _MSC_VER
# define snprintf sprintf_s
#endif

#include <cstdio>
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
      message += str;
      return *this;
    }
    knutMessage& operator<<(int i)
    {
      char buf[31+1];
      buf[31] = '\0';
      snprintf(buf, 31, "%d", i);
      message += std::string(buf);
      return *this;
    }
    knutMessage& operator<<(double i)
    {
      char buf[31+1];
      buf[31] = '\0';
      snprintf(buf, 31, "%lf", i);
      message += std::string(buf);
      return *this;
    }
    const std::string& str() const { return message; }
  private:
    std::string message;
};

class knutException
{
  public:
    knutException() : file(""), line(0), message("")
    { }
    knutException(const std::string& f, int l, const knutMessage& m) :
        file(f), line(l), message(m)
    { removePath(); }
    const std::string& getFile() const { return file; }
    int                getLine() const { return line; }
    const knutMessage& getMessage() const { return message; } 
  private:
    void removePath();
    std::string file;
    int         line;
    knutMessage message;
};

#ifdef DEBUG

#  define P_ASSERT(cond) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Critical assertion (" #cond ") has failed. ") )); }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Critical assertion (" #cond ") has failed. ")<<arg1 )); }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2 )); }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2) do{ if(!(cond) ) \
  throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3 )); }while(0)

#  define P_MESSAGE1(msg) throw(-1)
#  define P_MESSAGE2(msg,arg1) throw(-1)
#  define P_MESSAGE3(msg,arg1,arg2) throw(-1)
#  define P_MESSAGE4(msg,arg1,arg2,arg3) throw(-1)
#  define P_MESSAGE5(msg,arg1,arg2,arg3,arg4) throw(-1)

#  define P_ERROR(cond) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) throw(-1); }while(0)

#else

#  define P_ASSERT(cond) do{ }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2, arg3) do{ }while(0)

#  define P_MESSAGE1(msg) throw(knutException( __FILE__, __LINE__, knutMessage(msg) ))
#  define P_MESSAGE2(msg,arg1) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1 ))
#  define P_MESSAGE3(msg,arg1,arg2) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1<<arg2 ))
#  define P_MESSAGE4(msg,arg1,arg2,arg3) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1<<arg2<<arg3 ))
#  define P_MESSAGE5(msg,arg1,arg2,arg3,arg4) throw(knutException( __FILE__, __LINE__, knutMessage(msg)<<arg1<<arg2<<arg3<<arg4 ))

#  define P_ERROR(cond) do{ if(!(cond) ) \
    throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Assertion (" #cond ") has failed. ") )); }while(0)
#  define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) \
    throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Assertion (" #cond ") has failed. ")<<arg1 )); }while(0)
#  define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) \
    throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2 )); }while(0)
#  define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) \
    throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3 )); }while(0)
#  define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) \
    throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4 )); }while(0)
#  define P_ERROR_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) \
    throw(knutException( std::string(__FILE__), __LINE__, knutMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5 )); }while(0)
#endif

#endif

