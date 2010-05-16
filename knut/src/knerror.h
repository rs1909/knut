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

class KNMessage
{
  public:
    KNMessage(const char *str) : message(str)
    { }
    KNMessage(const std::string& str) : message(str)
    { }
    KNMessage(const KNMessage& m) : message(m.message)
    { }
    KNMessage& operator<<(const std::string& str)
    {
      message += str;
      return *this;
    }
    KNMessage& operator<<(int i)
    {
      char buf[31+1];
      buf[31] = '\0';
      snprintf(buf, 31, "%d", i);
      message += std::string(buf);
      return *this;
    }
    KNMessage& operator<<(double i)
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

class KNException
{
  public:
    KNException() : file(""), line(0), message("")
    { }
    KNException(const std::string& f, int l, const KNMessage& m) :
        file(f), line(l), message(m)
    { removePath(); }
    const std::string& getFile() const { return file; }
    int                getLine() const { return line; }
    const KNMessage& getMessage() const { return message; } 
  private:
    void removePath();
    std::string file;
    int         line;
    KNMessage message;
};

#ifdef DEBUG

#  define P_ASSERT(cond) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ") )); }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1 )); }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2 )); }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3 )); }while(0)
#  define P_ASSERT_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4 )); }while(0)
#  define P_ASSERT_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5 )); }while(0)
#  define P_ASSERT_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6 )); }while(0)
#  define P_ASSERT_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7 )); }while(0)
#  define P_ASSERT_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8 )); }while(0)
#  define P_ASSERT_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9 )); }while(0)
#  define P_ASSERT_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10 )); }while(0)
#  define P_ASSERT_X11(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ if(!(cond) ) \
  throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Critical assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<arg11 )); }while(0)

#  define P_MESSAGE1(msg) throw(-1)
#  define P_MESSAGE2(msg,arg1) throw(-1)
#  define P_MESSAGE3(msg,arg1,arg2) throw(-1)
#  define P_MESSAGE4(msg,arg1,arg2,arg3) throw(-1)
#  define P_MESSAGE5(msg,arg1,arg2,arg3,arg4) throw(-1)
#  define P_MESSAGE6(msg,arg1,arg2,arg3,arg4,arg5) throw(-1)
#  define P_MESSAGE7(msg,arg1,arg2,arg3,arg4,arg5,arg6) throw(-1)
#  define P_MESSAGE8(msg,arg1,arg2,arg3,arg4,arg5,arg6,arg7) throw(-1)

#  define P_ERROR(cond) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) throw(-1); }while(0)
#  define P_ERROR_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) throw(-1); }while(0)

#else

#  define P_ASSERT(cond) do{ }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2, arg3) do{ }while(0)
#  define P_ASSERT_X4(cond, arg1, arg2, arg3, arg4) do{ }while(0)
#  define P_ASSERT_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ }while(0)
#  define P_ASSERT_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ }while(0)
#  define P_ASSERT_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ }while(0)
#  define P_ASSERT_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ }while(0)
#  define P_ASSERT_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ }while(0)
#  define P_ASSERT_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ }while(0)
#  define P_ASSERT_X11(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ }while(0)

#  define P_MESSAGE1(msg) throw(KNException( __FILE__, __LINE__, KNMessage(msg) ))
#  define P_MESSAGE2(msg,arg1) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1 ))
#  define P_MESSAGE3(msg,arg1,arg2) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1<<arg2 ))
#  define P_MESSAGE4(msg,arg1,arg2,arg3) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3 ))
#  define P_MESSAGE5(msg,arg1,arg2,arg3,arg4) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4 ))
#  define P_MESSAGE6(msg,arg1,arg2,arg3,arg4,arg5) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5 ))
#  define P_MESSAGE7(msg,arg1,arg2,arg3,arg4,arg5,arg6) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6 ))
#  define P_MESSAGE8(msg,arg1,arg2,arg3,arg4,arg5,arg6,arg7) throw(KNException( __FILE__, __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7 ))

#  define P_ERROR(cond) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ") )); }while(0)
#  define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1 )); }while(0)
#  define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2 )); }while(0)
#  define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3 )); }while(0)
#  define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4 )); }while(0)
#  define P_ERROR_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5 )); }while(0)
#  define P_ERROR_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6 )); }while(0)
#  define P_ERROR_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7 )); }while(0)
#  define P_ERROR_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8 )); }while(0)
#  define P_ERROR_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9 )); }while(0)
#  define P_ERROR_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) \
    throw(KNException( std::string(__FILE__), __LINE__, KNMessage("Assertion (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10 )); }while(0)

#endif

#endif

