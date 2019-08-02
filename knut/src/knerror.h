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
#include "config.h"

constexpr const char * __KNRP__2 (const char * path, const char * orig)
{
  return (*orig == '\0') ? path : __KNRP__2 (path + 1, orig + 1);
}

constexpr const char * KNRP (const char * path)
{
  return __KNRP__2 (path, KNUT_SOURCE_DIR);
}

class KNMessage
{
  public:
    KNMessage() = default;
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
    KNMessage& operator<<(size_t i)
    {
      char buf[31+1];
      buf[31] = '\0';
      snprintf(buf, 31, "%zd", i);
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
    std::string message{""};
};

class KNException
{
  public:
    KNException () : file(""), line(0), ipos(0), message("")
    { }
    KNException (const char * f, int l, const KNMessage& m, size_t pos = 0) :
        file(f), line(l), ipos(pos), message(m)
    { removePath(); }
    const char* getFile () const { return file; }
    int         getLine () const { return line; }
    const KNMessage& getMessage () const { return message; }
    size_t     getPos () const { return ipos; }
    const std::string& str () { message << " at " << file << ":" << line; return message.str(); }
    const std::string& exprStr (const std::string& expression);
  private:
    void removePath () { }
    const char*  file;
    const int    line;
    const size_t ipos;
    KNMessage message{};
};

#ifdef DEBUG

#include <iostream>

#  define P_ASSERT(cond) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X1(cond, arg1) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<std::endl; std::abort(); } }while(0)
#  define P_ASSERT_X11(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical assertion (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<arg11<<std::endl; std::abort(); } }while(0)

#  define P_MESSAGE1(msg) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE2(msg,arg1) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE3(msg,arg1,arg2) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE4(msg,arg1,arg2,arg3) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE5(msg,arg1,arg2,arg3,arg4) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE6(msg,arg1,arg2,arg3,arg4,arg5) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<arg5<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE7(msg,arg1,arg2,arg3,arg4,arg5,arg6) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<std::endl; std::abort(); }while(0)
#  define P_MESSAGE8(msg,arg1,arg2,arg3,arg4,arg5,arg6,arg7) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<std::endl; std::abort(); }while(0)

#  define P_ERROR(cond) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) { \
  std::cerr<<KNRP(KNRP(__FILE__))<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<std::endl; std::abort(); } }while(0)
#  define P_ERROR_X11(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<arg11<<std::endl; std::abort(); } }while(0)

// with ...
#  define PE_MESSAGE1(pos, msg) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE2(pos, msg,arg1) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE3(pos, msg,arg1,arg2) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE4(pos, msg,arg1,arg2,arg3) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE5(pos, msg,arg1,arg2,arg3,arg4) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE6(pos, msg,arg1,arg2,arg3,arg4,arg5) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<arg5<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE7(pos, msg,arg1,arg2,arg3,arg4,arg5,arg6) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<std::endl; std::abort(); }while(0)
#  define PE_MESSAGE8(pos, msg,arg1,arg2,arg3,arg4,arg5,arg6,arg7) do{ std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" "<<msg<<" "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<std::endl; std::abort(); }while(0)

#  define PE_ERROR(cond, pos) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X1(cond, pos, arg1) do{ if(!(cond) ) { \
  std::cerr<<KNRP(KNRP(__FILE__))<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X2(cond, pos, arg1, arg2) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X3(cond, pos, arg1, arg2, arg3) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X4(cond, pos, arg1, arg2, arg3, arg4) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X5(cond, pos, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X6(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X7(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X8(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X9(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X10(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<std::endl; std::abort(); } }while(0)
#  define PE_ERROR_X11(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ if(!(cond) ) { \
  std::cerr<<KNRP(__FILE__)<<":"<<__LINE__<<" Critical error (" #cond ") has failed. "<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<arg11<<std::endl; std::abort(); } }while(0)

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

#  define P_MESSAGE1(msg) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg) ))
#  define P_MESSAGE2(msg,arg1) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1 ))
#  define P_MESSAGE3(msg,arg1,arg2) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2 ))
#  define P_MESSAGE4(msg,arg1,arg2,arg3) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3 ))
#  define P_MESSAGE5(msg,arg1,arg2,arg3,arg4) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4 ))
#  define P_MESSAGE6(msg,arg1,arg2,arg3,arg4,arg5) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5 ))
#  define P_MESSAGE7(msg,arg1,arg2,arg3,arg4,arg5,arg6) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6 ))
#  define P_MESSAGE8(msg,arg1,arg2,arg3,arg4,arg5,arg6,arg7) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7 ))

#  define P_ERROR(cond) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ") )); } }while(0)
#  define P_ERROR_X1(cond, arg1) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1 )); } }while(0)
#  define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2 )); } }while(0)
#  define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3 )); } }while(0)
#  define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4 )); } }while(0)
#  define P_ERROR_X5(cond, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5 )); } }while(0)
#  define P_ERROR_X6(cond, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6 )); } }while(0)
#  define P_ERROR_X7(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7 )); } }while(0)
#  define P_ERROR_X8(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8 )); } }while(0)
#  define P_ERROR_X9(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9 )); } }while(0)
#  define P_ERROR_X10(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10 )); } }while(0)
#  define P_ERROR_X11(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<arg11 )); } }while(0)

// with expression syntax error
#  define PE_MESSAGE1(pos, msg) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg), pos ))
#  define PE_MESSAGE2(pos, msg,arg1) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1, pos ))
#  define PE_MESSAGE3(pos, msg,arg1,arg2) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2, pos ))
#  define PE_MESSAGE4(pos, msg,arg1,arg2,arg3) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3, pos ))
#  define PE_MESSAGE5(pos, msg,arg1,arg2,arg3,arg4) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4, pos ))
#  define PE_MESSAGE6(pos, msg,arg1,arg2,arg3,arg4,arg5) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5, pos ))
#  define PE_MESSAGE7(pos, msg,arg1,arg2,arg3,arg4,arg5,arg6) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6, pos ))
#  define PE_MESSAGE8(pos, msg,arg1,arg2,arg3,arg4,arg5,arg6,arg7) throw(KNException( KNRP(__FILE__), __LINE__, KNMessage(msg)<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7, pos ))

#  define PE_ERROR(cond, pos) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. "), pos )); } }while(0)
#  define PE_ERROR_X1(cond, pos, arg1) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1, pos )); } }while(0)
#  define PE_ERROR_X2(cond, pos, arg1, arg2) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2, pos )); } }while(0)
#  define PE_ERROR_X3(cond, pos, arg1, arg2, arg3) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3, pos )); } }while(0)
#  define PE_ERROR_X4(cond, pos, arg1, arg2, arg3, arg4) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4, pos )); } }while(0)
#  define PE_ERROR_X5(cond, pos, arg1, arg2, arg3, arg4, arg5) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5, pos )); } }while(0)
#  define PE_ERROR_X6(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6, pos )); } }while(0)
#  define PE_ERROR_X7(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7, pos )); } }while(0)
#  define PE_ERROR_X8(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8, pos )); } }while(0)
#  define PE_ERROR_X9(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9, pos )); } }while(0)
#  define PE_ERROR_X10(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10, pos )); } }while(0)
#  define PE_ERROR_X11(cond, pos, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ if(!(cond) ) { \
  throw(KNException( KNRP(__FILE__), __LINE__, KNMessage("Critical error (" #cond ") has failed. ")<<arg1<<arg2<<arg3<<arg4<<arg5<<arg6<<arg7<<arg8<<arg9<<arg10<<arg11, pos )); } }while(0)
#endif

#endif

