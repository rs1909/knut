// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
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

class pddeMessage
{
 public:
	pddeMessage( const char *str ) : message(str) { }
	pddeMessage( const std::string& str ) : message(str) { }
	pddeMessage( const pddeMessage& m ) : message(m.message) { }
	pddeMessage& operator<<( const std::string& str ) { message += std::string(" ") + str; return *this; }
	pddeMessage& operator<<( int i ) { char buf[12+1]; buf[12]='\0'; snprintf( buf, 12, " %d", i ); message += std::string(buf); return *this; }
	pddeMessage& operator<<( double i ) { char buf[12+1]; buf[12]='\0'; snprintf( buf, 12, " %lf", i ); message += std::string(buf); return *this; }
	std::string message;
};

class pddeException
{
 public:
	pddeException() : file(""), line(0), message("") { }
	pddeException( const std::string& f, int l, const pddeMessage& m ) :
		file(f), line(l), message(m) { }
	std::string file;
	int         line;
	pddeMessage message;
};

#define P_MESSAGE(msg) throw(pddeException( __FILE__, __LINE__, pddeMessage(msg) ))

#ifdef DEBUG
#  define P_ASSERT(cond) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\ninternal error: " )); }while(0)
#  define P_ASSERT_X(cond, arg1) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\ninternal error: "<<"\n\t"<<arg1 )); }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\ninternal error: "<<"\n\t"<<arg1<<arg2 )); }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\ninternal error: "<<"\n\t"<<arg1<<arg2<<arg3 )); }while(0)
#else
#  define P_ASSERT(cond) do{ }while(0)
#  define P_ASSERT_X(cond, arg1) do{ }while(0)
#  define P_ASSERT_X2(cond, arg1, arg2) do{ }while(0)
#  define P_ASSERT_X3(cond, arg1, arg2, arg3) do{ }while(0)
#endif

#define P_ERROR(cond) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\nerror: " )); }while(0)
#define P_ERROR_X(cond, arg1) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\nerror: "<<"\n\t"<<arg1 )); }while(0)
#define P_ERROR_X2(cond, arg1, arg2) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\nerror: "<<arg1<<arg2 )); }while(0)
#define P_ERROR_X3(cond, arg1, arg2, arg3) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\nerror: "<<arg1<<arg2<<arg3 )); }while(0)
#define P_ERROR_X4(cond, arg1, arg2, arg3, arg4) do{ if(!(cond) ) \
		throw(pddeException( std::string(__FILE__), __LINE__, pddeMessage("(")<<#cond<<")"<<"\nerror: "<<arg1<<arg2<<arg3<<arg4 )); }while(0)
#endif
