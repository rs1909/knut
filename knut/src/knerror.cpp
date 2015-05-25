// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <string>
#include "knerror.h"

const std::string& KNException::exprStr (const std::string& expression)
{
  if (ipos != 0)
  {
    size_t start=ipos, end=ipos;
    while (expression [start] != '\n' && start > 0) start--;
    if (start != 0) start++; // start point to '\n' before this point
    size_t line = 1, k = 0;
    while (k < ipos) { if (expression [k] == '\n') { line++; }; k++; };
    while (end < expression.size () && expression [end] != '\n') end++;
    std::string msgUL(end-start, '-');
    msgUL [ipos - start] = '^';
    message << " at line " << line 
            << " column " << 1 + ipos - start 
            << '\n' + expression.substr (start, end - start) + '\n' + msgUL + '\n';
  }
  message << " at " << file << ":" << line;
  return message.str();
}
