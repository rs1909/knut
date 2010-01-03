// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cstring>
#include "knerror.h"
#include "config.h"


void KNException::removePath()
{
  std::string::size_type loc = file.find(KNUT_SOURCE_DIR);
  if (loc != std::string::npos)
  {
    file.erase(loc, loc + strlen(KNUT_SOURCE_DIR));
  }
  if (file[0] == '/') file.erase(0,1);
}
