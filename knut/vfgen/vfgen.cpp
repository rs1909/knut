//
//  vfgen.cpp -- Multi-format vector field file generator.
//
//  by Warren Weckesser
//
//
//  Copyright (C) 2008 Warren Weckesser
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License, Version 2, as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program.  The file LICENSE-GPL2 in the VFGEN source directory
//  contains a copy. You may also obtain a copy by writing to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <ginac/ginac.h>

#include "vf.h"

using namespace std;
using namespace GiNaC;

///////////////////////////////////////////////////////////
//  main
///////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  VectorField vf;

  if (argc != 2)
  {
    std::cerr << "Bad number of arguments.\n";
    exit(-1);
  }

  map<string, string> options;

  //
  //  Read the vector field file.  This just puts the strings into the
  //  appropriate fields.  It doesn't do any symbolic processing.
  //
  vf.ReadXML(argv[1]);

  //
  //  Process the strings to create the GiNaC symbolic expressions in the object.
  //
  int pserr = vf.ProcessSymbols();
  if (pserr == -1)
  {
    exit(-1);
  }

  //
  // Call the appropriate output function based on the first
  // command line argument.
  //
  if (vf.testHasNonconstantDelay()) cerr << "Nonconstant delays are not suported yet. ";
  else vf.PrintKnut(std::cout, options);
  return(0);
}  // end main()
