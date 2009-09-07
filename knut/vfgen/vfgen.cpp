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



//
// This function overrides the default print method for the Zlags_ function.
// Ginac thinks Zlags_ is a function, but in fact, in the code generated in
// a couple of the VFGEN commands that handle delay equations, Zlags_ is a
// two-dimensional array.  So we want Zlags_(1,2) printed like that, not as
// Zlags_(1.0,2.0).
//

static void Zlags_print(const ex& arg1, const ex& arg2, const print_context& c)
    {
    c.s << "Zlags_(";
    if (is_a<numeric>(arg1))
        c.s << ex_to<numeric>(arg1).to_int();
    else
        arg1.print(c);
    c.s << ",";
    if (is_a<numeric>(arg2))
        c.s << ex_to<numeric>(arg2).to_int();
    else
        arg2.print(c);
    c.s << ",idx)";
    }

//
// Add a function called "delay" to the functions known by ginac.
// (The corresponding macro DECLARE_FUNCTION_2P(delay) is in vf.h.)
//

REGISTER_FUNCTION(delay, dummy())
REGISTER_FUNCTION(Zlags_, print_func<print_csrc_float>(Zlags_print).
                          print_func<print_csrc_double>(Zlags_print).
                          print_func<print_python>(Zlags_print) )

static ex heaviside_evalf(const ex& x)
{
  if (is_exactly_a<numeric>(x))
  {
    if (ex_to<numeric>(x) > 0) return 1;
    else return 0;
  } else
    return heaviside(x).hold();
}

static ex heaviside_eval(const ex& x)
{
  return heaviside(x).hold();
}

static ex heaviside_deriv(const ex& x, unsigned deriv_param)
{
  return 0;
}

REGISTER_FUNCTION(heaviside, eval_func(heaviside_eval).
                             evalf_func(heaviside_evalf).
                             derivative_func(heaviside_deriv) )
                             
static ex ramp_evalf(const ex& x)
{
  if (is_exactly_a<numeric>(x))
  {
    if (ex_to<numeric>(x) > 0) return x;
    else return 0;
  } else
    return ramp(x).hold();
}

static ex ramp_eval(const ex& x)
{
  return ramp(x).hold();
}

static ex ramp_deriv(const ex& x, unsigned deriv_param)
{
  return heaviside(x);
}

REGISTER_FUNCTION(ramp, eval_func(ramp_eval).
                        evalf_func(ramp_evalf).
                        derivative_func(ramp_deriv) )

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

  map<string,string> options;
      
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
  if (vf.HasNonconstantDelay) cerr << "Nonconstant delays are not suported yet. ";
  else vf.PrintKnut(std::cout, options);
  return(0);
}  // end main()
