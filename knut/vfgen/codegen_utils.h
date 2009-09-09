//
//  codegen_utils.h
//
//
//  Prototypes for functions in codegen_utils.cpp
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

#ifndef CODEGEN_UTILS
#define CODEGEN_UTILS 1

#include <ginac/ginac.h>

char *DateTimeMsg();
void PrintVFGENComment(std::ostream &fout, const char *prefix);
void GetFromVector(std::ostream &fout, const char *skip, GiNaC::lst names, const char *vector,
                   const char *braces, size_t istart, const char *term);
void GetFromVector2(std::ostream &fout, const char *skip, GiNaC::lst names, const char *vector,
                    const char *bropen, const char *brclose, size_t istart, const char *term);
GiNaC::ex iterated_subs(GiNaC::ex f, GiNaC::lst e);

#endif
