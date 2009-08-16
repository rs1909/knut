//
//  codegen_utils.cpp
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

//  TO DO -- Add 'const' declarations where appropriate.


#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

static const char *datetimefmt = "Generated on %2d-%3s-%04d at %02d:%02d";
static char datetimebuf[40]; // Must be enough space to hold a copy of datetimefmt

char *DateTimeMsg()
    {
    time_t now_t;
    tm     now;
    const char *months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

    time(&now_t);
    now = *(localtime(&now_t));
    sprintf(datetimebuf,datetimefmt,
                  now.tm_mday,months[now.tm_mon],now.tm_year+1900,
                  now.tm_hour,now.tm_min);
    return(datetimebuf);
    }


void PrintVFGENComment(ofstream &fout, const char *prefix)
    {
    fout << prefix << "This file was generated by the program VFGEN (Version:" << VERSION << ")" << endl;
    fout << prefix << DateTimeMsg() << endl;
    }

void GetFromVector(ofstream &fout, const char *skip, lst names, const char *vector,
                      const char *braces, int istart, const char *term)
    {
    int n;

    n = names.nops();
    for (int i = 0; i < n; ++i)
        {
        fout << skip;
        fout.width(10);
        fout << left << names[i];
        fout.width(0);
        fout << " = " << vector << braces[0] << (i+istart) << braces[1] << term << endl;
        }
    }


void GetFromVector2(ofstream &fout, const char *skip, lst names, const char *vector,
                      const char *bropen, const char *brclose, int istart, const char *term)
    {
    int n;

    n = names.nops();
    for (int i = 0; i < n; ++i)
        {
        fout << skip;
        fout.width(10);
        fout << left << names[i];
        fout.width(0);
        fout << " = " << vector << bropen << (i+istart) << brclose << term << endl;
        }
    }

