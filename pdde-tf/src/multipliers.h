// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2007 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef MULTIPLIERS_H
#define MULTIPLIERS_H

#include "pointtype.h"

class Vector;

int unstableMultipliers(const Vector& mulRe, const Vector& mulIm, const int lp, const int pd, const int ns);
PtType bifurcationType(const Vector& mulRe, const Vector& mulIm, const int lp, const int pd, const int ns);

#endif
