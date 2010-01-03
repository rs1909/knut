// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2007 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef MULTIPLIERS_H
#define MULTIPLIERS_H

#include "pointtype.h"

class KNVector;

int unstableMultipliers(const KNVector& mulRe, const KNVector& mulIm, const int lp, const int pd, const int ns);
BifType bifurcationType(const KNVector& mulRe, const KNVector& mulIm, const int lp, const int pd, const int ns);

#endif
