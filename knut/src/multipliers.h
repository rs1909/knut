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

size_t unstableMultipliers(const KNVector& mulRe, const KNVector& mulIm, const size_t lp, const size_t pd, const size_t ns, size_t pt);
BifType bifurcationType(const KNVector& mulRe0, const KNVector& mulIm0, 
                        const KNVector& mulRe1, const KNVector& mulIm1,
                        const size_t lp, const size_t pd, const size_t ns, size_t pt0, size_t pt1);

#endif
