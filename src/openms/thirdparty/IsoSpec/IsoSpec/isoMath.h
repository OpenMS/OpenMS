/*
 *   Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#pragma once

#include <cmath>
#include <fenv.h>

#if !defined(ISOSPEC_G_FACT_TABLE_SIZE)
// 10M should be enough for anyone, right?
// Actually, yes. If anyone tries to input a molecule that has more than 10M atoms, 
// he deserves to get an exception thrown in his face.
#define ISOSPEC_G_FACT_TABLE_SIZE 1024 // *1024*10
#endif

namespace IsoSpec
{

extern double* g_lfact_table;

static inline double minuslogFactorial(int n) 
{ 
    // check if value is too small or too large to be found in the lookup table
    if (n < 2) 
        return 0.0;
    if (n >= ISOSPEC_G_FACT_TABLE_SIZE) return -lgamma(n+1);

    // store in lookup table if not present yet
    if (g_lfact_table[n] == 0.0)
        g_lfact_table[n] = -lgamma(n+1);

    return g_lfact_table[n];
}
double NormalCDFInverse(double p);
double NormalCDFInverse(double p, double mean, double stdev);
double NormalCDF(double x, double mean, double stdev);
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);

} // namespace IsoSpec

