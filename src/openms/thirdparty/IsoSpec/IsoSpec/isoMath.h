/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
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
#include <random>

#if !defined(ISOSPEC_G_FACT_TABLE_SIZE)
// 10M should be enough for anyone, right?
// Actually, yes. If anyone tries to input a molecule that has more than 10M atoms,
// he deserves to get an exception thrown in his face. OpenMS guys don't want to alloc
// a table of 10M to memoize the necessary values though, use something smaller for them.
  #if ISOSPEC_BUILDING_OPENMS
    #define ISOSPEC_G_FACT_TABLE_SIZE 1024
  #else
    #define ISOSPEC_G_FACT_TABLE_SIZE 1024*1024*10
  #endif
#endif

namespace IsoSpec
{

extern double* g_lfact_table;

static inline double minuslogFactorial(int n)
{
    if (n < 2)
        return 0.0;
    #if ISOSPEC_BUILDING_OPENMS
    if (n >= ISOSPEC_G_FACT_TABLE_SIZE)
        return -lgamma(n+1);
    #endif
    if (g_lfact_table[n] == 0.0)
        g_lfact_table[n] = -lgamma(n+1);

    return g_lfact_table[n];
}

const double pi = 3.14159265358979323846264338328;
const double logpi = 1.144729885849400174143427351353058711647294812915311571513623071472137769884826079783623270275489708;

double NormalCDFInverse(double p);
double NormalCDFInverse(double p, double mean, double stdev);
double NormalCDF(double x, double mean, double stdev);
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);

// Returns lower incomplete gamma function of a/2, x, where a is int and > 0.
double LowerIncompleteGamma2(int a, double x);

// Returns y such that LowerIncompleteGamma2(a, y) == x. Approximately.
double InverseLowerIncompleteGamma2(int a, double x);

// Computes the inverse Cumulative Distribution Funcion of the Chi-Square distribution with k degrees of freedom
inline double InverseChiSquareCDF2(int k, double x)
{
    return InverseLowerIncompleteGamma2(k, x*tgamma(static_cast<double>(k)/2.0)) * 2.0;
}

extern std::mt19937 random_gen;
extern std::uniform_real_distribution<double> stdunif;

inline double rdvariate_beta_1_b(double b, std::mt19937& rgen = random_gen)
{
    return 1.0 - pow(stdunif(rgen), 1.0/b);
}


size_t rdvariate_binom(size_t tries, double succ_prob, std::mt19937& rgen = random_gen);




}  // namespace IsoSpec
