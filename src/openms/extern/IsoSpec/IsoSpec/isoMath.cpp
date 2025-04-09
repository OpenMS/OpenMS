/*
 *   This file has been released into public domain by John D. Cook
 *   and is used here with some slight modifications (which are hereby
 *   also released into public domain),
 *
 *   This file is part of IsoSpec.
 */


// NOLINT(legal/copyright)


#include <cmath>
#include <cstdlib>
#include "isoMath.h"
#include "platform.h"
#include "btrd.h"

namespace IsoSpec
{


void release_g_lfact_table()
{
#if ISOSPEC_GOT_MMAN
    munmap(g_lfact_table, ISOSPEC_G_FACT_TABLE_SIZE*sizeof(double));
#else
    free(g_lfact_table);
#endif
}

double* alloc_lfact_table()
{
    double* ret;
# if ISOSPEC_GOT_MMAN
    ret = reinterpret_cast<double*>(mmap(nullptr, sizeof(double)*ISOSPEC_G_FACT_TABLE_SIZE, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
#else
    ret = reinterpret_cast<double*>(calloc(ISOSPEC_G_FACT_TABLE_SIZE, sizeof(double)));
#endif
    std::atexit(release_g_lfact_table);
    return ret;
}

double* g_lfact_table = alloc_lfact_table();


double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p < 0.5)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    else
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
}

double NormalCDFInverse(double p, double mean, double stdev)
{
    return mean + stdev * NormalCDFInverse(p);
}

double NormalCDF(double x, double mean, double stdev)
{
    x = (x-mean)/stdev * 0.7071067811865476;

    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double NormalPDF(double x, double mean, double stdev)
{
    double two_variance = stdev * stdev * 2.0;
    double delta = x-mean;
    return exp( -delta*delta / two_variance )      /     sqrt( two_variance * pi );
}

const double sqrt_pi = 1.772453850905516027298167483341145182798;

double LowerIncompleteGamma2(int a, double x)
{
    double base;
    double exp_minus_x = exp(-x);
    double current_s;
    if(a % 2 == 0)
    {
        base = 1 - exp_minus_x;
        current_s = 1.0;
        a--;
    }
    else
    {
        base = sqrt_pi * erf(sqrt(x));
        current_s = 0.5;
    }

    a = a/2;
    for(; a; a--)
    {
        base = base * current_s - pow(x, current_s) * exp_minus_x;
        current_s += 1.0;
    }

    return base;
}

double InverseLowerIncompleteGamma2(int a, double x)
{
    double l = 0.0;
    double p = tgamma(a);
    double s;

    do {
        s = (l+p) / 2.0;
        double v = LowerIncompleteGamma2(a, s);
        if (x < v)
            p = s;
        else
            l = s;
    } while((p-l)*1000.0 > p);

    return s;
}

std::random_device random_dev;
std::mt19937 random_gen(random_dev());
std::uniform_real_distribution<double> stdunif(0.0, 1.0);

size_t rdvariate_binom(size_t tries, double succ_prob, std::mt19937& rgen)
{
    if (succ_prob >= 1.0)
        return tries;
    return IsoSpec::boost_binomial_distribution_variate(tries, succ_prob, rgen);
}



}  // namespace IsoSpec

