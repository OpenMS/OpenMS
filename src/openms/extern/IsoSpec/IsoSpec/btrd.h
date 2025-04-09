/* This file was taken from Boost, as permitted by Boost licence,
 * with slight modifications. The reason is: we don't want to introduce
 * dependency on the whole Boost just for this one thing.
 *
 * Source: boost random/binomial_distribution.hpp header file, at version 1.71
 *
 * Copyright Steven Watanabe 2010
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org for most recent version including documentation.
 *
 */

#pragma once

#include "isoMath.h"
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <limits>


namespace IsoSpec {

typedef double RealType;
typedef int64_t IntType;


static const RealType btrd_binomial_table[10] = {
    0.08106146679532726,
    0.04134069595540929,
    0.02767792568499834,
    0.02079067210376509,
    0.01664469118982119,
    0.01387612882307075,
    0.01189670994589177,
    0.01041126526197209,
    0.009255462182712733,
    0.008330563433362871
};


/**
 * The binomial distribution is an integer valued distribution with
 * two parameters, @c t and @c p.  The values of the distribution
 * are within the range [0,t].
 *
 * The distribution function is
 * \f$\displaystyle P(k) = {t \choose k}p^k(1-p)^{t-k}\f$.
 *
 * The algorithm used is the BTRD algorithm described in
 *
 *  @blockquote
 *  "The generation of binomial random variates", Wolfgang Hormann,
 *  Journal of Statistical Computation and Simulation, Volume 46,
 *  Issue 1 & 2 April 1993 , pages 101 - 110
 *  @endblockquote
 */


// computes the correction factor for the Stirling approximation
// for log(k!)
static RealType fc(IntType k)
{
    if(k < 10) { return btrd_binomial_table[k]; }
    else
    {
        RealType ikp1 = RealType(1) / (k + 1);
        return (RealType(1)/12
             - (RealType(1)/360
             - (RealType(1)/1260)*(ikp1*ikp1))*(ikp1*ikp1))*ikp1;
    }
}

IntType btrd(IntType _t, RealType p, IntType m, std::mt19937& urng = random_gen)
{
    using std::floor;
    using std::abs;
    using std::log;

    RealType btrd_r = p/(1-p);
    RealType btrd_nr = (_t+1)*btrd_r;
    RealType btrd_npq = _t*p*(1-p);
    RealType sqrt_npq = sqrt(btrd_npq);
    RealType btrd_b = 1.15 + 2.53 * sqrt_npq;
    RealType btrd_a = -0.0873 + 0.0248*btrd_b + 0.01*p;
    RealType btrd_c = _t*p + 0.5;
    RealType btrd_alpha = (2.83 + 5.1/btrd_b) * sqrt_npq;
    RealType btrd_v_r = 0.92 - 4.2/btrd_b;
    RealType btrd_u_rv_r = 0.86*btrd_v_r;

    while(true) {
        RealType u;
        RealType v = stdunif(urng);
        if(v <= btrd_u_rv_r) {
            u = v/btrd_v_r - 0.43;
            return static_cast<IntType>(floor(
                (2*btrd_a/(0.5 - abs(u)) + btrd_b)*u + btrd_c));
        }

        if(v >= btrd_v_r) {
            u = stdunif(urng) - 0.5;
        } else {
            u = v/btrd_v_r - 0.93;
            u = ((u < 0)? -0.5 : 0.5) - u;
            v = stdunif(urng) * btrd_v_r;
        }

        RealType us = 0.5 - abs(u);
        IntType k = static_cast<IntType>(floor((2*btrd_a/us + btrd_b)*u + btrd_c));
        if(k < 0 || k > _t) continue;
        v = v*btrd_alpha/(btrd_a/(us*us) + btrd_b);
        RealType km = abs(k - m);
        if(km <= 15) {
            RealType f = 1;
            if(m < k) {
                IntType i = m;
                do {
                    ++i;
                    f = f*(btrd_nr/i - btrd_r);
                } while(i != k);
            } else if(m > k) {
                IntType i = k;
                do {
                    ++i;
                    v = v*(btrd_nr/i - btrd_r);
                } while(i != m);
            }
            if(v <= f) return k;
            else continue;
        } else {
            // final acceptance/rejection
            v = log(v);
            RealType rho =
                (km/btrd_npq)*(((km/3. + 0.625)*km + 1./6)/btrd_npq + 0.5);
            RealType t = -km*km/(2*btrd_npq);
            if(v < t - rho) return k;
            if(v > t + rho) continue;

            IntType nm = _t - m + 1;
            RealType h = (m + 0.5)*log((m + 1)/(btrd_r*nm))
                       + fc(m) + fc(_t - m);

            IntType nk = _t - k + 1;
            if(v <= h + (_t+1)*log(static_cast<RealType>(nm)/nk)
                      + (k + 0.5)*log(nk*btrd_r/(k+1))
                      - fc(k)
                      - fc(_t - k))
            {
                return k;
            } else {
                continue;
            }
        }
    }
}

IntType invert(IntType t, RealType p, std::mt19937& urng = random_gen)
{
    RealType q = 1 - p;
    RealType s = p / q;
    RealType a = (t + 1) * s;
    RealType r = pow((1 - p), static_cast<RealType>(t));
    RealType u = stdunif(urng);
    IntType x = 0;
    while(u > r) {
        u = u - r;
        ++x;
        RealType r1 = ((a/x) - s) * r;
        // If r gets too small then the round-off error
        // becomes a problem.  At this point, p(i) is
        // decreasing exponentially, so if we just call
        // it 0, it's close enough.  Note that the
        // minimum value of q_n is about 1e-7, so we
        // may need to be a little careful to make sure that
        // we don't terminate the first time through the loop
        // for float.  (Hence the test that r is decreasing)
        if(r1 < std::numeric_limits<RealType>::epsilon() && r1 < r) {
            break;
        }
        r = r1;
    }
    return x;
}


IntType boost_binomial_distribution_variate(IntType t_arg, RealType p_arg, std::mt19937& urng = random_gen)
{
    bool other_side = p_arg > 0.5;
    RealType fake_p = other_side ? 1.0 - p_arg : p_arg;
    IntType m = static_cast<IntType>((t_arg+1)*fake_p);
    IntType result;
    if(m < 11)
        result = invert(t_arg, fake_p, urng);
    else
        result = btrd(t_arg, fake_p, m, urng);

    if(other_side)
        return t_arg - result;
    else
        return result;
}

}  // namespace IsoSpec
