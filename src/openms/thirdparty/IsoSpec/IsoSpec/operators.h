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

#include <string.h>
#include "conf.h"
#include "isoMath.h"
#include "misc.h"

namespace IsoSpec
{

class KeyHasher
{
private:
    int dim;
public:
    KeyHasher(int dim);

    inline std::size_t operator()(const int* conf) const
    {
        // Following Boost...
        std::size_t seed = 0;
        for(int i = 0; i < dim; ++i )
            seed ^= conf[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    };
};


class ConfEqual
{
private:
    int size;
public:
    ConfEqual(int dim);

    inline bool operator()(const int* conf1, const int* conf2) const
    {
        // The memcmp() function returns zero if the two strings are identical, oth-
        // erwise returns the difference between the first two differing bytes
        // (treated as unsigned char values, so that `\200' is greater than `\0',
        // for example).  Zero-length strings are always identical.  This behavior
        // is not required by C and portable code should only depend on the sign of
        // the returned value.
        //                                          sacred man of memcmp.
        return memcmp(conf1, conf2, size) == 0;
    }
};


class ConfOrder
{
//configurations comparator
public:
    inline bool operator()(void* conf1,void* conf2) const
    {
        return *reinterpret_cast<double*>(conf1) < *reinterpret_cast<double*>(conf2);
    };
};



class ConfOrderMarginal
{
//configurations comparator
    const double*  logProbs;
    int dim;
public:
    ConfOrderMarginal(const double* logProbs, int dim);

    inline bool operator()(const Conf conf1, const Conf conf2)
    {// Return true if conf1 is less probable than conf2.
        return unnormalized_logProb(conf1,logProbs,dim) < unnormalized_logProb(conf2,logProbs,dim);
    };
};

class ConfOrderMarginalDescending
{
//configurations comparator
    const double*  logProbs;
    int dim;
public:
    ConfOrderMarginalDescending(const double* logProbs, int dim);

    inline bool operator()(const Conf conf1, const Conf conf2)
    {// Return true if conf1 is less probable than conf2.
        return unnormalized_logProb(conf1,logProbs,dim) > unnormalized_logProb(conf2,logProbs,dim);
    };
};

template<typename T> class ReverseOrder
{
public:
    inline ReverseOrder() {};
    inline bool operator()(const T a,const T b) const { return a > b; };
};

template<typename T> class TableOrder
{
	const T* tbl;
public:
	inline TableOrder(const T* _tbl) : tbl(_tbl) {};
	inline bool operator()(unsigned int i, unsigned int j) { return tbl[i] < tbl[j]; };
};

} // namespace IsoSpec

#include "marginalTrek++.h"

class PrecalculatedMarginal; // In case marginalTrek++.h us including us, and can't be included again...

namespace IsoSpec
{

class OrderMarginalsBySizeDecresing
{
    PrecalculatedMarginal const* const* const T;
public:
    OrderMarginalsBySizeDecresing(PrecalculatedMarginal const* const * const _T);
    bool operator()(int m1, int m2);
};


} // namespace IsoSpec



