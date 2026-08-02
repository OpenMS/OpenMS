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

#include "operators.h"
#include "marginalTrek++.h"

namespace IsoSpec
{

KeyHasher::KeyHasher(int _dim)
: dim(_dim-1)
{}

ConfEqual::ConfEqual(int dim)
: size( dim*sizeof(int) )
{}

ConfOrderMarginal::ConfOrderMarginal(const double* _logProbs, int _dim)
: logProbs(_logProbs), dim(_dim)
{}

ConfOrderMarginalDescending::ConfOrderMarginalDescending(const double* _logProbs, int _dim)
: logProbs(_logProbs), dim(_dim)
{}
}  // namespace IsoSpec
