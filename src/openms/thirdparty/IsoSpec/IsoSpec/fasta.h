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

namespace IsoSpec{

// We will work with C H N O S Se tuples */
extern const int aa_isotope_numbers[6];

extern const double aa_elem_masses[19];

extern const double aa_elem_nominal_masses[19];

extern const double aa_elem_probabilities[19];

extern const int aa_symbol_to_elem_counts[256*6];

}  // namespace IsoSpec
