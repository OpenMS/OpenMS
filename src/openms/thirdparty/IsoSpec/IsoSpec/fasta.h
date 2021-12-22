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

inline void parse_fasta(const char* fasta, int atomCounts[6])
{
    memset(atomCounts, 0, sizeof(decltype(atomCounts[0]))*6);

    for(size_t idx = 0; fasta[idx] != '\0'; ++idx)
    {
        const int* counts = &aa_symbol_to_elem_counts[fasta[idx]*6];
        for(int ii = 0; ii < 6; ++ii)
            atomCounts[ii] += counts[ii];
    }
}

}  // namespace IsoSpec
