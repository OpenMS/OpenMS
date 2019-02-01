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

namespace IsoSpec
{

#ifdef __cplusplus
extern "C" {
#endif


#define ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES 288

extern const int elem_table_atomicNo[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const double elem_table_probability[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const double elem_table_mass[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const int elem_table_massNo[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const int elem_table_extraNeutrons[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const char* elem_table_element[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const char* elem_table_symbol[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const bool elem_table_Radioactive[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];
extern const double elem_table_log_probability[ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES];


#ifdef __cplusplus
}
#endif

} // namespace IsoSpec

