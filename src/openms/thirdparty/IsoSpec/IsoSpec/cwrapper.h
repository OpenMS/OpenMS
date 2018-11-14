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

#define ISOSPEC_ALGO_LAYERED 0
#define ISOSPEC_ALGO_ORDERED 1
#define ISOSPEC_ALGO_THRESHOLD_ABSOLUTE 2
#define ISOSPEC_ALGO_THRESHOLD_RELATIVE 3
#define ISOSPEC_ALGO_LAYERED_ESTIMATE 4


#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

void * setupIso(int             dimNumber,
                const int*      isotopeNumbers,
                const int*      atomCounts,
                const double*   isotopeMasses,
                const double*   isotopeProbabilities);

void deleteIso(void* iso);

#define ISOSPEC_C_FN_HEADER(generatorType, dataType, method)\
dataType method##generatorType(void* generator);

#define ISOSPEC_C_FN_HEADER_GET_CONF_SIGNATURE(generatorType)\
void method##generatorType(void* generator);

#define ISOSPEC_C_FN_HEADERS(generatorType)\
ISOSPEC_C_FN_HEADER(generatorType, double, mass) \
ISOSPEC_C_FN_HEADER(generatorType, double, lprob) \
ISOSPEC_C_FN_HEADER(generatorType, double, prob) \
ISOSPEC_C_FN_HEADER_GET_CONF_SIGNATURE(generatorType) \
ISOSPEC_C_FN_HEADER(generatorType, bool, advanceToNextConfiguration) \
ISOSPEC_C_FN_HEADER(generatorType, void, delete)




//______________________________________________________THRESHOLD GENERATOR
void* setupIsoThresholdGenerator(void* iso,
                                 double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize);
ISOSPEC_C_FN_HEADERS(IsoThresholdGenerator)


//______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(void* iso,
                               double _target_coverage,
                               double _percentage_to_expand,
                               int _tabSize,
                               int _hashSize,
                               bool _do_trim);
ISOSPEC_C_FN_HEADERS(IsoLayeredGenerator)

//______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(void* iso,
                               int _tabSize,
                               int _hashSize);
ISOSPEC_C_FN_HEADERS(IsoOrderedGenerator)



void* setupThresholdTabulator(void* generator,
                              bool  get_masses,
                              bool  get_probs,
                              bool  get_lprobs,
                              bool  get_confs);

void deleteThresholdTabulator(void* tabulator);

const double* massesThresholdTabulator(void* tabulator);
const double* lprobsThresholdTabulator(void* tabulator);
const double* probsThresholdTabulator(void* tabulator);
const int*    confsThresholdTabulator(void* tabulator);
int confs_noThresholdTabulator(void* tabulator);



void* setupLayeredTabulator(void* generator,
                              bool  get_masses,
                              bool  get_probs,
                              bool  get_lprobs,
                              bool  get_confs);

void deleteLayeredTabulator(void* tabulator);

const double* massesLayeredTabulator(void* tabulator);
const double* lprobsLayeredTabulator(void* tabulator);
const double* probsLayeredTabulator(void* tabulator);
const int*    confsLayeredTabulator(void* tabulator);
int confs_noLayeredTabulator(void* tabulator);

void freeReleasedArray(void* array);

#ifdef __cplusplus
}
#endif

