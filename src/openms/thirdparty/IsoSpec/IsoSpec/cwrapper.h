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

void * isoFromFasta(const char* fasta, bool use_nominal_masses, bool add_water);

double getLightestPeakMassIso(void* iso);
double getHeaviestPeakMassIso(void* iso);
double getMonoisotopicPeakMassIso(void* iso);
double getModeLProbIso(void* iso);
double getModeMassIso(void* iso);
double getTheoreticalAverageMassIso(void* iso);
double getIsoVariance(void* iso);
double getIsoStddev(void* iso);
double* getMarginalLogSizeEstimates(void* iso, double target_total_prob);


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




// ______________________________________________________THRESHOLD GENERATOR
void* setupIsoThresholdGenerator(void* iso,
                                 double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize,
                                 bool reorder_marginals);
ISOSPEC_C_FN_HEADERS(IsoThresholdGenerator)


// ______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(void* iso,
                               int _tabSize,
                               int _hashSize,
                               bool reorder_marginals,
                               double t_prob_hint);
ISOSPEC_C_FN_HEADERS(IsoLayeredGenerator)

// ______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(void* iso,
                               int _tabSize,
                               int _hashSize);
ISOSPEC_C_FN_HEADERS(IsoOrderedGenerator)

void* setupIsoStrochasticGenerator(void* iso,
                                   size_t no_molecules,
                                   double precision,
                                   double beta_bias);
ISOSPEC_C_FN_HEADERS(IsoStochasticGenerator)


void* setupThresholdFixedEnvelope(void* iso,
                              double threshold,
                              bool absolute,
                              bool get_confs);

void* setupTotalProbFixedEnvelope(void* iso,
                              double taget_coverage,
                              bool optimize,
                              bool get_confs);

void freeReleasedArray(void* array);

void* setupFixedEnvelope(double* masses, double* probs, size_t size, bool mass_sorted, bool prob_sorted, double total_prob);
void deleteFixedEnvelope(void* tabulator, bool releaseEverything);

const double* massesFixedEnvelope(void* tabulator);
const double* probsFixedEnvelope(void* tabulator);
const int*    confsFixedEnvelope(void* tabulator);
int confs_noFixedEnvelope(void* tabulator);

double wassersteinDistance(void* tabulator1, void* tabulator2);
double orientedWassersteinDistance(void* tabulator1, void* tabulator2);
void* addEnvelopes(void* tabulator1, void* tabulator2);
void* convolveEnvelopes(void* tabulator1, void* tabulator2);

double getTotalProbOfEnvelope(void* envelope);
void scaleEnvelope(void* envelope, double factor);
void normalizeEnvelope(void* envelope);
void* binnedEnvelope(void* envelope, double width, double middle);
void* linearCombination(void* const * const envelopes, const double* intensities, size_t count);

void sortEnvelopeByMass(void* envelope);
void sortEnvelopeByProb(void* envelope);

void parse_fasta_c(const char* fasta, int atomCounts[6]);


#ifdef __cplusplus
}
#endif
