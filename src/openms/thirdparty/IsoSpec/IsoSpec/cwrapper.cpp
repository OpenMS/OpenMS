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


#include <tuple>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "cwrapper.h"
#include "misc.h"
#include "marginalTrek++.h"
#include "isoSpec++.h"
#include "tabulator.h"

using namespace IsoSpec;

extern "C"
{
void * setupIso(int             dimNumber,
                const int*      isotopeNumbers,
                const int*      atomCounts,
                const double*   isotopeMasses,
                const double*   isotopeProbabilities)
{
    const double** IM = new const double*[dimNumber];
    const double** IP = new const double*[dimNumber];
    int idx = 0;
    for(int i=0; i<dimNumber; i++)
    {
        IM[i] = &isotopeMasses[idx];
        IP[i] = &isotopeProbabilities[idx];
        idx += isotopeNumbers[i];
    }
    //TODO in place (maybe pass a numpy matrix??)

    Iso* iso = new Iso(dimNumber, isotopeNumbers, atomCounts, IM, IP);

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}

void deleteIso(void* iso)
{
    delete reinterpret_cast<Iso*>(iso);
}


#define ISOSPEC_C_FN_CODE(generatorType, dataType, method)\
dataType method##generatorType(void* generator){ return reinterpret_cast<generatorType*>(generator)->method(); }

#define ISOSPEC_C_FN_CODE_GET_CONF_SIGNATURE(generatorType)\
void get_conf_signature##generatorType(void* generator, int* space)\
{ reinterpret_cast<generatorType*>(generator)->get_conf_signature(space); }


#define ISOSPEC_C_FN_DELETE(generatorType) void delete##generatorType(void* generator){ delete reinterpret_cast<generatorType*>(generator); }

#define ISOSPEC_C_FN_CODES(generatorType)\
ISOSPEC_C_FN_CODE(generatorType, double, mass) \
ISOSPEC_C_FN_CODE(generatorType, double, lprob) \
ISOSPEC_C_FN_CODE(generatorType, double, prob) \
ISOSPEC_C_FN_CODE_GET_CONF_SIGNATURE(generatorType) \
ISOSPEC_C_FN_CODE(generatorType, bool, advanceToNextConfiguration) \
ISOSPEC_C_FN_DELETE(generatorType)



//______________________________________________________THRESHOLD GENERATOR
void* setupIsoThresholdGenerator(void* iso,
                                 double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize)
{
    IsoThresholdGenerator* iso_tmp = new IsoThresholdGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        threshold,
        _absolute,
        _tabSize,
        _hashSize);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoThresholdGenerator)


//______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(void* iso,
                     double _target_coverage,
                     double _percentage_to_expand,
                     int _tabSize,
                     int _hashSize,
                     bool _do_trim)
{
    IsoLayeredGenerator* iso_tmp = new IsoLayeredGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        _target_coverage,
        _percentage_to_expand,
        _tabSize,
        _hashSize,
        _do_trim);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoLayeredGenerator)


//______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(void* iso,
                               int _tabSize,
                               int _hashSize)
{
    IsoOrderedGenerator* iso_tmp = new IsoOrderedGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        _tabSize,
        _hashSize);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoOrderedGenerator)

//______________________________________________________ Threshold Tabulator

void* setupThresholdTabulator(void* generator,
                     bool  get_masses,
                     bool  get_probs,
                     bool  get_lprobs,
                     bool  get_confs)
{
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(reinterpret_cast<IsoThresholdGenerator*>(generator),
                                         get_masses,
                                         get_probs,
                                         get_lprobs,
                                         get_confs);

    return reinterpret_cast<void*>(tabulator);
}

void deleteThresholdTabulator(void* t)
{
    delete reinterpret_cast<Tabulator<IsoThresholdGenerator>*>(t);
}

const double* massesThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoThresholdGenerator>*>(tabulator)->masses(true);
}

const double* lprobsThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoThresholdGenerator>*>(tabulator)->lprobs(true);
}

const double* probsThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoThresholdGenerator>*>(tabulator)->probs(true);
}

const int*    confsThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoThresholdGenerator>*>(tabulator)->confs(true);
}

int confs_noThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoThresholdGenerator>*>(tabulator)->confs_no();
}


//______________________________________________________ Layered Tabulator

void* setupLayeredTabulator(void* generator,
                     bool  get_masses,
                     bool  get_probs,
                     bool  get_lprobs,
                     bool  get_confs)
{
    Tabulator<IsoLayeredGenerator>* tabulator = new Tabulator<IsoLayeredGenerator>(reinterpret_cast<IsoLayeredGenerator*>(generator),
                                         get_masses,
                                         get_probs,
                                         get_lprobs,
                                         get_confs);

    return reinterpret_cast<void*>(tabulator);
}

void deleteLayeredTabulator(void* t)
{
    delete reinterpret_cast<Tabulator<IsoLayeredGenerator>*>(t);
}

const double* massesLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoLayeredGenerator>*>(tabulator)->masses(true);
}

const double* lprobsLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoLayeredGenerator>*>(tabulator)->lprobs(true);
}

const double* probsLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoLayeredGenerator>*>(tabulator)->probs(true);
}

const int*    confsLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoLayeredGenerator>*>(tabulator)->confs(true);
}

int confs_noLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<Tabulator<IsoLayeredGenerator>*>(tabulator)->confs_no();
}

void freeReleasedArray(void* array)
{
    free(array);
}

}  //extern "C" ends here
