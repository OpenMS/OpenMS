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


#include "isoSpec++.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <queue>
#include <utility>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <limits>
#include <memory>
#include <cassert>
#include <cctype>
#include "platform.h"
#include "conf.h"
#include "dirtyAllocator.h"
#include "operators.h"
#include "summator.h"
#include "marginalTrek++.h"
#include "misc.h"
#include "element_tables.h"
#include "fasta.h"



namespace IsoSpec
{

Iso::Iso() :
disowned(false),
dimNumber(0),
isotopeNumbers(new int[0]),
atomCounts(new int[0]),
confSize(0),
allDim(0),
marginals(new Marginal*[0])
{}


Iso::Iso(
    int             _dimNumber,
    const int*      _isotopeNumbers,
    const int*      _atomCounts,
    const double* const *  _isotopeMasses,
    const double* const *  _isotopeProbabilities
) :
disowned(false),
dimNumber(_dimNumber),
isotopeNumbers(array_copy<int>(_isotopeNumbers, _dimNumber)),
atomCounts(array_copy<int>(_atomCounts, _dimNumber)),
confSize(_dimNumber * sizeof(int)),
allDim(0),
marginals(nullptr)
{
    for(int ii = 0; ii < dimNumber; ++ii)
        allDim += isotopeNumbers[ii];

    std::unique_ptr<double[]> masses(new double[allDim]);
    std::unique_ptr<double[]> probs(new double[allDim]);
    size_t idx = 0;

    for(int ii = 0; ii < dimNumber; ++ii)
        for(int jj = 0; jj < isotopeNumbers[ii]; ++jj)
        {
            masses[idx] = _isotopeMasses[ii][jj];
            probs[idx]  = _isotopeProbabilities[ii][jj];
            ++idx;
        }

    allDim = 0;  // setupMarginals will recalculate it, assuming it's set to 0

    try{
        setupMarginals(masses.get(), probs.get());
    }
    catch(...)
    {
        delete[] isotopeNumbers;
        delete[] atomCounts;
	// Since we're throwing in a constructor, the destructor won't run, and we don't need to NULL these.
	// However, this is not the fast code path and we can afford two unneeded instructions to keep
	// some static analysis tools happy.
	isotopeNumbers = nullptr;
	atomCounts = nullptr;
        throw;
    }
}

Iso::Iso(
    int             _dimNumber,
    const int*      _isotopeNumbers,
    const int*      _atomCounts,
    const double*   _isotopeMasses,
    const double*   _isotopeProbabilities
) :
disowned(false),
dimNumber(_dimNumber),
isotopeNumbers(array_copy<int>(_isotopeNumbers, _dimNumber)),
atomCounts(array_copy<int>(_atomCounts, _dimNumber)),
confSize(_dimNumber * sizeof(int)),
allDim(0),
marginals(nullptr)
{
    try{
        setupMarginals(_isotopeMasses, _isotopeProbabilities);
    }
    catch(...)
    {
        delete[] isotopeNumbers;
        delete[] atomCounts;
	// Since we're throwing in a constructor, the destructor won't run, and we don't need to NULL these.
	// However, this is not the fast code path and we can afford two unneeded instructions to keep
	// some static analysis tools happy.
	isotopeNumbers = nullptr;
	atomCounts = nullptr;
        throw;
    }
}

Iso::Iso(Iso&& other) :
disowned(other.disowned),
dimNumber(other.dimNumber),
isotopeNumbers(other.isotopeNumbers),
atomCounts(other.atomCounts),
confSize(other.confSize),
allDim(other.allDim),
marginals(other.marginals)
{
    other.disowned = true;
}


Iso::Iso(const Iso& other, bool fullcopy) :
disowned(!fullcopy),
dimNumber(other.dimNumber),
isotopeNumbers(fullcopy ? array_copy<int>(other.isotopeNumbers, dimNumber) : other.isotopeNumbers),
atomCounts(fullcopy ? array_copy<int>(other.atomCounts, dimNumber) : other.atomCounts),
confSize(other.confSize),
allDim(other.allDim),
marginals(fullcopy ? new Marginal*[dimNumber] : other.marginals)
{
    if(fullcopy)
    {
        for(int ii = 0; ii < dimNumber; ii++)
            marginals[ii] = new Marginal(*other.marginals[ii]);
    }
}

Iso Iso::FromFASTA(const char* fasta, bool use_nominal_masses, bool add_water)
{
    int atomCounts[6];

    parse_fasta(fasta, atomCounts);

    if(add_water)
    {
        atomCounts[1] += 2;
        atomCounts[3] += 1;
    }

    const int dimNr = atomCounts[5] > 0 ? 6 : 5;

    return Iso(dimNr, aa_isotope_numbers, atomCounts, use_nominal_masses ? aa_elem_nominal_masses : aa_elem_masses, aa_elem_probabilities);
}

inline void Iso::setupMarginals(const double* _isotopeMasses, const double* _isotopeProbabilities)
{
    if (marginals == nullptr)
    {
        int ii = 0;
        marginals = new Marginal*[dimNumber];
        try
        {
            while(ii < dimNumber)
            {
                marginals[ii] = new Marginal(
                        &_isotopeMasses[allDim],
                        &_isotopeProbabilities[allDim],
                        isotopeNumbers[ii],
                        atomCounts[ii]
                    );
                allDim += isotopeNumbers[ii];
                ii++;
            }
        }
        catch(...)
        {
            ii--;
            while(ii >= 0)
            {
                delete marginals[ii];
                ii--;
            }
            delete[] marginals;
            marginals = nullptr;
            throw;
        }
    }
}

Iso::~Iso()
{
    if(!disowned)
    {
    if (marginals != nullptr)
        dealloc_table(marginals, dimNumber);
    delete[] isotopeNumbers;
    delete[] atomCounts;
    }
}

bool Iso::doMarginalsNeedSorting() const
{
    int nontrivial_marginals = 0;
    for(int ii = 0; ii < dimNumber; ii++)
    {
        if(marginals[ii]->get_isotopeNo() > 1)
            nontrivial_marginals++;
        if(nontrivial_marginals > 1)
            return true;
    }
    return false;
}

double Iso::getLightestPeakMass() const
{
    double mass = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        mass += marginals[ii]->getLightestConfMass();
    return mass;
}

double Iso::getHeaviestPeakMass() const
{
    double mass = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        mass += marginals[ii]->getHeaviestConfMass();
    return mass;
}

double Iso::getMonoisotopicPeakMass() const
{
    double mass = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        mass += marginals[ii]->getMonoisotopicConfMass();
    return mass;
}

double Iso::getUnlikeliestPeakLProb() const
{
    double lprob = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        lprob += marginals[ii]->getSmallestLProb();
    return lprob;
}

double Iso::getModeMass() const
{
    double mass = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        mass += marginals[ii]->getModeMass();
    return mass;
}

double Iso::getModeLProb() const
{
    double ret = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        ret += marginals[ii]->getModeLProb();
    return ret;
}

double Iso::getTheoreticalAverageMass() const
{
    double mass = 0.0;
    for (int ii = 0; ii < dimNumber; ii++)
        mass += marginals[ii]->getTheoreticalAverageMass();
    return mass;
}

double Iso::variance() const
{
    double ret = 0.0;
    for(int ii = 0; ii < dimNumber; ii++)
        ret += marginals[ii]->variance();
    return ret;
}


Iso::Iso(const char* formula, bool use_nominal_masses) :
disowned(false),
allDim(0),
marginals(nullptr)
{
    std::vector<double> isotope_masses;
    std::vector<double> isotope_probabilities;

    dimNumber = parse_formula(formula, isotope_masses, isotope_probabilities, &isotopeNumbers, &atomCounts, &confSize, use_nominal_masses);

    setupMarginals(isotope_masses.data(), isotope_probabilities.data());
}


void Iso::addElement(int atomCount, int noIsotopes, const double* isotopeMasses, const double* isotopeProbabilities)
{
    Marginal* m = new Marginal(isotopeMasses, isotopeProbabilities, noIsotopes, atomCount);
    realloc_append<int>(&isotopeNumbers, noIsotopes, dimNumber);
    realloc_append<int>(&atomCounts, atomCount, dimNumber);
    realloc_append<Marginal*>(&marginals, m, dimNumber);
    dimNumber++;
    confSize += sizeof(int);
    allDim += noIsotopes;
}

void Iso::saveMarginalLogSizeEstimates(double* priorities, double target_total_prob) const
{
    /*
     * We shall now use Gaussian approximations of the marginal multinomial distributions to estimate
     * how many configurations we shall need to visit from each marginal. This should be approximately
     * proportional to the volume of the optimal P-ellipsoid of the marginal, which, in turn is defined
     * by the quantile function of the chi-square distribution plus some modifications.
     *
     * We're dropping the constant factor and the (monotonic) exp() transform - these will be used as keys
     * for sorting, so only the ordering is important.
     */

    double K = allDim - dimNumber;

    double log_R2 = log(InverseChiSquareCDF2(K, target_total_prob));

    for(int ii = 0; ii < dimNumber; ii++)
        priorities[ii] = marginals[ii]->getLogSizeEstimate(log_R2);
}

unsigned int parse_formula(const char* formula, std::vector<double>& isotope_masses, std::vector<double>& isotope_probabilities, int** isotopeNumbers, int** atomCounts, unsigned int* confSize, bool use_nominal_masses)
{
    // This function is NOT guaranteed to be secure against malicious input. It should be used only for debugging.
    size_t slen = strlen(formula);
    // Yes, it would be more elegant to use std::string here, but it's the only promiment place where it would be used in IsoSpec, and avoiding it here
    // means we can run the whole thing through Clang's memory sanitizer without the need for instrumented libc++/libstdc++. That's worth messing with char pointers a
    // little bit.
    std::vector<std::pair<const char*, size_t> > elements;
    std::vector<int> numbers;

    if(slen == 0)
        throw std::invalid_argument("Invalid formula: can't be empty");

    if(!isdigit(formula[slen-1]))
        throw std::invalid_argument("Invalid formula: every element must be followed by a number - write H2O1 and not H2O for water");

    for(size_t ii = 0; ii < slen; ii++)
        if(!isdigit(formula[ii]) && !isalpha(formula[ii]))
            throw std::invalid_argument("Invalid formula: contains invalid (non-digit, non-alpha) character");

    size_t position = 0;

    while(position < slen)
    {
        size_t elem_end = position;
        while(isalpha(formula[elem_end]))
            elem_end++;
        size_t digit_end = elem_end;
        while(isdigit(formula[digit_end]))
            digit_end++;
        elements.emplace_back(&formula[position], elem_end-position);
        numbers.push_back(std::stoi(&formula[elem_end]));
        position = digit_end;
    }

    std::vector<int> element_indexes;

    for (unsigned int i = 0; i < elements.size(); i++)
    {
        int idx = -1;
        for(int j = 0; j < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; j++)
        {
            if ((strlen(elem_table_symbol[j]) == elements[i].second) && (strncmp(elements[i].first, elem_table_symbol[j], elements[i].second) == 0))
            {
                idx = j;
                break;
            }
        }
        if(idx < 0)
            throw std::invalid_argument("Invalid formula");
        element_indexes.push_back(idx);
    }

    std::vector<int> _isotope_numbers;
    const double* masses = use_nominal_masses ? elem_table_massNo : elem_table_mass;

    for(std::vector<int>::iterator it = element_indexes.begin(); it != element_indexes.end(); ++it)
    {
        int num = 0;
        int at_idx = *it;
        int elem_ID = elem_table_ID[at_idx];
        while(at_idx < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_ID[at_idx] == elem_ID)
        {
            isotope_masses.push_back(masses[at_idx]);
            isotope_probabilities.push_back(elem_table_probability[at_idx]);
            at_idx++;
            num++;
        }
        _isotope_numbers.push_back(num);
    }

    const unsigned int dimNumber = elements.size();

    *isotopeNumbers = array_copy<int>(_isotope_numbers.data(), dimNumber);
    *atomCounts = array_copy<int>(numbers.data(), dimNumber);
    *confSize = dimNumber * sizeof(int);

    return dimNumber;
}


/*
 * ----------------------------------------------------------------------------------------------------------
 */



IsoGenerator::IsoGenerator(Iso&& iso, bool alloc_partials) :
    Iso(std::move(iso)),
    mode_lprob(getModeLProb()),
    partialLProbs(alloc_partials ? new double[dimNumber+1] : nullptr),
    partialMasses(alloc_partials ? new double[dimNumber+1] : nullptr),
    partialProbs(alloc_partials ? new double[dimNumber+1] : nullptr)
{
    for(int ii = 0; ii < dimNumber; ++ii)
        marginals[ii]->ensureModeConf();
    if(alloc_partials)
    {
        partialLProbs[dimNumber] = 0.0;
        partialMasses[dimNumber] = 0.0;
        partialProbs[dimNumber] = 1.0;
    }
}


IsoGenerator::~IsoGenerator()
{
    if(partialLProbs != nullptr)
        delete[] partialLProbs;
    if(partialMasses != nullptr)
        delete[] partialMasses;
    if(partialProbs != nullptr)
        delete[] partialProbs;
}



/*
 * ----------------------------------------------------------------------------------------------------------
 */

static const double minsqrt = -1.3407796239501852e+154;  // == constexpr(-sqrt(std::numeric_limits<double>::max()));

IsoThresholdGenerator::IsoThresholdGenerator(Iso&& iso, double _threshold, bool _absolute, int tabSize, int hashSize, bool reorder_marginals)
: IsoGenerator(std::move(iso)),
Lcutoff(_threshold <= 0.0 ? minsqrt : (_absolute ? log(_threshold) : log(_threshold) + mode_lprob))
{
    counter = new int[dimNumber];
    maxConfsLPSum = new double[dimNumber-1];
    marginalResultsUnsorted = new PrecalculatedMarginal*[dimNumber];

    empty = false;

    const bool marginalsNeedSorting = doMarginalsNeedSorting();

    for(int ii = 0; ii < dimNumber; ii++)
    {
        counter[ii] = 0;
        marginalResultsUnsorted[ii] = new PrecalculatedMarginal(std::move(*(marginals[ii])),
                                                        Lcutoff - mode_lprob + marginals[ii]->fastGetModeLProb(),
                                                        marginalsNeedSorting,
                                                        tabSize,
                                                        hashSize);

        if(!marginalResultsUnsorted[ii]->inRange(0))
            empty = true;
    }

    if(reorder_marginals && dimNumber > 1)
    {
        OrderMarginalsBySizeDecresing<PrecalculatedMarginal> comparator(marginalResultsUnsorted);
        int* tmpMarginalOrder = new int[dimNumber];

        for(int ii = 0; ii < dimNumber; ii++)
            tmpMarginalOrder[ii] = ii;

        std::sort(tmpMarginalOrder, tmpMarginalOrder + dimNumber, comparator);
        marginalResults = new PrecalculatedMarginal*[dimNumber];

        for(int ii = 0; ii < dimNumber; ii++)
            marginalResults[ii] = marginalResultsUnsorted[tmpMarginalOrder[ii]];

        marginalOrder = new int[dimNumber];
        for(int ii = 0; ii < dimNumber; ii++)
            marginalOrder[tmpMarginalOrder[ii]] = ii;

        delete[] tmpMarginalOrder;
    }
    else
    {
        marginalResults = marginalResultsUnsorted;
        marginalOrder = nullptr;
    }

    lProbs_ptr_start = marginalResults[0]->get_lProbs_ptr();

    if(dimNumber > 1)
        maxConfsLPSum[0] = marginalResults[0]->fastGetModeLProb();

    for(int ii = 1; ii < dimNumber-1; ii++)
        maxConfsLPSum[ii] = maxConfsLPSum[ii-1] + marginalResults[ii]->fastGetModeLProb();

    lProbs_ptr = lProbs_ptr_start;

    partialLProbs_second = partialLProbs;
    partialLProbs_second++;

    if(!empty)
    {
        recalc(dimNumber-1);
        counter[0]--;
        lProbs_ptr--;
    }
    else
    {
        terminate_search();
        lcfmsv = std::numeric_limits<double>::infinity();
    }
}

void IsoThresholdGenerator::terminate_search()
{
    for(int ii = 0; ii < dimNumber; ii++)
    {
        counter[ii] = marginalResults[ii]->get_no_confs()-1;
        partialLProbs[ii] = -std::numeric_limits<double>::infinity();
    }
    partialLProbs[dimNumber] = -std::numeric_limits<double>::infinity();
    lProbs_ptr = lProbs_ptr_start + marginalResults[0]->get_no_confs()-1;
}

size_t IsoThresholdGenerator::count_confs()
{
    if(empty)
        return 0;

    if(dimNumber == 1)
        return marginalResults[0]->get_no_confs();

    const double* lProbs_ptr_l = marginalResults[0]->get_lProbs_ptr() + marginalResults[0]->get_no_confs();

    std::unique_ptr<const double* []> lProbs_restarts(new const double*[dimNumber]);

    for(int ii = 0; ii < dimNumber; ii++)
        lProbs_restarts[ii] = lProbs_ptr_l;

    size_t count = 0;

    while(*lProbs_ptr_l < lcfmsv)
        lProbs_ptr_l--;

    while(true)
    {
        count += lProbs_ptr_l - lProbs_ptr_start + 1;

        int idx = 0;
        int * cntr_ptr = counter;

        while(idx < dimNumber - 1)
        {
            *cntr_ptr = 0;
            idx++;
            cntr_ptr++;
            (*cntr_ptr)++;

            partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
            if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= Lcutoff)
            {
                short_recalc(idx-1);
                lProbs_ptr_l = lProbs_restarts[idx];
                while(*lProbs_ptr_l < lcfmsv)
                    lProbs_ptr_l--;
                for(idx--; idx > 0; idx--)
                    lProbs_restarts[idx] = lProbs_ptr_l;
                break;
            }
        }
        if(idx == dimNumber - 1)
        {
            reset();
            return count;
        }
    }
}

void IsoThresholdGenerator::reset()
{
    if(empty)
    {
        terminate_search();
        return;
    }

    partialLProbs[dimNumber] = 0.0;

    memset(counter, 0, sizeof(int)*dimNumber);
    recalc(dimNumber-1);
    counter[0]--;

    lProbs_ptr = lProbs_ptr_start - 1;
}

IsoThresholdGenerator::~IsoThresholdGenerator()
{
    delete[] counter;
    delete[] maxConfsLPSum;
    if (marginalResultsUnsorted != marginalResults)
        delete[] marginalResultsUnsorted;
    dealloc_table(marginalResults, dimNumber);
    if(marginalOrder != nullptr)
        delete[] marginalOrder;
}


/*
 * ------------------------------------------------------------------------------------------------------------------------
 */


IsoLayeredGenerator::IsoLayeredGenerator(Iso&& iso, int tabSize, int hashSize, bool reorder_marginals, double t_prob_hint)
: IsoGenerator(std::move(iso))
{
    counter = new int[dimNumber];
    maxConfsLPSum = new double[dimNumber-1];
    currentLThreshold = nextafter(mode_lprob, -std::numeric_limits<double>::infinity());
    lastLThreshold = (std::numeric_limits<double>::min)();
    marginalResultsUnsorted = new LayeredMarginal*[dimNumber];
    resetPositions = new const double*[dimNumber];
    marginalsNeedSorting = doMarginalsNeedSorting();

    memset(counter, 0, sizeof(int)*dimNumber);

    for(int ii = 0; ii < dimNumber; ii++)
        marginalResultsUnsorted[ii] = new LayeredMarginal(std::move(*(marginals[ii])), tabSize, hashSize);

    if(reorder_marginals && dimNumber > 1)
    {
        double* marginal_priorities = new double[dimNumber];

        saveMarginalLogSizeEstimates(marginal_priorities, t_prob_hint);

        int* tmpMarginalOrder = new int[dimNumber];

        for(int ii = 0; ii < dimNumber; ii++)
            tmpMarginalOrder[ii] = ii;

        TableOrder<double> TO(marginal_priorities);

        std::sort(tmpMarginalOrder, tmpMarginalOrder + dimNumber, TO);
        marginalResults = new LayeredMarginal*[dimNumber];

        for(int ii = 0; ii < dimNumber; ii++)
            marginalResults[ii] = marginalResultsUnsorted[tmpMarginalOrder[ii]];

        marginalOrder = new int[dimNumber];
        for(int ii = 0; ii < dimNumber; ii++)
            marginalOrder[tmpMarginalOrder[ii]] = ii;

        delete[] tmpMarginalOrder;
        delete[] marginal_priorities;
    }
    else
    {
        marginalResults = marginalResultsUnsorted;
        marginalOrder = nullptr;
    }

    lProbs_ptr_start = marginalResults[0]->get_lProbs_ptr();

    if(dimNumber > 1)
        maxConfsLPSum[0] = marginalResults[0]->fastGetModeLProb();

    for(int ii = 1; ii < dimNumber-1; ii++)
        maxConfsLPSum[ii] = maxConfsLPSum[ii-1] + marginalResults[ii]->fastGetModeLProb();

    lProbs_ptr = lProbs_ptr_start;

    partialLProbs_second = partialLProbs;
    partialLProbs_second++;

    counter[0]--;
    lProbs_ptr--;
    lastLThreshold = 10.0;
    IsoLayeredGenerator::nextLayer(-0.00001);
}

bool IsoLayeredGenerator::nextLayer(double offset)
{
    size_t first_mrg_size = marginalResults[0]->get_no_confs();

    if(lastLThreshold < getUnlikeliestPeakLProb())
        return false;

    lastLThreshold = currentLThreshold;
    currentLThreshold += offset;

    for(int ii = 0; ii < dimNumber; ii++)
    {
        marginalResults[ii]->extend(currentLThreshold - mode_lprob + marginalResults[ii]->fastGetModeLProb(), marginalsNeedSorting);
        counter[ii] = 0;
    }

    lProbs_ptr_start = marginalResults[0]->get_lProbs_ptr();  // vector relocation might have happened

    lProbs_ptr = lProbs_ptr_start + first_mrg_size - 1;

    for(int ii = 0; ii < dimNumber; ii++)
        resetPositions[ii] = lProbs_ptr;

    recalc(dimNumber-1);

    return true;
}

bool IsoLayeredGenerator::carry()
{
    // If we reached this point, a carry is needed

    int idx = 0;

    int * cntr_ptr = counter;

    while(idx < dimNumber-1)
    {
        *cntr_ptr = 0;
        idx++;
        cntr_ptr++;
        (*cntr_ptr)++;
        partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
        if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= currentLThreshold)
        {
            partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
            partialProbs[idx] = partialProbs[idx+1] * marginalResults[idx]->get_prob(counter[idx]);
            recalc(idx-1);
            lProbs_ptr = resetPositions[idx];

            while(*lProbs_ptr <= last_lcfmsv)
                lProbs_ptr--;

            for(int ii = 0; ii < idx; ii++)
                resetPositions[ii] = lProbs_ptr;

            return true;
        }
    }

    return false;
}


void IsoLayeredGenerator::terminate_search()
{
    for(int ii = 0; ii < dimNumber; ii++)
    {
        counter[ii] = marginalResults[ii]->get_no_confs()-1;
        partialLProbs[ii] = -std::numeric_limits<double>::infinity();
    }
    partialLProbs[dimNumber] = -std::numeric_limits<double>::infinity();
    lProbs_ptr = lProbs_ptr_start + marginalResults[0]->get_no_confs()-1;
}

IsoLayeredGenerator::~IsoLayeredGenerator()
{
    delete[] counter;
    delete[] maxConfsLPSum;
    delete[] resetPositions;
    if (marginalResultsUnsorted != marginalResults)
        delete[] marginalResultsUnsorted;
    dealloc_table(marginalResults, dimNumber);
    if(marginalOrder != nullptr)
      delete[] marginalOrder;
}


/*
 * ------------------------------------------------------------------------------------------------------------------------
 */


IsoOrderedGenerator::IsoOrderedGenerator(Iso&& iso, int _tabSize, int _hashSize) :
IsoGenerator(std::move(iso), false), allocator(dimNumber, _tabSize)
{
    partialLProbs = &currentLProb;
    partialMasses = &currentMass;
    partialProbs = &currentProb;

    marginalResults = new MarginalTrek*[dimNumber];

    for(int i = 0; i < dimNumber; i++)
        marginalResults[i] = new MarginalTrek(std::move(*(marginals[i])), _tabSize, _hashSize);

    logProbs        = new const std::vector<double>*[dimNumber];
    masses          = new const std::vector<double>*[dimNumber];
    marginalConfs   = new const std::vector<int*>*[dimNumber];

    for(int i = 0; i < dimNumber; i++)
    {
        masses[i] = &marginalResults[i]->conf_masses();
        logProbs[i] = &marginalResults[i]->conf_lprobs();
        marginalConfs[i] = &marginalResults[i]->confs();
    }

    topConf = allocator.newConf();
    memset(
            reinterpret_cast<char*>(topConf) + sizeof(double),
            0,
            sizeof(int)*dimNumber
    );

    *(reinterpret_cast<double*>(topConf)) =
    combinedSum(
                getConf(topConf),
                logProbs,
                dimNumber
    );

    pq.push(topConf);
}


IsoOrderedGenerator::~IsoOrderedGenerator()
{
    dealloc_table<MarginalTrek*>(marginalResults, dimNumber);
    delete[] logProbs;
    delete[] masses;
    delete[] marginalConfs;
    partialLProbs = nullptr;
    partialMasses = nullptr;
    partialProbs = nullptr;
}


bool IsoOrderedGenerator::advanceToNextConfiguration()
{
    if(pq.size() < 1)
        return false;


    topConf = pq.top();
    pq.pop();

    int* topConfIsoCounts = getConf(topConf);

    currentLProb = *(reinterpret_cast<double*>(topConf));
    currentMass = combinedSum( topConfIsoCounts, masses, dimNumber );
    currentProb = exp(currentLProb);

    ccount = -1;
    for(int j = 0; j < dimNumber; ++j)
    {
        if(marginalResults[j]->probeConfigurationIdx(topConfIsoCounts[j] + 1))
        {
            if(ccount == -1)
            {
                topConfIsoCounts[j]++;
                *(reinterpret_cast<double*>(topConf)) = combinedSum(topConfIsoCounts, logProbs, dimNumber);
                pq.push(topConf);
                topConfIsoCounts[j]--;
                ccount = j;
            }
            else
            {
                void* acceptedCandidate = allocator.newConf();
                int* acceptedCandidateIsoCounts = getConf(acceptedCandidate);
                memcpy(acceptedCandidateIsoCounts, topConfIsoCounts, confSize);

                acceptedCandidateIsoCounts[j]++;

                *(reinterpret_cast<double*>(acceptedCandidate)) = combinedSum(acceptedCandidateIsoCounts, logProbs, dimNumber);

                pq.push(acceptedCandidate);
            }
        }
        if(topConfIsoCounts[j] > 0)
            break;
    }
    if(ccount >=0)
        topConfIsoCounts[ccount]++;

    return true;
}


/*
 * ---------------------------------------------------------------------------------------------------
 */


IsoStochasticGenerator::IsoStochasticGenerator(Iso&& iso, size_t no_molecules, double _precision, double _beta_bias) :
IsoGenerator(std::move(iso)),
ILG(std::move(*this)),
to_sample_left(no_molecules),
precision(_precision),
beta_bias(_beta_bias),
confs_prob(0.0),
chasing_prob(0.0)
{}

/*
 * ---------------------------------------------------------------------------------------------------
 */





}  // namespace IsoSpec

