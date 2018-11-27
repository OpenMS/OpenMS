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


#include <cmath>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <tuple>
#include <unordered_map>
#include <queue>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <stdexcept>
#include <string>
#include <limits>
#include <assert.h>
#include <ctype.h>
#include "platform.h"
#include "conf.h"
#include "dirtyAllocator.h"
#include "operators.h"
#include "summator.h"
#include "marginalTrek++.h"
#include "isoSpec++.h"
#include "misc.h"
#include "element_tables.h"


using namespace std;

namespace IsoSpec
{

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
marginals(nullptr),
modeLProb(0.0)
{
    try{
        setupMarginals(_isotopeMasses, _isotopeProbabilities);
    }
    catch(...)
    {
        delete[] isotopeNumbers;
        delete[] atomCounts;
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
marginals(other.marginals),
modeLProb(other.modeLProb)
{
    other.disowned = true;
}


Iso::Iso(const Iso& other, bool fullcopy) :
disowned(fullcopy ? throw std::logic_error("Not implemented") : true),
dimNumber(other.dimNumber),
isotopeNumbers(fullcopy ? array_copy<int>(other.isotopeNumbers, dimNumber) : other.isotopeNumbers),
atomCounts(fullcopy ? array_copy<int>(other.atomCounts, dimNumber) : other.atomCounts),
confSize(other.confSize),
allDim(other.allDim),
marginals(fullcopy ? throw std::logic_error("Not implemented") : other.marginals),
modeLProb(other.modeLProb)
{}


inline void Iso::setupMarginals(const double* const * _isotopeMasses, const double* const * _isotopeProbabilities)
{
    if (marginals == nullptr)
    {
        int ii = 0;
        try
        {
            marginals = new Marginal*[dimNumber];
            while(ii < dimNumber)
            {
                allDim += isotopeNumbers[ii];
                marginals[ii] = new Marginal(
                        _isotopeMasses[ii],
                        _isotopeProbabilities[ii],
                        isotopeNumbers[ii],
                        atomCounts[ii]
                    );
                modeLProb += marginals[ii]->getModeLProb();
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


double Iso::getLightestPeakMass() const
{
    double mass = 0.0;
    for (int ii=0; ii<dimNumber; ii++)
        mass += marginals[ii]->getLightestConfMass();
    return mass;
}

double Iso::getHeaviestPeakMass() const
{
    double mass = 0.0;
    for (int ii=0; ii<dimNumber; ii++)
        mass += marginals[ii]->getHeaviestConfMass();
    return mass;
}



Iso::Iso(const char* formula) :
disowned(false),
allDim(0),
marginals(nullptr),
modeLProb(0.0)
{
    std::vector<const double*> isotope_masses;
    std::vector<const double*> isotope_probabilities;

    dimNumber = parse_formula(formula, isotope_masses, isotope_probabilities, &isotopeNumbers, &atomCounts, &confSize);

    setupMarginals(isotope_masses.data(), isotope_probabilities.data());
}

unsigned int parse_formula(const char* formula, std::vector<const double*>& isotope_masses, std::vector<const double*>& isotope_probabilities, int** isotopeNumbers, int** atomCounts, unsigned int* confSize)
{
    // This function is NOT guaranteed to be secure against malicious input. It should be used only for debugging.
    size_t slen = strlen(formula);
    // Yes, it would be more elegant to use std::string here, but it's the only promiment place where it would be used in IsoSpec, and avoiding it here
    // means we can run the whole thing through Clang's memory sanitizer without the need for instrumented libc++/libstdc++. That's worth messing with char pointers a
    // little bit.
    std::vector<std::pair<const char*, size_t> > elements;
    std::vector<int> numbers;

    if(slen == 0)
        throw invalid_argument("Invalid formula: can't be empty");

    if(!isdigit(formula[slen-1]))
        throw invalid_argument("Invalid formula: every element must be followed by a number - write H2O1 and not H2O for water");

    for(size_t ii=0; ii<slen; ii++)
        if(!isdigit(formula[ii]) && !isalpha(formula[ii]))
            throw invalid_argument("Ivalid formula: contains invalid (non-digit, non-alpha) character");

    size_t position = 0;
    size_t elem_end = 0;
    size_t digit_end = 0;

    while(position < slen)
    {
        elem_end = position;
        while(isalpha(formula[elem_end]))
            elem_end++;
        digit_end = elem_end;
        while(isdigit(formula[digit_end]))
            digit_end++;
        elements.emplace_back(&formula[position], elem_end-position);
        numbers.push_back(atoi(&formula[elem_end]));
        position = digit_end;
    }

    std::vector<int> element_indexes;

    for (unsigned int i=0; i<elements.size(); i++)
    {
        int idx = -1;
        for(int j=0; j<ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; j++)
        {
            if ((strlen(elem_table_symbol[j]) == elements[i].second) && (strncmp(elements[i].first, elem_table_symbol[j], elements[i].second) == 0))
            {
                idx = j;
                break;
            }
        }
        if(idx < 0)
            throw invalid_argument("Invalid formula");
        element_indexes.push_back(idx);
    }

    vector<int> _isotope_numbers;

    for(vector<int>::iterator it = element_indexes.begin(); it != element_indexes.end(); ++it)
    {
        int num = 0;
        int at_idx = *it;
        int atomicNo = elem_table_atomicNo[at_idx];
        while(at_idx < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_atomicNo[at_idx] == atomicNo)
        {
            at_idx++;
            num++;
        }
        _isotope_numbers.push_back(num);
    }

    for(vector<int>::iterator it = element_indexes.begin(); it != element_indexes.end(); ++it)
    {
        isotope_masses.push_back(&elem_table_mass[*it]);
        isotope_probabilities.push_back(&elem_table_probability[*it]);
    };

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
    partialLProbs(alloc_partials ? new double[dimNumber+1] : nullptr),
    partialMasses(alloc_partials ? new double[dimNumber+1] : nullptr),
    partialProbs(alloc_partials ? new double[dimNumber+1] : nullptr)
{
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



IsoThresholdGenerator::IsoThresholdGenerator(Iso&& iso, double _threshold, bool _absolute, int tabSize, int hashSize, bool reorder_marginals)
: IsoGenerator(std::move(iso)),
Lcutoff(_threshold <= 0.0 ? std::numeric_limits<double>::lowest() : (_absolute ? log(_threshold) : log(_threshold) + modeLProb))
{
    counter = new int[dimNumber];
    maxConfsLPSum = new double[dimNumber-1];
    marginalResultsUnsorted = new PrecalculatedMarginal*[dimNumber];

    empty = false;

    for(int ii=0; ii<dimNumber; ii++)
    {
        counter[ii] = 0;
        marginalResultsUnsorted[ii] = new PrecalculatedMarginal(std::move(*(marginals[ii])),
                                                        Lcutoff - modeLProb + marginals[ii]->getModeLProb(),
                                                        true,
                                                        tabSize,
                                                        hashSize);

        if(!marginalResultsUnsorted[ii]->inRange(0))
            empty = true;
    }

    if(reorder_marginals)
    {
        OrderMarginalsBySizeDecresing comparator(marginalResultsUnsorted);
        int* tmpMarginalOrder = new int[dimNumber];

        for(int ii=0; ii<dimNumber; ii++)
            tmpMarginalOrder[ii] = ii;

        std::sort(tmpMarginalOrder, tmpMarginalOrder + dimNumber, comparator);
        marginalResults = new PrecalculatedMarginal*[dimNumber];
        
        for(int ii=0; ii<dimNumber; ii++)
            marginalResults[ii] = marginalResultsUnsorted[tmpMarginalOrder[ii]];

        marginalOrder = new int[dimNumber];
        for(int ii = 0; ii<dimNumber; ii++)
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
        maxConfsLPSum[0] = marginalResults[0]->getModeLProb();

    for(int ii=1; ii<dimNumber-1; ii++)
        maxConfsLPSum[ii] = maxConfsLPSum[ii-1] + marginalResults[ii]->getModeLProb();

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
    for(int ii=0; ii<dimNumber; ii++)
    {
        counter[ii] = marginalResults[ii]->get_no_confs()-1;
        partialLProbs[ii] = -std::numeric_limits<double>::infinity();
    }
    partialLProbs[dimNumber] = -std::numeric_limits<double>::infinity();
    lProbs_ptr = lProbs_ptr_start + marginalResults[0]->get_no_confs()-1;
}

size_t IsoThresholdGenerator::count_confs()
{
    // Smarter algorithm forthcoming in 2.0
    size_t ret = 0;
    while(advanceToNextConfiguration())
        ret++;
    reset();
    return ret;
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

    for(int i = 0; i<dimNumber; i++)
        marginalResults[i] = new MarginalTrek(std::move(*(marginals[i])), _tabSize, _hashSize);

    logProbs        = new const vector<double>*[dimNumber];
    masses          = new const vector<double>*[dimNumber];
    marginalConfs   = new const vector<int*>*[dimNumber];

    for(int i = 0; i<dimNumber; i++)
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




#if !ISOSPEC_BUILDING_R

void printConfigurations(
    const   std::tuple<double*,double*,int*,int>& results,
    int     dimNumber,
    int*    isotopeNumbers
){
    int m = 0;

    for(int i=0; i<std::get<3>(results); i++){

        std::cout << "Mass = "  << std::get<0>(results)[i] <<
        "\tand log-prob = "     << std::get<1>(results)[i] <<
        "\tand prob = "                 << exp(std::get<1>(results)[i]) <<
        "\tand configuration =\t";


        for(int j=0; j<dimNumber; j++){
            for(int k=0; k<isotopeNumbers[j]; k++ )
            {
                std::cout << std::get<2>(results)[m] << " ";
                m++;
            }
            std::cout << '\t';
        }


        std::cout << std::endl;
    }
}

#endif /* !ISOSPEC_BUILDING_R */



IsoLayeredGenerator::IsoLayeredGenerator( Iso&&     iso,
                        double    _targetCoverage,
                        double    _percentageToExpand,
                        int       _tabSize,
                        int       _hashSize,
                        bool      trim
) : IsoGenerator(std::move(iso)),
allocator(dimNumber, _tabSize),
candidate(new int[dimNumber]),
targetCoverage(_targetCoverage >= 1.0 ? 10000.0 : _targetCoverage), // If the user wants the entire spectrum,
                                                                    // give it to him - and make sure we don't terminate
                                                                    // early because of rounding errors
percentageToExpand(_percentageToExpand),
do_trim(trim),
layers(0),
generator_position(-1)
{
    marginalResults = new MarginalTrek*[dimNumber];

    for(int i = 0; i<dimNumber; i++)
        marginalResults[i] = new MarginalTrek(std::move(*(marginals[i])), _tabSize, _hashSize);

    logProbs        = new const vector<double>*[dimNumber];
    masses          = new const vector<double>*[dimNumber];
    marginalConfs   = new const vector<int*>*[dimNumber];

    for(int i = 0; i<dimNumber; i++)
    {
        masses[i] = &marginalResults[i]->conf_masses();
        logProbs[i] = &marginalResults[i]->conf_lprobs();
        marginalConfs[i] = &marginalResults[i]->confs();
    }

    void* topConf = allocator.newConf();
    memset(reinterpret_cast<char*>(topConf) + sizeof(double), 0, sizeof(int)*dimNumber);
    *(reinterpret_cast<double*>(topConf)) = combinedSum(getConf(topConf), logProbs, dimNumber);

    current = new std::vector<void*>();
    next    = new std::vector<void*>();

    current->push_back(topConf);

    lprobThr = (*reinterpret_cast<double*>(topConf));

    if(targetCoverage > 0.0)
        while(advanceToNextLayer()) {};
}


IsoLayeredGenerator::~IsoLayeredGenerator()
{
    if(current != nullptr)
        delete current;
    if(next != nullptr)
        delete next;
    delete[] logProbs;
    delete[] masses;
    delete[] marginalConfs;
    delete[] candidate;
    dealloc_table(marginalResults, dimNumber);
}

bool IsoLayeredGenerator::advanceToNextLayer()
{
    layers += 1;
    double maxFringeLprob = -std::numeric_limits<double>::infinity();

    if(current == nullptr)
        return false;
    int accepted_in_this_layer = 0;
    Summator prob_in_this_layer(totalProb);

    void* topConf;

    while(current->size() > 0)
    {
        topConf = current->back();
        current->pop_back();

        double top_lprob = getLProb(topConf);

        if(top_lprob >= lprobThr)
        {
            newaccepted.push_back(topConf);
            accepted_in_this_layer++;
            prob_in_this_layer.add(exp(top_lprob));
        }
        else
        {
            next->push_back(topConf);
            continue;
        }

        int* topConfIsoCounts = getConf(topConf);

        for(int j = 0; j < dimNumber; ++j)
        {
            // candidate cannot refer to a position that is
            // out of range of the stored marginal distribution.
            if(marginalResults[j]->probeConfigurationIdx(topConfIsoCounts[j] + 1))
            {
                memcpy(candidate, topConfIsoCounts, confSize);
                candidate[j]++;

                void*       acceptedCandidate          = allocator.newConf();
                int*        acceptedCandidateIsoCounts = getConf(acceptedCandidate);
                memcpy(     acceptedCandidateIsoCounts, candidate, confSize);

                double newConfProb = combinedSum(
                    candidate,
                    logProbs,
                    dimNumber
                );



                *(reinterpret_cast<double*>(acceptedCandidate)) = newConfProb;

                if(newConfProb >= lprobThr)
                    current->push_back(acceptedCandidate);
                else
        {
                    next->push_back(acceptedCandidate);
            if(newConfProb > maxFringeLprob)
                maxFringeLprob = top_lprob;
        }
            }
            if(topConfIsoCounts[j] > 0)
                break;
        }
    }

    if(next == nullptr || next->size() < 1)
        return false;
    else
    {
        if(prob_in_this_layer.get() < targetCoverage)
        {
            std::vector<void*>* nnew = current;
            nnew->clear();
            current = next;
            next = nnew;
            int howmany = floor(current->size()*percentageToExpand);
            lprobThr = getLProb(quickselect(current->data(), howmany, 0, current->size()));
            totalProb = prob_in_this_layer;
        }
        else
        {
            delete next;
            next = nullptr;
            delete current;
            current = nullptr;
            int start = 0;
            int end = accepted_in_this_layer - 1;
            void* swapspace;

            if(do_trim)
            {
                void** lastLayer = &(newaccepted.data()[newaccepted.size()-accepted_in_this_layer]);

                Summator qsprob(totalProb);
                while(totalProb.get() < targetCoverage)
                {
                    if(start == end)
                        break;

                    // Partition part

                    int len = end - start;
#if ISOSPEC_BUILDING_R
            int pivot = len/2 + start;  // We're very definitely NOT switching to R to use a RNG, and if R sees us use C RNG it complains...
#else
            int pivot = rand() % len + start;
#endif
                    void* pval = lastLayer[pivot];
                    double pprob = getLProb(pval);
                    mswap(lastLayer[pivot], lastLayer[end-1]);
                    int loweridx = start;
                    for(int i=start; i<end-1; i++)
                    {
                        if(getLProb(lastLayer[i]) > pprob)
                        {
                            mswap(lastLayer[i], lastLayer[loweridx]);
                            loweridx++;
                        }
                    }
                    mswap(lastLayer[end-1], lastLayer[loweridx]);

                    // Selection part

                    Summator leftProb(qsprob);
                    for(int i=start; i<=loweridx; i++)
                    {
                        leftProb.add(exp(getLProb(lastLayer[i])));
                    }
                    if(leftProb.get() < targetCoverage)
                    {
                        start = loweridx+1;
                        qsprob = leftProb;
                    }
                    else
                        end = loweridx;
                }
                int accend = newaccepted.size()-accepted_in_this_layer+start+1;

                totalProb = qsprob;
                newaccepted.resize(accend);
                return true;
            }
            else // No trimming
            {
                totalProb = prob_in_this_layer;
                return true;
            }
        }
    }
    return true;

}

bool IsoLayeredGenerator::advanceToNextConfiguration()
{
    generator_position++;
    if(generator_position < newaccepted.size())
    {
        partialLProbs[0] = getLProb(newaccepted[generator_position]);
        partialMasses[0] = combinedSum(getConf(newaccepted[generator_position]), masses, dimNumber);
        partialProbs[0] = exp(partialLProbs[0]);
        return true;
    }
    else
        return false;
}




} // namespace IsoSpec

