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
#include "lang.h"
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
	setupMarginals(_isotopeMasses, _isotopeProbabilities);
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
        marginals = new Marginal*[dimNumber];
        for(int i=0; i<dimNumber;i++)
        {
        allDim += isotopeNumbers[i];
        marginals[i] = new Marginal(
                _isotopeMasses[i],
                _isotopeProbabilities[i],
                isotopeNumbers[i],
                atomCounts[i]
            );
            modeLProb += marginals[i]->getModeLProb();
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



inline int str_to_int(const string& s)
{
    char* endptr[1];
    const char* c_s = s.c_str();
    int ret = (int) strtol(c_s, endptr, 10);
    if (c_s == endptr[0])
        throw invalid_argument("Invalid formula");
    return ret;
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

    string cpp_formula(formula);
    int last_modeswitch = 0;
    int mode = 0;
    int pos = 0;
    std::vector<string> elements;
    std::vector<int> numbers;
    while(formula[pos] != '\0')
    {
        if(isdigit(formula[pos]) && mode == 0)
        {
            elements.push_back(cpp_formula.substr(last_modeswitch, pos-last_modeswitch));
            last_modeswitch = pos;
            mode = 1;
        }
        else if(isalpha(formula[pos]) && mode == 1)
        {
            numbers.push_back(str_to_int(cpp_formula.substr(last_modeswitch, pos-last_modeswitch)));
            last_modeswitch = pos;
            mode = 0;
        }
        pos++;
    }

    numbers.push_back(str_to_int(cpp_formula.substr(last_modeswitch, pos)));


    if(elements.size() != numbers.size())
        throw invalid_argument("Invalid formula");

    std::vector<int> element_indexes;

    for (unsigned int i=0; i<elements.size(); i++)
    {
        int idx = -1;
        for(int j=0; j<ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; j++)
        {
            if (elements[i].compare(elem_table_symbol[j]) == 0)
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



IsoGenerator::IsoGenerator(Iso&& iso) :
    Iso(std::move(iso)),
    partialLProbs(new double[dimNumber+1+ISOSPEC_PADDING]),
    partialMasses(new double[dimNumber+1+ISOSPEC_PADDING]),
    partialExpProbs(new double[dimNumber+1+ISOSPEC_PADDING])
{
    partialLProbs[dimNumber] = 0.0;
    partialMasses[dimNumber] = 0.0;
    partialExpProbs[dimNumber] = 1.0;
}


IsoGenerator::~IsoGenerator() 
{
    if(partialLProbs != nullptr)
        delete[] partialLProbs; 
    if(partialMasses != nullptr)
        delete[] partialMasses; 
    if(partialExpProbs != nullptr)
        delete[] partialExpProbs; 
}



/*
 * ----------------------------------------------------------------------------------------------------------
 */




PrecalculatedMarginal** Iso::get_MT_marginal_set(double Lcutoff, bool absolute, int tabSize, int hashSize)
{
    PrecalculatedMarginal** ret = new PrecalculatedMarginal*[dimNumber];

    if(absolute)
        Lcutoff -= modeLProb;

    for(int ii = 0; ii<dimNumber - 1; ii++)
        ret[ii] = new PrecalculatedMarginal(std::move(*(marginals[ii])),
                                            Lcutoff + marginals[ii]->getModeLProb(),
                                            true,
                                            tabSize,
                                            hashSize);


    const unsigned int ii = dimNumber - 1;
    ret[ii] = new SyncMarginal(std::move(*(marginals[ii])),
                            Lcutoff + marginals[ii]->getModeLProb(),
                            tabSize,
                            hashSize);
    return ret;
}


IsoThresholdGeneratorMT::IsoThresholdGeneratorMT(Iso&& iso, double _threshold, PrecalculatedMarginal** PMs, bool _absolute)
: IsoGenerator(Iso(iso, false)),
Lcutoff(_threshold <= 0.0 ? std::numeric_limits<double>::lowest() : (_absolute ? log(_threshold) : log(_threshold) + modeLProb)),
last_marginal(static_cast<SyncMarginal*>(PMs[dimNumber-1]))
{
    counter = new unsigned int[dimNumber+ISOSPEC_PADDING];
    maxConfsLPSum = new double[dimNumber-1];

    marginalResults = PMs;

    bool empty = false;
    for(int ii=0; ii<dimNumber-1; ii++)
    {
        counter[ii] = 0;

        if(!marginalResults[ii]->inRange(0))
            empty = true;
    }

    marginalResults[dimNumber-1] = last_marginal;
    counter[dimNumber-1] = last_marginal->getNextConfIdx();
    if(!last_marginal->inRange(counter[dimNumber-1]))
        empty = true;


    if(dimNumber > 1)
        maxConfsLPSum[0] = marginalResults[0]->getModeLProb();

    for(int ii=1; ii<dimNumber-1; ii++)
        maxConfsLPSum[ii] = maxConfsLPSum[ii-1] + marginalResults[ii]->getModeLProb();


    if(!empty)
    {
        recalc(dimNumber-1);
        counter[0]--;
    }
    else
        terminate_search();
}

bool IsoThresholdGeneratorMT::advanceToNextConfiguration()
{
    counter[0]++;
    partialLProbs[0] = partialLProbs[1] + marginalResults[0]->get_lProb(counter[0]);
    if(partialLProbs[0] >= Lcutoff)
    {
        partialMasses[0] = partialMasses[1] + marginalResults[0]->get_mass(counter[0]);
        partialExpProbs[0] = partialExpProbs[1] * marginalResults[0]->get_eProb(counter[0]);
        return true;
    }

    // If we reached this point, a carry is needed

    int idx = 0;

    while(idx<dimNumber-2)
    {
        counter[idx] = 0;
        idx++;
        counter[idx]++;
        partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
        if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= Lcutoff)
        {
            partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
            partialExpProbs[idx] = partialExpProbs[idx+1] * marginalResults[idx]->get_eProb(counter[idx]);
            recalc(idx-1);
            return true;
        }
    }

    counter[idx] = 0;
    idx++;
    counter[idx] = last_marginal->getNextConfIdx();
    if(last_marginal->inRange(counter[idx]))
    {
        partialLProbs[idx] = partialLProbs[idx+1] + last_marginal->get_lProb(counter[idx]);
        if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= Lcutoff)
        {
            partialMasses[idx] = partialMasses[idx+1] + last_marginal->get_mass(counter[idx]);
            partialExpProbs[idx] = partialExpProbs[idx+1] * last_marginal->get_eProb(counter[idx]);
            recalc(idx-1);
            return true;
        }
    }
    terminate_search();
    return false;
}

void IsoThresholdGeneratorMT::terminate_search()
{
    for(int ii=0; ii<dimNumber; ii++)
        counter[ii] = marginalResults[ii]->get_no_confs();
}

/*
 * ----------------------------------------------------------------------------------------------------------
 */










IsoThresholdGenerator::IsoThresholdGenerator(Iso&& iso, double _threshold, bool _absolute, int tabSize, int hashSize)
: IsoGenerator(std::move(iso)),
Lcutoff(_threshold <= 0.0 ? std::numeric_limits<double>::lowest() : (_absolute ? log(_threshold) : log(_threshold) + modeLProb))
{
    counter = new int[dimNumber];
    maxConfsLPSum = new double[dimNumber-1];
    marginalResults = new PrecalculatedMarginal*[dimNumber];

    bool empty = false;

    for(int ii=0; ii<dimNumber; ii++)
    {
        counter[ii] = 0;
        marginalResults[ii] = new PrecalculatedMarginal(std::move(*(marginals[ii])),
                                                        Lcutoff - modeLProb + marginals[ii]->getModeLProb(),
                                                        true,
                                                        tabSize,
                                                        hashSize);

        if(!marginalResults[ii]->inRange(0))
            empty = true;
    }

    if(dimNumber > 1)
        maxConfsLPSum[0] = marginalResults[0]->getModeLProb();

    for(int ii=1; ii<dimNumber-1; ii++)
        maxConfsLPSum[ii] = maxConfsLPSum[ii-1] + marginalResults[ii]->getModeLProb();

    if(!empty)
    {
        recalc(dimNumber-1);
        counter[0]--;
    }
    else
        terminate_search();


}

bool IsoThresholdGenerator::advanceToNextConfiguration()
{
    counter[0]++;
    partialLProbs[0] = partialLProbs[1] + marginalResults[0]->get_lProb(counter[0]);
    if(partialLProbs[0] >= Lcutoff)
    {
        partialMasses[0] = partialMasses[1] + marginalResults[0]->get_mass(counter[0]);
        partialExpProbs[0] = partialExpProbs[1] * marginalResults[0]->get_eProb(counter[0]);
        return true;
    }

    // If we reached this point, a carry is needed

    int idx = 0;

    while(idx<dimNumber-1)
    {
        counter[idx] = 0;
        idx++;
        counter[idx]++;
        partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
        if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= Lcutoff)
        {
            partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
            partialExpProbs[idx] = partialExpProbs[idx+1] * marginalResults[idx]->get_eProb(counter[idx]);
            recalc(idx-1);
            return true;
        }
    }

    terminate_search();
    return false;
}

void IsoThresholdGenerator::terminate_search()
{
    for(int ii=0; ii<dimNumber; ii++)
        counter[ii] = marginalResults[ii]->get_no_confs();
}

/*
 * ------------------------------------------------------------------------------------------------------------------------
 */

IsoOrderedGenerator::IsoOrderedGenerator(Iso&& iso, int _tabSize, int _hashSize) :
IsoGenerator(std::move(iso)), allocator(dimNumber, _tabSize)
{
    delete[] partialLProbs;
    delete[] partialMasses;
    delete[] partialExpProbs;

    partialLProbs = &currentLProb;
    partialMasses = &currentMass;
    partialExpProbs = &currentEProb;

    marginalResults = new MarginalTrek*[dimNumber];

    for(int i = 0; i<dimNumber; i++)
        marginalResults[i] = new MarginalTrek(std::move(*(marginals[i])), _tabSize, _hashSize);

    logProbs        = new const vector<double>*[dimNumber];
    masses          = new const vector<double>*[dimNumber];
    marginalConfs   = new const vector<int*>*[dimNumber];
    candidate	    = new int[dimNumber];

    for(int i = 0; i<dimNumber; i++)
    {
        masses[i] = &marginalResults[i]->conf_masses();
        logProbs[i] = &marginalResults[i]->conf_probs();
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
    delete[] candidate;
    partialLProbs = nullptr;
    partialMasses = nullptr;
    partialExpProbs = nullptr;
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
    currentEProb = exp(currentLProb);

    ccount = -1;
    for(int j = 0; j < dimNumber; ++j)
    {
        // candidate cannot refer to a position that is
        // out of range of the stored marginal distribution.
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

IsoLayeredGenerator::IsoLayeredGenerator(Iso&& iso, double _delta, int tabSize, int hashSize)
: IsoGenerator(std::move(iso)),
current_layer_lcutoff(nextafter(modeLProb, std::numeric_limits<double>::infinity())),
delta(_delta)
{

    counter = new int[dimNumber];
    maxConfsLPSum = new double[dimNumber-1];

    marginalResults = new LayeredMarginal*[dimNumber];

    final_cutoff = -100000.0; // FIXME

    for(int ii=0; ii<dimNumber; ii++)
    {
        counter[ii] = 0;

        marginalResults[ii] = new LayeredMarginal(std::move(*(marginals[ii])),
                                                            tabSize,
                                                            hashSize);

        final_cutoff += marginalResults[ii]->getSmallestLProb();
    }

    if(dimNumber > 1)
        maxConfsLPSum[0] = marginalResults[0]->getModeLProb();

    for(int ii=1; ii<dimNumber-1; ii++)
        maxConfsLPSum[ii] = maxConfsLPSum[ii-1] + marginalResults[ii]->getModeLProb();


    
    for(int ii=dimNumber-1; ii>0; ii--)
    {
    	partialLProbs[ii] = partialLProbs[ii+1] + marginalResults[ii]->getModeLProb();
	partialMasses[ii] = partialMasses[ii+1] + marginalResults[ii]->getModeMass();
	partialExpProbs[ii] = partialExpProbs[ii+1] + marginalResults[ii]->getModeEProb();
    }

    last_counters = new int[dimNumber-1];
    nextLayer(_delta);
}


bool IsoLayeredGenerator::nextLayer(double logCutoff_delta)
{
    last_layer_lcutoff = current_layer_lcutoff;
    current_layer_lcutoff += logCutoff_delta;

    std::cout << "=============================================================================================\nNextLayer\n=============================================================================================\n\n\n";

    std::cout << "nextLayer LLC: " << last_layer_lcutoff << " CCC: " << current_layer_lcutoff << std::endl;

    if(last_layer_lcutoff<final_cutoff)
        return false;

    for(int ii=0; ii<dimNumber; ii++)
    {
        std::cout << "Extending marginal " << ii << std::endl;
        marginalResults[ii]->extend(current_layer_lcutoff - modeLProb + marginals[ii]->getModeLProb());
        std::cout << "Done extending marginal" << std::endl;
    }

	memset(counter, 0, dimNumber * sizeof(unsigned int));

    recalc(dimNumber-1);

    counter[0] = marginalResults[0]->get_no_confs();

    for(int ii=0; ii<dimNumber-1; ii++)
        last_counters[ii] = counter[0];

    return true;
}


bool IsoLayeredGenerator::advanceToNextConfiguration_internal()
{
    std::cout << "cntr: " << counter[0] << " " << counter[1] << std::endl;
    counter[0]--;
    double cprob = partialLProbs[1] + marginalResults[0]->get_lProb(counter[0]);

    if(cprob <= last_layer_lcutoff)
    {
        partialLProbs[0] = cprob;
        partialMasses[0] = partialMasses[1] + marginalResults[0]->get_mass(counter[0]);
        partialExpProbs[0] = partialExpProbs[1] * marginalResults[0]->get_eProb(counter[0]);
        std::cout << "carry not needed" << std::endl;
        return true;
    }

    // If we reached this point, a carry is needed
    std::cout << "CARRY, cntr:" << counter[0] << " " << counter[1] << std::endl;

    int idx = 0;
    while(idx < dimNumber-1)
    {
        std::cout << "carry loop" << std::endl;
	std::cout << "last_counters: ";
	printArray(last_counters, 1);
        counter[idx] = 0;
	counter[0] = last_counters[idx];
        idx++;
        counter[idx]++;
        std::cout << "carry loop, cntr after carry: " << counter[0] << " " << counter[1] << std::endl;
        partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
        std::cout << "maxConfsLPSum[idx-1]: " << maxConfsLPSum[idx-1] << std::endl;
        std::cout << "partialLProbs[idx]: " << partialLProbs[idx] << std::endl;
        if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= current_layer_lcutoff)
        {
            counter[0] = last_counters[idx];

            for(int ii=idx-1; ii >=0; ii--)
            {
                std::cout << "PLP[i] = PLP[i+1] + MR[i]->LP[c[i]]: " << partialLProbs[ii+1] << " " << marginalResults[ii]->get_lProb(counter[ii]) << std::endl;
                std::cout << "i: " << ii << std::endl;
                std::cout << "cntr: " << counter[0] << " " << counter[1] << std::endl;
                partialLProbs[ii] = partialLProbs[ii+1] + marginalResults[ii]->get_lProb(counter[ii]);
            }
            std::cout << "Prob: " << partialLProbs[0] << std::endl;
            while(partialLProbs[0] < current_layer_lcutoff)
            {
                counter[0]--;
                partialLProbs[0] = partialLProbs[1] + marginalResults[0]->get_lProb(counter[0]);
                std::cout << "Prob. too low, decreasing to: " << counter[0] << " " << counter[1] << " prob: " << partialLProbs[0] << std::endl;
            }

            last_counters[idx] = counter[0];

            for(int ii=0; ii<idx; ii++)
                last_counters[ii] = last_counters[idx];

            if(partialLProbs[0] <= last_layer_lcutoff)
                return true;

            partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
            partialExpProbs[idx] = partialExpProbs[idx+1] * marginalResults[idx]->get_eProb(counter[idx]);
            recalc(idx-1);
            std::cout << "XXX" << std::endl;
//            return true;

        }
    }

    return false; // layer switch is needed
}

IsoLayeredGenerator::~IsoLayeredGenerator()
{
    dealloc_table(marginalResults, dimNumber);
    delete[] counter;
    delete[] maxConfsLPSum;
    delete[] last_counters;
}




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

} // namespace IsoSpec
