/*!
    Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.

    This file is part of IsoSpec.

    IsoSpec is free software: you can redistribute it and/or modify
    it under the terms of the Simplified ("2-clause") BSD licence.

    IsoSpec is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the Simplified BSD Licence
    along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
*/

#pragma once

#include <tuple>
#include <unordered_map>
#include <queue>
#include <limits>
#include "platform.h"
#include "dirtyAllocator.h"
#include "summator.h"
#include "operators.h"
#include "marginalTrek++.h"


#if ISOSPEC_BUILDING_R
#include <Rcpp.h>
using namespace Rcpp;
#endif /* ISOSPEC_BUILDING_R */


namespace IsoSpec
{

// This function is NOT guaranteed to be secure against malicious input. It should be used only for debugging.
unsigned int parse_formula(const char* formula,
                           std::vector<const double*>& isotope_masses,
                           std::vector<const double*>& isotope_probabilities,
                           int** isotopeNumbers,
                           int** atomCounts,
                           unsigned int* confSize);


//! The Iso class for the calculation of the isotopic distribution.
/*!
    It contains full description of the molecule for which one would like to calculate the isotopic distribution.
*/
class ISOSPEC_EXPORT_SYMBOL Iso {
private:

    //! Set up the marginal isotopic envelopes, corresponding to subisotopologues.
    /*!
        \param _isotopeMasses A table of masses of isotopes of the elements in the chemical formula,
                              e.g. {12.0, 13.003355, 1.007825, 2.014102} for C100H202.
        \param _isotopeProbabilities A table of isotope frequencies of the elements in the chemical formula,
                                     e.g. {.989212, .010788, .999885, .000115} for C100H202.
    */
    void setupMarginals(const double* const * _isotopeMasses,
                        const double* const * _isotopeProbabilities);
public:
    bool            disowned;       /*!< A variable showing if the Iso class was specialized by its child-class. If so, then the description of the molecules has been transfered there and Iso is a carcass class, dead as a dodo, an ex-class if you will. */
protected:
    int             dimNumber;      /*!< The number of elements in the chemical formula of the molecule. */
    int*            isotopeNumbers; /*!< A table with numbers of isotopes for each element. */
    int*            atomCounts;     /*!< A table with numbers of isotopes for each element. */
    unsigned int    confSize;       /*!< The number of bytes needed to represent the counts of isotopes present in the extended chemical formula. */
    int             allDim;         /*!< The total number of isotopes of elements present in a chemical formula, e.g. for H20 it is 2+3=5. */
    Marginal**      marginals;      /*!< The table of pointers to the distributions of individual subisotopologues. */
    double          modeLProb;      /*!< The log-probability of the mode of the isotopic distribution. */

public:
    //! General constructror.
    /*!
        \param _dimNumber The number of elements in the formula, e.g. for C100H202 it would be 2, as there are only carbon and hydrogen atoms.
        \param _isotopeNumbers A table with numbers of isotopes for each element, e.g. for C100H202 it would be {2, 2}, because both C and H have two stable isotopes.
        \param _atomCounts Number of atoms of each element in the formula, e.g. for C100H202 corresponds to {100, 202}.
        \param _isotopeMasses A table of masses of isotopes of the elements in the chemical formula, e.g. {12.0, 13.003355, 1.007825, 2.014102} for C100H202.
        \param _isotopeProbabilities A table of isotope frequencies of the elements in the chemical formula, e.g. {.989212, .010788, .999885, .000115} for C100H202.
    */
    Iso(
        int             _dimNumber,
        const int*      _isotopeNumbers,
        const int*      _atomCounts,
        const double* const *  _isotopeMasses,
        const double* const *  _isotopeProbabilities
    );

    //! Constructor from the formula object.
    Iso(const char* formula);

    //! The move constructor.
    Iso(Iso&& other);

    //! The copy constructor.
    /*!
        \param other The other instance of the Iso class.
        \param fullcopy If false, copy only the number of atoms in the formula, the size of the configuration, the total number of isotopes, and the probability of the mode isotopologue.
    */
    Iso(const Iso& other, bool fullcopy);

    //! Destructor.
    virtual ~Iso();

    //! Get the mass of the lightest peak in the isotopic distribution.
    double getLightestPeakMass() const;

    //! Get the mass of the heaviest peak in the isotopic distribution.
    double getHeaviestPeakMass() const;

    //! Get the log-probability of the mode-configuration (if there are many modes, they share this value).
    inline double getModeLProb() const { return modeLProb; };

    //! Get the number of elements in the chemical formula of the molecule.
    inline int getDimNumber() const { return dimNumber; };

    //! Get the total number of isotopes of elements present in a chemical formula.
    inline int getAllDim() const { return allDim; };
};


//! The generator of isotopologues.
/*!
    This class provides the common interface for all isotopic generators.
*/
class ISOSPEC_EXPORT_SYMBOL IsoGenerator : public Iso
{
protected:
    double* partialLProbs;  /*!< The prefix sum of the log-probabilities of the current isotopologue. */
    double* partialMasses;  /*!< The prefix sum of the masses of the current isotopologue. */
    double* partialProbs;/*!< The prefix product of the probabilities of the current isotopologue. */

public:
    //! Advance to the next, not yet visited, most probable isotopologue.
    /*!
        \return Return false if it is not possible to advance.
    */
    virtual bool advanceToNextConfiguration() = 0;

    //! Get the log-probability of the current isotopologue.
    /*!
        \return The log-probability of the current isotopologue.
    */
    virtual double lprob() const { return partialLProbs[0]; };

    //! Get the mass of the current isotopologue.
    /*!
        \return The mass of the current isotopologue.
    */
    virtual double mass()  const { return partialMasses[0]; };

    //! Get the probability of the current isotopologue.
    /*!
        \return The probability of the current isotopologue.
    */
    virtual double prob() const { return partialProbs[0]; };

    //TODO: what is this???
    virtual void get_conf_signature(int* space) const = 0;

    //! Move constructor.
    IsoGenerator(Iso&& iso, bool alloc_partials = true);

    //! Destructor.
    virtual ~IsoGenerator();
};



//! The generator of isotopologues sorted by their probability of occurrence.
/*!
    The subsequent isotopologues are generated with diminishing probability, starting from the mode.
    This algorithm take O(N*log(N)) to compute the N isotopologues because of using the Priority Queue data structure.
    Obtaining the N isotopologues can be achieved in O(N) if they are not required to be spit out in the descending order.
*/
class ISOSPEC_EXPORT_SYMBOL IsoOrderedGenerator: public IsoGenerator
{
private:
    MarginalTrek**              marginalResults;            /*!< Table of pointers to marginal distributions of subisotopologues. */
    std::priority_queue<void*,std::vector<void*>,ConfOrder> pq; /*!< The priority queue used to generate isotopologues ordered by descending probability. */
    void*                       topConf;                    /*!< Most probable configuration. */
    DirtyAllocator              allocator;                  /*!< Structure used for alocating memory for isotopologues. */
    const std::vector<double>** logProbs;                   /*!< Obtained log-probabilities. */
    const std::vector<double>** masses;                     /*!< Obtained masses. */
    const std::vector<int*>**   marginalConfs;              /*!< Obtained counts of isotopes. */
    double                      currentLProb;               /*!< The log-probability of the current isotopologue. */
    double                      currentMass;                /*!< The mass of the current isotopologue. */
    double                      currentProb;                /*!< The probability of the current isotopologue. */
    int                         ccount;

public:
    bool advanceToNextConfiguration() override final;

    //! Save the counts of isotopes in the space.
    /*!
        \param space An array where counts of isotopes shall be written. 
                     Must be as big as the overall number of isotopes.
    */
    inline void get_conf_signature(int* space) const override final
    {
        int* c = getConf(topConf);

        if (ccount >= 0)
            c[ccount]--;

        for(int ii=0; ii<dimNumber; ii++)
        {
            memcpy(space, marginalResults[ii]->confs()[c[ii]], isotopeNumbers[ii]*sizeof(int));
            space += isotopeNumbers[ii];
        }

        if (ccount >= 0)
            c[ccount]++;
    };

    //! The move-contstructor.
    IsoOrderedGenerator(Iso&& iso, int _tabSize  = 1000, int _hashSize = 1000);

    //! Destructor.
    virtual ~IsoOrderedGenerator();
};




//! The generator of isotopologues above a given threshold value.
/*!
    Attention: the calculated configurations are only partially ordeded and the user should not assume they will be ordered.
    This algorithm computes N isotopologues in O(N).
    It is a considerable advantage w.r.t. the IsoOrderedGenerator.
*/
class ISOSPEC_EXPORT_SYMBOL IsoThresholdGenerator: public IsoGenerator
{
private:

    int*                    counter;            /*!< An array storing the position of an isotopologue in terms of the subisotopologues ordered by decreasing probability. */
    double*                 maxConfsLPSum;      
    const double            Lcutoff;            /*!< The logarithm of the lower bound on the calculated probabilities. */
    PrecalculatedMarginal** marginalResults;
    PrecalculatedMarginal** marginalResultsUnsorted;
    int* marginalOrder;

    const double* lProbs_ptr;
    const double* lProbs_ptr_start;
    double* partialLProbs_second;
    double partialLProbs_second_val, lcfmsv;
    bool empty;

public:
    inline void get_conf_signature(int* space) const override final
    {
        counter[0] = lProbs_ptr - lProbs_ptr_start;
        if(marginalOrder != nullptr)
            for(int ii=0; ii<dimNumber; ii++)
            {
                int jj = marginalOrder[ii];
                memcpy(space, marginalResultsUnsorted[ii]->get_conf(counter[jj]), isotopeNumbers[ii]*sizeof(int));
                space += isotopeNumbers[ii];
            }
        else
            for(int ii=0; ii<dimNumber; ii++)
            {
                memcpy(space, marginalResultsUnsorted[ii]->get_conf(counter[ii]), isotopeNumbers[ii]*sizeof(int));
                space += isotopeNumbers[ii];
            }

    };

    //! The move-constructor.
    /*!
        \param iso An instance of the Iso class.
        \param _threshold The threshold value.
        \param _absolute If true, the _threshold is interpreted as the absolute minimal peak height for the isotopologues.
                         If false, the _threshold is the fraction of the heighest peak's probability.
        \param tabSize The size of the extension of the table with configurations.
        \param hashSize The size of the hash-table used to store subisotopologues and check if they have been already calculated.
    */
    IsoThresholdGenerator(Iso&& iso, double _threshold, bool _absolute=true, int _tabSize=1000, int _hashSize=1000, bool reorder_marginals = true);

    inline ~IsoThresholdGenerator()
    {
        delete[] counter;
        delete[] maxConfsLPSum;
        if (marginalResultsUnsorted != marginalResults)
            delete[] marginalResultsUnsorted;
        dealloc_table(marginalResults, dimNumber); 
        if(marginalOrder != nullptr)
            delete[] marginalOrder;
    };

    // Perform highly aggressive inling as this function is often called as while(advanceToNextConfiguration()) {}
    // which leads to an extremely tight loop and some compilers miss this (potentially due to the length of the function). 
    ISOSPEC_FORCE_INLINE bool advanceToNextConfiguration() override final
    {
        lProbs_ptr++;

        if(ISOSPEC_LIKELY(*lProbs_ptr >= lcfmsv))
        {
            return true;
        }

        // If we reached this point, a carry is needed

        int idx = 0;
        lProbs_ptr = lProbs_ptr_start;

        int * cntr_ptr = counter;

        while(idx<dimNumber-1)
        {
            // counter[idx] = 0;
            *cntr_ptr = 0;
            idx++;
            cntr_ptr++;
            // counter[idx]++;
            (*cntr_ptr)++;
            partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
            if(partialLProbs[idx] + maxConfsLPSum[idx-1] >= Lcutoff)
            {
                partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
                partialProbs[idx] = partialProbs[idx+1] * marginalResults[idx]->get_prob(counter[idx]);
                recalc(idx-1);
                return true;
            }
        }

        terminate_search();
        return false;
    }


    ISOSPEC_FORCE_INLINE double lprob() const override final { return partialLProbs_second_val + (*(lProbs_ptr)); };
    ISOSPEC_FORCE_INLINE double mass()  const override final { return partialMasses[1] + marginalResults[0]->get_mass(lProbs_ptr - lProbs_ptr_start); };
    ISOSPEC_FORCE_INLINE double prob()  const override final { return partialProbs[1] * marginalResults[0]->get_prob(lProbs_ptr - lProbs_ptr_start); };

    //! Block the subsequent search of isotopologues.
    void terminate_search();

    /*! Reset the generator to the beginning of the sequence. Allows it to be reused, eg. to go through the conf space once, calculate
        the amount of space needed to store configurations, then to allocate that memory, and go through it again, this time saving 
        configurations (and *is* in fact faster than allocating a std::vector and depending on it to grow as needed. This is cheaper
        than throwing away the generator and making a new one too: marginal distributions don't need to be recalculated. */
    void reset();

    /*! Count the number of configurations in the distribution. This can be used to pre-allocate enough memory to store it (e.g.
     * std::vector's reserve() method - this is faster than depending on the vector's dynamic resizing, even though it means that
     * the configuration space is walked through twice. This method has to be called before the first call to advanceToNextConfiguration
     * and has undefined results (incl. segfaults) otherwise. */
    size_t count_confs();

private:
    //! Recalculate the current partial log-probabilities, masses, and probabilities.
    ISOSPEC_FORCE_INLINE void recalc(int idx)
    {
        for(; idx > 0; idx--)
        {
            partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
            partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
            partialProbs[idx] = partialProbs[idx+1] * marginalResults[idx]->get_prob(counter[idx]);
        }
        partialLProbs_second_val = *partialLProbs_second;
        partialLProbs[0] = *partialLProbs_second + marginalResults[0]->get_lProb(counter[0]);
        lcfmsv = Lcutoff - partialLProbs_second_val;
    }
};



//! The class that represents isotopologues above a given joint probability value.
/*!
    This class generates subsequent isotopologues that ARE NOT GUARANTEED TO BE ORDERED BY probability.
    The overal set of isotopologues is guaranteed to surpass a given threshold of probability contained in the
    isotopic distribution.
    This calculations are performed in O(N) operations, where N is the total number of the output isotopologues.

    This class is not a true generator yet - the generator methods have been implemented for compatibility, but
    the class actually performs all computations during the initialization and stores them, and the generator methods
    only walk through the array of precomputed values. . It will be reimplemented as a true generator in 2.0.
*/
class ISOSPEC_EXPORT_SYMBOL IsoLayeredGenerator : public IsoGenerator
{
private:
    Summator                totalProb;
    std::vector<void*>      newaccepted;
    DirtyAllocator allocator;
    int* candidate;
    const std::vector<double>** logProbs;                   /*!< Obtained log-probabilities. */
    const std::vector<double>** masses;                     /*!< Obtained masses. */
    const std::vector<int*>**   marginalConfs;              /*!< Obtained counts of isotopes. */
    MarginalTrek** marginalResults;
    std::vector<void*>*         current;
    std::vector<void*>*         next;
    double                      lprobThr;
    double                      targetCoverage;
    double                      percentageToExpand;
    bool                        do_trim;
    int layers;
    size_t generator_position;

    bool advanceToNextLayer(); 

public:
    bool advanceToNextConfiguration() override final;

    inline void get_conf_signature(int* space) const override final
    {
        int* conf = getConf(newaccepted[generator_position]);
        for(int ii=0; ii<dimNumber; ii++)
        {
            memcpy(space, marginalResults[ii]->confs()[conf[ii]], isotopeNumbers[ii]*sizeof(int));
            space += isotopeNumbers[ii];
        }
    };


    IsoLayeredGenerator(Iso&& iso, double _targetCoverage, double _percentageToExpand = 0.3, int _tabSize  = 1000, int _hashSize = 1000, bool trim = false);
    virtual ~IsoLayeredGenerator();

    void terminate_search();

};




#if !ISOSPEC_BUILDING_R

void printConfigurations(
    const   std::tuple<double*,double*,int*,int>& results,
    int     dimNumber,
    int*    isotopeNumbers
);
#endif /* !ISOSPEC_BUILDING_R */

} // namespace IsoSpec

