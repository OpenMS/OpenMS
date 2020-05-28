/*!
    Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.

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

#include <unordered_map>
#include <queue>
#include <limits>
#include <string>
#include <vector>
#include "platform.h"
#include "dirtyAllocator.h"
#include "summator.h"
#include "operators.h"
#include "marginalTrek++.h"



namespace IsoSpec
{

// This function is NOT guaranteed to be secure against malicious input. It should be used only for debugging.
unsigned int parse_formula(const char* formula,
                           std::vector<double>& isotope_masses,
                           std::vector<double>& isotope_probabilities,
                           int** isotopeNumbers,
                           int** atomCounts,
                           unsigned int* confSize,
                           bool use_nominal_masses = false);


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
    void setupMarginals(const double* _isotopeMasses,
                        const double* _isotopeProbabilities);
    bool            disowned;       /*!< A variable showing if the Iso class was specialized by its child-class. If so, then the description of the molecules has been transfered there and Iso is a carcass class, dead as a dodo, an ex-class if you will. */

 protected:
    int             dimNumber;      /*!< The number of elements in the chemical formula of the molecule. */
    int*            isotopeNumbers; /*!< A table with numbers of isotopes for each element. */
    int*            atomCounts;     /*!< A table with numbers of isotopes for each element. */
    unsigned int    confSize;       /*!< The number of bytes needed to represent the counts of isotopes present in the extended chemical formula. */
    int             allDim;         /*!< The total number of isotopes of elements present in a chemical formula, e.g. for H20 it is 2+3=5. */
    Marginal**      marginals;      /*!< The table of pointers to the distributions of individual subisotopologues. */

    bool doMarginalsNeedSorting() const;

 public:
    Iso();

    //! General constructror.
    /*!
        \param _dimNumber The number of elements in the formula, e.g. for C100H202 it would be 2, as there are only carbon and hydrogen atoms.
        \param _isotopeNumbers A table with numbers of isotopes for each element, e.g. for C100H202 it would be {2, 2}, because both C and H have two stable isotopes.
        \param _atomCounts Number of atoms of each element in the formula, e.g. for C100H202 corresponds to {100, 202}.
        \param _isotopeMasses A table of tables of masses of isotopes of the elements in the chemical formula, e.g. {{12.0, 13.003355}, {1.007825, 2.014102}} for C100H202.
        \param _isotopeProbabilities A table of tables of isotope frequencies of the elements in the chemical formula, e.g. {{.989212, .010788}, {.999885, .000115}} for C100H202.
    */
    Iso(
        int             _dimNumber,
        const int*      _isotopeNumbers,
        const int*      _atomCounts,
        const double*   _isotopeMasses,
        const double*   _isotopeProbabilities
    );
    Iso(
        int             _dimNumber,
        const int*      _isotopeNumbers,
        const int*      _atomCounts,
        const double* const *  _isotopeMasses,
        const double* const *  _isotopeProbabilities
    );

    //! Constructor from the formula object.
    Iso(const char* formula, bool use_nominal_masses = false);  // NOLINT(runtime/explicit) - constructor deliberately left to be used as a conversion

    //! Constructor from C++ std::string chemical formula.
    inline Iso(const std::string& formula, bool use_nominal_masses = false) : Iso(formula.c_str(), use_nominal_masses) {}  // NOLINT(runtime/explicit) - constructor deliberately left to be used as a conversion

    //! Constructor (named) from aminoacid FASTA sequence as C string.
    /*!
        \param fasta An aminoacid FASTA sequence. May be upper/lower/mixed case, may contain selenocystein (U). Subisotopologues will be in order: CHNOS, possibly with Se added at an end if present.
        \use_nominal_masses Whether to use nucleon number instead of the real mass of each isotope during calculations.
        \add_water Whether the chain should have the terminating -H and -OH groups at the N and C terminus, respectively.
    */
    static Iso FromFASTA(const char* fasta, bool use_nominal_masses = false, bool add_water = true);

    //! Constructor (named) from aminoacid FASTA sequence as C++ std::string. See above for details.
    static inline Iso FromFASTA(const std::string& fasta, bool use_nominal_masses = false, bool add_water = true) { return FromFASTA(fasta.c_str(), use_nominal_masses, add_water); }

    //! The move constructor.
    Iso(Iso&& other);

    /* We're not exactly following standard copy and assign semantics with Iso objects, so delete the default assign constructor just in case, so noone tries to use it. Copy ctor declared below. */
    Iso& operator=(const Iso& other) = delete;

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

    /*!
        Get the mass of the monoisotopic peak in the isotopic distribution. Monoisotopc molecule is defined as
        consisting only of the most frequent isotopes of each element. These are often (but not always) the
        lightest ones, making this often (but again, not always) equal to getLightestPeakMass()
    */
    double getMonoisotopicPeakMass() const;

    //! Get the log-probability of the mode-configuration (if there are many modes, they share this value).
    double getModeLProb() const;

    //! Get the logprobability of the least probable subisotopologue.
    double getUnlikeliestPeakLProb() const;

    //! Get the mass of the mode-configuration (if there are many modes, it is undefined which one will be selected).
    double getModeMass() const;

    //! Get the theoretical average mass of the molecule.
    double getTheoreticalAverageMass() const;

    //! Get the theoretical variance of the distribution.
    double variance() const;

    //! Get the standard deviation of the theoretical distribution.
    double stddev() const { return sqrt(variance()); }

    //! Get the number of elements in the chemical formula of the molecule.
    inline int getDimNumber() const { return dimNumber; }

    //! Get the total number of isotopes of elements present in a chemical formula.
    inline int getAllDim() const { return allDim; }

    //! Add an element to the molecule. Note: this method can only be used BEFORE Iso is used to construct an IsoGenerator instance.
    void addElement(int atomCount, int noIsotopes, const double* isotopeMasses, const double* isotopeProbabilities);

    //! Save estimates of logarithms of target sizes of marginals using Gaussian approximation into argument array. Array priorities must have length equal to dimNumber.
    void saveMarginalLogSizeEstimates(double* priorities, double target_total_prob) const;
};


//! The generator of isotopologues.
/*!
    This class provides the common interface for all isotopic generators.
*/
class ISOSPEC_EXPORT_SYMBOL IsoGenerator : public Iso
{
 public:
    const double mode_lprob;

 protected:
    double* partialLProbs;  /*!< The prefix sum of the log-probabilities of the current isotopologue. */
    double* partialMasses;  /*!< The prefix sum of the masses of the current isotopologue. */
    double* partialProbs;   /*!< The prefix product of the probabilities of the current isotopologue. */

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
    virtual double lprob() const { return partialLProbs[0]; }

    //! Get the mass of the current isotopologue.
    /*!
        \return The mass of the current isotopologue.
    */
    virtual double mass()  const { return partialMasses[0]; }

    //! Get the probability of the current isotopologue.
    /*!
        \return The probability of the current isotopologue.
    */
    virtual double prob() const { return partialProbs[0]; }

    //! Write the signature of configuration into target memory location. It must be large enough to accomodate it.
    virtual void get_conf_signature(int* space) const = 0;

    //! Move constructor.
    IsoGenerator(Iso&& iso, bool alloc_partials = true);  // NOLINT(runtime/explicit) - constructor deliberately left to be used as a conversion

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
    MarginalTrek**              marginalResults;                    /*!< Table of pointers to marginal distributions of subisotopologues. */
    std::priority_queue<void*, std::vector<void*>, ConfOrder> pq;   /*!< The priority queue used to generate isotopologues ordered by descending probability. */
    void*                       topConf;                            /*!< Most probable configuration. */
    DirtyAllocator              allocator;                          /*!< Structure used for alocating memory for isotopologues. */
    const std::vector<double>** logProbs;                           /*!< Obtained log-probabilities. */
    const std::vector<double>** masses;                             /*!< Obtained masses. */
    const std::vector<int*>**   marginalConfs;                      /*!< Obtained counts of isotopes. */
    double                      currentLProb;                       /*!< The log-probability of the current isotopologue. */
    double                      currentMass;                        /*!< The mass of the current isotopologue. */
    double                      currentProb;                        /*!< The probability of the current isotopologue. */
    int                         ccount;

 public:
    IsoOrderedGenerator(const IsoOrderedGenerator& other) = delete;
    IsoOrderedGenerator& operator=(const IsoOrderedGenerator& other) = delete;

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

        for(int ii = 0; ii < dimNumber; ii++)
        {
            memcpy(space, marginalResults[ii]->confs()[c[ii]], isotopeNumbers[ii]*sizeof(int));
            space += isotopeNumbers[ii];
        }

        if (ccount >= 0)
            c[ccount]++;
    };

    //! The move-contstructor.
    IsoOrderedGenerator(Iso&& iso, int _tabSize  = 1000, int _hashSize = 1000);  // NOLINT(runtime/explicit) - constructor deliberately left to be used as a conversion

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
    IsoThresholdGenerator(const IsoThresholdGenerator& other) = delete;
    IsoThresholdGenerator& operator=(const IsoThresholdGenerator& other) = delete;

    inline void get_conf_signature(int* space) const override final
    {
        counter[0] = lProbs_ptr - lProbs_ptr_start;
        if(marginalOrder != nullptr)
        {
            for(int ii = 0; ii < dimNumber; ii++)
            {
                int jj = marginalOrder[ii];
                memcpy(space, marginalResultsUnsorted[ii]->get_conf(counter[jj]), isotopeNumbers[ii]*sizeof(int));
                space += isotopeNumbers[ii];
            }
        }
        else
        {
            for(int ii = 0; ii < dimNumber; ii++)
            {
                memcpy(space, marginalResultsUnsorted[ii]->get_conf(counter[ii]), isotopeNumbers[ii]*sizeof(int));
                space += isotopeNumbers[ii];
            }
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
    IsoThresholdGenerator(Iso&& iso, double _threshold, bool _absolute = true, int _tabSize = 1000, int _hashSize = 1000, bool reorder_marginals = true);

    ~IsoThresholdGenerator();

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

        while(idx < dimNumber-1)
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


    ISOSPEC_FORCE_INLINE double lprob() const override final { return partialLProbs_second_val + (*(lProbs_ptr)); }
    ISOSPEC_FORCE_INLINE double mass()  const override final { return partialMasses[1] + marginalResults[0]->get_mass(lProbs_ptr - lProbs_ptr_start); }
    ISOSPEC_FORCE_INLINE double prob()  const override final { return partialProbs[1] * marginalResults[0]->get_prob(lProbs_ptr - lProbs_ptr_start); }

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

    ISOSPEC_FORCE_INLINE void short_recalc(int idx)
    {
        for(; idx > 0; idx--)
            partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
        partialLProbs_second_val = *partialLProbs_second;
        partialLProbs[0] = *partialLProbs_second + marginalResults[0]->get_lProb(counter[0]);
        lcfmsv = Lcutoff - partialLProbs_second_val;
    }
};





class ISOSPEC_EXPORT_SYMBOL IsoLayeredGenerator : public IsoGenerator
{
 private:
    int*                    counter;            /*!< An array storing the position of an isotopologue in terms of the subisotopologues ordered by decreasing probability. */
    double*                 maxConfsLPSum;
    double currentLThreshold, lastLThreshold;
    LayeredMarginal** marginalResults;
    LayeredMarginal** marginalResultsUnsorted;
    int* marginalOrder;

    const double* lProbs_ptr;
    const double* lProbs_ptr_start;
    const double** resetPositions;
    double* partialLProbs_second;
    double partialLProbs_second_val, lcfmsv, last_lcfmsv;
    bool marginalsNeedSorting;


 public:
    IsoLayeredGenerator(const IsoLayeredGenerator& other) = delete;
    IsoLayeredGenerator& operator=(const IsoLayeredGenerator& other) = delete;

    inline void get_conf_signature(int* space) const override final
    {
        counter[0] = lProbs_ptr - lProbs_ptr_start;
        if(marginalOrder != nullptr)
        {
            for(int ii = 0; ii < dimNumber; ii++)
            {
                int jj = marginalOrder[ii];
                memcpy(space, marginalResultsUnsorted[ii]->get_conf(counter[jj]), isotopeNumbers[ii]*sizeof(int));
                space += isotopeNumbers[ii];
            }
        }
        else
        {
            for(int ii = 0; ii < dimNumber; ii++)
            {
                memcpy(space, marginalResultsUnsorted[ii]->get_conf(counter[ii]), isotopeNumbers[ii]*sizeof(int));
                space += isotopeNumbers[ii];
            }
        }
    };

    inline double get_currentLThreshold() const { return currentLThreshold; }

    IsoLayeredGenerator(Iso&& iso, int _tabSize = 1000, int _hashSize = 1000, bool reorder_marginals = true, double t_prob_hint = 0.99);  // NOLINT(runtime/explicit) - constructor deliberately left to be used as a conversion

    ~IsoLayeredGenerator();

    ISOSPEC_FORCE_INLINE bool advanceToNextConfiguration() override final
    {
        do
        {
            if(advanceToNextConfigurationWithinLayer())
                return true;
        } while(IsoLayeredGenerator::nextLayer(-2.0));
        return false;
    }

    ISOSPEC_FORCE_INLINE bool advanceToNextConfigurationWithinLayer()
    {
        do{
            lProbs_ptr++;

            if(ISOSPEC_LIKELY(*lProbs_ptr >= lcfmsv))
                return true;
        }
        while(carry());  // NOLINT(whitespace/empty_loop_body) - cpplint bug, that's not an empty loop body, that's a do{...}while(...) construct
        return false;
    }

    ISOSPEC_FORCE_INLINE double lprob() const override final { return partialLProbs_second_val + (*(lProbs_ptr)); };
    ISOSPEC_FORCE_INLINE double mass()  const override final { return partialMasses[1] + marginalResults[0]->get_mass(lProbs_ptr - lProbs_ptr_start); };
    ISOSPEC_FORCE_INLINE double prob()  const override final { return partialProbs[1] * marginalResults[0]->get_prob(lProbs_ptr - lProbs_ptr_start); };

    //! Block the subsequent search of isotopologues.
    void terminate_search();


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
        partialLProbs[0] = partialLProbs_second_val + marginalResults[0]->get_lProb(counter[0]);
        lcfmsv = currentLThreshold - partialLProbs_second_val;
        last_lcfmsv = lastLThreshold - partialLProbs_second_val;
    }

    bool nextLayer(double offset);

 private:
    bool carry();
};



class IsoStochasticGenerator : IsoGenerator
{
    IsoLayeredGenerator ILG;
    size_t to_sample_left;
    const double precision;
    const double beta_bias;
    double confs_prob;
    double chasing_prob;
    size_t current_count;

 public:
    IsoStochasticGenerator(Iso&& iso, size_t no_molecules, double precision = 0.9999, double beta_bias = 5.0);

    ISOSPEC_FORCE_INLINE size_t count() const { return current_count; }

    ISOSPEC_FORCE_INLINE double mass() const override final { return ILG.mass(); }

    ISOSPEC_FORCE_INLINE double prob() const override final { return static_cast<double>(count()); }

    ISOSPEC_FORCE_INLINE double lprob() const override final { return log(prob()); }

    ISOSPEC_FORCE_INLINE void get_conf_signature(int* space) const override final { ILG.get_conf_signature(space); }

    ISOSPEC_FORCE_INLINE bool advanceToNextConfiguration() override final
    {
        /* This function will be used mainly in very small, tight loops, therefore it makes sense to
         * aggressively inline it, despite its seemingly large body.
         */
        while(true)
        {
            double curr_conf_prob_left, current_prob;

            if(to_sample_left <= 0)
                return false;

            if(confs_prob < chasing_prob)
            {
                // Beta was last
                current_count = 1;
                to_sample_left--;
                ILG.advanceToNextConfiguration();
                current_prob = ILG.prob();
                confs_prob += current_prob;
                while(confs_prob <= chasing_prob)
                {
                    ILG.advanceToNextConfiguration();
                    current_prob = ILG.prob();
                    confs_prob += current_prob;
                }
                if(to_sample_left <= 0)
                    return true;
                curr_conf_prob_left = confs_prob - chasing_prob;
            }
            else
            {
                // Binomial was last
                current_count = 0;
                ILG.advanceToNextConfiguration();
                current_prob = ILG.prob();
                confs_prob += current_prob;
                curr_conf_prob_left = current_prob;
            }

            double prob_left_to_1 = precision - chasing_prob;
            double expected_confs = curr_conf_prob_left * to_sample_left / prob_left_to_1;

            if(expected_confs <= beta_bias)
            {
                // Beta mode: we keep making beta jumps until we leave the current configuration
                chasing_prob += rdvariate_beta_1_b(to_sample_left) * prob_left_to_1;
                while(chasing_prob <= confs_prob)
                {
                    current_count++;
                    to_sample_left--;
                    if(to_sample_left == 0)
                        return true;
                    prob_left_to_1 = precision - chasing_prob;
                    chasing_prob += rdvariate_beta_1_b(to_sample_left) * prob_left_to_1;
                }
                if(current_count > 0)
                    return true;
            }
            else
            {
                // Binomial mode: a single binomial step
                size_t rbin = rdvariate_binom(to_sample_left, curr_conf_prob_left/prob_left_to_1);
                current_count += rbin;
                to_sample_left -= rbin;
                chasing_prob = confs_prob;
                if(current_count > 0)
                    return true;
            }
        };
    }
};


}  // namespace IsoSpec
