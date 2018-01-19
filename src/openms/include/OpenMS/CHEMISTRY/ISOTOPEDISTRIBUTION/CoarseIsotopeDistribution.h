// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_COARSEID_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_COARSEID_H

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>


namespace OpenMS
{

  /**
    *         @ingroup Chemistry
    *         @brief Isotope pattern generator for coarse (nominal mass) isotope distributions.
    *
    *         The most important value which should be set is the max isotope value.
    *         It limits the number of calculated isotopes.
    *         If e.g., set to 3, only the first three isotopes, Monoisotopic mass, +1 and +2 are
    *         calculated.
    *         By default all possible isotopes are calculated, which leads to a large
    *         number of values, if the mass value is large!
    **/

  class OPENMS_DLLAPI CoarseIsotopeDistribution 
    : public IsotopePatternGenerator
  {

 public:
    CoarseIsotopeDistribution();

    CoarseIsotopeDistribution(Size max_isotope);

    CoarseIsotopeDistribution(const IsotopeDistribution& isotope_distribution);

    /// @name Accessors
    //@{
    /** @brief sets the maximal isotope with @p max_isotope

            sets the maximal isotope which is included in the distribution
            and used to limit the calculations. This is useful as distributions
            with numerous isotopes tend to have a lot of numerical zeros at the end
    */
    void setMaxIsotope(Size max_isotope);

    /// returns the currently set maximum isotope
    Size getMaxIsotope() const;
    
    Size getMin() const;
    Size getMax() const;
    
    void clear();
    //@}
    /// equality operator, returns true if the @p isotope_distribution is identical to this, false else
    // bool operator==(const CoarseIsotopeDistribution& isotope_distribution) const;

    /// inequality operator, returns true if the @p isotope_distribution differs from this, false else
    bool operator!=(const CoarseIsotopeDistribution& isotope_distribution) const;


    void run(const EmpiricalFormula&);

    /**
       @brief Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported

       Implementation using the averagine model proposed by Senko et al. in
       "Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions"
    */
    void estimateFromPeptideWeight(double average_weight);

    /**
       @brief Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs

       @param average_weight: Average weight to estimate an EmpiricalFormula for
       @param S: The exact number of Sulfurs in this molecule

       @pre S <= average_weight / average_weight(sulfur)
       @pre average_weight >= 0
    */
    void estimateFromPeptideWeightAndS(double average_weight, UInt S);

    /**
       @brief Estimate Nucleotide Isotopedistribution from weight and number of isotopes that should be reported

       averagine model from Zubarev, R. A.; Demirev, P. A. in
       "Isotope  depletion  of  large biomolecules: Implications for molecular mass measurements."
    */
    void estimateFromRNAWeight(double average_weight);

    /**
       @brief Estimate Nucleotide Isotopedistribution from weight and number of isotopes that should be reported
       averagine model from Zubarev, R. A.; Demirev, P. A. in
       "Isotope  depletion  of  large biomolecules: Implications for molecular mass measurements."
    */
    void estimateFromDNAWeight(double average_weight);

    /**

       @brief Estimate Isotopedistribution from weight, average composition, and number of isotopes that should be reported

    */
    void estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P);

    /**
       @brief Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition

       @param average_weight: Average weight to estimate an IsotopeDistribution for
       @param S: The exact numbers of Sulfurs in this molecule
       @param C: The approximate relative stoichiometry of Carbons to other elements (excluding Sulfur) in this molecule
       @param H: The approximate relative stoichiometry of Hydrogens to other elements (excluding Sulfur) in this molecule
       @param N: The approximate relative stoichiometry of Nitrogens to other elements (excluding Sulfur) in this molecule
       @param O: The approximate relative stoichiometry of Oxygens to other elements (excluding Sulfur) in this molecule
       @param P: The approximate relative stoichiometry of Phosphoruses to other elements (excluding Sulfur) in this molecule

       @pre S, C, H, N, O, P >= 0
       @pre average_weight >= 0
    */
    void estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P);

    /**
       @brief Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
       fragment's average weight, and a list of isolated precursor isotopes.

       The max_depth of the isotopic distribution is set to max(precursor_isotopes)+1.
       @param average_weight_precursor: average weight of the precursor peptide
       @param average_weight_fragment: average weight of the fragment
       @param precursor_isotopes: the precursor isotopes that were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.

       @pre average_weight_precursor >= average_weight_fragment
       @pre average_weight_fragment > 0
       @pre average_weight_precursor > 0
       @pre precursor_isotopes.size() > 0
    */
    void estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes);

    /**
       @brief Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
       number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
       and a list of isolated precursor isotopes.

       The max_depth of the isotopic distribution is set to max(precursor_isotopes)+1.
       @param average_weight_precursor: average weight of the precursor peptide
       @param S_precursor: The exact number of Sulfurs in the precursor peptide
       @param average_weight_fragment: average weight of the fragment
       @param S_fragment: The exact number of Sulfurs in the fragment
       @param precursor_isotopes: the precursor isotopes that were isolated

       @pre S_fragment <= average_weight_fragment / average_weight(sulfur)
       @pre S_precursor - S_fragment <= (average_weight_precursor - average_weight_fragment) / average_weight(sulfur)
       @pre average_weight_precursor >= average_weight_fragment
       @pre average_weight_precursor > 0
       @pre average_weight_fragment > 0
       @pre precursor_isotopes.size() > 0
    */
    void estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, const std::set<UInt>& precursor_isotopes);

    /**
       @brief Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
       fragment's average weight, and a list of isolated precursor isotopes.

       The max_depth of the isotopic distribution is set to max(precursor_isotopes)+1.
       @param average_weight_precursor: average weight of the precursor nucleotide
       @param average_weight_fragment: average weight of the fragment
       @param precursor_isotopes: the precursor isotopes that were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.

       @pre average_weight_precursor >= average_weight_fragment
       @pre average_weight_precursor > 0
       @pre average_weight_fragment > 0
       @pre precursor_isotopes.size() > 0
    */
    void estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes);

    /**
       @brief Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
       fragment's average weight, and a list of isolated precursor isotopes.

       The max_depth of the isotopic distribution is set to max(precursor_isotopes)+1.
       @param average_weight_precursor: average weight of the precursor nucleotide
       @param average_weight_fragment: average weight of the fragment
       @param precursor_isotopes: the precursor isotopes that were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.

       @pre average_weight_precursor >= average_weight_fragment
       @pre average_weight_precursor > 0
       @pre average_weight_fragment > 0
       @pre precursor_isotopes.size() > 0
    */
    void estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes);

    /**
       @brief Estimate fragment IsotopeDistribution from the precursor's average weight,
       fragment's average weight, a list of isolated precursor isotopes, and average composition

       The max_depth of the isotopic distribution is set to max(precursor_isotopes)+1.
       @param average_weight_precursor: average weight of the precursor molecule
       @param average_weight_fragment: average weight of the fragment molecule
       @param precursor_isotopes: the precursor isotopes that were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.
       @param C: The approximate relative stoichiometry of Carbons to other elements in this molecule
       @param H: The approximate relative stoichiometry of Hydrogens to other elements in this molecule
       @param N: The approximate relative stoichiometry of Nitrogens to other elements in this molecule
       @param O: The approximate relative stoichiometry of Oxygens to other elements in this molecule
       @param S: The approximate relative stoichiometry of Sulfurs to other elements in this molecule
       @param P: The approximate relative stoichiometry of Phosphoruses to other elements in this molecule

       @pre S, C, H, N, O, P >= 0
       @pre average_weight_precursor >= average_weight_fragment
       @pre average_weight_precursor > 0
       @pre average_weight_fragment > 0
       @pre precursor_isotopes.size() > 0
    */
    void estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes, double C, double H, double N, double O, double S, double P);

    /**
       @brief Calculate isotopic distribution for a fragment molecule

       This calculates the isotopic distribution for a fragment molecule given
       the isotopic distribution of the fragment and complementary fragment
       (as if they were precursors), and which precursor isotopes were isolated.
       Do consider normalising the distribution afterwards to get conditional probabilities.
       Equations come from Rockwood, AL; Kushnir, MA; Nelson, GJ. in
       "Dissociation of Individual Isotopic Peaks: Predicting Isotopic Distributions of Product Ions in MSn"
       @param fragment_isotope_dist the isotopic distribution of the fragment (as if it was a precursor).
       @param comp_fragment_isotope_dist the isotopic distribution of the complementary fragment (as if it was a precursor).
       @param precursor_isotopes a list of which precursor isotopes were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.
    */
    void calcFragmentIsotopeDist(const IsotopeDistribution& fragment_isotope_dist, const IsotopeDistribution& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes);

    CoarseIsotopeDistribution& operator=(const CoarseIsotopeDistribution& iso);

    /// equality operator, returns true if the @p isotope_distribution is identical to this, false else
    bool operator==(const CoarseIsotopeDistribution& isotope_distribution) const;

    /// operator which adds this distribution and the @p isotope_distribution to return IsotopeDisribution (similar to convolve distributions)
    CoarseIsotopeDistribution operator+(const CoarseIsotopeDistribution& isotope_distribution) const;

    /// operator which adds @p isotope_distribution to this (similar to convolve distributions)
    CoarseIsotopeDistribution& operator+=(const CoarseIsotopeDistribution& isotope_distribution);

    /// operator which multiplies this distribution by @p factor (similar to @p factor times applying operator '+')
    CoarseIsotopeDistribution operator*(Size factor) const;

    /// operator which multiplies this distribution by @p factor (similar to @p factor times applying operator '+=')
    CoarseIsotopeDistribution& operator*=(Size factor);
    
 protected:
    /// convolves the distributions @p left and @p right and stores the result in @p result
    void convolve_(ContainerType & result, const ContainerType & left, const ContainerType & right) const;

    /// convolves the distribution @p input @p factor times and stores the result in @p result
    void convolvePow_(ContainerType & result, const ContainerType & input, Size factor) const;

    /// convolves the distribution @p input with itself and stores the result in @p result
    void convolveSquare_(ContainerType & result, const ContainerType & input) const;


    /** @brief calculates the fragment distribution for a fragment molecule and stores it in @p result.

        @param fragment_isotope_dist the isotopic distribution of the fragment (as if it was a precursor).
        @param comp_fragment_isotope_dist the isotopic distribution of the complementary fragment (as if it was a precursor).
        @param precursor_isotopes which precursor isotopes were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.
    */
    void calcFragmentIsotopeDist_(ContainerType& result, const ContainerType& fragment_isotope_dist, const ContainerType& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes);

    /// fill a gapped isotope pattern (i.e. certain masses are missing), with zero probability masses
    ContainerType fillGaps_(const ContainerType& id) const;

 protected:
    /// maximal isotopes which is used to calculate the distribution
    Size max_isotope_;

  };

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_COARSEID_H
