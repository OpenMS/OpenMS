// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <set>

namespace OpenMS
{
  /**
    * @ingroup Chemistry
    * @brief Isotope pattern generator for coarse isotope distributions.
    * 
    * This algorithm generates theoretical pattern distributions for empirical
    * formulas with resolution of 1Da. It assumes that every isotope has atomic
    * mass that is rounded to the closest integer in Daltons, therefore it
    * produces coarse distributions (it does not discriminate between 13C, N15
    * and O18 peaks).  For example for (13)Carbon it assumes that the mass of
    * the isotope is 13Da instead of 13.0033548378. In case you need fine
    * resolution, please consider using FineIsotopePatternGenerator.
    *
    * The output is a list of pairs containing nominal isotope probabilities
    * paired with a number that is either an accurate or rounded (integer)
    * mass. The accurate masses assume the nominal isotopes are mostly due to
    * (13)Carbon.  To return accurate vs rounded masses, use setRoundMasses
    * accordingly.  The default is to return accurate masses (note that setting
    * this option will not influence the probabilities and still produce a
    * coarse distributions spaced at ca 1Da). For example, using rounded mass,
    * for a C100 molecule, you will get:
    *
    * @code
    *     1200 : 0.341036528
    *     1201 : 0.368855864
    *     1202 : 0.197477505
    *     1203 : 0.0697715357
    * @endcode
    *
    * while accurate mass will produce:
    *
    * @code
    *     1200       : 0.341036528
    *     1201.00335 : 0.368855864
    *     1202.00671 : 0.197477505
    *     1203.01006 : 0.0697715357
    * @endcode
    *
    * The other important value which needs to be set is the max isotope value.
    * This value can be set using the setMaxIsotope method. It is an upper
    * bound for the number of isotopes which are calculated If e.g., set to 3,
    * only the first three isotopes, Monoisotopic mass, +1 and +2 are
    * calculated.
    *
    * @note By default all possible isotopes are calculated, which leads to a large
    * number of values, if the mass value is large!
    *
    * @note If you need fine isotope distributions, consider using the
    * FineIsotopePatternGenerator.
    *
    * See also method run()
    **/

  class OPENMS_DLLAPI CoarseIsotopePatternGenerator 
    : public IsotopePatternGenerator
  {

 public:
    CoarseIsotopePatternGenerator();

    CoarseIsotopePatternGenerator(const Size& max_isotope);

    CoarseIsotopePatternGenerator(const Size& max_isotope, const bool round_masses);

    virtual ~CoarseIsotopePatternGenerator();

    /// @name Accessors
    //@{
    /** @brief sets the maximal isotope with @p max_isotope

            sets the maximal isotope which is included in the distribution
            and used to limit the calculations. This is useful as distributions
            with numerous isotopes tend to have a lot of numerical zeros at the end
    */
    void setMaxIsotope(const Size& max_isotope);

    /// sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)
    void setRoundMasses(const bool round_masses);

    /// returns the currently set maximum isotope
    Size getMaxIsotope() const;

    /// returns the current value of the flag to return expected masses (true) or atomic numbers (false).
    bool getRoundMasses() const;
    //@}

    /**
      * @brief Creates an isotope distribution from an empirical sum formula
      *
      * Iterates through all elements, convolves them according to the number
      * of atoms from that element and sums up the result.
      *
      **/
    IsotopeDistribution run(const EmpiricalFormula&) const override;

    /**
       @brief Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported

       Implementation using the averagine model proposed by Senko et al. in
       "Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions"
    */
    IsotopeDistribution estimateFromPeptideWeight(double average_weight);

    /**
       @brief Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs

       @param average_weight: Average weight to estimate an EmpiricalFormula for
       @param S: The exact number of Sulfurs in this molecule

       @pre S <= average_weight / average_weight(sulfur)
       @pre average_weight >= 0
    */
    IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, UInt S);

    /**
       @brief Estimate Nucleotide Isotopedistribution from weight and number of isotopes that should be reported

       averagine model from Zubarev, R. A.; Demirev, P. A. in
       "Isotope  depletion  of  large biomolecules: Implications for molecular mass measurements."
    */
    IsotopeDistribution estimateFromRNAWeight(double average_weight);

    /**
       @brief Estimate Nucleotide Isotopedistribution from weight and number of isotopes that should be reported
       averagine model from Zubarev, R. A.; Demirev, P. A. in
       "Isotope  depletion  of  large biomolecules: Implications for molecular mass measurements."
    */
    IsotopeDistribution estimateFromDNAWeight(double average_weight);

    /**

       @brief Estimate Isotopedistribution from weight, average composition, and number of isotopes that should be reported

    */
    IsotopeDistribution estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P);

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
    IsotopeDistribution estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P);

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
    IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes);

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
    IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, const std::set<UInt>& precursor_isotopes);

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
    IsotopeDistribution estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes);

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
    IsotopeDistribution estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes);

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
    IsotopeDistribution estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes, double C, double H, double N, double O, double S, double P);

    /**
       @brief Calculate isotopic distribution for a fragment molecule

       This calculates the isotopic distribution for a fragment molecule given
       the isotopic distribution of the fragment and complementary fragment (as
       if they were precursors), and which precursor isotopes were isolated.

       @note Do consider normalising the distribution afterwards to get conditional probabilities.

       Equations come from Rockwood, AL; Kushnir, MA; Nelson, GJ. in
       "Dissociation of Individual Isotopic Peaks: Predicting Isotopic Distributions of Product Ions in MSn"

       @param fragment_isotope_dist the isotopic distribution of the fragment (as if it was a precursor).
       @param comp_fragment_isotope_dist the isotopic distribution of the complementary fragment (as if it was a precursor).
       @param precursor_isotopes a list of which precursor isotopes were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.
       @param fragment_mono_mass the monoisotopic mass of the fragment.
       @pre fragment_isotope_dist and comp_fragment_isotope_dist are gapless (no missing isotopes between the min/max isotopes of the dist)
    */
    IsotopeDistribution calcFragmentIsotopeDist(const IsotopeDistribution& fragment_isotope_dist, const IsotopeDistribution& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes, const double fragment_mono_mass) const;

    CoarseIsotopePatternGenerator& operator=(const CoarseIsotopePatternGenerator& iso);

    /// convolves the distributions @p left and @p right and stores the result in @p result
    IsotopeDistribution::ContainerType convolve_(const IsotopeDistribution::ContainerType & left, const IsotopeDistribution::ContainerType & right) const;

    /// convolves the distribution @p input @p factor times and stores the result in @p result
    IsotopeDistribution::ContainerType convolvePow_(const IsotopeDistribution::ContainerType & input, Size factor) const;

    /// convolves the distribution @p input with itself and stores the result in @p result
    IsotopeDistribution::ContainerType convolveSquare_(const IsotopeDistribution::ContainerType & input) const;

    /// converts the masses of distribution @p input from atomic numbers to accurate masses
    IsotopeDistribution::ContainerType correctMass_(const IsotopeDistribution::ContainerType & input, const double mono_weight) const;

  protected:

    /** @brief calculates the fragment distribution for a fragment molecule and stores it in @p result.

        @param fragment_isotope_dist the isotopic distribution of the fragment (as if it was a precursor).
        @param comp_fragment_isotope_dist the isotopic distribution of the complementary fragment (as if it was a precursor).
        @param precursor_isotopes which precursor isotopes were isolated. 0 corresponds to the mono-isotopic molecule (M0), 1->M1, etc.
    */
    IsotopeDistribution calcFragmentIsotopeDist_(const IsotopeDistribution::ContainerType& fragment_isotope_dist, const IsotopeDistribution::ContainerType& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes) const;

    /// fill a gapped isotope pattern (i.e. certain masses are missing), with zero probability masses
    IsotopeDistribution::ContainerType fillGaps_(const IsotopeDistribution::ContainerType& id) const;

 protected:
    /// maximal isotopes which is used to calculate the distribution
    Size max_isotope_;
    /// flag to determine whether masses should be rounded or not
    bool round_masses_;

  };

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_COARSEID_H
