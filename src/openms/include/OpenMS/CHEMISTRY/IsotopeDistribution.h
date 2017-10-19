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

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_H

#include <OpenMS/CONCEPT/Types.h>

#include <utility>
#include <vector>
#include <set>

namespace OpenMS
{
  /**
        @ingroup Chemistry

        @brief Isotope distribution class

        Holds an isotope distribution with the weight value and according
        probability. Distribution can be add using the '+' or '+=' operators.

        The most important value which should be set is the max isotope value.
        This value can be set using the setMaxIsotope method. It is an upper
        bound for the number of isotopes which are calculated. E.g. if it is set
        to 3, only the first three isotopes, Monoisotopic mass, +1 and +2 are
        calculated.
        By default all possible isotopes are calculated, which leads to a large
        number of values, if the mass value is large!
    */
  class OPENMS_DLLAPI IsotopeDistribution
  {
public:

    /// @name typedefs
    //@{
    /// container type, first holds the weight of the isotope, second the probability
    typedef std::vector<std::pair<Size, double> > ContainerType;
    typedef ContainerType::iterator iterator;
    typedef ContainerType::iterator Iterator;
    typedef ContainerType::const_iterator const_iterator;
    typedef ContainerType::const_iterator ConstIterator;

    typedef ContainerType::reverse_iterator reverse_iterator;
    typedef ContainerType::reverse_iterator ReverseIterator;
    typedef ContainerType::const_reverse_iterator const_reverse_iterator;
    typedef ContainerType::const_reverse_iterator ConstReverseIterator;
    //@}

    /// @name Constructors and Destructors
    //@{
    /** Default constructor, note max_isotope must be set later
            @see setMaxIsotope(Size max_isotope)
    */
    IsotopeDistribution();

    /// Detailed constructor which sets the @p max_isotope
    explicit IsotopeDistribution(Size max_isotope);

    /// Copy constructor
    IsotopeDistribution(const IsotopeDistribution & isotope_distribution);

    /// Destructor
    virtual ~IsotopeDistribution();
    //@}

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

    /// overwrites the container which holds the distribution using @p distribution
    void set(const ContainerType & distribution);

    /// returns the container which holds the distribution
    const ContainerType & getContainer() const;

    /// returns the maximal weight isotope which is stored in the distribution
    Size getMax() const;

    /// returns the minimal weight isotope which is stored in the distribution
    Size getMin() const;

    /// returns the size of the distribution which is the number of isotopes in the distribution
    Size size() const;

    /// clears the distribution and resets max isotope to 0
    void clear();

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

    /** @brief re-normalizes the sum of the probabilities of the isotopes to 1

            The re-normalisation is needed as in distributions with a lot of isotopes (and with high max isotope)
            the calculations tend to be inexact.
    */
    void renormalize();

    /** @brief Trims the right side of the isotope distribution to isotopes with a significant contribution.

            If the isotope distribution is calculated for large masses (and with high max isotope)
            it might happen that many entries contain only small numbers. This function can be
            used to remove these entries.

            Do consider normalising the distribution afterwards.
    */
    void trimRight(double cutoff);

    /** @brief Trims the left side of the isotope distribution to isotopes with a significant contribution.

            If the isotope distribution is calculated for large masses (and with high max isotope)
            it might happen that many entries contain only small numbers. This function can be
            used to remove these entries.

            Do consider normalising the distribution afterwards.
    */
    void trimLeft(double cutoff);
    //@}

    /// @name Operators
    //@{
    /// Assignment operator
    IsotopeDistribution & operator=(const IsotopeDistribution & isotope_distribution);

    /// operator which adds this distribution and the @p isotope_distribution to return IsotopeDisribution (similar to convolve distributions)
    IsotopeDistribution operator+(const IsotopeDistribution & isotope_distribution) const;

    /// operator which adds @p isotope_distribution to this (similar to convolve distributions)
    IsotopeDistribution & operator+=(const IsotopeDistribution & isotope_distribution);

    /// operator which multiplies this distribution by @p factor (similar to @p factor times applying operator '+')
    IsotopeDistribution operator*(Size factor) const;

    /// operator which multiplies this distribution by @p factor (similar to @p factor times applying operator '+=')
    IsotopeDistribution & operator*=(Size factor);

    /// equality operator, returns true if the @p isotope_distribution is identical to this, false else
    bool operator==(const IsotopeDistribution & isotope_distribution) const;

    /// inequality operator, returns true if the @p isotope_distribution differs from this, false else
    bool operator!=(const IsotopeDistribution & isotope_distribution) const;
    //@}

    /// @name Iterators
    //@{
    inline Iterator begin() { return distribution_.begin(); }

    inline Iterator end()   { return distribution_.end(); }

    inline ConstIterator begin() const { return distribution_.begin(); }

    inline ConstIterator end() const { return distribution_.end(); }

    inline ReverseIterator rbegin() { return distribution_.rbegin(); }

    inline ReverseIterator rend()   { return distribution_.rend(); }

    inline ConstReverseIterator rbegin() const { return distribution_.rbegin(); }

    inline ConstReverseIterator rend() const { return distribution_.rend(); }
    //@}

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

    /// maximal isotopes which is used to calculate the distribution
    Size max_isotope_;

    /// stores the isotope distribution
    ContainerType distribution_;
  };

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_H
