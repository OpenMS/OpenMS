// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Ahmed Khalil $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//
#pragma once

#include <iosfwd>
#include <map>
#include <set>
#include <string>

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  class String;
  class Element;
  class ElementDB;
  class IsotopeDistribution;
  class IsotopePatternGenerator;
  class CoarseIsotopePatternGenerator;

  /**
    @ingroup Chemistry

    @brief Representation of an empirical formula

    The formula can be used as follows: elements are represented through its symbol or full name.
    The symbol or name is followed by a number. If not, the frequency is set to one. Examples
    are CH3OH or CarbonHydrogen3OH. The names must start with an capital letter (symbols always have
    an upper-case letter at the beginning). Additionally charges can be used with '+' followed by
    a number, if no number follows the charge of +1 is set. Negative charges can be added using a '-'
    sign. However, negative charges are only set if the last element in the string also has a number.
    E.g. H4C-1, does not set a negative charge, only -1 Carbon atoms, correctly it should be stated
    H4C-1-.

    This class also supports the usage of specific isotopes. By default "C" describes not one isotope
    but a natural distribution or different isotopes. This distribution can be accessed via the
    member getIsotopeDistribution().

    If one wants only use a specific isotope, it can be specified using "(",")"
    brackets. For example, to specify 13C a heavy isotope of carbon it is
    expressed as "(13)C". The isotope distribution of that instance contains
    only one isotope, 13C itself with a frequency of 100%. It is possible to
    mix isotopes, for example "(13)C1CH6O" specifies an ethanol molecule with
    one 12C and one 13C isotope.

    Instances EmpiricalFormula support a (limited) set of mathematical operations. Additions and subtractions
    are supported in different flavors. However, one must be careful, because this can lead to negative
    frequencies. In most cases this might be misleading, however, the class therefore supports difference
    formulae. E.g. formula differences of reactions from post-translational modifications.
  */

  class OPENMS_DLLAPI EmpiricalFormula
  {

protected:
	  /// Internal typedef for the used map type
	  typedef std::map<const Element*, SignedSize> MapType_;

public:
    /** @name Typedefs
    */
    //@{
    /// Iterators
    typedef MapType_::const_iterator ConstIterator;
    typedef MapType_::const_iterator const_iterator;
    typedef MapType_::iterator Iterator;
    typedef MapType_::iterator iterator;
    //@}

    /** @name Constructors and Destructors
    */
    //@{
    /// Default constructor
    EmpiricalFormula();

    /// Copy constructor
    EmpiricalFormula(const EmpiricalFormula&) = default;

    /// Move constructor
    EmpiricalFormula(EmpiricalFormula&&) = default;

    /**
      Constructor from an OpenMS String

      @throw throws ParseError if the formula cannot be parsed
    */
    explicit EmpiricalFormula(const String& rhs);

    /// Constructor with element pointer and number
    EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge = 0);

    /// Destructor
    virtual ~EmpiricalFormula();
    //@}

     /**
     @brief Create EmpiricalFormula object from a String

     @param rhs Input string

     @throws Exception::ParseError if the formula cannot be parsed
   */
    static EmpiricalFormula fromString(const String& rhs)
    {
      EmpiricalFormula ef(rhs);
      return ef;
    }

    /** @name Accessors
    */
    //@{
    /// returns the monoisotopic (most abundant isotope per element) weight of the formula (includes proton charges)
    double getMonoWeight() const;

    /// returns the sum of the lightest isotopes per element in the formula (includes proton charges)
    double getLightestIsotopeWeight() const;

    /// returns the average weight of the formula (includes proton charges)
    double getAverageWeight() const;

    /// returns the total number of discrete isotopes
    double calculateTheoreticalIsotopesNumber() const;

    /**
      @brief Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry

      @param average_weight: Average weight to estimate an EmpiricalFormula for
      @param C: The approximate relative stoichiometry of Carbons to other elements in this molecule
      @param H: The approximate relative stoichiometry of Hydrogens to other elements in this molecule
      @param N: The approximate relative stoichiometry of Nitrogens to other elements in this molecule
      @param O: The approximate relative stoichiometry of Oxygens to other elements in this molecule
      @param S: The approximate relative stoichiometry of Sulfurs to other elements in this molecule
      @param P: The approximate relative stoichiometry of Phosphoruses to other elements in this molecule

      @return bool flag for whether the approximation succeeded without requesting negative hydrogens. true = no problems, 1 = negative hydrogens requested.
    */
    bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P);

    /**
      @brief Fills this EmpiricalFormula with an approximate elemental composition for a given monoisotopic weight and approximate elemental stoichiometry

      @param mono_weight: Monoisotopic weight to estimate an EmpiricalFormula for
      @param C: The approximate relative stoichiometry of Carbons to other elements in this molecule
      @param H: The approximate relative stoichiometry of Hydrogens to other elements in this molecule
      @param N: The approximate relative stoichiometry of Nitrogens to other elements in this molecule
      @param O: The approximate relative stoichiometry of Oxygens to other elements in this molecule
      @param S: The approximate relative stoichiometry of Sulfurs to other elements in this molecule
      @param P: The approximate relative stoichiometry of Phosphoruses to other elements in this molecule

      @return bool flag for whether the approximation succeeded without requesting negative hydrogens. true = no problems, 1 = negative hydrogens requested.
    */
    bool estimateFromMonoWeightAndComp(double mono_weight, double C, double H, double N, double O, double S, double P);

    /**
      @brief Fills this EmpiricalFormula with an approximate elemental composition for a given average weight,
      exact number of sulfurs, and approximate elemental stoichiometry

      @param average_weight: Average weight to estimate an EmpiricalFormula for
      @param S: The exact number of Sulfurs in this molecule
      @param C: The approximate relative stoichiometry of Carbons to other elements (excluding Sulfur) in this molecule
      @param H: The approximate relative stoichiometry of Hydrogens to other elements (excluding Sulfur) in this molecule
      @param N: The approximate relative stoichiometry of Nitrogens to other elements (excluding Sulfur) in this molecule
      @param O: The approximate relative stoichiometry of Oxygens to other elements (excluding Sulfur) in this molecule
      @param P: The approximate relative stoichiometry of Phosphoruses to other elements (excluding Sulfur) in this molecule

      @return bool flag for whether the approximation succeeded without requesting negative hydrogens. true = no problems, false = negative hydrogens requested.
   */
    bool estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P);


    /**
      @brief returns the isotope distribution of the formula
      The details of the calculation of the isotope distribution
      are described in the doc to the CoarseIsotopePatternGenerator class.

      @param method: the method that will be used for the calculation of the IsotopeDistribution
    */
    IsotopeDistribution getIsotopeDistribution(const IsotopePatternGenerator& method) const;

    /**
      @brief returns the fragment isotope distribution of this given a precursor formula
      and conditioned on a set of isolated precursor isotopes.

      The max_depth of the isotopic distribution is set to max(precursor_isotopes)+1.
      @param precursor: the empirical formula of the precursor
      @param precursor_isotopes: the precursor isotopes that were isolated
      @param method: the method that will be used for the calculation of the IsotopeDistribution
      @return the conditional IsotopeDistribution of the fragment
    */
    IsotopeDistribution getConditionalFragmentIsotopeDist(const EmpiricalFormula& precursor,
                                                          const std::set<UInt>& precursor_isotopes,
                                                          const CoarseIsotopePatternGenerator& method) const;

    /// returns the number of atoms for a certain @p element (can be negative)
    SignedSize getNumberOf(const Element* element) const;

    /// returns the atoms total (not absolute: negative counts for certain elements will reduce the overall count). Negative result is possible.
    SignedSize getNumberOfAtoms() const;

    /// returns the charge
    Int getCharge() const;

    /// sets the charge
    void setCharge(Int charge);

    /// returns the formula as a string (charges are not included)
    String toString() const;

    /// returns the formula as a map (charges are not included)
    std::map<std::string, int> toMap() const;
    //@}

    /** Assignment
    */
    //@{

    /// Assignment operator
    EmpiricalFormula& operator=(const EmpiricalFormula&) = default;

    /// Move assignment operator
    EmpiricalFormula& operator=(EmpiricalFormula&&) & = default;

    /// adds the elements of the given formula
    EmpiricalFormula& operator+=(const EmpiricalFormula& rhs);

    /// multiplies the elements and charge with a factor
    EmpiricalFormula operator*(const SignedSize& times) const;

    /// adds the elements of the given formula and returns a new formula
    EmpiricalFormula operator+(const EmpiricalFormula& rhs) const;

    /// subtracts the elements of a formula
    EmpiricalFormula& operator-=(const EmpiricalFormula& rhs);

    /// subtracts the elements of a formula an returns a new formula
    EmpiricalFormula operator-(const EmpiricalFormula& rhs) const;

    //@}

    /**@name Predicates
    */
    //@{
    /// returns true if the formula does not contain a element
    bool isEmpty() const;

    /// returns true if charge != 0
    bool isCharged() const;

    /// returns true if the formula contains the element
    bool hasElement(const Element* element) const;

    /// returns true if all elements from @p ef are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula
    bool contains(const EmpiricalFormula& ef) const;

    /// returns true if the formulas contain equal elements in equal quantities
    bool operator==(const EmpiricalFormula& rhs) const;

    /// returns true if the formulas differ in elements composition
    bool operator!=(const EmpiricalFormula& rhs) const;

    /// less operator
    bool operator<(const EmpiricalFormula& rhs) const;

    //@}

    /// writes the formula to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const EmpiricalFormula& formula);

    /** @name Iterators
    */
    //@{
    inline ConstIterator begin() const { return formula_.begin(); }

    inline ConstIterator end() const { return formula_.end(); }

    inline Iterator begin() { return formula_.begin(); }

    inline Iterator end() { return formula_.end(); }
    //@}

    /** @name Static member functions
     */
    // @TODO: make these static member variables instead?
    //@{
    /// Efficiently generates a formula for hydrogen
    static EmpiricalFormula hydrogen(int n_atoms = 1);

    /// Efficiently generates a formula for water
    static EmpiricalFormula water(int n_molecules = 1);
    //@}

protected:

    /// remove elements with count 0
    void removeZeroedElements_();

    MapType_ formula_;

    Int charge_;

    Int parseFormula_(std::map<const Element*, SignedSize>& ef, const String& formula) const;

  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const EmpiricalFormula& formula);

} // namespace OpenMS
