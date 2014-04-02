// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_CHEMISTRY_EMPIRICALFORMULA_H
#define OPENMS_CHEMISTRY_EMPIRICALFORMULA_H

#include <iosfwd>
#include <vector>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  class Element;
  class ElementDB;
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

          If one wants only use a specific isotope, it can be specified using "(",")" brackets. For example,
          to specify 14C a heavy isotope of carbon it is expressed as "(14)C". The isotope distribution
          of that instance contains only one isotope, 14C itself with a frequency of 100%.

          Instances EmpiricalFormula support a (limited) set of mathematical operations. Additions and subtractions
          are supported in different flavors. However, one must be careful, because this can lead to negative
          frequencies. In most cases this might be misleading, however, the class therefore supports difference
          formulae. E.g. formula differences of reactions from post-translational modifications.
  */

  class OPENMS_DLLAPI EmpiricalFormula
  {

public:

    /** @name Typedefs
    */
    //@{
    /// Iterators
    typedef Map<const Element *, SignedSize>::ConstIterator ConstIterator;
    typedef Map<const Element *, SignedSize>::ConstIterator const_iterator;

    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    EmpiricalFormula();

    /// copy constructor
    EmpiricalFormula(const EmpiricalFormula & rhs);

    /** constructor from an OpenMS String

            @throw throws ParseError if the formula cannot be parsed
    */
    explicit EmpiricalFormula(const String & rhs);

    /// constructor with element pointer and number
    EmpiricalFormula(SignedSize number, const Element * element, SignedSize charge = 0);

    /// destructor
    virtual ~EmpiricalFormula();
    //@}

    /** @name Accessors
    */
    //@{
    /// returns the mono isotopic weight of the formula (includes proton charges)
    double getMonoWeight() const;

    /// returns the average weight of the formula (includes proton charges)
    double getAverageWeight() const;

    /** @brief returns the isotope distribution of the formula
        *	The details of the calculation of the isotope distribution
        * are described in the doc to the IsotopeDistribution class.
        *	@param max_depth: this parameter gives the max isotope which is considered, if 0 all are reported
        */
    IsotopeDistribution getIsotopeDistribution(UInt max_depth) const;

    /// returns a pointer to the element with name or symbol or 0 if no such element is fount
    const Element * getElement(const String & name) const;

    /// returns a pointer to the element with given atomic number or 0 if none if found
    const Element * getElement(UInt atomic_number) const;

    /// returns a pointer to the element db which is used with this class
    const ElementDB * getElementDB() const;

    /// returns the number of atoms with the given atomic_number
    Size getNumberOf(UInt atomic_number) const;

    /// returns the number of atoms with the given name
    Size getNumberOf(const String & name) const;

    /// returns the number of atoms
    Size getNumberOf(const Element * element) const;

    /// returns the atoms total
    Size getNumberOfAtoms() const;

    /// returns the charge
    SignedSize getCharge() const;

    /// sets the charge
    void setCharge(SignedSize charge);

    /// returns the formula as a string (charges are not included)
    String toString() const;
    //@}

    /** Assignment
    */
    //@{
    /// assignment operator
    EmpiricalFormula & operator=(const EmpiricalFormula & rhs);

    /** assignment operator which assigns an string to the formula

            @throw throws ParseError if the formula cannot be parsed
    */
    EmpiricalFormula & operator=(const String & rhs);

    /// adds the elements of the given formula
    EmpiricalFormula & operator+=(const EmpiricalFormula & rhs);

    /** adds the elements from the given formula, which is given as a OpenMS String

            @throw throws ParseError if the formula cannot be parsed
    */
    EmpiricalFormula & operator+=(const String & rhs);

    /// multiplies the elements and charge with a factor
    EmpiricalFormula operator*(const SignedSize & times) const;

    /// adds the elements of the given formula and returns a new formula
    EmpiricalFormula operator+(const EmpiricalFormula & rhs) const;

    /** adds the elements of the given formula (given as a String) and returns a new formula

            @throw throws ParseError if the formula cannot be parsed
    */
    EmpiricalFormula operator+(const String & rhs) const;

    /// subtracts the elements of a formula
    EmpiricalFormula & operator-=(const EmpiricalFormula & rhs);

    /** subtracts the elements of a formula given as string

            @throw throws ParseError if the formula cannot be parsed
    */
    EmpiricalFormula & operator-=(const String & rhs);

    /// subtracts the elements of a formula an returns a new formula
    EmpiricalFormula operator-(const EmpiricalFormula & rhs) const;

    /** subtracts the elements of a formula given as a String and returns a new formula

            @throw throws ParseError if the formula cannot be parsed
    */
    EmpiricalFormula operator-(const String & rhs) const;
    //@}

    /**@name Predicates
    */
    //@{
    /// returns true if the formula does not contain a element
    bool isEmpty() const;

    /// returns true if charge != 0
    bool isCharged() const;

    /// returns true if the formula contains the element
    bool hasElement(const Element * element) const;

    /// returns true if the formula contains the element, given with its name or symbol
    bool hasElement(const String & name) const;

    /// returns true if the formula contains the element with the given atomic number
    bool hasElement(UInt atomic_number) const;

    /// returns true if the formulas contain equal elements in equal quantities
    bool operator==(const EmpiricalFormula & rhs) const;

    /** returns true if the formulas contain equal elements in equal quantities

            @throw throws ParseError if the formula cannot be parsed
    */
    bool operator==(const String & rhs) const;

    /// returns true if the formulas differ in elements composition
    bool operator!=(const EmpiricalFormula & rhs) const;

    /** returns true if the formulas differ in elements composition

            @throw throws ParseError if the formula cannot be parsed
    */
    bool operator!=(const String & rhs) const;
    //@}

    /// writes the formula to a stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const EmpiricalFormula & formula);

    /** @name Iterators
    */
    //@{
    inline ConstIterator begin() const { return formula_.begin(); }

    inline ConstIterator end() const { return formula_.end(); }
    //@}

protected:

    /// remove elements with count 0
    void removeZeroedElements_();

    Map<const Element *, SignedSize> formula_;

    SignedSize charge_;

    void readElementsFromFile_(const String & file_name);

    SignedSize parseFormula_(Map<const Element *, SignedSize> & ef, const String & formula) const;

    const ElementDB * element_db_;
  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const EmpiricalFormula & formula);

} // namespace OpenMS
#endif
