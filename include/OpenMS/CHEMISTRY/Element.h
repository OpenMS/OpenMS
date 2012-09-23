// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_CHEMISTRY_ELEMENT_H
#define OPENMS_CHEMISTRY_ELEMENT_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#define OPENMS_CHEMISTRY_ELEMENT_NAME_DEFAULT "unknown"
#define OPENMS_CHEMISTRY_ELEMENT_SYMBOL_DEFAULT "??"
#define OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT 0.0
#define OPENMS_CHEMISTRY_ELEMENT_ATOMICNUMBER_DEFAULT 0

namespace OpenMS
{
  /** @ingroup Chemistry

          @brief Representation of an element
  */
  class OPENMS_DLLAPI Element
  {
public:

    /** @name Constructor and Destructors
    */
    //@{
    /// default constructor
    Element();

    /// copy constructor
    Element(const Element & element);

    /// detailed constructor
    Element(const String & name,
            const String & symbol,
            UInt atomic_number,
            DoubleReal average_weight,
            DoubleReal mono_weight,
            const IsotopeDistribution & isotopes);

    /// destructor
    virtual ~Element();
    //@}

    /** @name Accessors
    */
    //@{
    /// sets unique atomic number
    void setAtomicNumber(UInt atomic_number);

    /// returns the unique atomic number
    UInt getAtomicNumber() const;

    /// sets the average weight of the element
    void setAverageWeight(DoubleReal weight);

    /// returns the average weight of the element
    DoubleReal getAverageWeight() const;

    /// sets the mono isotopic weight of the element
    void setMonoWeight(DoubleReal weight);

    /// returns the mono isotopic weight of the element
    DoubleReal getMonoWeight() const;

    /// sets the isotope distribution of the element
    void setIsotopeDistribution(const IsotopeDistribution & isotopes);

    /// returns the isotope distribution of the element
    const IsotopeDistribution & getIsotopeDistribution() const;

    /// set the name of the element
    void setName(const String & name);

    /// returns the name of the element
    const String & getName() const;

    /// sets symbol of the element
    void setSymbol(const String & symbol);

    /// returns symbol of the element
    const String & getSymbol() const;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    Element & operator=(const Element & element);
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const Element & element) const;

    /// inequality operator
    bool operator!=(const Element & element) const;
    //@}

    /// writes the element to an output stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Element & element);

protected:

    /// name of the element
    String name_;

    /// symbol of the element
    String symbol_;

    /// atomic number of the element
    UInt atomic_number_;

    /// average weight over all isotopes
    DoubleReal average_weight_;

    /// mono isotopic weight of the most frequent isotope
    DoubleReal mono_weight_;

    /// distribution of the isotopes
    IsotopeDistribution isotopes_;
  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream &, const Element &);

} // namespace OpenMS

#endif
