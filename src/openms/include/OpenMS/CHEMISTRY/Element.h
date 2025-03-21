// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <string>

#define OPENMS_CHEMISTRY_ELEMENT_NAME_DEFAULT "unknown"
#define OPENMS_CHEMISTRY_ELEMENT_SYMBOL_DEFAULT "??"
#define OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT 0.0
#define OPENMS_CHEMISTRY_ELEMENT_ATOMICNUMBER_DEFAULT 0

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Representation of an element

      This contains information on an element and its isotopes, including a
      common name, atomic symbol and mass/abundance of its isotopes.
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
    Element(const std::string & name,
            const std::string & symbol,
            unsigned int atomic_number,
            double average_weight,
            double mono_weight,
            const IsotopeDistribution & isotopes);

    /// destructor
    virtual ~Element();
    //@}

    /** @name Accessors
    */
    //@{
    /// sets unique atomic number
    void setAtomicNumber(unsigned int atomic_number);

    /// returns the unique atomic number
    unsigned int getAtomicNumber() const;

    /// sets the average weight of the element
    void setAverageWeight(double weight);

    /// returns the average weight of the element
    double getAverageWeight() const;

    /// sets the mono isotopic weight of the element
    void setMonoWeight(double weight);

    /// returns the mono isotopic weight of the element
    double getMonoWeight() const;

    /// sets the isotope distribution of the element
    void setIsotopeDistribution(const IsotopeDistribution & isotopes);

    /// returns the isotope distribution of the element
    const IsotopeDistribution & getIsotopeDistribution() const;

    /// set the name of the element
    void setName(const std::string & name);

    /// returns the name of the element
    const std::string & getName() const;

    /// sets symbol of the element
    void setSymbol(const std::string & symbol);

    /// returns symbol of the element
    const std::string & getSymbol() const;
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

    /// less operator
    bool operator<(const Element & element) const;
    //@}

    /// writes the element to an output stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Element & element);

protected:

    /// name of the element
    std::string name_;

    /// symbol of the element
    std::string symbol_;

    /// atomic number of the element
    unsigned int atomic_number_;

    /// average weight over all isotopes
    double average_weight_;

    /// mono isotopic weight of the most frequent isotope
    double mono_weight_;

    /// distribution of the isotopes (mass and natural frequency)
    IsotopeDistribution isotopes_;
  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream &, const Element &);

} // namespace OpenMS

