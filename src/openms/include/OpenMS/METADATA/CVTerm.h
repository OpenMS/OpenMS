// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

namespace OpenMS
{

  /**
      @brief Representation of controlled vocabulary term

      This class simply stores a CV term, its value and unit if necessary.

      Representation of a CV term used by CVMappings

      @ingroup Metadata
  */
  class OPENMS_DLLAPI CVTerm
  {
public:

    struct Unit
    {

      /// Default constructor
      Unit() = default;

      Unit(const String& p_accession, const String& p_name, const String& p_cv_ref) :
        accession(p_accession),
        name(p_name),
        cv_ref(p_cv_ref)
      {
      }

      /// Copy constructor
      Unit(const Unit &) = default;

      /// Move constructor
      Unit(Unit&&) = default;

      /// Destructor
      virtual ~Unit()
      {
      }

      /// Assignment operator
      Unit& operator=(const Unit&) = default;

      /// Move assignment operator
      Unit& operator=(Unit&&)& = default;

      bool operator==(const Unit& rhs) const
      {
        return accession == rhs.accession &&
               name == rhs.name &&
               cv_ref == rhs.cv_ref;
      }

      bool operator!=(const Unit& rhs) const
      {
        return !(*this == rhs);
      }

      String accession;
      String name;
      String cv_ref;
    };

    /// Default constructor
    CVTerm() = default;

    /// Detailed constructor
    CVTerm(const String& accession, const String& name = "", const String& cv_identifier_ref = "", const String& value = "", const Unit& unit = Unit());

    /// Copy constructor
    CVTerm(const CVTerm&) = default;

    /// Move constructor
    CVTerm(CVTerm&&) = default;

    /// Destructor
    virtual ~CVTerm();

    /// Assignment operator
    CVTerm& operator=(const CVTerm&) = default;

    /// Move assignment operator
    CVTerm& operator=(CVTerm&&)& = default;

    /** @name Accessors
    */
    //@{
    /// sets the accession string of the term
    void setAccession(const String& accession);

    /// returns the accession string of the term
    const String& getAccession() const;

    /// sets the name of the term
    void setName(const String& name);

    /// returns the name of the term
    const String& getName() const;

    /// sets the cv identifier reference string, e.g. UO for unit obo
    void setCVIdentifierRef(const String& cv_identifier_ref);

    /// returns the cv identifier reference string
    const String& getCVIdentifierRef() const;

    /// set the value of the term
    void setValue(const DataValue& value);

    /// returns the value of the term
    const DataValue& getValue() const;

    /// sets the unit of the term
    void setUnit(const Unit& unit);

    /// returns the unit
    const Unit& getUnit() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVTerm& rhs) const;

    /// inequality operator
    bool operator!=(const CVTerm& rhs) const;

    /// checks whether the term has a value
    bool hasValue() const;

    /// checks whether the term has a unit
    bool hasUnit() const;
    //}

protected:

    String accession_;

    String name_;

    String cv_identifier_ref_;

    Unit unit_;

    DataValue value_;
  };

} // namespace OpenMS

