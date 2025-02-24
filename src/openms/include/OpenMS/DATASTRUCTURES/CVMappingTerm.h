// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  /**
      @brief Representation of controlled vocabulary term

      This class simply stores CV terms read from e.g. an OBO-file

      @ingroup Datastructures
  */
  ///Representation of a CV term used by CVMappings
  class OPENMS_DLLAPI CVMappingTerm
  {
public:

    /// Defaults constructor
    CVMappingTerm();

    /// Copy constructor
    CVMappingTerm(const CVMappingTerm& rhs);

    /// Destructor
    virtual ~CVMappingTerm();

    /// Assignment operator
    CVMappingTerm& operator=(const CVMappingTerm& rhs);

    /** @name Accessors
    */
    //@{
    /// sets the accession string of the term
    void setAccession(const String& accession);

    /// returns the accession string of the term
    const String& getAccession() const;

    /// sets whether the term name should be used, instead of the accession
    void setUseTermName(bool use_term_name);

    /// returns whether the term name should be used, instead of the accession
    bool getUseTermName() const;

    /// sets whether the term itself can be used (or only its children)
    void setUseTerm(bool use_term);

    /// returns true if the term can be used, false if only children are allowed
    bool getUseTerm() const;

    /// sets the name of the term
    void setTermName(const String& term_name);

    /// returns the name of the term
    const String& getTermName() const;

    /// sets whether this term can be repeated
    void setIsRepeatable(bool is_repeatable);

    /// returns true if this term can be repeated, false otherwise
    bool getIsRepeatable() const;

    /// sets whether children of this term are allowed
    void setAllowChildren(bool allow_children);

    /// returns true if the children of this term are allowed to be used
    bool getAllowChildren() const;

    /// sets the cv identifier reference string, e.g. UO for unit obo
    void setCVIdentifierRef(const String& cv_identifier_ref);

    /// returns the cv identifier reference string
    const String& getCVIdentifierRef() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVMappingTerm& rhs) const;

    /// inequality operator
    bool operator!=(const CVMappingTerm& rhs) const;
    //}

protected:

    String accession_;

    bool use_term_name_;

    bool use_term_;

    String term_name_;

    bool is_repeatable_;

    bool allow_children_;

    String cv_identifier_ref_;
  };

} // namespace OpenMS

