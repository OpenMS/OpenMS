// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <map>

namespace OpenMS
{
  /**
      @brief Representation of controlled vocabulary term list

      This class should be used to inherit from, to allow to add
      an arbitrary number of CV terms to the inherited class

      @ingroup Metadata
  */
  ///Representation of a CV term used by CVMappings
  class OPENMS_DLLAPI CVTermList :
    public MetaInfoInterface
  {
public:

    /// Defaults constructor
    CVTermList() = default;

    /// Copy constructor
    CVTermList(const CVTermList&) = default;

    // note: we implement the move constructor ourselves due to a bug in MSVS
    // 2015/2017 which cannot produce a default move constructor for classes
    // that contain STL containers (other than vector).

    /// Move constructor
    CVTermList(CVTermList&&) noexcept;

    /// Destructor
    virtual ~CVTermList();

    /// Assignment operator
    CVTermList& operator=(const CVTermList& rhs) & = default;

    /// Move assignment operator
    CVTermList& operator=(CVTermList&&) & = default;

    /** @name Accessors
    */
    //@{
    /// sets the CV terms
    void setCVTerms(const std::vector<CVTerm>& terms);

    /// replaces the specified CV term
    void replaceCVTerm(const CVTerm& cv_term);

    /// replaces the specified CV terms using the given accession number
    void replaceCVTerms(const std::vector<CVTerm>& cv_terms, const String& accession);

    /// replaces all cv terms with a map (can be obtained via getCVTerms)
    void replaceCVTerms(const std::map<String, std::vector<CVTerm> >& cv_term_map);

    /// merges the given map into the member map, no duplicate checking
    void consumeCVTerms(const std::map<String, std::vector<CVTerm> >& cv_term_map);

    /// returns the accession string of the term
    const std::map<String, std::vector<CVTerm> >& getCVTerms() const;

    /// adds a CV term
    void addCVTerm(const CVTerm& term);

    /// checks whether the spellings of the CV terms stored are correct
    //bool checkCVTerms(const ControlledVocabulary& cv) const;

    /// corrects the CVTerm names, according to the loaded CV
    //void correctCVTermNames();
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVTermList& cv_term_list) const;

    /// inequality operator
    bool operator!=(const CVTermList& cv_term_list) const;

    /// checks whether the term has a value
    bool hasCVTerm(const String& accession) const;

    /// checks whether the stored terms fulfill a given CVMappingRule
    /// TODO : implement
    //bool checkCVTerms(const CVMappingRule & rule, const ControlledVocabulary & cv) const;

    /// return true if no terms are available
    bool empty() const;
    //}

protected:

  std::map<String, std::vector<CVTerm> > cv_terms_;

  };

} // namespace OpenMS

