// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <map>

namespace OpenMS
{

  class CVTermList;

  /**
      @brief Interface to the controlled vocabulary term list

      This class is an interface to CVTermList, providing the same interface
      and it can be used instead (direct plug-in). It can be used to inherit
      from instead of CVTermList. The advantage is that this class only stores
      a single pointer to a CVTermList which may be more efficient for an empty
      list than storing the whole object.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI CVTermListInterface :
    public MetaInfoInterface
  {

  public:

    /** @name Constructors and Assignment
    */
    //@{
    // Constructor
    CVTermListInterface();
    /// Copy constructor
    CVTermListInterface(const CVTermListInterface & rhs);
    /// Move constructor
    CVTermListInterface(CVTermListInterface&&) noexcept;
    // Destructor (non virtual)
    ~CVTermListInterface();

    /// Assignment operator
    CVTermListInterface & operator=(const CVTermListInterface & rhs);
    /// Move assignment operator
    CVTermListInterface& operator=(CVTermListInterface&&) noexcept;
    //@}

    /// equality operator
    bool operator==(const CVTermListInterface& rhs) const;

    /// inequality operator
    bool operator!=(const CVTermListInterface& rhs) const;

    void replaceCVTerms(std::map<String, std::vector<CVTerm> > & cv_terms);

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

    /// checks whether the term has a value
    bool hasCVTerm(const String& accession) const;

    bool empty() const;

  private:

    void createIfNotExists_(); 

    CVTermList* cvt_ptr_;
  };

} // namespace OpenMS

