// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/ExposedVector.h>

namespace OpenMS
{
  /**
    @brief Container for peptide identifications from multiple spectra.
    
    This class represents a collection of PeptideIdentification objects, typically 
    originating from multiple spectra within a single identification run or experiment.
    
    PeptideIdentificationList is commonly used in scenarios where multiple spectra 
    have been searched against a database, resulting in a list of identifications that 
    need to be processed together. Examples include:
    
    - Storing all peptide identifications from an LC-MS/MS run
    - Collecting identifications for protein inference algorithms
    - Aggregating results from multiple search engines
    - Managing identifications in feature mapping workflows
    
    Note: The class uses composition over inheritance (via ExposedVector) to avoid the 
    well-known pitfalls of inheriting directly from STL containers
    
    Usage example:
    @code
    PeptideIdentificationList peptide_ids;
    
    // Add identifications from search results
    PeptideIdentification id1, id2;
    id1.setRT(1234.5);
    id1.setMZ(567.8);
    peptide_ids.push_back(id1);
    peptide_ids.push_back(id2);
    
    // Process all identifications
    for (auto& id : peptide_ids)
    {
      // Apply filtering, scoring, etc.
      IDFilter::keepNBestHits(id, 5);
    }
    
    // Use with algorithms that operate on collections
    IDFilter::keepNBestHits(peptide_ids, 1);
    BasicProteinInferenceAlgorithm algo;
    algo.run(peptide_ids, protein_ids);
    @endcode
    
    The container provides full vector-like functionality including element access,
    iteration, size management, and all standard algorithms.

    @see PeptideIdentification, ProteinIdentification, IDFilter, BasicProteinInferenceAlgorithm, ExposedVector
    
    @ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideIdentificationList : public ExposedVector<PeptideIdentification>
  {
  public:
    EXPOSED_VECTOR_INTERFACE(PeptideIdentification)
    
    /// @name Copy semantics
    //@{
    /// Copy constructor
    PeptideIdentificationList(const PeptideIdentificationList& other) = default;
    
    /// Copy assignment operator
    PeptideIdentificationList& operator=(const PeptideIdentificationList& other) = default;

    /// Move constructor
    PeptideIdentificationList(PeptideIdentificationList&& other) noexcept = default;
    //@}
    
    /// @name Additional constructors for std::vector compatibility
    //@{
    /// Constructor from std::vector
    PeptideIdentificationList(const std::vector<PeptideIdentification>& vec) 
      : ExposedVector<PeptideIdentification>(vec.begin(), vec.end()) {}
    
    /// Move constructor from std::vector
    PeptideIdentificationList(std::vector<PeptideIdentification>&& vec) 
      : ExposedVector<PeptideIdentification>(std::make_move_iterator(vec.begin()), 
                                            std::make_move_iterator(vec.end())) {}
    
    /// Constructor from initializer list
    PeptideIdentificationList(std::initializer_list<PeptideIdentification> init) 
      : ExposedVector<PeptideIdentification>(init.begin(), init.end()) {}
    //@}
    
    /// @name Assignment operators for std::vector compatibility
    //@{
    /// Assignment from std::vector
    PeptideIdentificationList& operator=(const std::vector<PeptideIdentification>& vec)
    {
      assign(vec.begin(), vec.end());
      return *this;
    }
    
    /// Move assignment from std::vector
    PeptideIdentificationList& operator=(std::vector<PeptideIdentification>&& vec)
    {
      data_ = std::move(vec);
      return *this;
    }
    
    /// Assignment from initializer list
    PeptideIdentificationList& operator=(std::initializer_list<PeptideIdentification> init)
    {
      assign(init.begin(), init.end());
      return *this;
    }
    //@}
    
    /// @name List-level score metadata methods
    //@{
    /// Returns the list-level score type
    const String& getScoreType() const;
    /// Sets the list-level score type
    void setScoreType(const String& type);
    
    /// Returns whether higher scores are better for the list
    bool isHigherScoreBetter() const;
    /// Sets whether higher scores are better for the list
    void setHigherScoreBetter(bool value);
    
    //@}
    
  private:
    String score_type_;
    bool higher_score_better_ = false;    
  };

} //namespace OpenMS
