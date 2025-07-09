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
    
    This class represents a collection of PeptideIdentification objects,
    typically from multiple spectra in a single identification run.
    It provides all the functionality of std::vector<PeptideIdentification>
    while maintaining type safety and allowing for future extensions.
    
    Uses composition over inheritance to avoid the pitfalls of inheriting
    from STL containers.
    
    @ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideIdentificationList : public ExposedVector<PeptideIdentification>
  {
  public:
    EXPOSED_VECTOR_INTERFACE(PeptideIdentification)
    
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
    
    /// @name Conversion operators for compatibility with existing code
    //@{
    /// Implicit conversion to std::vector<PeptideIdentification>&
    operator std::vector<PeptideIdentification>&() 
    {
      return data_;
    }
    
    /// Implicit conversion to const std::vector<PeptideIdentification>&
    operator const std::vector<PeptideIdentification>&() const
    {
      return data_;
    }
    
    /// Explicit access to underlying vector for template compatibility
    std::vector<PeptideIdentification>& asVector()
    {
      return data_;
    }
    
    /// Explicit access to underlying vector for template compatibility (const)
    const std::vector<PeptideIdentification>& asVector() const
    {
      return data_;
    }
    //@}
  };

} //namespace OpenMS