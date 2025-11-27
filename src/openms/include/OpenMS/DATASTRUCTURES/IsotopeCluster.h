// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <set>

namespace OpenMS
{
  /**
    @brief Stores information about an isotopic cluster (i.e. potential peptide charge variants)
    
    An isotopic cluster represents a group of related peaks that likely originate from the same
    peptide but with different isotopic compositions. This structure stores the indices of these
    peaks and the scans they appear in, along with charge state information when available.
    
    The structure is typically used in mass spectrometry data analysis to group related peaks
    and track their charge states for further processing.
    
    @ingroup Datastructures
  */
  struct OPENMS_DLLAPI IsotopeCluster
  {
    /**
      @brief An index pair typically representing (scan_index, peak_index) in an MSExperiment
      
      The first value usually refers to the scan/spectrum index, while the second value
      refers to the peak index within that scan/spectrum.
    */
    typedef std::pair<Size, Size> IndexPair;
    
    /**
      @brief A set of index pairs, usually referring to peaks in an MSExperiment
      
      This collection stores unique pairs of indices that point to specific peaks
      in specific scans of a mass spectrometry experiment.
    */
    typedef std::set<IndexPair> IndexSet;

    /**
      @brief Index set with associated charge estimate
      
      Extends the basic IndexSet with charge state information for the peaks.
      This allows tracking which peaks belong to the same isotopic pattern
      and what charge state they represent.
    */
    struct ChargedIndexSet :
      public IndexSet
    {
      /**
        @brief Default constructor
        
        Initializes the charge to 0, which by convention means "no charge estimate"
      */
      ChargedIndexSet() :
        charge(0)
      {
      }

      /// Charge estimate (convention: zero means "no charge estimate")
      Int charge;
    };

    /**
      @brief Default constructor
      
      Initializes an empty isotope cluster with no peaks and no scans
    */
    IsotopeCluster() :
      peaks(),
      scans()
    {
    }

    /// Peaks in this cluster, with their charge state information
    ChargedIndexSet peaks;

    /// The scan indices where this cluster appears
    std::vector<Size> scans;
  };

} // namespace OPENMS

