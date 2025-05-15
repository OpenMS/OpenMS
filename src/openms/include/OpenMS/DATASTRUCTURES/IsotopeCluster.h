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
  ///Stores information about an isotopic cluster (i.e. potential peptide charge variants)
  struct OPENMS_DLLAPI IsotopeCluster
  {
    /// An index e.g. in an MSExperiment
    typedef std::pair<Size, Size> IndexPair;
    /// A set of index pairs, usually referring to an MSExperiment.
    typedef std::set<IndexPair> IndexSet;

    ///index set with associated charge estimate
    struct ChargedIndexSet :
      public IndexSet
    {
      ChargedIndexSet() :
        charge(0)
      {
      }

      /// charge estimate (convention: zero means "no charge estimate")
      Int charge;
    };

    IsotopeCluster() :
      peaks(),
      scans()
    {
    }

    /// peaks in this cluster
    ChargedIndexSet peaks;

    /// the scans of this cluster
    std::vector<Size> scans;
  };

} // namespace OPENMS

