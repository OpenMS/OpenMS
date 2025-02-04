// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer$
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Utilities for Feature and ConsensusMap

    @ingroup Datastructures
  */
  template <class MapType>
  class MapUtilities
  {
  public:
    /// applies a function on all PeptideHits or only assigned ones
    template <class T>
    void applyFunctionOnPeptideHits(T&& f, bool include_unassigned = true)
    {
      for (auto& feat : static_cast<MapType&>(*this))
      {
        applyFunctionOnPeptideHits_(feat.getPeptideIdentifications(), f);
      }
      if (include_unassigned)
      {
        applyFunctionOnPeptideHits_(static_cast<MapType&>(*this).getUnassignedPeptideIdentifications(), f);
      }
    }

    /// applies a function on all PeptideIDs or only assigned ones
    template <class T>
    void applyFunctionOnPeptideIDs(T&& f, bool include_unassigned = true)
    {
      for (auto& feat : static_cast<MapType&>(*this))
      {
        applyFunctionOnPeptideIDs_(feat.getPeptideIdentifications(), f);
      }
      if (include_unassigned)
      {
        applyFunctionOnPeptideIDs_(static_cast<MapType&>(*this).getUnassignedPeptideIdentifications(), f);
      }
    }

    /// applies a const function on all PeptideHits or only assigned ones
    template <class T>
    void applyFunctionOnPeptideHits(T&& f, bool include_unassigned = true) const
    {
      for (const auto& feat : static_cast<MapType const&>(*this))
      {
        applyFunctionOnPeptideHits_(feat.getPeptideIdentifications(), f);
      }
      if (include_unassigned)
      {
        applyFunctionOnPeptideHits_(static_cast<MapType const&>(*this).getUnassignedPeptideIdentifications(), f);
      }
    }

    /// applies a const function on all PeptideIDs or only assigned ones
    template <class T>
    void applyFunctionOnPeptideIDs(T&& f, bool include_unassigned = true) const
    {
      for (const auto& feat : static_cast<MapType const&>(*this))
      {
        applyFunctionOnPeptideIDs_(feat.getPeptideIdentifications(), f);
      }
      if (include_unassigned)
      {
        applyFunctionOnPeptideIDs_(static_cast<MapType const&>(*this).getUnassignedPeptideIdentifications(), f);
      }
    }

  private:
    template <class T>
    void applyFunctionOnPeptideIDs_(std::vector<PeptideIdentification>& idvec, T&& f)
    {
      for (auto& id : idvec)
      {
        f(id);
      }
    }

    template <class T>
    void applyFunctionOnPeptideHits_(std::vector<PeptideIdentification>& idvec, T&& f)
    {
      for (auto& id : idvec)
      {
        for (auto& hit : id.getHits())
        {
          f(hit);
        }
      }
    }

    template <class T>
    void applyFunctionOnPeptideIDs_(const std::vector<PeptideIdentification>& idvec, T&& f) const
    {
      for (const auto& id : idvec)
      {
        f(id);
      }
    }

    template <class T>
    void applyFunctionOnPeptideHits_(const std::vector<PeptideIdentification>& idvec, T&& f) const
    {
      for (const auto& id : idvec)
      {
        for (const auto& hit : id.getHits())
        {
          f(hit);
        }
      }
    }
  };

} // namespace OpenMS

