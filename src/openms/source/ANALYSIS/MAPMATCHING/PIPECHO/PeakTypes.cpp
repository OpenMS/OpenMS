// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "PeakTypes.h"

namespace OpenMS::PipEcho
{

/******************************************************************************/
bool Acceptor::is_donor_compatible(const Donor& donor,
                                   const Window& window) const
{
  return donor.feature.getCharge() == this->feature.getCharge()
         && std::fabs(donor.feature.getRT() - this->feature.getRT())
              <= window.rt_tol
         && std::fabs(donor.feature.getMZ() - this->feature.getMZ())
              <= window.mz_tol;
}

/******************************************************************************/
std::optional<const Donor*> Acceptor::fetch_and_mark_best_donor()
{
  std::optional<scored_t> match;

  auto check = [&](scored_t& scored) {
    if (! match.has_value() || match->first.mbr_score < scored.first.mbr_score)
    {
      match = scored;
    }
  };

  std::for_each(targets.begin(), targets.end(), check);
  std::for_each(decoys.begin(), decoys.end(), check);

  if (match.has_value())
  {
    // FIXME: Set feature peak decoy flag!
    return match->second;
  }
  else { return {}; }
}

} // namespace OpenMS::PipEcho
