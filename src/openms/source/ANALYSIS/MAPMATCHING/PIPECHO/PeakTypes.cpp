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
void Acceptor::push_back(const Score& score,
                         std::shared_ptr<const Donor> donor,
                         DonorType dtype)
{
  auto better = [&score](auto local) -> bool {
    return ! local.has_value() || local->first.mbr_score < score.mbr_score;
  };

  scored_t entry = std::make_pair(score, donor);

  switch (dtype)
  {
    case DonorType::Target:
      if (better(target)) target = entry;
      break;

    case DonorType::Decoy:
      if (better(decoy)) decoy = entry;
      break;
  }
}
} // namespace OpenMS::PipEcho
