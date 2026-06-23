// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "Util.h"

#include <OpenMS/MATH/MathFunctions.h>

namespace OpenMS::PipEcho::Util
{

/**************************************************************************/
std::optional<PeptideHit> feature_hit(const Feature& feature)
{
  const PeptideIdentificationList& peps = feature.getPeptideIdentifications();

  // FIXME: Is this the correct way to know if a feature has an
  // ID?  There seems to be two ways to store an ID on a feature,
  // but most code I see uses the "old way".
  if (! peps.empty() && ! peps[0].getHits().empty())
  {
    return peps[0].getHits()[0];
  }

  return {};
}

/******************************************************************************/
bool feature_is_decoy(const Feature& feature)
{
  auto hit = feature_hit(feature);

  return hit.has_value()
         && hit->getTargetDecoyType() == PeptideHit::TargetDecoyType::DECOY;
}

/**************************************************************************/
std::optional<double> feature_mass_error(const Feature& feature)
{
  auto hit = Util::feature_hit(feature);
  if (! hit.has_value()) return {};
  double experimental = feature.getMZ();
  double theoretical = hit->getSequence().getMZ(hit->getCharge());
  return Math::getPPM(experimental, theoretical);
}

} // namespace OpenMS::PipEcho::Util
