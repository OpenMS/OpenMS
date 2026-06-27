// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "Util.h"

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include <algorithm>
#include <cmath>
#include <vector>

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
  // Prefer FeatureFinderIdentification's observed mass-trace deviation:
  // "masserror_ppm" is a per-isotope getPPM(observed, theoretical) from the DIA
  // mass-difference score (-1 = "no signal" sentinel). It is the identified
  // feature's REAL MS1 mass accuracy -- the quantity FlashLFQ calibrates its ppm
  // distribution from. feature.getMZ() is FFID's THEORETICAL assay m/z (FFID seeds
  // the position m/z), so getPPM(feature.getMZ(), theoretical) is identically 0
  // and the calibration distribution collapses, leaving the mass feature dead.
  // (Requires masserror_ppm to survive ProteomicsLFQ's feature-meta strip, where
  // it is whitelisted in keep_meta.)
  // A real mass error can be more negative than -1 ppm, so the "no signal"
  // sentinel must be matched EXACTLY (-1), not as "<= -1". (DIAScoring pushes a
  // literal -1.0 for an isotope with no observed signal.)
  auto valid = [](double e) { return std::isfinite(e) && e != -1.0; };
  if (feature.metaValueExists("masserror_ppm"))
  {
    const DataValue dv = feature.getMetaValue("masserror_ppm");
    if (dv.valueType() == DataValue::DOUBLE_LIST)
    {
      std::vector<double> errs = dv.toDoubleList();
      if (! errs.empty() && valid(errs.front())) { return errs.front(); } // monoisotopic
      // monoisotopic missing -> robust median of the remaining valid isotopes
      std::vector<double> rest;
      for (double e : errs) { if (valid(e)) { rest.push_back(e); } }
      if (! rest.empty())
      {
        std::sort(rest.begin(), rest.end());
        return rest[rest.size() / 2];
      }
    }
    else if (dv.valueType() == DataValue::DOUBLE_VALUE)
    {
      const double e = dv;
      if (valid(e)) { return e; }
    }
  }

  // Fallbacks when masserror_ppm is unavailable (e.g. a non-FFID feature path):
  // the identification's OBSERVED precursor m/z vs theoretical, else feature.getMZ()
  // (which is theoretical for FFID, giving a 0 error -- the dead-feature case).
  const PeptideIdentificationList& peps = feature.getPeptideIdentifications();
  if (peps.empty() || peps[0].getHits().empty()) { return {}; }
  const PeptideHit& hit = peps[0].getHits()[0];
  const double theoretical = hit.getSequence().getMZ(hit.getCharge());
  double observed = peps[0].hasMZ() ? peps[0].getMZ() : feature.getMZ();
  if (! std::isfinite(observed) || observed <= 0.0) { observed = feature.getMZ(); }
  return Math::getPPM(observed, theoretical);
}

/**************************************************************************/
std::optional<std::vector<double>> feature_obs_envelope(const Feature& feature)
{
  if (! feature.metaValueExists("pipecho_obs_envelope")) { return {}; }
  const DataValue dv = feature.getMetaValue("pipecho_obs_envelope");
  if (dv.valueType() != DataValue::DOUBLE_LIST) { return {}; }
  std::vector<double> env = dv.toDoubleList();
  if (env.empty()) { return {}; }
  return env;
}

/**************************************************************************/
std::optional<double> feature_ion_mobility(const Feature& feature)
{
  if (feature.metaValueExists("IM_median"))
  {
    return (double)feature.getMetaValue("IM_median");
  }
  if (feature.metaValueExists(Constants::UserParam::IM))
  {
    return (double)feature.getMetaValue(Constants::UserParam::IM);
  }
  return {};
}

/**************************************************************************/
std::optional<double> feature_im_width(const Feature& feature)
{
  if (! feature.metaValueExists("IM_min") || ! feature.metaValueExists("IM_max"))
  {
    return {};
  }
  double width = (double)feature.getMetaValue("IM_max")
                 - (double)feature.getMetaValue("IM_min");
  if (! (width > 0)) return {}; // NaN-safe: only accept a strictly positive width
  return width;
}

} // namespace OpenMS::PipEcho::Util
