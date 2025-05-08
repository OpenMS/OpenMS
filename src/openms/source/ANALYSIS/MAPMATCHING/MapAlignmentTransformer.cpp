// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

using std::vector;

namespace OpenMS
{

  bool MapAlignmentTransformer::storeOriginalRT_(MetaInfoInterface& meta_info,
                                                 double original_rt)
  {
    if (meta_info.metaValueExists("original_RT")) return false;
    meta_info.setMetaValue("original_RT", original_rt);
    return true;
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    PeakMap& msexp, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    msexp.clearRanges();

    // Transform spectra
    for (PeakMap::iterator mse_iter = msexp.begin();
         mse_iter != msexp.end(); ++mse_iter)
    {
      double rt = mse_iter->getRT();
      if (store_original_rt) storeOriginalRT_(*mse_iter, rt);
      mse_iter->setRT(trafo.apply(rt));
    }

    // Also transform chromatograms
    for (Size i = 0; i < msexp.getNrChromatograms(); ++i)
    {
      MSChromatogram& chromatogram = msexp.getChromatogram(i);
      vector<double> original_rts;
      if (store_original_rt) original_rts.reserve(chromatogram.size());
      for (Size j = 0; j < chromatogram.size(); j++)
      {
        double rt = chromatogram[j].getRT();
        if (store_original_rt) original_rts.push_back(rt);
        chromatogram[j].setRT(trafo.apply(rt));
      }
      if (store_original_rt && !chromatogram.metaValueExists("original_rt"))
      {
        chromatogram.setMetaValue("original_rt", original_rts);
      }
    }

    msexp.updateRanges();
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    FeatureMap& fmap, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    for (vector<Feature>::iterator fmit = fmap.begin(); fmit != fmap.end();
         ++fmit)
    {
      applyToFeature_(*fmit, trafo, store_original_rt);
    }

    // adapt RT values of unassigned peptides:
    if (!fmap.getUnassignedPeptideIdentifications().empty())
    {
      transformRetentionTimes(fmap.getUnassignedPeptideIdentifications(), trafo,
                              store_original_rt);
    }
  }


  void MapAlignmentTransformer::applyToBaseFeature_(
    BaseFeature& feature, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    // transform feature position:
    double rt = feature.getRT();
    if (store_original_rt) storeOriginalRT_(feature, rt);
    feature.setRT(trafo.apply(rt));

    // adapt RT values of annotated peptides:
    if (!feature.getPeptideIdentifications().empty())
    {
      transformRetentionTimes(feature.getPeptideIdentifications(), trafo,
                              store_original_rt);
    }
  }


  void MapAlignmentTransformer::applyToFeature_(
    Feature& feature, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    applyToBaseFeature_(feature, trafo, store_original_rt);

    // loop over all convex hulls
    vector<ConvexHull2D>& convex_hulls = feature.getConvexHulls();
    for (vector<ConvexHull2D>::iterator chiter = convex_hulls.begin();
         chiter != convex_hulls.end(); ++chiter)
    {
      // transform all hull point positions within convex hull
      ConvexHull2D::PointArrayType points = chiter->getHullPoints();
      chiter->clear();
      for (ConvexHull2D::PointArrayType::iterator points_iter = points.begin();
           points_iter != points.end(); ++points_iter)
      {
        double rt = (*points_iter)[Feature::RT];
        (*points_iter)[Feature::RT] = trafo.apply(rt);
      }
      chiter->setHullPoints(points);
    }

    // recurse into subordinates
    for (vector<Feature>::iterator subiter = feature.getSubordinates().begin();
         subiter != feature.getSubordinates().end(); ++subiter)
    {
      applyToFeature_(*subiter, trafo, store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    ConsensusMap& cmap, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    for (ConsensusMap::Iterator cmit = cmap.begin(); cmit != cmap.end(); ++cmit)
    {
      applyToConsensusFeature_(*cmit, trafo, store_original_rt);
    }

    // adapt RT values of unassigned peptides:
    if (!cmap.getUnassignedPeptideIdentifications().empty())
    {
      transformRetentionTimes(cmap.getUnassignedPeptideIdentifications(), trafo,
                              store_original_rt);
    }
  }


  void MapAlignmentTransformer::applyToConsensusFeature_(
    ConsensusFeature& feature, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    applyToBaseFeature_(feature, trafo, store_original_rt);

    // apply to grouped features (feature handles):
    for (ConsensusFeature::HandleSetType::const_iterator it =
           feature.getFeatures().begin(); it != feature.getFeatures().end();
         ++it)
    {
      double rt = it->getRT();
      it->asMutable().setRT(trafo.apply(rt));
    }
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    vector<PeptideIdentification>& pep_ids,
    const TransformationDescription& trafo, bool store_original_rt)
  {
    for (vector<PeptideIdentification>::iterator pep_it = pep_ids.begin();
         pep_it != pep_ids.end(); ++pep_it)
    {
      if (pep_it->hasRT())
      {
        double rt = pep_it->getRT();
        if (store_original_rt) storeOriginalRT_(*pep_it, rt);
        pep_it->setRT(trafo.apply(rt));
      }
    }
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    IdentificationData& id_data, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    // update RTs in-place:
    id_data.applyToObservations([&](IdentificationData::Observation& obs)
      {
        if (store_original_rt)
        {
          storeOriginalRT_(obs, obs.rt);
        }
        obs.rt = trafo.apply(obs.rt);
      });
  }

}
