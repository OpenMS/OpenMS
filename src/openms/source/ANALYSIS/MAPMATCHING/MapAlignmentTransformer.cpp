// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    for (auto& mse : msexp)
    {
      double rt = mse.getRT();
      if (store_original_rt) storeOriginalRT_(mse, rt);
      mse.setRT(trafo.apply(rt));
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
    for (Feature& fm : fmap)
    {
      applyToFeature_(fm, trafo, store_original_rt);
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
    for (ConvexHull2D& ch : convex_hulls)
    {
      // transform all hull point positions within convex hull
      ConvexHull2D::PointArrayType points = ch.getHullPoints();
      ch.clear();
      for (ConvexHull2D::PointType& point : points)
      {
        double rt = point[Feature::RT];
        point[Feature::RT] = trafo.apply(rt);
      }
      ch.setHullPoints(points);
    }

    // recurse into subordinates
    for (Feature& sub : feature.getSubordinates())
    {
      applyToFeature_(sub, trafo, store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    ConsensusMap& cmap, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    for (ConsensusFeature& cm : cmap)
    {
      applyToConsensusFeature_(cm, trafo, store_original_rt);
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
    for (const FeatureHandle& handle : feature.getFeatures())
    {
      double rt = handle.getRT();
      handle.asMutable().setRT(trafo.apply(rt));
    }
  }


  void MapAlignmentTransformer::transformRetentionTimes(
    vector<PeptideIdentification>& pep_ids,
    const TransformationDescription& trafo, bool store_original_rt)
  {
    for (PeptideIdentification& pepid : pep_ids)
    {
      if (pepid.hasRT())
      {
        double rt = pepid.getRT();
        if (store_original_rt) storeOriginalRT_(pepid, rt);
        pepid.setRT(trafo.apply(rt));
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
