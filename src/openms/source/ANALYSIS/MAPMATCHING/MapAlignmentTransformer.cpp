// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using std::vector;

namespace OpenMS
{

  void MapAlignmentTransformer::checkInputSizes_(Size maps_size,
                                                 Size trafos_size)
  {
    if (trafos_size != maps_size)
    {
      String msg = "MapAlignmentTransformer expects one transformation per "
        "input map (got " + String(trafos_size) + " transformations and " + 
        String(maps_size) + " input maps)";
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       msg);
    }
  }


  bool MapAlignmentTransformer::storeOriginalRT_(MetaInfoInterface& interface,
                                                 double original_rt)
  {
    if (interface.metaValueExists("original_RT")) return false;
    interface.setMetaValue("original_RT", original_rt);
    return true;
  }


  void MapAlignmentTransformer::transformPeakMaps(
    vector<MSExperiment<> >& maps,
    const vector<TransformationDescription>& trafos, bool store_original_rt)
  {
    checkInputSizes_(maps.size(), trafos.size());

    vector<TransformationDescription>::const_iterator trafo_it = trafos.begin();
    for (vector<MSExperiment<> >::iterator map_it = maps.begin();
         map_it != maps.end(); ++map_it, ++trafo_it)
    {
      transformSinglePeakMap(*map_it, *trafo_it, store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformSinglePeakMap(
    MSExperiment<>& msexp, const TransformationDescription& trafo,
    bool store_original_rt)
  {
    msexp.clearRanges();

    // Transform spectra
    for (MSExperiment<>::iterator mse_iter = msexp.begin();
         mse_iter != msexp.end(); ++mse_iter)
    {
      double rt = mse_iter->getRT();
      if (store_original_rt) storeOriginalRT_(*mse_iter, rt);
      mse_iter->setRT(trafo.apply(rt));
    }

    // Also transform chromatograms
    for (Size i = 0; i < msexp.getNrChromatograms(); ++i)
    {
      MSChromatogram<ChromatogramPeak>& chromatogram = msexp.getChromatogram(i);
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


  void MapAlignmentTransformer::transformFeatureMaps(
    vector<FeatureMap>& maps, const vector<TransformationDescription>& trafos,
    bool store_original_rt)
  {
    checkInputSizes_(maps.size(), trafos.size());

    vector<TransformationDescription>::const_iterator trafo_it = trafos.begin();
    for (vector<FeatureMap>::iterator map_it = maps.begin(); 
         map_it != maps.end(); ++map_it, ++trafo_it)
    {
      transformSingleFeatureMap(*map_it, *trafo_it, store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformSingleFeatureMap(
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
      transformSinglePeptideIdentification(
        fmap.getUnassignedPeptideIdentifications(), trafo, store_original_rt);
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
      transformSinglePeptideIdentification(feature.getPeptideIdentifications(),
                                           trafo, store_original_rt);
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


  void MapAlignmentTransformer::transformConsensusMaps(
    vector<ConsensusMap>& maps,
    const vector<TransformationDescription>& trafos, bool store_original_rt)
  {
    checkInputSizes_(maps.size(), trafos.size());

    vector<TransformationDescription>::const_iterator trafo_it = trafos.begin();
    for (vector<ConsensusMap>::iterator map_it = maps.begin();
         map_it != maps.end(); ++map_it, ++trafo_it)
    {
      transformSingleConsensusMap(*map_it, *trafo_it, store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformSingleConsensusMap(
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
      transformSinglePeptideIdentification(
        cmap.getUnassignedPeptideIdentifications(), trafo, store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformPeptideIdentifications(
    vector<vector<PeptideIdentification> >& maps,
    const vector<TransformationDescription>& trafos, bool store_original_rt)
  {
    checkInputSizes_(maps.size(), trafos.size());

    vector<TransformationDescription>::const_iterator trafo_it = trafos.begin();
    for (vector<vector<PeptideIdentification> >::iterator pep_it = 
           maps.begin(); pep_it != maps.end(); ++pep_it, ++trafo_it)
    {
      transformSinglePeptideIdentification(*pep_it, *trafo_it,
                                           store_original_rt);
    }
  }


  void MapAlignmentTransformer::transformSinglePeptideIdentification(
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

}
