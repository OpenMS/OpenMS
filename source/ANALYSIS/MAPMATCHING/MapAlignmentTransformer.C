// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using std::vector;


namespace OpenMS
{
  void MapAlignmentTransformer::transformPeakMaps(vector<MSExperiment<> > & maps,
                                                  const vector<TransformationDescription> & given_trafos)
  {
    typedef vector<MSExperiment<> >::iterator TMapIterator;
    typedef vector<TransformationDescription>::const_iterator TTrafoIterator;

    if (given_trafos.size() != maps.size())
    {
      throw Exception::IllegalArgument(__FILE__,
                                       __LINE__,
                                       __PRETTY_FUNCTION__,
                                       String("MapAlignmentTransformer expects one given transformation (got: ")
                                       + given_trafos.size()
                                       + ") per input map (got: "
                                       + maps.size()
                                       + "), these numbers are not equal");
    }

    TMapIterator mapIter = maps.begin();
    TTrafoIterator trafoIter = given_trafos.begin();

    for (; mapIter != maps.end() && trafoIter != given_trafos.end();
         ++mapIter,
         ++trafoIter)
    {
      transformSinglePeakMap(*mapIter, *trafoIter);
    }
  }

  void MapAlignmentTransformer::transformSinglePeakMap(MSExperiment<> & msexp,
                                                       const TransformationDescription & trafo)
  {
    msexp.clearRanges();

    // Transform spectra
    for (MSExperiment<>::iterator mse_iter = msexp.begin(); mse_iter != msexp.end(); ++mse_iter)
    {
      DoubleReal rt = mse_iter->getRT();
      mse_iter->setRT(trafo.apply(rt));
    }

    // Also transform chromatograms
    DoubleReal rt;
    std::vector<MSChromatogram<ChromatogramPeak> > chromatograms;
    for (Size i = 0; i < msexp.getChromatograms().size(); i++)
    {
      MSChromatogram<ChromatogramPeak> chromatogram = msexp.getChromatograms()[i];
      for (Size j = 0; j < chromatogram.size(); j++)
      {
        rt = chromatogram[j].getRT();
        chromatogram[j].setRT(trafo.apply(rt));
      }
      chromatograms.push_back(chromatogram);
    }
    msexp.setChromatograms(chromatograms);

    msexp.updateRanges();
  }

  void MapAlignmentTransformer::transformFeatureMaps(vector<FeatureMap<> > & maps,
                                                     const vector<TransformationDescription> & given_trafos)
  {
    typedef vector<FeatureMap<> >::iterator TFeatureMapsIterator;
    typedef vector<TransformationDescription>::const_iterator TTrafoIterator;

    if (given_trafos.size() != maps.size())
    {
      throw Exception::IllegalArgument(__FILE__,
                                       __LINE__,
                                       __PRETTY_FUNCTION__,
                                       String("MapAlignmentTransformer expects one given transformation (got: ")
                                       + given_trafos.size()
                                       + ") per input map (got: "
                                       + maps.size()
                                       + "), these numbers are not equal");
    }

    TFeatureMapsIterator fmapIter = maps.begin();
    TTrafoIterator trafoIter = given_trafos.begin();

    for (; fmapIter != maps.end() && trafoIter != given_trafos.end();
         ++fmapIter,
         ++trafoIter)
    {
      transformSingleFeatureMap(*fmapIter, *trafoIter);
    }
  }

  void MapAlignmentTransformer::transformSingleFeatureMap(FeatureMap<> & fmap,
                                                          const TransformationDescription & trafo)
  {
    for (vector<Feature>::iterator fmit = fmap.begin(); fmit != fmap.end(); ++fmit)
    {
      applyToFeature_(*fmit, trafo);
    }

    // adapt RT values of unassigned peptides:
    if (!fmap.getUnassignedPeptideIdentifications().empty())
    {
      transformSinglePeptideIdentification(
        fmap.getUnassignedPeptideIdentifications(), trafo);
    }
  }

  void MapAlignmentTransformer::applyToBaseFeature_(BaseFeature & feature,
                                                    const TransformationDescription & trafo)
  {
    // transform feature position:
    DoubleReal rt = feature.getRT();
    feature.setRT(trafo.apply(rt));

    // adapt RT values of annotated peptides:
    if (!feature.getPeptideIdentifications().empty())
    {
      transformSinglePeptideIdentification(feature.getPeptideIdentifications(),
                                           trafo);
    }
  }

  void MapAlignmentTransformer::applyToFeature_(Feature & feature,
                                                const TransformationDescription & trafo)
  {
    applyToBaseFeature_(feature, trafo);

    // loop over all convex hulls
    vector<ConvexHull2D> & convex_hulls = feature.getConvexHulls();
    for (vector<ConvexHull2D>::iterator chiter = convex_hulls.begin();
         chiter != convex_hulls.end(); ++chiter)
    {
      // transform all hull point positions within convex hull
      ConvexHull2D::PointArrayType points = chiter->getHullPoints();
      chiter->clear();
      for (ConvexHull2D::PointArrayType::iterator points_iter = points.begin();
           points_iter != points.end();
           ++points_iter
           )
      {
        DoubleReal rt = (*points_iter)[Feature::RT];
        (*points_iter)[Feature::RT] = trafo.apply(rt);
      }
      chiter->setHullPoints(points);
    }

    // recurse into subordinates
    for (vector<Feature>::iterator subiter = feature.getSubordinates().begin();
         subiter != feature.getSubordinates().end();
         ++subiter)
    {
      applyToFeature_(*subiter, trafo);
    }
  }

  void MapAlignmentTransformer::transformConsensusMaps(vector<ConsensusMap> & maps,
                                                       const vector<TransformationDescription> & given_trafos)
  {
    typedef vector<ConsensusMap>::iterator TMapIterator;
    typedef vector<TransformationDescription>::const_iterator TTrafoIterator;

    if (given_trafos.size() != maps.size())
    {
      throw Exception::IllegalArgument(__FILE__,
                                       __LINE__,
                                       __PRETTY_FUNCTION__,
                                       String("MapAlignmentTransformer expects one given transformation (got: ")
                                       + given_trafos.size()
                                       + ") per input map (got: "
                                       + maps.size()
                                       + "), these numbers are not equal");
    }

    TMapIterator mapIter = maps.begin();
    TTrafoIterator trafoIter = given_trafos.begin();

    for (; mapIter != maps.end() && trafoIter != given_trafos.end();
         ++mapIter,
         ++trafoIter)
    {
      transformSingleConsensusMap(*mapIter, *trafoIter);
    }
  }

  void MapAlignmentTransformer::transformSingleConsensusMap(ConsensusMap & cmap,
                                                            const TransformationDescription & trafo)
  {
    for (ConsensusMap::Iterator cmit = cmap.begin(); cmit != cmap.end();
         ++cmit)
    {
      applyToConsensusFeature_(*cmit, trafo);
    }

    // adapt RT values of unassigned peptides:
    if (!cmap.getUnassignedPeptideIdentifications().empty())
    {
      transformSinglePeptideIdentification(
        cmap.getUnassignedPeptideIdentifications(), trafo);
    }
  }

  void MapAlignmentTransformer::transformPeptideIdentifications(vector<vector<PeptideIdentification> > & maps,
                                                                const vector<TransformationDescription> & given_trafos)
  {
    typedef vector<vector<PeptideIdentification> >::iterator TPeptideIdentificationsIterator;
    typedef vector<TransformationDescription>::const_iterator TTrafoIterator;

    if (given_trafos.size() != maps.size())
    {
      throw Exception::IllegalArgument(__FILE__,
                                       __LINE__,
                                       __PRETTY_FUNCTION__,
                                       String("MapAlignmentTransformer expects one given transformation (got: ")
                                       + given_trafos.size()
                                       + ") per input map (got: "
                                       + maps.size()
                                       + "), these numbers are not equal");
    }

    TPeptideIdentificationsIterator pepIdentIter = maps.begin();
    TTrafoIterator trafoIter = given_trafos.begin();

    for (; pepIdentIter != maps.end() && trafoIter != given_trafos.end();
         ++pepIdentIter,
         ++trafoIter)
    {
      transformSinglePeptideIdentification(*pepIdentIter, *trafoIter);
    }
  }

  void MapAlignmentTransformer::transformSinglePeptideIdentification(vector<PeptideIdentification> & pepids,
                                                                     const TransformationDescription & trafo)
  {
    const UInt meta_index_RT = MetaInfo::registry().getIndex("RT");
    for (UInt pepid_index = 0; pepid_index < pepids.size(); ++pepid_index)
    {
      PeptideIdentification & pepid = pepids[pepid_index];
      DataValue dv = pepid.getMetaValue(meta_index_RT);
      if (dv != DataValue::EMPTY)
      {
        DoubleReal rt(dv);
        rt = trafo.apply(rt);
        pepid.setMetaValue(meta_index_RT, rt);
      }
    }

  }

  void MapAlignmentTransformer::applyToConsensusFeature_(ConsensusFeature & feature,
                                                         const TransformationDescription & trafo)
  {
    typedef ConsensusFeature::HandleSetType::const_iterator TConstHandleSetIterator;

    applyToBaseFeature_(feature, trafo);

    // apply to grouped features (feature handles):
    for (TConstHandleSetIterator it = feature.getFeatures().begin();
         it != feature.getFeatures().end();
         ++it)
    {
      DoubleReal rt = it->getRT();
      it->asMutable().setRT(trafo.apply(rt));
    }
  }

}
