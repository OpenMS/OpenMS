// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

#include <map>
#include <vector>

namespace OpenMS
{
class IsotopeDistribution;
  
class OPENMS_DLLAPI FeatureFinderAlgorithmMetaboIdent :
  public DefaultParamHandler
{
public:
  struct Row 
  {
    String name;
    String formula;
    double mass;
    std::vector<Int> charges;
    std::vector<double> rts;
    std::vector<double> rt_ranges;
    std::vector<double> iso_distrib;
  };
  
  using MetaboIdentTable = std::vector<Row>;

  /// default constructor
  FeatureFinderAlgorithmMetaboIdent();

  PeakMap& getMSData() { return ms_data_; }
  const PeakMap& getMSData() const { return ms_data_; }

  const PeakMap& getChromatograms() const { return chrom_data_; }

  const TargetedExperiment& getLibrary() const { return library_; }

  TransformationDescription extractTransformations() const
  {
    TransformationDescription trafo;
    TransformationDescription::DataPoints points;
    for (FeatureMap::ConstIterator it = features.begin();
          it != features.end(); ++it)
    {
      TransformationDescription::DataPoint point;
      point.first = it->getMetaValue("expected_rt");
      point.second = it->getRT();
      point.note = it->getMetaValue("PeptideRef");
      points.push_back(point);
    }
    trafo.setDataPoints(points);
  }

protected:
  double rt_window_; ///< RT window width
  double mz_window_; ///< m/z window width
  bool mz_window_ppm_; ///< m/z window width is given in PPM (not Da)?

  double isotope_pmin_; ///< min. isotope probability for peptide assay
  Size n_isotopes_; ///< number of isotopes for peptide assay

  double peak_width_;
  double min_peak_width_;
  double signal_to_noise_;

  String elution_model_;

  // output file (before filtering)
  String candidates_out_;

  Size debug_level_;

  void updateMembers_() override;

  MSExperiment ms_data_;
  PeakMap chrom_data_; ///< accumulated chromatograms (XICs)

  MRMFeatureFinderScoring feat_finder_; ///< OpenSWATH feature finder

  TargetedExperiment library_; ///< accumulated assays for targets
  
  CoarseIsotopePatternGenerator iso_gen_; ///< isotope pattern generator
  std::map<String, double> isotope_probs_; ///< isotope probabilities of transitions
  std::map<String, double> target_rts_; ///< RTs of targets (assays)
  
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  typedef std::vector<Feature*> FeatureGroup; ///< group of (overlapping) features

  /// Boundaries for a mass trace in a feature
  struct MassTraceBounds
  {
    Size sub_index;
    double rt_min, rt_max, mz_min, mz_max;
  };

  /// Boundaries for all mass traces per feature
  typedef std::map<UInt64, vector<MassTraceBounds> > FeatureBoundsMap;

  /// Predicate for filtering features by overall quality
  struct FeatureFilterQuality
  {
    bool operator()(const Feature& feature)
    {
      return feature.metaValueExists("FFMetId_remove");
    }
  } feature_filter_;

  /// Comparison functor for features
  struct FeatureCompare
  {
    bool operator()(const Feature& f1, const Feature& f2)
    {
      const String& ref1 = f1.getMetaValue("PeptideRef");
      const String& ref2 = f2.getMetaValue("PeptideRef");
      if (ref1 == ref2)
      {
        return f1.getRT() < f2.getRT();
      }
      return ref1 < ref2;
    }
  } feature_compare_;
};

} // namespace OpenMS
