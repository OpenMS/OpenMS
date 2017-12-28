// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_CONFIDENCESCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_CONFIDENCESCORING_H

#include <cmath> // for "exp"
#include <ctime> // for "time" (random number seed)
#include <limits> // for "infinity"
#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

namespace OpenMS
{

class OPENMS_DLLAPI ConfidenceScoring :
    public ProgressLogger
{
public:

    /// Constructor
    explicit ConfidenceScoring(bool test_mode_=false) :
      generator_(), rand_gen_(generator_, boost::uniform_int<>())
    {
      if (!test_mode_) rand_gen_.engine().seed(time(nullptr)); // seed with current time
    }

    virtual ~ConfidenceScoring() {}

protected:

    /// Mapping: Q3 m/z <-> transition intensity (maybe not unique!)
    typedef boost::bimap<double, boost::bimaps::multiset_of<double> > 
    BimapType;

    /// Binomial GLM
    struct
    {
      double intercept;
      double rt_coef;
      double int_coef;

      double operator()(double diff_rt, double dist_int)
      {
        double lm = intercept + rt_coef * diff_rt * diff_rt + 
          int_coef * dist_int;
        return 1.0 / (1.0 + exp(-lm));
      }
    } glm_;

    /// Helper for RT normalization (range 0-100)
    struct
    {
      double min_rt;
      double max_rt;
      
      double operator()(double rt)
      {
        return (rt - min_rt) / (max_rt - min_rt) * 100;
      }
    } rt_norm_;

    TargetedExperiment library_; // assay library

    IntList decoy_index_; // indexes of assays to use as decoys

    Size n_decoys_; // number of decoys to use (per feature/true assay)

    Map<String, IntList> transition_map_; // assay (ID) -> transitions (indexes)

    Size n_transitions_; // number of transitions to consider

    /// RT transformation to map measured RTs to assay RTs
    TransformationDescription rt_trafo_;

    boost::mt19937 generator_; // random number generation engine

    /// Random number generator (must be initialized in init. list of c'tor!)
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_gen_;

    /// Randomize the list of decoy indexes
    void chooseDecoys_();

    /// Manhattan distance
    double manhattanDist_(DoubleList x, DoubleList y);

    /// Get the retention time of an assay
    double getAssayRT_(const TargetedExperiment::Peptide& assay);

    /// Extract the @p n_transitions highest intensities from @p intensity_map,
    /// store them in @p intensities
    void extractIntensities_(BimapType& intensity_map, Size n_transitions,
                             DoubleList& intensities);

    /// Score the assay @p assay against feature data (@p feature_rt,
    /// @p feature_intensities), optionally using only the specified transitions
    /// (@p transition_ids)
    double scoreAssay_(const TargetedExperiment::Peptide& assay, 
                           double feature_rt, DoubleList& feature_intensities,
                           const std::set<String>& transition_ids = std::set<String>());

    /// Score a feature
    void scoreFeature_(Feature& feature);

public:

    void initialize(TargetedExperiment library, Size n_decoys, Size n_transitions, TransformationDescription rt_trafo)
    {
      library_ = TargetedExperiment(library); 
      n_decoys_ = n_decoys;
      n_transitions_ = n_transitions;
      rt_trafo_ = rt_trafo;
    }

    void initializeGlm(double intercept, double rt_coef, double int_coef)
    {
      glm_.intercept = intercept;
      glm_.rt_coef = rt_coef;
      glm_.int_coef = int_coef;
    }

    /**
      @brief Score a feature map -> make sure the class is properly initialized

      both functions initializeGlm and initialize need to be called first.

      The input to the program is 
      - a transition library which contains peptides with corresponding assays.
      - a feature map where each feature corresponds to an assay (mapped with
        MetaValue "PeptideRef") and each feature has as many subordinates as the
        assay has transitions (mapped with MetaValue "native_id").

    */
    void scoreMap(FeatureMap & features)
    {
      // are there enough assays in the library?
      Size n_assays = library_.getPeptides().size();
      if (n_assays < 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "There need to be at least 2 assays in the library for ConfidenceScoring.");

      }
      if (n_assays - 1 < n_decoys_)
      {
        LOG_WARN << "Warning: Parameter 'decoys' (" << n_decoys_ 
                 << ") is higher than the number of unrelated assays in the "
                 << "library (" << n_assays - 1 << "). "
                 << "Using all unrelated assays as decoys." << std::endl;
      }
      if (n_assays - 1 <= n_decoys_) n_decoys_ = 0; // use all available assays

      decoy_index_.resize(n_assays);
      for (Size i = 0; i < n_assays; ++i) decoy_index_[i] = boost::numeric_cast<Int>(i);

      // build mapping between assays and transitions:
      LOG_DEBUG << "Building transition map..." << std::endl;
      for (Size i = 0; i < library_.getTransitions().size(); ++i)
      {
        const String& ref = library_.getTransitions()[i].getPeptideRef();
        transition_map_[ref].push_back(boost::numeric_cast<Int>(i));
      }
      // find min./max. RT in the library:
      LOG_DEBUG << "Determining retention time range..." << std::endl;
      rt_norm_.min_rt = std::numeric_limits<double>::infinity();
      rt_norm_.max_rt = -std::numeric_limits<double>::infinity();
      for (std::vector<TargetedExperiment::Peptide>::const_iterator it = 
             library_.getPeptides().begin(); it != library_.getPeptides().end();
           ++it)
      {
        double current_rt = getAssayRT_(*it);
        if (current_rt == -1.0) continue; // indicates a missing value
        rt_norm_.min_rt = std::min(rt_norm_.min_rt, current_rt);
        rt_norm_.max_rt = std::max(rt_norm_.max_rt, current_rt);
      }

      // log scoring progress:
      LOG_DEBUG << "Scoring features..." << std::endl;
      startProgress(0, features.size(), "scoring features");

      for (FeatureMap::Iterator feat_it = features.begin(); 
           feat_it != features.end(); ++feat_it)
      {
        LOG_DEBUG << "Feature " << feat_it - features.begin() + 1 
                  << " (ID '" << feat_it->getUniqueId() << "')"<< std::endl;
        scoreFeature_(*feat_it);
        setProgress(feat_it - features.begin());
      }
      endProgress();

    }

};

}

#endif // OPENMS_ANALYSIS_OPENSWATH_CONFIDENCESCORING

