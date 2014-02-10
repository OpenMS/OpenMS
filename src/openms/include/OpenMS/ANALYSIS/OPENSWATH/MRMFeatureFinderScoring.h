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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFINDERSCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFINDERSCORING_H

#define USE_SP_INTERFACE

// Actual scoring
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgScoring.h>

// Kernel classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left, const std::pair<double, double>& right);

namespace OpenMS
{

  /**
  @brief The MRMFeatureFinder finds and scores peaks of transitions that co-elute.

  It does so using an internal peakpicker for each chromatogram and then
  creating consensus / meta-peaks (MRMFeatures) that contain the information of
  all corresponding chromatograms at the peak-position. It then goes on to
  score those MRMFeatures using different criteria described in the
  MRMScoring class.

  @htmlinclude OpenMS_MRMFeatureFinderScoring.parameters

  */
  class OPENMS_DLLAPI MRMFeatureFinderScoring :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:
    ///Type definitions
    //@{

    // All the filters expect MSSpectrum<PeakT>, thus we give it an "MSSpectrum"
    // but filled with Chromatogram Peaks.
    
    // this is the type in which we store the chromatograms for this analysis
    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightTargetedExperiment TargetedExpType;
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightProtein ProteinType;
    typedef OpenSwath::LightModification ModificationType;
    // a transition group holds the MSSpectra with the Chromatogram peaks from above
    typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; 
    typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;
    //@}

    /// Constructor
    MRMFeatureFinderScoring();

    /// Destructor
    ~MRMFeatureFinderScoring();

    /// Picker and prepare functions
    //@{
    /** @brief Pick features in one experiment containing chromatogram
     *
     * Function for for wrapping in Python, only uses OpenMS datastructures and
     * does not return the map.
     *
     * @param chromatograms The input chromatograms
     * @param output The output features with corresponding scores
     * @param transition_exp The transition list describing the experiment
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the transition list
     * @param swath_map Optional SWATH-MS (DIA) map corresponding from which the chromatograms were extracted
     *
    */
    void pickExperiment(MSExperiment<Peak1D> & chromatograms, FeatureMap<Feature>& output, TargetedExperiment& transition_exp,
                        TransformationDescription trafo, MSExperiment<Peak1D>& swath_map);

    /** @brief Pick features in one experiment containing chromatogram
     *
     * @param input The input chromatograms
     * @param output The output features with corresponding scores
     * @param transition_exp The transition list describing the experiment
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param swath_map Optional SWATH-MS (DIA) map corresponding from which
     *                  the chromatograms were extracted. Use empty map if no
     *                  data is available.
     * @param transition_group_map Output mapping of transition groups
     *
    */
    void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap<Feature>& output, OpenSwath::LightTargetedExperiment& transition_exp,
                        TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType& transition_group_map);

    /** @brief Prepares the internal mappings of peptides and proteins.
     *
     * Calling this method _is_ required before calling scorePeakgroups.
     *
     * @param transition_exp The transition list describing the experiment
     *
    */
    void prepareProteinPeptideMaps_(OpenSwath::LightTargetedExperiment& transition_exp);
    //@}

    /** @brief Score all peak groups of a transition group
     *
     * Iterate through all features found along the chromatograms of the
     * transition group and score each one individually.
     *
     * @param transition_group The MRMTransitionGroup to be scored (input)
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param swath_map Optional SWATH-MS (DIA) map corresponding from which
     *                  the chromatograms were extracted. Use empty map if no
     *                  data is available.
     * @param output The output features with corresponding scores (the found
     *               features will be added to this FeatureMap).
     *
    */
    void scorePeakgroups(MRMTransitionGroupType& transition_group, TransformationDescription & trafo,
                         OpenSwath::SpectrumAccessPtr swath_map, FeatureMap<Feature>& output);

    /** @brief Set the flag for strict mapping
    */
    void setStrictFlag(bool f)
    {
      strict_ = f;
    }

    /** @brief Add an MS1 map containing spectra
     *
     * For DIA (SWATH-MS), an optional MS1 map can be supplied which can be
     * used to extract precursor ion signal and provides additional scores. If
     * no MS1 map is provided, the respective scores are not calculated.
     *
     * @param ms1_map The raw mass spectrometric MS1 data
     *
    */
    void setMS1Map(OpenSwath::SpectrumAccessPtr ms1_map)
    {
      ms1_map_ = ms1_map;
    }

    /** @brief Map the chromatograms to the transitions.
     *
     * Map an input chromatogram experiment (mzML) and transition list (TraML)
     * onto each other when they share identifiers, e.g. if the transition id
     * is the same as the chromatogram native id.
     *
     * @param input The input chromatograms
     * @param transition_exp The transition list describing the experiment
     * @param transition_group_map Mapping of transition groups
     * @param trafo Optional transformation of the experimental retention time
     *              to the normalized retention time space used in the
     *              transition list.
     * @param rt_extraction_window The used retention time extraction window
     *
    */
    void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, OpenSwath::LightTargetedExperiment& transition_exp,
                                       TransitionGroupMapType& transition_group_map, TransformationDescription trafo, double rt_extraction_window);
private:

    /// Synchronize members with param class
    void updateMembers_();

    // parameters
    double rt_extraction_window_;
    double quantification_cutoff_;
    int stop_report_after_feature_;
    bool write_convex_hull_;
    bool strict_;

    // scoring parameters
    double rt_normalization_factor_;
    int add_up_spectra_;
    double spacing_for_spectra_resampling_;

    // members
    std::map<OpenMS::String, const PeptideType*> PeptideRefMap_;
    std::map<OpenMS::String, const ProteinType*> ProteinRefMap_;
    OpenSwath_Scores_Usage su_;
    OpenMS::DIAScoring diascoring_;
    OpenMS::EmgScoring emgscoring_;

    // data 
    OpenSwath::SpectrumAccessPtr ms1_map_;

  };
}

#undef run_identifier
#endif
