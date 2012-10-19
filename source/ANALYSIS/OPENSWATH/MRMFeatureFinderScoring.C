// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left, const std::pair<double, double>& right)
{
  return left.first < right.first;
}

namespace OpenMS
{

  MRMFeatureFinderScoring::MRMFeatureFinderScoring() :
    DefaultParamHandler("MRMFeatureFinderScoring"),
    ProgressLogger()
  {
    defaults_.setValue("stop_report_after_feature", -1, "Stop reporting after feature (ordered by quality; -1 means do not stop).");
    defaults_.setValue("rt_extraction_window", -1.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 500 means to extract around +/- 500 s of the expected elution). For this to work, the TraML input file needs to contain normalized RT values.");
    defaults_.setValue("rt_normalization_factor", 1.0, "The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)");
    defaults_.setValue("quantification_cutoff", 0.0, "Cutoff below which peaks should not be used for quantification any more", StringList::create("advanced"));
    defaults_.setMinFloat("quantification_cutoff", 0.0);
    defaults_.setValue("write_convex_hull", "false", "Whether to write out all points of all features into the featureXML", StringList::create("advanced"));
    defaults_.setValidStrings("write_convex_hull", StringList::create("true,false"));
    defaults_.setValue("add_up_spectra", 1, "Add up spectra around the peak apex (needs to be a non-even integer)", StringList::create("advanced"));
    defaults_.setMinInt("add_up_spectra", 1);
    defaults_.setValue("spacing_for_spectra_resampling", 0.005, "If spectra are to be added, use this spacing to add them up", StringList::create("advanced"));
    defaults_.setMinFloat("spacing_for_spectra_resampling", 0.0);

    defaults_.insert("TransitionGroupPicker:", MRMTransitionGroupPicker().getDefaults());

    defaults_.insert("DIAScoring:", DIAScoring().getDefaults());

    defaults_.insert("EMGScoring:", EmgScoring().getDefaults());

    // One can turn on / off each score individually
    Param scores_to_use;
    scores_to_use.setValue("use_shape_score", "true", "Use the shape score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_shape_score", StringList::create("true,false"));
    scores_to_use.setValue("use_coelution_score", "true", "Use the coelution score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_coelution_score", StringList::create("true,false"));
    scores_to_use.setValue("use_rt_score", "true", "Use the retention time score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_rt_score", StringList::create("true,false"));
    scores_to_use.setValue("use_library_score", "true", "Use the library score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_library_score", StringList::create("true,false"));
    scores_to_use.setValue("use_elution_model_score", "true", "Use the elution model (EMG) score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_elution_model_score", StringList::create("true,false"));
    scores_to_use.setValue("use_intensity_score", "true", "Use the intensity score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_intensity_score", StringList::create("true,false"));
    scores_to_use.setValue("use_nr_peaks_score", "true", "Use the number of peaks score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_nr_peaks_score", StringList::create("true,false"));
    scores_to_use.setValue("use_total_xic_score", "true", "Use the total XIC score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_total_xic_score", StringList::create("true,false"));
    scores_to_use.setValue("use_sn_score", "true", "Use the SN (signal to noise) score", StringList::create("advanced"));
    scores_to_use.setValidStrings("use_sn_score", StringList::create("true,false"));
    defaults_.insert("Scores:", scores_to_use);

    // write defaults into Param object param_
    defaultsToParam_();

    strict_ = true;
  }

  MRMFeatureFinderScoring::~MRMFeatureFinderScoring()
  {
  }

  void MRMFeatureFinderScoring::updateMembers_()
  {
    /*
    handle_params();
  }

  void MRMFeatureFinderScoring::handle_params()
  {
    */
    stop_report_after_feature_ = (int)param_.getValue("stop_report_after_feature");
    rt_extraction_window_ = (DoubleReal)param_.getValue("rt_extraction_window");
    rt_normalization_factor_ = (DoubleReal)param_.getValue("rt_normalization_factor");
    quantification_cutoff_ = (DoubleReal)param_.getValue("quantification_cutoff");
    write_convex_hull_ = param_.getValue("write_convex_hull").toBool();
    add_up_spectra_ = param_.getValue("add_up_spectra");
    spacing_for_spectra_resampling_ = param_.getValue("spacing_for_spectra_resampling");

    diascoring_.setParameters(param_.copy("DIAScoring:", true));
    emgscoring_.setFitterParam(param_.copy("EmgScoring:", true));

    use_coelution_score_     = param_.getValue("Scores:use_coelution_score").toBool();
    use_shape_score_         = param_.getValue("Scores:use_shape_score").toBool();
    use_rt_score_            = param_.getValue("Scores:use_rt_score").toBool();
    use_library_score_       = param_.getValue("Scores:use_library_score").toBool();
    use_elution_model_score_ = param_.getValue("Scores:use_elution_model_score").toBool();
    use_intensity_score_     = param_.getValue("Scores:use_intensity_score").toBool();
    use_total_xic_score_     = param_.getValue("Scores:use_total_xic_score").toBool();
    use_nr_peaks_score_      = param_.getValue("Scores:use_nr_peaks_score").toBool();
    use_sn_score_            = param_.getValue("Scores:use_sn_score").toBool();
  }

  void MRMFeatureFinderScoring::mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input,
                                                              TargetedExpType& transition_exp, TransitionGroupMapType& transition_group_map,
                                                              TransformationDescription trafo, double rt_extraction_window)
  {
    double rt_min, rt_max, expected_rt;
    trafo.invert();

    std::map<String, int> chromatogram_map;
    Size nr_chromatograms = input->getNrChromatograms();
    for (Size i = 0; i < input->getNrChromatograms(); i++)
    {
      chromatogram_map[input->getChromatogramNativeID(i)] = i;
    }

    // Iterate thorugh all transitions and store the transition with the
    // corresponding chromatogram in the corresponding transition group
    Size progress = 0;
    startProgress(0, nr_chromatograms, "Mapping transitions to chromatograms ");
    for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
    {
      // get the current transition and try to find the corresponding chromatogram
      const TransitionType* transition = &transition_exp.getTransitions()[i];
      if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end())
      {
        std::cerr << "Error: Transition " + transition->getNativeID() + " from group " +
        transition->getPeptideRef() + " does not have a corresponding chromatogram" << std::endl;
        if (strict_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                           "Error: Transition " + transition->getNativeID() + " from group " +
                                           transition->getPeptideRef() + " does not have a corresponding chromatogram");
        }
        continue;
      }
      MSChromatogram<ChromatogramPeak> chromatogram_old;
      OpenSwath::ChromatogramPtr cptr = input->getChromatogramById(chromatogram_map[transition->getNativeID()]);
      OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram_old, cptr);
      RichPeakChromatogram chromatogram;

      // Create the chromatogram information
      // Get the expected retention time, apply the RT-transformation
      // (which describes the normalization) and then take the difference.
      // Note that we inverted the transformation in the beginning because
      // we want to transform from normalized to real RTs here and not the
      // other way round.
      rt_max = rt_min = 0;
      expected_rt = PeptideRTMap_[transition->getPeptideRef()];
      double de_normalized_experimental_rt = trafo.apply(expected_rt);
      rt_max = de_normalized_experimental_rt + rt_extraction_window;
      rt_min = de_normalized_experimental_rt - rt_extraction_window;
      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
      {
        if (rt_extraction_window >= 0 && (it->getRT() < rt_min || it->getRT() > rt_max))
        {
          continue;
        }
        ChromatogramPeak peak;
        peak.setMZ(it->getRT());
        peak.setIntensity(it->getIntensity());
        chromatogram.push_back(peak);
      }
      if (chromatogram.empty())
      {
        std::cerr << "Error: Could not find any points for chromatogram " + transition->getNativeID() + \
        ". Maybe your retention time transformation is off?" << std::endl;
        if (strict_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                           "Error: Could not find any points for chromatogram " + transition->getNativeID() + \
                                           ". Maybe your retention time transformation is off?");
        }
      }
      chromatogram.setMetaValue("product_mz", transition->getProductMZ());
      chromatogram.setMetaValue("precursor_mz", transition->getPrecursorMZ());
      chromatogram.setNativeID(transition->getNativeID());

      // Create new transition group if there is none for this peptide
      if (transition_group_map.find(transition->getPeptideRef()) == transition_group_map.end())
      {
        MRMTransitionGroupType transition_group;
        transition_group.setTransitionGroupID(transition->getPeptideRef());
        transition_group_map[transition->getPeptideRef()] = transition_group;
      }

      // Now add the transition and the chromatogram to the group
      MRMTransitionGroupType& transition_group = transition_group_map[transition->getPeptideRef()];
      transition_group.addTransition(*transition, transition->getNativeID());
      transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());

      setProgress(++progress);
    }
    endProgress();

    // The assumption is that for each transition that is in the TargetedExperiment we have exactly one chromatogram
    for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); trgroup_it++)
    {
      if (trgroup_it->second.getChromatograms().size() > 0 && (trgroup_it->second.getChromatograms().size() != trgroup_it->second.getTransitions().size()))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error: Could not match all transition to all chromatograms:\nFor chromatogram " + \
                                         trgroup_it->second.getTransitionGroupID() + " I found " + String(trgroup_it->second.getChromatograms().size()) + \
                                         " chromatograms but " + String(trgroup_it->second.getTransitions().size()) + " transitions.");
      }
    }
  }

}
