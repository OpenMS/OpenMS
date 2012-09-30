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

bool SortDoubleDoublePairFirst(const std::pair<double, double> & left, const std::pair<double, double> & right)
{
  return left.first < right.first;
}

namespace OpenMS
{

  MRMFeatureFinderScoring::MRMFeatureFinderScoring() :
    DefaultParamHandler("MRMFeatureFinderScoring"),
    ProgressLogger()
  {

    defaults_.setValue("sgolay_frame_length", 15, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("sgolay_polynomial_order", 3, "Order or the polynomial that is fitted.");
    defaults_.setValue("gauss_width", 50.0, "Gaussian width in seconds, estimated peak size.");

    defaults_.setValue("peak_width", 40.0, "Estimated peak width in seconds.");
    defaults_.setMinFloat("peak_width", 0.0);
    defaults_.setValue("signal_to_noise", 1.0, "Signal to noise.");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("sn_win_len", 1000.0, "Signal to noise window length.");
    defaults_.setValue("sn_bin_count", 30, "Signal to noise bin count.");

    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("stop_report_after_feature", -1, "Stop reporting after feature (ordered by quality; 1 means do not stop).");

    defaults_.setValue("rt_extraction_window", -1.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 500 means to extract around +/- 500 s of the expected elution).");
    defaults_.setValue("rt_normalization_factor", 1.0, "The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)");

    defaults_.setValue("emgfitter_maxiterations", 100, "Maximal number of iterators the EMG fitter is allowed to perform (influences performance)");

    defaults_.setValue("dia_byseries_ppm_diff", 10.0, "DIA b/y series minimal difference in ppm to consider.");
    defaults_.setMinFloat("dia_byseries_ppm_diff", 0.0);
    defaults_.setValue("dia_byseries_intensity_min", 300.0, "DIA b/y series minimum intensity to consider.");
    defaults_.setMinFloat("dia_byseries_intensity_min", 0.0);
    defaults_.setValue("dia_extraction_window", 0.05, "DIA extraction window in Th.");
    defaults_.setMinFloat("dia_extraction_window", 0.0);
    defaults_.setValue("dia_centroided", "false", "Use centroded DIA data.");
    defaults_.setValidStrings("dia_centroided", StringList::create("true,false"));
    defaults_.setValue("quantification_cutoff", 0.0, "Cutoff below which peaks should not be used for quantification any more");
    defaults_.setMinFloat("quantification_cutoff", 0.0);

    defaults_.setValue("use_shape_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_shape_score", StringList::create("true,false"));
    defaults_.setValue("use_coelution_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_coelution_score", StringList::create("true,false"));
    defaults_.setValue("use_rt_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_rt_score", StringList::create("true,false"));
    defaults_.setValue("use_library_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_library_score", StringList::create("true,false"));
    defaults_.setValue("use_elution_model_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_elution_model_score", StringList::create("true,false"));
    defaults_.setValue("use_intensity_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_intensity_score", StringList::create("true,false"));
    defaults_.setValue("use_nr_peaks_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_nr_peaks_score", StringList::create("true,false"));
    defaults_.setValue("use_total_xic_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_total_xic_score", StringList::create("true,false"));
    defaults_.setValue("use_sn_score", "true", "Use the retention time score", StringList::create("advanced"));
    defaults_.setValidStrings("use_sn_score", StringList::create("true,false"));

    defaults_.setValue("do_local_fdr", "false", "Use the local FDR score", StringList::create("advanced"));
    defaults_.setValidStrings("do_local_fdr", StringList::create("true,false"));
    defaults_.setValue("local_fdr_lib", "", "Library for the local FDR score", StringList::create("advanced"));

    defaults_.setValue("write_convex_hull", "false", "Whether to write out all points of all features into the featureXML", StringList::create("advanced"));
    defaults_.setValidStrings("write_convex_hull", StringList::create("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();

    strict = true;
  }

  MRMFeatureFinderScoring::~MRMFeatureFinderScoring()
  {
  }

  void MRMFeatureFinderScoring::handle_params()
  {
    sgolay_frame_length_ = (UInt)param_.getValue("sgolay_frame_length");
    sgolay_polynomial_order_ = (UInt)param_.getValue("sgolay_polynomial_order");
    gauss_width_ = (DoubleReal)param_.getValue("gauss_width");
    peak_width_ = (DoubleReal)param_.getValue("peak_width");
    signal_to_noise_ = (DoubleReal)param_.getValue("signal_to_noise");
    sn_win_len_ = (DoubleReal)param_.getValue("sn_win_len");
    sn_bin_count_ = (UInt)param_.getValue("sn_bin_count");

    stop_after_feature_ = (int)param_.getValue("stop_after_feature");
    stop_report_after_feature_ = (int)param_.getValue("stop_report_after_feature");
    stop_after_intensity_ratio_ = (DoubleReal)param_.getValue("stop_after_intensity_ratio");
    rt_extraction_window_ = (DoubleReal)param_.getValue("rt_extraction_window");

    emgfitter_maxiterations_ = (int)param_.getValue("emgfitter_maxiterations");

    dia_byseries_ppm_diff_ = (DoubleReal)param_.getValue("dia_byseries_ppm_diff");
    dia_byseries_intensity_min_ = (DoubleReal)param_.getValue("dia_byseries_intensity_min");
    dia_extract_window_ = (DoubleReal)param_.getValue("dia_extraction_window");
    dia_centroided_ = param_.getValue("dia_centroided").toBool();

    rt_normalization_factor_ = (DoubleReal)param_.getValue("rt_normalization_factor");
    quantification_cutoff_ = (DoubleReal)param_.getValue("quantification_cutoff");

    use_coelution_score_ = param_.getValue("use_coelution_score").toBool();
    use_shape_score_ = param_.getValue("use_shape_score").toBool();
    use_rt_score_ = param_.getValue("use_rt_score").toBool();
    use_library_score_ = param_.getValue("use_library_score").toBool();
    use_elution_model_score_ = param_.getValue("use_elution_model_score").toBool();
    use_intensity_score_ = param_.getValue("use_intensity_score").toBool();
    use_total_xic_score_ = param_.getValue("use_total_xic_score").toBool();
    use_nr_peaks_score_ = param_.getValue("use_nr_peaks_score").toBool();
    use_sn_score_ = param_.getValue("use_sn_score").toBool();

    do_local_fdr_ = param_.getValue("do_local_fdr").toBool();
    write_convex_hull_ = param_.getValue("write_convex_hull").toBool();

    Param fitter_param;
    fitter_param.setValue("max_iteration", emgfitter_maxiterations_);
    emgscoring.setFitterParam(fitter_param);
  }

  void MRMFeatureFinderScoring::mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input,
                                                              TargetedExpType & transition_exp, TransitionGroupMapType & transition_group_map,
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
      const TransitionType * transition = &transition_exp.getTransitions()[i];
      if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end())
      {
        std::cerr << "Error: Transition " + transition->getNativeID() + " from group " +
        transition->getPeptideRef() + " does not have a corresponding chromatogram" << std::endl;
        if (strict)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                           "Error: Transition " + transition->getNativeID() + " from group " +
                                           transition->getPeptideRef() + " does not have a corresponding chromatogram");
        }
        continue;
      }
      MSChromatogram<ChromatogramPeak> * chromatogram_old = new MSChromatogram<ChromatogramPeak>;
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
      expected_rt = PeptideRTMap[transition->getPeptideRef()];
      double de_normalized_experimental_rt = trafo.apply(expected_rt);
      rt_max = de_normalized_experimental_rt + rt_extraction_window;
      rt_min = de_normalized_experimental_rt - rt_extraction_window;
      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old->begin(); it != chromatogram_old->end(); ++it)
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
        if (strict)
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
      MRMTransitionGroupType & transition_group = transition_group_map[transition->getPeptideRef()];
      transition_group.addTransition(*transition, transition->getNativeID());
      transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());

      setProgress(++progress);
      delete chromatogram_old;
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
