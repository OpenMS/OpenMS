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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/QUANTITATION/FeatureFindingIntact.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------
/**
    @page TOPP_FeatureFinderIntact TOPP_FeatureFinderIntact

    @brief TOPP_FeatureFinderIntact The intact protein feature detection for quantification (centroided).
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderIntact :
    public TOPPBase,
    public ProgressLogger
{
public:
  TOPPFeatureFinderIntact():
    TOPPBase("FeatureFinderIntact", "The intact protein feature detection for quantification", false, {}, false), ProgressLogger()
  {
    this->setLogType(CMD);
  }

private:
  double local_rt_range_;
  double local_mz_range_;
  double charge_lower_bound_;
  double charge_upper_bound_;
  bool use_smoothed_intensities_;

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

//    registerOutputFile_("out", "<file>", "", "FeatureXML file with metabolite features");
//    setValidFormats_("out", ListUtils::create<String>("featureXML"));
  }

public:
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    // TODO: need to update this value (came from FeatureFindingMetabo)
    local_mz_range_ = 6.5 ; // MZ range where to look for isotopic mass traces (-> decides size of isotopes =(local_mz_range_ * charge))
    local_rt_range_ = 15.0 ; // RT range where to look for coeluting mass traces
    charge_lower_bound_ = 7;
    charge_upper_bound_ = 30;
    use_smoothed_intensities_ = true; // for intensity of a mass trace

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    String in = getStringOption_("in");

    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    mz_data_file.getOptions().setMSLevels(ms_level);
    /// for test purpose : reduce in_ex
    mz_data_file.getOptions().setRTRange(DRange<1>(DPosition<1>(120.0), DPosition<1>(210)));
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                         " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = ms_peakmap[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__,
                                           "Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    // make sure the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    //-------------------------------------------------------------
    // Mass traces detection
    //-------------------------------------------------------------
    vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtdet.run(ms_peakmap, m_traces);

    //-------------------------------------------------------------
    // Elution peak detection
    //-------------------------------------------------------------
    std::vector<MassTrace> m_traces_final;
    std::vector<MassTrace> splitted_mtraces;
    ElutionPeakDetection epdet;
    // fill mass traces with smoothed data as well .. bad design..
    epdet.detectPeaks(m_traces, splitted_mtraces);
//    epdet.filterByPeakWidth(splitted_mtraces, m_traces_final);
    m_traces_final = splitted_mtraces;

    //-------------------------------------------------------------
    // building feature hypotheses
    //-------------------------------------------------------------
//    std::vector<FeatureHypothesis> feat_hypos;
//    build_feature_hypotheses_(m_traces_final, feat_hypos);
    FeatureFindingIntact ffi;
    ffi.updateMembers_(); // TODO: change this with param handler

    std::vector<FeatureFindingIntact::FeatureHypothesis> feat_hypos;
    ffi.buildFeatureHypotheses_(m_traces_final, feat_hypos);

    //-------------------------------------------------------------
    // clustering feature hypotheses
    //-------------------------------------------------------------

    //-------------------------------------------------------------
    // resolve conflicts in cluster
    //-------------------------------------------------------------

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
  }

  void build_feature_hypotheses_(vector<MassTrace>& input_mtraces, vector<FeatureHypothesis> output_hypotheses)
  {
    output_hypotheses.clear();
    if (input_mtraces.empty())
    {
      return;
    }
    this->startProgress(0, input_mtraces.size(), "assembling mass traces to features");

    // *********************************************************** //
    // Step 1 Preparation
    // *********************************************************** //

    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    // TODO: build isotope model for isotope ratio filtering

    double total_intensity(0.0);
    for (Size i = 0; i < input_mtraces.size(); ++i)
    {
      total_intensity += input_mtraces[i].getIntensity(use_smoothed_intensities_);
    }

    // *********************************************************** //
    // Step 2 Iterate through all mass traces to find likely matches
    // and generate isotopic / charge hypotheses
    // *********************************************************** //
    Size progress(0);
    for (SignedSize i = 0; i < (SignedSize)input_mtraces.size(); ++i)
    {
      this->setProgress(progress);

      ++progress;

      std::vector<const MassTrace*> local_traces;
      double ref_trace_mz(input_mtraces[i].getCentroidMZ());
      double ref_trace_rt(input_mtraces[i].getCentroidRT());

      local_traces.push_back(&input_mtraces[i]);

      for (Size ext_idx = i + 1; ext_idx < input_mtraces.size(); ++ext_idx)
      {
        // traces are sorted by m/z, so we can break when we leave the allowed window
        double diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
        if (diff_mz > local_mz_range_) break;

        double diff_rt = std::fabs(input_mtraces[ext_idx].getCentroidRT() - ref_trace_rt);
        if (diff_rt <= local_rt_range_)
        {
          // std::cout << " accepted!" << std::endl;
          local_traces.push_back(&input_mtraces[ext_idx]);
        }
      }
      findLocalFeatures_(local_traces, total_intensity, output_hypotheses);
    }
    this->endProgress();
    OPENMS_LOG_INFO << "feature _ hypos size:" << output_hypotheses.size() << std::endl;

    // sort feature candidates by their score
    std::sort(output_hypotheses.begin(), output_hypotheses.end(), CmpHypothesesByScore());

    // *********************************************************** //
    // Step 3 Iterate through all hypotheses, starting with the highest
    // scoring one. Remove hypotheses sharing the same mono-iso traces
    // *********************************************************** //

  }

  void findLocalFeatures_(const std::vector<const MassTrace*>& candidates, const double total_intensity, std::vector<FeatureHypothesis>& output_hypotheses) const
  {
    // not storing hypothesis with only one mass trace (with only mono), while FeatureFindingMetabo does

    for (Size charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
      FeatureHypothesis fh_tmp;
      fh_tmp.addMassTrace(*candidates[0]); // ref_mtrace (which is mono here)
      fh_tmp.setScore((candidates[0]->getIntensity(use_smoothed_intensities_)) / total_intensity);

      Size last_iso_idx(0); // largest index of found iso index
      Size iso_pos_max(static_cast<Size>(std::floor(charge * local_mz_range_))); // TODO: change this w/ mass_upper_bound ?
      for (Size iso_pos = 1; iso_pos <= iso_pos_max; ++iso_pos)
      {
        // expected m/z window for iso_pos -> 13C isotope peak position

        // Find mass trace that best agrees with current hypothesis of charge & isotopic position
        double best_so_far(0.0);
        Size best_idx(0);
        for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
        {
          // Score current mass trace candidates against hypothesis
          double rt_score(scoreRT_(*candidates[0], *candidates[mt_idx]));
          double mz_score(scoreMZ_(*candidates[0], *candidates[mt_idx], iso_pos, charge));

          // TODO : change this with stored model intensities
          double int_score(1.0);
          std::vector<double> tmp_ints(fh_tmp.getAllIntensities());
          tmp_ints.push_back(candidates[mt_idx]->getIntensity(use_smoothed_intensities_));
          int_score = computeAveragineSimScore_(tmp_ints, candidates[mt_idx]->getCentroidMZ() * charge);

          double total_pair_score(0.0);
          if (rt_score > 0.0 && mz_score > 0.0 && int_score > 0.0)
          {
            total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(int_score));
          }
          if (total_pair_score > best_so_far)
          {
            best_so_far = total_pair_score;
            best_idx = mt_idx;
          }
        } // end of mt_idx

        // Store mass trace that best agrees with current hypothesis of charge and isotopic position
        if (best_so_far > 0.0)
        {
          fh_tmp.addMassTrace(*candidates[best_idx]);
          double weighted_score(((candidates[best_idx]->getIntensity(use_smoothed_intensities_)) * best_so_far) / total_intensity);

          fh_tmp.setScore(fh_tmp.getScore() + weighted_score);
          fh_tmp.setCharge(charge);
          last_iso_idx = best_idx;

          output_hypotheses.push_back(fh_tmp);
        }
        else {break;}
      } // end for iso_pos
    } // end of charge

  } // end of findLocalFeatures_(...)

  double scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    std::map<double, std::vector<double> > coinciding_rts;

    std::pair<Size, Size> tr1_fwhm_idx(tr1.getFWHMborders());
    std::pair<Size, Size> tr2_fwhm_idx(tr2.getFWHMborders());

    double tr1_length(tr1.getFWHM());
    double tr2_length(tr2.getFWHM());
    double max_length = (tr1_length > tr2_length) ? tr1_length : tr2_length;

    // Extract peak shape between FWHM borders for both peaks
    for (Size i = tr1_fwhm_idx.first; i <= tr1_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr1[i].getRT()].push_back(tr1[i].getIntensity());
    }
    for (Size i = tr2_fwhm_idx.first; i <= tr2_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr2[i].getRT()].push_back(tr2[i].getIntensity());
    }

    // Look at peaks at the same RT
    std::vector<double> x, y, overlap_rts;
    for (std::map<double, std::vector<double> >::const_iterator m_it = coinciding_rts.begin(); m_it != coinciding_rts.end(); ++m_it)
    {
      if (m_it->second.size() == 2)
      {
        x.push_back(m_it->second[0]);
        y.push_back(m_it->second[1]);
        overlap_rts.push_back(m_it->first);
      }
    }

    double overlap(0.0);
    if (overlap_rts.size() > 0)
    {
      double start_rt(*(overlap_rts.begin())), end_rt(*(overlap_rts.rbegin()));
      overlap = std::fabs(end_rt - start_rt);
    }

    double proportion(overlap / max_length);
    if (proportion < 0.7)
    {
      return 0.0;
    }
    return computeCosineSim_(x, y);
  }

  double scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const
  {
    double diff_mz(std::fabs(tr2.getCentroidMZ() - tr1.getCentroidMZ()));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));

    double mz_score(0.0);
    /// mz scoring by expected mean w/ C13
    double mu = (Constants::C13C12_MASSDIFF_U * iso_pos) / charge; // using '1.0033548378'
    double sd = (0.0016633 * iso_pos - 0.0004751) / charge;
    double sigma_mult(3.0);

    //standard deviation including the estimated isotope deviation
    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));

    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);
    }

    return mz_score;
  }

  double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);
    double y_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
      y_squared_sum += y[i] * y[i];
    }

    double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));
    return (denom > 0.0) ? mixed_sum / denom : 0.0;
  }

  double computeAveragineSimScore_(const std::vector<double>& hypo_ints, const double& mol_weight) const
  {
    CoarseIsotopePatternGenerator solver(hypo_ints.size());
    auto isodist = solver.estimateFromPeptideWeight(mol_weight);

    IsotopeDistribution::ContainerType averagine_dist = isodist.getContainer();
    double max_int(0.0), theo_max_int(0.0);
    for (Size i = 0; i < hypo_ints.size(); ++i)
    {
      if (hypo_ints[i] > max_int)
      {
        max_int = hypo_ints[i];
      }

      if (averagine_dist[i].getIntensity() > theo_max_int)
      {
        theo_max_int = averagine_dist[i].getIntensity();
      }
    }

    // compute normalized intensities
    std::vector<double> averagine_ratios, hypo_isos;
    for (Size i = 0; i < hypo_ints.size(); ++i)
    {
      averagine_ratios.push_back(averagine_dist[i].getIntensity() / theo_max_int);
      hypo_isos.push_back(hypo_ints[i] / max_int);
    }

    double iso_score = computeCosineSim_(averagine_ratios, hypo_isos);
    return iso_score;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderIntact tool;
  return tool.main(argc, argv);
}

/// @endcond