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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumExtractor.h>

namespace OpenMS
{
  SpectrumExtractor::SpectrumExtractor() :
    DefaultParamHandler("SpectrumExtractor")
  {
    getDefaultParameters(defaults_);

    // write defaults into Param object param_
    defaultsToParam_();
  }

  SpectrumExtractor::~SpectrumExtractor() {}

  void SpectrumExtractor::setRTWindow(const double& rt_window)
  {
    rt_window_ = rt_window;
  }

  double SpectrumExtractor::getRTWindow() const
  {
    return rt_window_;
  }

  void SpectrumExtractor::setMinScore(const double& min_score)
  {
    min_score_ = min_score;
  }

  double SpectrumExtractor::getMinScore() const
  {
    return min_score_;
  }

  void SpectrumExtractor::setMinForwardMatch(const double& min_forward_match)
  {
    min_forward_match_ = min_forward_match;
  }

  double SpectrumExtractor::getMinForwardMatch() const
  {
    return min_forward_match_;
  }

  void SpectrumExtractor::setMinReverseMatch(const double& min_reverse_match)
  {
    min_reverse_match_ = min_reverse_match;
  }

  double SpectrumExtractor::getMinReverseMatch() const
  {
    return min_reverse_match_;
  }

  void SpectrumExtractor::setMZTolerance(const double& mz_tolerance)
  {
    mz_tolerance_ = mz_tolerance;
  }

  double SpectrumExtractor::getMZTolerance() const
  {
    return mz_tolerance_;
  }

  void SpectrumExtractor::setMZToleranceUnits(const String& mz_tolerance_units)
  {
    mz_tolerance_units_ = mz_tolerance_units;
  }

  String SpectrumExtractor::getMZToleranceUnits() const
  {
    return mz_tolerance_units_;
  }

  void SpectrumExtractor::setSGolayFrameLength(const UInt& sgolay_frame_length)
  {
    sgolay_frame_length_ = sgolay_frame_length;
  }

  UInt SpectrumExtractor::getSGolayFrameLength() const
  {
    return sgolay_frame_length_;
  }

  void SpectrumExtractor::setSGolayPolynomialOrder(const UInt& sgolay_polynomial_order)
  {
    sgolay_polynomial_order_ = sgolay_polynomial_order;
  }

  UInt SpectrumExtractor::getSGolayPolynomialOrder() const
  {
    return sgolay_polynomial_order_;
  }

  void SpectrumExtractor::setGaussWidth(const double& gauss_width)
  {
    gauss_width_ = gauss_width;
  }

  double SpectrumExtractor::getGaussWidth() const
  {
    return gauss_width_;
  }

  void SpectrumExtractor::setUseGauss(const bool& use_gauss)
  {
    use_gauss_ = use_gauss;
  }

  bool SpectrumExtractor::getUseGauss() const
  {
    return use_gauss_;
  }

  void SpectrumExtractor::setSignalToNoise(const double& signal_to_noise)
  {
    signal_to_noise_ = signal_to_noise;
  }

  double SpectrumExtractor::getSignalToNoise() const
  {
    return signal_to_noise_;
  }

  void SpectrumExtractor::setTICWeight(const double& tic_weight)
  {
    tic_weight_ = tic_weight;
  }

  double SpectrumExtractor::getTICWeight() const
  {
    return tic_weight_;
  }

  void SpectrumExtractor::setFWHMWeight(const double& fwhm_weight)
  {
    fwhm_weight_ = fwhm_weight;
  }

  double SpectrumExtractor::getFWHMWeight() const
  {
    return fwhm_weight_;
  }

  void SpectrumExtractor::setSNRWeight(const double& snr_weight)
  {
    snr_weight_ = snr_weight;
  }

  double SpectrumExtractor::getSNRWeight() const
  {
    return snr_weight_;
  }

  void SpectrumExtractor::updateMembers_()
  {
    rt_window_ = (double)param_.getValue("rt_window");
    min_score_ = (double)param_.getValue("min_score");
    min_forward_match_ = (double)param_.getValue("min_forward_match");
    min_reverse_match_ = (double)param_.getValue("min_reverse_match");
    mz_tolerance_ = (double)param_.getValue("mz_tolerance");
    mz_tolerance_units_ = (String)param_.getValue("mz_tolerance_units");

    sgolay_frame_length_ = (UInt)param_.getValue("sgolay_frame_length");
    sgolay_polynomial_order_ = (UInt)param_.getValue("sgolay_polynomial_order");
    gauss_width_ = (double)param_.getValue("gauss_width");
    use_gauss_ = (bool)param_.getValue("use_gauss").toBool();
    signal_to_noise_ = (double)param_.getValue("signal_to_noise");

    tic_weight_ = (double)param_.getValue("tic_weight");
    fwhm_weight_ = (double)param_.getValue("fwhm_weight");
    snr_weight_ = (double)param_.getValue("snr_weight");
  }

  void SpectrumExtractor::getDefaultParameters(Param& params)
  {
    params.clear();

    params.setValue("rt_window", 30, "Retention time window in seconds.");

    params.setValue("min_score", 0.7, "Minimum score.");
    params.setMinFloat("min_score", 0.0);
    params.setMaxFloat("min_score", 1.0);

    params.setValue("min_forward_match", 0.9, "Minimum forward match.");
    params.setMinFloat("min_forward_match", 0.0);
    params.setMaxFloat("min_forward_match", 1.0);

    params.setValue("min_reverse_match", 0.9, "Minimum reverse match.");
    params.setMinFloat("min_reverse_match", 0.0);
    params.setMaxFloat("min_reverse_match", 1.0);

    params.setValue("mz_tolerance", 0.1, "Mass to Charge tolerance.");

    params.setValue("mz_tolerance_units", "Da", "Mass to Charge tolerance units.");
    params.setValidStrings("mz_tolerance_units", ListUtils::create<String>("ppm,Da"));

    params.setValue("sgolay_frame_length", 15, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    params.setValue("sgolay_polynomial_order", 3, "Order of the polynomial that is fitted.");
    params.setValue("gauss_width", 0.2, "Gaussian width in Da or ppm, estimated peak size.");
    params.setValue("use_gauss", "true", "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)");
    params.setValidStrings("use_gauss", ListUtils::create<String>("false,true"));
    params.setValue("signal_to_noise", 1.0, "Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.");
    params.setMinFloat("signal_to_noise", 0.0);

    params.setValue("tic_weight", 1.0, "TIC weight when scoring spectra.");
    params.setMinFloat("tic_weight", 0.0);
    params.setValue("fwhm_weight", 1.0, "FWHM weight when scoring spectra.");
    params.setMinFloat("fwhm_weight", 0.0);
    params.setValue("snr_weight", 1.0, "SNR weight when scoring spectra.");
    params.setMinFloat("snr_weight", 0.0);
  }

  void SpectrumExtractor::pickSpectrum(const MSSpectrum& spectrum, MSSpectrum& picked_spectrum)
  {
    if (!spectrum.isSorted())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Spectrum must be sorted by position");
    }

    LOG_DEBUG << " ====  Picking spectrum " << spectrum.getNativeID() << " with " << spectrum.size() << " peaks ";
    if (spectrum.empty())
    {
        LOG_DEBUG << std::endl;
        LOG_DEBUG << " - Error: spectrum is empty, abort picking."  << std::endl;
        return;
    }
    LOG_DEBUG << "(start at RT " << spectrum[0].getMZ() << " to RT " << spectrum[spectrum.size() - 1].getMZ() << ") " << std::endl;

    // Smooth the spectrum
    MSSpectrum smoothed_spectrum = spectrum;
    if (!use_gauss_)
    {
      SavitzkyGolayFilter sgolay;
      Param filter_parameters = sgolay.getParameters();
      filter_parameters.setValue("frame_length", sgolay_frame_length_);
      filter_parameters.setValue("polynomial_order", sgolay_polynomial_order_);
      sgolay.setParameters(filter_parameters);
      sgolay.filter(smoothed_spectrum);
    }
    else
    {
      GaussFilter gauss;
      Param filter_parameters = gauss.getParameters();
      filter_parameters.setValue("gaussian_width", gauss_width_);
      gauss.setParameters(filter_parameters);
      gauss.filter(smoothed_spectrum);
    }

    // Find initial seeds (peak picking)
    Param pepi_param = PeakPickerHiRes().getDefaults();
    pepi_param.setValue("signal_to_noise", signal_to_noise_);
    // disable spacing constraints, since we're dealing with spectrum
    pepi_param.setValue("spacing_difference", 0.0);
    pepi_param.setValue("spacing_difference_gap", 0.0);
    pepi_param.setValue("report_FWHM", "true");
    pepi_param.setValue("report_FWHM_unit", "absolute");
    picked_spectrum.clear(true);
    PeakPickerHiRes pp;
    pp.setParameters(pepi_param);
    pp.pick(smoothed_spectrum, picked_spectrum);
    LOG_DEBUG << "Found " << picked_spectrum.size() << " spectrum peaks." << std::endl;
  }

  void SpectrumExtractor::annotateSpectrum(
    const std::vector<MSSpectrum>& spectra,
    const TargetedExperiment& targeted_exp,
    std::vector<MSSpectrum>& annotated_spectra
  )
  {
    std::vector<ReactionMonitoringTransition> transitions = targeted_exp.getTransitions();
    for (UInt i=0; i<spectra.size(); ++i)
    {
      MSSpectrum spectrum = spectra[i];
      // It is supposed to have RT in minutes in target list file, therefore we divide by 60.0
      const double spectrum_rt = spectrum.getRT() / 60.0;
      const double rt_left_lim = spectrum_rt - getRTWindow() / 60.0 / 2.0;
      const double rt_right_lim = spectrum_rt + getRTWindow() / 60.0 / 2.0;
      std::vector<Precursor> precursors = spectrum.getPrecursors();
      if (precursors.size() < 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Spectrum does not contain precursor info.");
      }
      const double spectrum_mz = precursors[0].getMZ();
      const double mz_left_lim = spectrum_mz - getMZTolerance();
      const double mz_right_lim = spectrum_mz + getMZTolerance();

      LOG_DEBUG << "[" << i << "]\trt: " << spectrum_rt << "\tmz: " << spectrum_mz << std::endl;

      for (UInt j=0; j<transitions.size(); ++j)
      {
        const double target_rt = targeted_exp.getPeptideByRef(transitions[j].getPeptideRef()).getRetentionTime();
        const double target_mz = transitions[j].getPrecursorMZ();

        if (target_rt >= rt_left_lim && target_rt <= rt_right_lim && target_mz >= mz_left_lim && target_mz <= mz_right_lim)
        {
          LOG_DEBUG << "target_rt: " << target_rt << "\ttarget_mz: " << target_mz << std::endl;
          LOG_DEBUG << "pushed thanks to transition: " << j << " with name: " << transitions[j].getPeptideRef() << std::endl << std::endl;
          spectrum.setName(transitions[j].getPeptideRef());
          annotated_spectra.push_back(spectrum);
          break;
        }
      }
    }
    LOG_DEBUG << "the annotated variable has " << annotated_spectra.size() << " elements instead of " << spectra.size() << std::endl;
  }

  void SpectrumExtractor::scoreSpectrum(
    std::vector<MSSpectrum>& annotated,
    std::vector<MSSpectrum>& picked,
    std::vector<MSSpectrum>& scored
  )
  {
    std::vector<double> scores(annotated.size());
    for (UInt i=0; i<annotated.size(); ++i)
    {
      double total_tic = 0;
      for (UInt j=0; j<annotated[i].size(); ++j)
      {
        total_tic += annotated[i][j].getIntensity();
      }

      double avgFWHM = 0;
      for (UInt j=0; j<picked[i].getFloatDataArrays()[0].size(); ++j)
      {
        avgFWHM += picked[i].getFloatDataArrays()[0][j];
      }
      avgFWHM /= picked[i].getFloatDataArrays()[0].size();

      SignalToNoiseEstimatorMedian<MSSpectrum> sne;
      Param p;
      p.setValue("win_len", 40.0);
      p.setValue("noise_for_empty_window", 2.0);
      p.setValue("min_required_elements", 10);
      sne.setParameters(p);
      MSSpectrum::const_iterator it;
      sne.init(annotated[i].begin(),annotated[i].end());
      double avgSNR = 0.0;
      for (it=annotated[i].begin(); it!=annotated[i].end(); ++it)
      {
        avgSNR += sne.getSignalToNoise(it);
      }
      avgSNR /= annotated[i].size();

      double log10_total_tic = log10(total_tic);
      double inverse_avgFWHM = 1.0 / avgFWHM;
      double score = log10_total_tic * getTICWeight() + inverse_avgFWHM * getFWHMWeight() + avgSNR * getSNRWeight();

      MSSpectrum spectrum = annotated[i];
      spectrum.getFloatDataArrays().resize(5);
      spectrum.getFloatDataArrays()[1].setName("score");
      spectrum.getFloatDataArrays()[1].push_back(score);
      spectrum.getFloatDataArrays()[2].setName("log10_total_tic");
      spectrum.getFloatDataArrays()[2].push_back(log10_total_tic);
      spectrum.getFloatDataArrays()[3].setName("inverse_avgFWHM");
      spectrum.getFloatDataArrays()[3].push_back(inverse_avgFWHM);
      spectrum.getFloatDataArrays()[4].setName("avgSNR");
      spectrum.getFloatDataArrays()[4].push_back(avgSNR);
      scored.push_back(spectrum);
    }
  }
}
