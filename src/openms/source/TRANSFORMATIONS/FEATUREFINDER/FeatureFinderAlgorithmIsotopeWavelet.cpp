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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <iostream>
#include <algorithm>

namespace OpenMS
{
  FeatureFinderAlgorithmIsotopeWavelet::FeatureFinderAlgorithmIsotopeWavelet()
  {
    this->defaults_.setValue("max_charge", 3, "The maximal charge state to be considered.");
    this->defaults_.setMinInt("max_charge", 1);

    this->defaults_.setValue("intensity_threshold", -1., "The final threshold t' is build upon the formula: t' = av+t*sd, "
                                                         "where t is the intensity_threshold, av the average intensity within the wavelet transformed signal "
                                                         "and sd the standard deviation of the transform. "
                                                         "If you set intensity_threshold=-1, t' will be zero.\n"
                                                         "As the 'optimal' value for this parameter is highly data dependent, we would recommend to start "
                                                         "with -1, which will also extract features with very low signal-to-noise ratio. Subsequently, one "
                                                         "might increase the threshold to find an optimized trade-off between false positives and true positives. "
                                                         "Depending on the dynamic range of your spectra, suitable value ranges include: -1, [0:10], and if your data "
                                                         "features even very high intensity values, t can also adopt values up to around 30. "
                                                         "Please note that this parameter is not of an integer type, s.t. you can also use t:=0.1, e.g.");
    this->defaults_.setValue("intensity_type", "ref", "Determines the intensity type returned for the identified features. 'ref' (default) returns the sum of the intensities of each isotopic peak within an isotope pattern. 'trans' refers to the intensity of the monoisotopic peak within the wavelet transform. 'corrected' refers also to the transformed intensity with an attempt to remove the effects of the convolution. While the latter ones might be preferable for qualitative analyses, 'ref' might be the best option to obtain quantitative results. Please note that intensity values might be spoiled (in particular for the option 'ref'), as soon as patterns overlap (see also the explanations given in the class documentation of FeatureFinderAlgorihtmIsotopeWavelet).", ListUtils::create<String>("advanced"));
    this->defaults_.setValidStrings("intensity_type", ListUtils::create<String>("ref,trans,corrected"));

    this->defaults_.setValue("check_ppm", "false", "Enables/disables a ppm test vs. the averagine model, i.e. "
                                                   "potential peptide masses are checked for plausibility. In addition, a heuristic correcting potential mass shifts induced by the wavelet is applied.", ListUtils::create<String>("advanced"));
    this->defaults_.setValidStrings("check_ppm", ListUtils::create<String>("true,false"));

    this->defaults_.setValue("hr_data", "false", "Must be true in case of high-resolution data, i.e. "
                                                 "for spectra featuring large m/z-gaps (present in FTICR and Orbitrap data, e.g.). Please check "
                                                 "a single MS scan out of your recording, if you are unsure.");
    this->defaults_.setValidStrings("hr_data", ListUtils::create<String>("true,false"));

    this->defaults_.setValue("sweep_line:rt_votes_cutoff", 5, "Defines the minimum number of "
                                                              "subsequent scans where a pattern must occur to be considered as a feature.", ListUtils::create<String>("advanced"));
    this->defaults_.setMinInt("sweep_line:rt_votes_cutoff", 0);
    this->defaults_.setValue("sweep_line:rt_interleave", 1, "Defines the maximum number of "
                                                            "scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing. There is usually no reason to change the default value.", ListUtils::create<String>("advanced"));
    this->defaults_.setMinInt("sweep_line:rt_interleave", 0);

    this->defaultsToParam_();
  }

  FeatureFinderAlgorithmIsotopeWavelet::~FeatureFinderAlgorithmIsotopeWavelet()
  {
  }

  MSSpectrum* FeatureFinderAlgorithmIsotopeWavelet::createHRData(const UInt i)
  {
    MSSpectrum spec((*this->map_)[i]);

    const MSSpectrum& specr((*this->map_)[i]);

    for (UInt j = 0; j < spec.size() - 1; ++j)
    {
      spec[j].setMZ(-1 * (specr[j + 1].getMZ() - specr[j].getMZ()));
      spec[j].setIntensity((specr[j].getIntensity() + specr[j + 1].getIntensity()));
    }
    spec[spec.size() - 1].setMZ(-1); spec[spec.size() - 1].setIntensity(-1);

    ConstRefVector<MSSpectrum> c_sorted_spec(spec.begin(), spec.end());
    //Sort in ascending order according to the intensities present in the transform
    c_sorted_spec.sortByPosition();

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
    std::ofstream ofilex("spacings.trans");
    for (UInt j = 0; j < spec.size() - 1; ++j)
    {
      ofilex << ::std::setprecision(12) << std::fixed << spec[j].getMZ() << "\t" << spec[j].getIntensity() << std::endl;
    }
    ofilex.close();
#endif

    UInt pos = 0;
    while (c_sorted_spec[pos].getIntensity() <= 0)
    {
      if (++pos >= c_sorted_spec.size())
      {
        std::cout << "Detected empty scan or a scan that cannot be interpolated with zeros in HR mode. " << std::endl;
        std::cout << "Please check scan # " << i << " of your data set." << std::endl;
        exit(-1);
      }
    }
    double bound = -1 * c_sorted_spec[pos].getMZ();

    if (bound > (1. / max_charge_) / 2.)
    {
      //that might be case for simulated spectra,
      //which might show a very artificial spacing
      bound = (1. / max_charge_) / 2. / 4.;
    }

    MSSpectrum* new_spec = new MSSpectrum;
    new_spec->reserve(200000);
    new_spec->setRT(((*this->map_)[i]).getRT());
    PeakType p; p.setMZ(specr[0].getMZ()); p.setIntensity(specr[0].getIntensity());
    new_spec->push_back(p);

    UInt count;
    for (UInt j = 0; j < spec.size() - 1; ++j)
    {
      count = 0;
      while (-spec[j].getMZ() - count * bound > bound)
      {
        ++count;
        p.setMZ(specr[j].getMZ() + count * bound); p.setIntensity(0);
        new_spec->push_back(p);
      }
      p.setMZ(specr[j + 1].getMZ()); p.setIntensity(specr[j + 1].getIntensity());
      new_spec->push_back(p);
    }

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
    std::ofstream ofiley("new_spec.trans");
    for (UInt j = 0; j < new_spec->size(); ++j)
    {
      ofiley << ::std::setprecision(12) << std::fixed << (*new_spec)[j].getMZ() << "\t" << (*new_spec)[j].getIntensity() << std::endl;
    }
    ofiley.close();
#endif

    return new_spec;
  }

  void FeatureFinderAlgorithmIsotopeWavelet::run()
  {
    double max_mz = this->map_->getMax()[1];
    double min_mz = this->map_->getMin()[1];

    Size max_size = 0;

    //Check for useless RT_votes_cutoff_ parameter
    if (RT_votes_cutoff_ > this->map_->size())
    {
      real_RT_votes_cutoff_ = 0;
    }
    else
    {
      real_RT_votes_cutoff_ = RT_votes_cutoff_;
    }

    this->ff_->setLogType(ProgressLogger::CMD);
    progress_counter_ = 0;
    this->ff_->startProgress(0, 2 * this->map_->size() * max_charge_, "analyzing spectra");

    IsotopeWaveletTransform<PeakType>* iwt = new IsotopeWaveletTransform<PeakType>(min_mz, max_mz, max_charge_, max_size, hr_data_, intensity_type_);
    for (UInt i = 0; i < this->map_->size(); ++i)
    {
      const MSSpectrum& c_ref((*this->map_)[i]);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
      std::cout << ::std::fixed << ::std::setprecision(6) << "Spectrum " << i + 1 << " (" << (*this->map_)[i].getRT() << ") of " << this->map_->size() << " ... ";
      std::cout.flush();
#endif

      if (c_ref.size() <= 1)                 //unable to do transform anything
      {
#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
        std::cout << "scan empty or consisting of a single data point. Skipping." << std::endl;
#endif
        this->ff_->setProgress(progress_counter_ += 2);
        continue;
      }

      if (!hr_data_)                   //LowRes data
      {
        iwt->initializeScan((*this->map_)[i]);
        for (UInt c = 0; c < max_charge_; ++c)
        {
          MSSpectrum c_trans(c_ref);

          iwt->getTransform(c_trans, c_ref, c);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::stringstream stream;
          stream << "cpu_lowres_" << c_ref.getRT() << "_" << c + 1 << ".trans\0";
          std::ofstream ofile(stream.str().c_str());
          for (UInt k = 0; k < c_ref.size(); ++k)
          {
            ofile << ::std::setprecision(8) << std::fixed << c_trans[k].getMZ() << "\t" << c_trans[k].getIntensity() << "\t" << c_ref[k].getIntensity() << std::endl;
          }
          ofile.close();
#endif

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::cout << "transform O.K. ... "; std::cout.flush();
#endif
          this->ff_->setProgress(++progress_counter_);

          iwt->identifyCharge(c_trans, c_ref, i, c, intensity_threshold_, check_PPMs_);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::cout << "charge recognition O.K. ... "; std::cout.flush();
#endif
          this->ff_->setProgress(++progress_counter_);
        }
      }
      else                   //HighRes data
      {
        for (UInt c = 0; c < max_charge_; ++c)
        {
          MSSpectrum* new_spec = createHRData(i);
          iwt->initializeScan(*new_spec, c);
          MSSpectrum c_trans(*new_spec);

          iwt->getTransformHighRes(c_trans, *new_spec, c);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::stringstream stream;
          stream << "cpu_highres_" << new_spec->getRT() << "_" << c + 1 << ".trans\0";
          std::ofstream ofile(stream.str().c_str());
          for (UInt k = 0; k < new_spec->size(); ++k)
          {
            ofile << ::std::setprecision(8) << std::fixed << c_trans[k].getMZ() << "\t" << c_trans[k].getIntensity() << "\t" << (*new_spec)[k].getIntensity() << std::endl;
          }
          ofile.close();
#endif

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::cout << "transform O.K. ... "; std::cout.flush();
#endif
          this->ff_->setProgress(++progress_counter_);

          iwt->identifyCharge(c_trans, *new_spec, i, c, intensity_threshold_, check_PPMs_);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::cout << "charge recognition O.K. ... "; std::cout.flush();
#endif
          this->ff_->setProgress(++progress_counter_);

          delete (new_spec); new_spec = nullptr;
        }
      }


      iwt->updateBoxStates(*this->map_, i, RT_interleave_, real_RT_votes_cutoff_);
#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
      std::cout << "updated box states." << std::endl;
#endif

      std::cout.flush();
    }

    this->ff_->endProgress();

    //Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_
    iwt->updateBoxStates(*this->map_, INT_MAX, RT_interleave_, real_RT_votes_cutoff_);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
    std::cout << "Final mapping."; std::cout.flush();
#endif
    *this->features_ = iwt->mapSeeds2Features(*this->map_, real_RT_votes_cutoff_);

    delete (iwt);
  }

  const String FeatureFinderAlgorithmIsotopeWavelet::getProductName()
  {
    return "isotope_wavelet";
  }

  FeatureFinderAlgorithm* FeatureFinderAlgorithmIsotopeWavelet::create()
  {
    return new FeatureFinderAlgorithmIsotopeWavelet();
  }

  void FeatureFinderAlgorithmIsotopeWavelet::updateMembers_()
  {
    max_charge_ = this->param_.getValue("max_charge");
    intensity_threshold_ = this->param_.getValue("intensity_threshold");
    RT_votes_cutoff_ = this->param_.getValue("sweep_line:rt_votes_cutoff");
    RT_interleave_ = this->param_.getValue("sweep_line:rt_interleave");
    IsotopeWavelet::setMaxCharge(max_charge_);
    check_PPMs_ = ((String)(this->param_.getValue("check_ppm")) == "true");
    hr_data_ = ((String)(this->param_.getValue("hr_data")) == "true");
    intensity_type_ = ((String)(this->param_.getValue("intensity_type")));
  }

}
