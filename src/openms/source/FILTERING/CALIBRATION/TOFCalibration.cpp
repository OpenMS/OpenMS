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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>
#include <OpenMS/MATH/STATISTICS/QuadraticRegression.h>

namespace OpenMS
{

  TOFCalibration::TOFCalibration() :
    DefaultParamHandler("TOFCalibration"), ProgressLogger()
  {
    subsections_.push_back("PeakPicker");
    check_defaults_ = false;   // class has no own parameters
  }

  TOFCalibration::~TOFCalibration()
  {
  }

  void TOFCalibration::calculateCalibCoeffs_(MSExperiment<> & calib_spectra)
  {
    // flight times are needed later
    calib_peaks_ft_ = calib_spectra;


    // convert flight times of peaks into m/z values
    applyTOFConversion_(calib_spectra);
    std::vector<std::vector<unsigned int> > monoiso_peaks;
    getMonoisotopicPeaks_(calib_spectra, monoiso_peaks);

    startProgress(0, calib_spectra.size(), "quadratic fitting of calibrant spectra");
    // do the quadratic fitting for each calibration spectra separately
    for (unsigned int spec = 0; spec < calib_spectra.size(); ++spec)
    {
      std::vector<unsigned int> monoiso_peaks_scan;
      std::vector<double> exp_masses;
      // match the m/z-values to the expected masses
      matchMasses_(calib_spectra, monoiso_peaks, monoiso_peaks_scan, exp_masses, spec);

      // the actual quadratic fitting part
      Size n = exp_masses.size();
      if (n < 3)
      {
        continue;
      }

      // matrix containing the observations
      std::vector<double> x;
      // vector containing the expected masses
      std::vector<double> y;

      for (Size i = 0; i < n; i++)
      {
        // get the flight time
        double xi = ((calib_peaks_ft_.begin() + spec)->begin() + monoiso_peaks_scan[i])->getMZ();
        x.push_back(xi);
        y.push_back(exp_masses[i]);
      }

      Math::QuadraticRegression qr;
      qr.computeRegression(x.begin(), x.end(), y.begin());

#ifdef DEBUG_CALIBRATION
      std::cout << "chi^2: " << qr.getChiSquared() << std::endl;//DEBUG
      std::cout << "a: " << qr.getA() << "b: " << qr.getB()
            << "c: " << qr.getC() << std::endl;//DEBUG
#endif
      // store the coefficients
      coeff_quad_fit_.push_back(qr.getA());
      coeff_quad_fit_.push_back(qr.getB());
      coeff_quad_fit_.push_back(qr.getC());

      // determine the errors in ppm
      for (Size p = 0; p < n; ++p)
      {
#ifdef DEBUG_CALIBRATION
        std::cout << exp_masses[p]
                  << "\t" << mQ_(calib_peaks_ft_[spec][monoiso_peaks_scan[p]].getMZ(), spec) - exp_masses[p] << std::endl;
#endif
        errors_[exp_masses[p]].push_back((mQ_(calib_peaks_ft_[spec][monoiso_peaks_scan[p]].getMZ(), spec) - exp_masses[p]));
      }
      setProgress(spec);
    }
    endProgress();

    if (coeff_quad_fit_.empty())
    {
      String mess = String("Data can't be calibrated, not enough reference masses found: ") + coeff_quad_fit_.size() / 3;
      throw Exception::UnableToCalibrate(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToCalibrate", mess.c_str());
    }
    averageErrors_();
    averageCoefficients_();
  }

  void TOFCalibration::averageCoefficients_()
  {
    a_ = 0;
    b_ = 0;
    c_ = 0;
    for (unsigned int i = 0; i < coeff_quad_fit_.size(); i += 3)
    {
      a_ += coeff_quad_fit_[i];
      b_ += coeff_quad_fit_[i + 1];
      c_ += coeff_quad_fit_[i + 2];
    }
    a_ /= (coeff_quad_fit_.size() / 3);
    b_ /= (coeff_quad_fit_.size() / 3);
    c_ /= (coeff_quad_fit_.size() / 3);
  }

  void TOFCalibration::averageErrors_()
  {
    for (unsigned int p = 0; p < exp_masses_.size(); ++p)
    {
      // mean
      if (errors_[exp_masses_[p]].size() > 0)
      {
        double sum = 0;
        for (unsigned int i = 0; i < errors_[exp_masses_[p]].size(); ++i)
        {
          sum += errors_[exp_masses_[p]][i];

        }
        error_medians_.push_back(sum / errors_[exp_masses_[p]].size());
        calib_masses_.push_back(exp_masses_[p]);
      }
    }
  }

  void TOFCalibration::matchMasses_(MSExperiment<> & calib_peaks,
                                    std::vector<std::vector<unsigned int> > & monoiso_peaks,
                                    std::vector<unsigned int> & obs_masses,
                                    std::vector<double> & exp_masses, unsigned int idx)
  {
    for (unsigned int i = 0; i < monoiso_peaks[idx].size(); ++i)
    {
      for (unsigned int j = 0; j < exp_masses_.size(); ++j)
      {
        if (fabs(((calib_peaks.begin() + idx)->begin() + (monoiso_peaks[idx])[i])->getMZ() - exp_masses_[j]) < 1)
        {
          obs_masses.push_back((monoiso_peaks[idx])[i]);
          exp_masses.push_back(exp_masses_[j]);
          break;
        }
      }
    }
#ifdef DEBUG_CALIBRATION

    std::cout << "\n\n---------\nmatching monoisotopic peaks\n";

    for (unsigned int i = 0; i < obs_masses.size(); ++i)
    {
      std::cout << ((calib_peaks_ft_.begin() + idx)->begin() + obs_masses[i])->getMZ()
                << "\t" << exp_masses[i]
                << std::endl;

    }

#endif
  }

  void TOFCalibration::getMonoisotopicPeaks_(MSExperiment<> & calib_peaks, std::vector<std::vector<unsigned int> > & monoiso_peaks)
  {

    MSExperiment<>::iterator spec_iter = calib_peaks.begin();
    MSExperiment<>::SpectrumType::iterator peak_iter, help_iter;
#ifdef DEBUG_CALIBRATION
    spec_iter = calib_peaks.begin();
    std::cout << "\n\nbefore---------\n\n";
    // iterate through all spectra
    for (; spec_iter != calib_peaks.end(); ++spec_iter)
    {
      peak_iter = spec_iter->begin();
      // go through current scan
      for (; peak_iter != spec_iter->end(); ++peak_iter)
      {
        std::cout << peak_iter->getMZ() << std::endl;
      }
    }

#endif
    spec_iter = calib_peaks.begin();
    // iterate through all spectra
    for (; spec_iter != calib_peaks.end(); ++spec_iter)
    {
      peak_iter = spec_iter->begin();
      help_iter = peak_iter;
      std::vector<unsigned int> vec;
      // go through current scan
      while (peak_iter < spec_iter->end())
      {
        while (peak_iter + 1 < spec_iter->end() && ((peak_iter + 1)->getMZ() - peak_iter->getMZ() < 1.2))
        {
          ++peak_iter;
        }

        vec.push_back(distance(spec_iter->begin(), help_iter));

        help_iter = peak_iter + 1;
        ++peak_iter;

      }
      monoiso_peaks.push_back(vec);

    }

#ifdef DEBUG_CALIBRATION


    std::cout << "\n\nafter---------\n\n";

    for (unsigned int i = 0; i < monoiso_peaks.size(); ++i)
    {
      for (unsigned int j = 0; j < monoiso_peaks[i].size(); ++j)
      {
        std::cout << i << "\t" << ((calib_peaks.begin() + i)->begin() + (monoiso_peaks[i])[j])->getMZ() << std::endl;
      }
      std::cout << "--------------\n";

    }
    std::cout << "--------------\n\n\n";
#endif
  }

  void TOFCalibration::applyTOFConversion_(MSExperiment<> & calib_spectra)
  {
    MSExperiment<>::iterator spec_iter = calib_spectra.begin();
    MSExperiment<>::SpectrumType::iterator peak_iter;
    unsigned int idx = 0;

    //two point conversion
    if (ml3s_.empty())
    {
      for (; spec_iter != calib_spectra.end(); ++spec_iter)
      {
        peak_iter = spec_iter->begin();
        double ml1, ml2;
        if (ml1s_.size() == 1)
        {
          ml1 = ml1s_[0];
          ml2 = ml2s_[0];
        }
        else
        {
          ml1 = ml1s_[idx];
          ml2 = ml2s_[idx];
        }

        // go through current scan
        for (; peak_iter != spec_iter->end(); ++peak_iter)
        {
          double time = peak_iter->getMZ();
          peak_iter->setPos(ml1 / 1E12 * (time * 1000 - ml2));
        }
        ++idx;
      }
    }
    else
    {
      // three point conversion
      for (; spec_iter != calib_spectra.end(); ++spec_iter)
      {
        peak_iter = spec_iter->begin();
        double ml1, ml2, ml3;
        if (ml1s_.size() == 1)
        {
          ml1 = ml1s_[0];
          ml2 = ml2s_[0];
          ml3 = ml3s_[0];
        }
        else
        {
          ml1 = ml1s_[idx];
          ml2 = ml2s_[idx];
          ml3 = ml3s_[idx];
        }

        // go through current scan
        for (; peak_iter != spec_iter->end(); ++peak_iter)
        {
          double time = peak_iter->getMZ();
          peak_iter->setPos((-ml2 - (0.1E7 * (-5E5 + sqrt(0.25E12 - ml1 * ml2 * ml3 + ml1 * ml3 * time))) / (ml1 * ml3) + time) / ml3);
        }
        ++idx;
      }
    }

  }

} //namespace openms
