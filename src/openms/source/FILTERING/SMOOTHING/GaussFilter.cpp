// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <cmath>

namespace OpenMS
{

  GaussFilter::GaussFilter() :
    ProgressLogger(),
    DefaultParamHandler("GaussFilter"),
    spacing_(0.01)
  {
    //Parameter settings
    defaults_.setValue("gaussian_width", 0.2, "Use a gaussian filter width which has approximately the same width as your mass peaks (FWHM in m/z).");
    defaults_.setValue("ppm_tolerance", 10.0, "Gaussian width, depending on the m/z position.\nThe higher the value, the wider the peak and therefore the wider the gaussian.");
    defaults_.setValue("use_ppm_tolerance", "false", "If true, instead of the gaussian_width value, the ppm_tolerance is used. The gaussian is calculated in each step anew, so this is much slower.");
    defaults_.setValidStrings("use_ppm_tolerance", {"true","false"});
    defaults_.setValue("write_log_messages", "false", "true: Warn if no signal was found by the Gauss filter algorithm.");
    defaultsToParam_();
  }

  void GaussFilter::updateMembers_()
  {
    gauss_algo_.initialize(
      (double)param_.getValue("gaussian_width"), 
      spacing_,
      (double)param_.getValue("ppm_tolerance"), 
      param_.getValue("use_ppm_tolerance").toBool());

    write_log_messages_ = param_.getValue("write_log_messages").toBool();
  }

  void GaussFilter::filter(MSSpectrum & spectrum)
  {
    // make sure the right data type is set
    spectrum.setType(SpectrumSettings::PROFILE);
    bool found_signal = false;
    const Size data_size = spectrum.size();
    std::vector<double> mz_in(data_size), int_in(data_size), mz_out(data_size), int_out(data_size);

    // copy spectrum to container
    for (Size p = 0; p < spectrum.size(); ++p)
    {
      mz_in[p] = spectrum[p].getMZ();
      int_in[p] = static_cast<double>(spectrum[p].getIntensity());
    }

    // apply filter
    auto mz_out_it = mz_out.begin();
    auto int_out_it = int_out.begin();
    found_signal = gauss_algo_.filter(mz_in.begin(), mz_in.end(), int_in.begin(), mz_out_it, int_out_it);

    // If all intensities are zero in the scan and the scan has a reasonable size, throw an exception.
    // This is the case if the Gaussian filter is smaller than the spacing of raw data
    if (!found_signal && spectrum.size() >= 3)
    {
      if (write_log_messages_)
      {
        String error_message = "Found no signal. The Gaussian width is probably smaller than the spacing in your profile data. Try to use a bigger width.";
        if (spectrum.getRT() > 0.0)
        {
          error_message += String(" The error occurred in the spectrum with retention time ") + spectrum.getRT() + ".";
        }
        OPENMS_LOG_WARN << error_message << std::endl;
      }
    }
    else
    {
      // copy the new data into the spectrum
      auto mz_it = mz_out.begin();
      auto int_it = int_out.begin();
      for (Size p = 0; mz_it != mz_out.end(); mz_it++, int_it++, p++)
      {
        spectrum[p].setIntensity(*int_it);
        spectrum[p].setMZ(*mz_it);
      }
    }
  }

  void GaussFilter::filter(MSChromatogram & chromatogram)
  {
    if (param_.getValue("use_ppm_tolerance").toBool())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "GaussFilter: Cannot use ppm tolerance on chromatograms");
    }

    bool found_signal = false;
    const Size data_size = chromatogram.size();
    std::vector<double> rt_in(data_size), int_in(data_size), rt_out(data_size), int_out(data_size);

    // copy spectrum to container
    for (Size p = 0; p < chromatogram.size(); ++p)
    {
      rt_in[p] = chromatogram[p].getRT();
      int_in[p] = chromatogram[p].getIntensity();
    }

    // apply filter
    auto mz_out_it = rt_out.begin();
    auto int_out_it = int_out.begin();
    found_signal = gauss_algo_.filter(rt_in.begin(), rt_in.end(), int_in.begin(), mz_out_it, int_out_it);

    // If all intensities are zero in the scan and the scan has a reasonable size, throw an exception.
    // This is the case if the Gaussian filter is smaller than the spacing of raw data
    if (!found_signal && chromatogram.size() >= 3)
    {
      if (write_log_messages_)
      {
        String error_message = "Found no signal. The Gaussian width is probably smaller than the spacing in your chromatogram data. Try to use a bigger width.";
        if (chromatogram.getMZ() > 0.0)
        {
          error_message += String(" The error occurred in the chromatogram with m/z time ") + chromatogram.getMZ() + ".";
        }
        OPENMS_LOG_ERROR << error_message << std::endl;
      }
    }
    else
    {
      // copy the new data into the spectrum
      auto mz_it = rt_out.begin();
      auto int_it = int_out.begin();
      for (Size p = 0; mz_it != rt_out.end(); mz_it++, int_it++, p++)
      {
        chromatogram[p].setIntensity(*int_it);
        chromatogram[p].setMZ(*mz_it);
      }
    }
  }

  void GaussFilter::filterExperiment(PeakMap & map)
  {
    Size progress = 0;
    startProgress(0, map.size() + map.getChromatograms().size(), "smoothing data");
    for (Size i = 0; i < map.size(); ++i)
    {
      filter(map[i]);
      setProgress(++progress);
    }

    for (Size i = 0; i < map.getChromatograms().size(); ++i)
    {
      filter(map.getChromatogram(i));
      setProgress(++progress);
    }
    endProgress();
  }

}
