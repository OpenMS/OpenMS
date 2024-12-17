// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>

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
    defaults_.setValidStrings("write_log_messages", {"true","false"});
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
      for (Size p = 0; p < spectrum.size(); ++p)
      {
        spectrum[p].setIntensity(int_out[p]);
        spectrum[p].setMZ(mz_out[p]);
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
          error_message += String(" The error occurred in the chromatogram with m/z ratio ") + chromatogram.getMZ() + ".";
        }
        OPENMS_LOG_ERROR << error_message << std::endl;
      }
    }
    else
    {
      // copy the new data into the chromatogram
      for (Size p = 0; p < chromatogram.size(); ++p)
      {
        chromatogram[p].setIntensity(int_out[p]);
        chromatogram[p].setRT(rt_out[p]);
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
