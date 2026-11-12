// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Mobilogram.h>
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
    defaults_.setMinFloat("gaussian_width", 0.0);
    defaults_.setValue("ppm_tolerance", 10.0, "Gaussian width, depending on the m/z position.\nThe higher the value, the wider the peak and therefore the wider the gaussian.");
    defaults_.setMinFloat("ppm_tolerance", 0.0);
    defaults_.setValue("use_ppm_tolerance", "false",
                       "If true, instead of the gaussian_width value, the ppm_tolerance is used. The gaussian is calculated in each step anew, so "
                       "this is much slower.");
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

  namespace
  {
    /**
      @brief Smooths the 1D peak container @p container (shared by the filter() overloads)

      @param[in] report_no_signal report_no_signal(): called if no signal was found although the container has at least three peaks;
                 the container is left untouched then
    */
    template <typename Container, typename ReportNoSignal>
    void filterContainer_(GaussFilterAlgorithm& gauss_algo, Container& container, ReportNoSignal&& report_no_signal)
    {
      bool found_signal = false;
      const Size data_size = container.size();
      std::vector<double> pos_in(data_size), int_in(data_size), pos_out(data_size), int_out(data_size);

      // copy container to vectors
      for (Size p = 0; p < container.size(); ++p)
      {
        pos_in[p] = container[p].getPos();
        int_in[p] = static_cast<double>(container[p].getIntensity());
      }

      // apply filter
      auto pos_out_it = pos_out.begin();
      auto int_out_it = int_out.begin();
      found_signal = gauss_algo.filter(pos_in.begin(), pos_in.end(), int_in.begin(), pos_out_it, int_out_it);

      // If all intensities are zero in the scan and the scan has a reasonable size, throw an exception.
      // This is the case if the Gaussian filter is smaller than the spacing of raw data
      if (!found_signal && container.size() >= 3)
      {
        report_no_signal();
      }
      else
      {
        // copy the new data into the container
        for (Size p = 0; p < container.size(); ++p)
        {
          container[p].setIntensity(int_out[p]);
          container[p].setPos(pos_out[p]);
        }
      }
    }
  } // namespace

  void GaussFilter::filter(MSSpectrum & spectrum)
  {
    // make sure the right data type is set
    spectrum.setType(SpectrumSettings::SpectrumType::PROFILE);
    filterContainer_(gauss_algo_, spectrum, [&]()
    {
      if (write_log_messages_)
      {
        std::string error_message = "Found no signal. The Gaussian width is probably smaller than the spacing in your profile data. Try to use a bigger width.";
        if (spectrum.getRT() > 0.0)
        {
          error_message +=std::string(" The error occurred in the spectrum with retention time ") + spectrum.getRT() + ".";
        }
        OPENMS_LOG_WARN << error_message << std::endl;
      }
    });
  }

  void GaussFilter::filter(MSChromatogram & chromatogram)
  {
    if (param_.getValue("use_ppm_tolerance").toBool())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "GaussFilter: Cannot use ppm tolerance on chromatograms");
    }

    filterContainer_(gauss_algo_, chromatogram, [&]()
    {
      if (write_log_messages_)
      {
        std::string error_message = "Found no signal. The Gaussian width is probably smaller than the spacing in your chromatogram data. Try to use a bigger width.";
        if (chromatogram.getMZ() > 0.0)
        {
          error_message +=std::string(" The error occurred in the chromatogram with m/z ratio ") + chromatogram.getMZ() + ".";
        }
        OPENMS_LOG_ERROR << error_message << std::endl;
      }
    });
  }

  void GaussFilter::filter(Mobilogram & mobilogram)
  {
    if (param_.getValue("use_ppm_tolerance").toBool())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "GaussFilter: Cannot use ppm tolerance on mobilogram");
    }

    filterContainer_(gauss_algo_, mobilogram, [&]()
    {
      if (write_log_messages_)
      {
        std::string error_message = "Found no signal. The Gaussian width is probably smaller than the spacing in your mobilogram data. Try to use a bigger width.";
        if (mobilogram.getRT() > 0.0)
        {
          error_message +=std::string(" The error occurred in the mobilogram with RT ") + mobilogram.getRT() + ".";
        }
        OPENMS_LOG_ERROR << error_message << std::endl;
      }
    });
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
