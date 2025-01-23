// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <limits>
#include <cmath>

namespace OpenMS
{
  /**
      @brief Linear Resampling of raw data.

      This class can be used to generate uniform data from non-uniform raw data (e.g. ESI-TOF or MALDI-TOF experiments).
      Therefore the intensity at every position x in the input raw data is spread to the two
      adjacent resampling points.
      This method preserves the area of the input signal and also the centroid position of a peak.
      Therefore it is recommended for quantitation as well as for ProteinIdentification experiments.

      @note Use this method only for high resolution data (< 0.1 Th between two adjacent raw data points).
           The resampling rate should be >= the precision.

      @htmlinclude OpenMS_LinearResampler.parameters
  */
  class OPENMS_DLLAPI LinearResampler :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:

    /// Constructor
    LinearResampler() :
      DefaultParamHandler("LinearResampler")
    {
      defaults_.setValue("spacing", 0.05, "Spacing of the resampled output peaks.");
      defaultsToParam_();
    }

    /// Destructor.
    ~LinearResampler() override
    {
    }

    /**
        @brief Applies the resampling algorithm to an MSSpectrum, without alignment between spectra.
    */
    void raster(MSSpectrum& spectrum) const
    {
      //return if nothing to do
      if (spectrum.empty()) return;

      typename MSSpectrum::iterator first = spectrum.begin();
      typename MSSpectrum::iterator last = spectrum.end();

      double end_pos = (last - 1)->getMZ();
      double start_pos = first->getMZ();
      int number_raw_points = static_cast<int>(spectrum.size());
      int number_resampled_points = static_cast<int>(ceil((end_pos - start_pos) / spacing_ + 1));

      std::vector<Peak1D> resampled_peak_container;
      resampled_peak_container.resize(number_resampled_points);

      // generate the resampled peaks at positions origin+i*spacing_
      std::vector<Peak1D>::iterator it = resampled_peak_container.begin();
      for (int i = 0; i < number_resampled_points; ++i)
      {
        it->setMZ(start_pos + i * spacing_);
        ++it;
      }

      // spread the intensity h of the data point at position x to the left and right
      // adjacent resampled peaks
      double distance_left = 0.;
      double distance_right = 0.;
      int left_index = 0;
      int right_index = 0;

      it = resampled_peak_container.begin();
      for (int i = 0; i < number_raw_points; ++i)
      {
        int help = static_cast<int>(floor(((first + i)->getMZ() - start_pos) / spacing_));
        left_index = (help < 0) ? 0 : help;
        help = distance(first, last) - 1;
        right_index = (left_index >= help) ? help : left_index + 1;

        // compute the distance between x and the left adjacent resampled peak
        distance_left = fabs((first + i)->getMZ() - (it + left_index)->getMZ()) / spacing_;
        //std::cout << "Distance left " << distance_left << std::endl;
        // compute the distance between x and the right adjacent resampled peak
        distance_right = fabs((first + i)->getMZ() - (it + right_index)->getMZ());
        //std::cout << "Distance right " << distance_right << std::endl;


        // add the distance_right*h to the left resampled peak and distance_left*h to the right resampled peak
        double intensity = static_cast<double>((it + left_index)->getIntensity());
        intensity += static_cast<double>((first + i)->getIntensity()) * distance_right / spacing_;
        (it + left_index)->setIntensity(intensity);
        intensity = static_cast<double>((it + right_index)->getIntensity());
        intensity += static_cast<double>((first + i)->getIntensity()) * distance_left;
        (it + right_index)->setIntensity(intensity);
      }

      spectrum.swap(resampled_peak_container);
    }

    /**
        @brief Resamples the data in an MSExperiment, without alignment between spectra.
    */
    void rasterExperiment(PeakMap& exp)
    {
      startProgress(0, exp.size(), "resampling of data");
      for (Size i = 0; i < exp.size(); ++i)
      {
        raster(exp[i]);
        setProgress(i);
      }
      endProgress();
    }

protected:

    /// Spacing of the resampled data
    double spacing_;

    void updateMembers_() override
    {
      spacing_ =  param_.getValue("spacing");
    }

  };


} // namespace OpenMS

