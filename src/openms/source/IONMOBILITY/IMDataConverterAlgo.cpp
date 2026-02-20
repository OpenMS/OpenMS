// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//
// This file contains the IMDataConverter::splitExperimentByIonMobility method
// which depends on SpectraMerger (Algo tier). It is compiled into OpenMS_Algo
// to avoid a Core→Algo back-dependency.

#include <OpenMS/IONMOBILITY/IMDataConverter.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>

#include <cassert>
#include <climits>

namespace OpenMS
{
  std::tuple<std::vector<MSExperiment>, Math::BinContainer> IMDataConverter::splitExperimentByIonMobility(MSExperiment&& in,
                                                                                                          UInt number_of_bins,
                                                                                                          double bin_extension_abs,
                                                                                                          double mz_binning_width,
                                                                                                          MZ_UNITS mz_binning_width_unit)
  {
    if (number_of_bins == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot split into 0 bins.", String(number_of_bins));
    }
    if (bin_extension_abs < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Overlap must not be negative.", String(bin_extension_abs));
    }
    std::vector<MSExperiment> results(number_of_bins);
    in.updateRanges();
    // find the IM range
    const auto range_IM = RangeMobility(in.spectrumRanges());
    if (range_IM.getSpan() / number_of_bins < bin_extension_abs * 2)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Bin size (") + String(range_IM.getSpan() / number_of_bins) + ") is smaller than the overlap.", String(bin_extension_abs*2));
    }

    // compute the bins
    const auto bins = Math::createBins(range_IM.getMin(), range_IM.getMax(), number_of_bins, bin_extension_abs);

    // results for each IM-frame: all spectra per bin, to get merged
    MSExperiment binned_spectra;

    SpectraMerger merger;
    auto p = merger.getParameters();
    const auto ms_levels = in.getMSLevels();
    p.setValue("block_method:ms_levels", IntList(ms_levels.begin(), ms_levels.end())); // merge all MS levels
    p.setValue("mz_binning_width", mz_binning_width);
    p.setValue("mz_binning_width_unit", String(MZ_UNIT_NAMES[(int)mz_binning_width_unit]));
    p.setValue("block_method:rt_block_size", INT_MAX);
    p.setValue("block_method:rt_max_length", 10e10);


    for (auto& frame : in)
    {
      // For data without ion mobility, simply append the result (only
      // collapse for scans that actually have a float data array).
      if (! frame.containsIMData())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Spectrum does not contain 'wide' IM data.", frame.getNativeID());
      }

      MSExperiment frame_melt = IMDataConverter::reshapeIMFrameToMany(std::move(frame));
      for (size_t i = 0; i < bins.size(); ++i)
      {
        binned_spectra.clear(false);
        // check if spectrum goes into this bin
        for (auto&& spec : frame_melt)
        {
          if (bins[i].contains(spec.getDriftTime()))
          { // spectrum goes into this bin
            binned_spectra.addSpectrum(std::move(spec));
          }
        }
        // collapse spectra in this bin
        if (!binned_spectra.empty())
        {
          merger.setParameters(p);
          merger.mergeSpectraBlockWise(binned_spectra);
          assert(binned_spectra.size() == 1);
          results[i].addSpectrum(std::move(binned_spectra.getSpectra().back()));
          results[i].getSpectra().back().setDriftTime(bins[i].center());
        }
      }
    }
    for (auto& result : results)
    {
      result.ExperimentalSettings::operator=(in);
      result.updateRanges();
    }
    return {std::move(results), std::move(bins)};
  }

} // namespace OpenMS
