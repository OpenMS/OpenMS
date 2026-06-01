// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataAggregatingConsumer.h>

#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>

namespace OpenMS
{

  MSDataAggregatingConsumer::~MSDataAggregatingConsumer()
  {
    // flush remaining spectra
    if (!s_list.empty())
    {
      MSSpectrum tmps = SpectrumAddition::addUpSpectra(s_list, -1, true);
      copySpectrumMeta(s_list[0], tmps, false);
      next_consumer_->consumeSpectrum(tmps);
    }
  }

  void MSDataAggregatingConsumer::consumeSpectrum(SpectrumType & s)
  {
    // aggregate by RT
    double RT = s.getRT();

    if (rt_initialized_ && std::fabs(RT - previous_rt_) < 1e-5)
    {
      // need to aggregate spectrum
      s_list.push_back(s);
    }
    else
    {
      // consume the previous list
      if (rt_initialized_ && !s_list.empty())
      {
        MSSpectrum tmps = SpectrumAddition::addUpSpectra(s_list, -1, true);
        copySpectrumMeta(s_list[0], tmps, false);
        next_consumer_->consumeSpectrum(tmps);
      }

      // start new spectrum list
      int expected_size = s_list.size();
      s_list.clear();
      s_list.reserve(expected_size);
      s_list.push_back(s);
    }

    previous_rt_ = RT;
    rt_initialized_ = true;
  }

  void MSDataAggregatingConsumer::consumeChromatogram(ChromatogramType & c)
  {
    // NOP
    next_consumer_->consumeChromatogram(c);
  }

} // namespace OpenMS

