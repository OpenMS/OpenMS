// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------


#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/QC/TIC.h>

using namespace std;

namespace OpenMS
{

  TIC::Result TIC::compute(const MSExperiment& exp, float bin_size, UInt ms_level)
  {
    TIC::Result result;
    MSChromatogram tic = exp.calculateTIC(bin_size, ms_level);
    if (!tic.empty())
    {
      for (const auto& p : tic)
      {
        result.intensities.push_back(p.getIntensity());
        result.retention_times.push_back(p.getRT());
      }

      UInt max_int = *max_element(result.intensities.begin(), result.intensities.end());

      for (const auto& i : result.intensities)
      {
        if (max_int != 0)
        {
          result.relative_intensities.push_back((double)i / max_int * 100);
        }
        else
        {
          result.relative_intensities.push_back(0.0);
        }
      }

      result.area = result.intensities[0];

      for (size_t i = 1; i < result.intensities.size(); ++i)
      {
        result.area += result.intensities[i];
        if (result.intensities[i] > result.intensities[i - 1] * 10) // detect 10x jumps between two subsequent scans
        {
          ++result.jump;
        }
        if (result.intensities[i] < result.intensities[i - 1] / 10) // detect 10x falls between two subsequent scans
        {
          ++result.fall;
        }
      }
    }
    return result;
  }

  bool TIC::Result::operator==(const Result& rhs) const
  {
    return intensities == rhs.intensities && retention_times == rhs.retention_times && area == rhs.area && fall == rhs.fall && jump == rhs.jump;
  }

  /// Returns the name of the metric
  const String& TIC::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status TIC::requirements() const
  {
    return QCBase::Status(QCBase::Requires::RAWMZML);
  }

  void TIC::addMetaDataMetricsToMzTab(OpenMS::MzTabMetaData& meta, vector<TIC::Result>& tics)
  {
    // Adding TIC information to meta data
    for (Size i = 0; i < tics.size(); ++i)
    {
      if (tics[i].intensities.empty())
      {
        continue; // no MS1 spectra
      }
      MzTabParameter tic {};
      tic.setCVLabel("total ion current");
      tic.setAccession("MS:1000285");
      tic.setName("TIC_" + String(i + 1));
      String value("[");
      value += String(tics[i].retention_times[0], false) + ", " + String((UInt64)tics[i].intensities[0]);
      for (Size j = 1; j < tics[i].intensities.size(); ++j)
      {
        value += ", " + String(tics[i].retention_times[j], false) + ", " + String((UInt64)tics[i].intensities[j]);
      }
      value += "]";
      tic.setValue(value);
      meta.custom[meta.custom.size()] = tic;
    }
  }
} // namespace OpenMS
