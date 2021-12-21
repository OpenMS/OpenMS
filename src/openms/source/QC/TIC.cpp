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
// $Maintainer: Chris Bielow $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------


#include <OpenMS/QC/TIC.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>
#include <OpenMS/FORMAT/MzTab.h>

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

      for (const auto& i: result.intensities)
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
        if (result.intensities[i] > result.intensities[i-1] * 10) // detect 10x jumps between two subsequent scans
        {
          ++result.jump;
        }
        if (result.intensities[i] < result.intensities[i-1] / 10) // detect 10x falls between two subsequent scans
        {
          ++result.fall;
        }
      }
    }
    return result;
  }

  bool TIC::Result::operator==(const Result& rhs) const
  {
    return intensities == rhs.intensities
          && retention_times == rhs.retention_times
          && area == rhs.area
          && fall == rhs.fall
          && jump == rhs.jump;
  }

  /// Returns the name of the metric
  const String& TIC::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status TIC::requires() const
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
      MzTabParameter tic{};
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
}
