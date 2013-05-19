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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

namespace OpenMS
{

  IsobaricQuantifierStatistics::IsobaricQuantifierStatistics() :
    channel_count(0),
    iso_number_ms2_negative(0),
    iso_number_reporter_negative(0),
    iso_number_reporter_different(0),
    iso_solution_different_intensity(0),
    iso_total_intensity_negative(0),
    number_ms2_total(0),
    number_ms2_empty(0),
    empty_channels()
  {
  }

  void IsobaricQuantifierStatistics::reset()
  {
    channel_count = 0;
    iso_number_ms2_negative = 0;
    iso_number_reporter_negative = 0;
    iso_number_reporter_different = 0;
    iso_solution_different_intensity = 0;
    iso_total_intensity_negative = 0;
    number_ms2_total = 0;
    number_ms2_empty = 0;
    empty_channels.clear();
  }

  IsobaricQuantifierStatistics::IsobaricQuantifierStatistics(const IsobaricQuantifierStatistics& other)
  {
    channel_count = other.channel_count;
    iso_number_ms2_negative = other.iso_number_ms2_negative;
    iso_number_reporter_negative = other.iso_number_reporter_negative;
    iso_number_reporter_different = other.iso_number_reporter_different;
    iso_solution_different_intensity = other.iso_solution_different_intensity;
    iso_total_intensity_negative = other.iso_total_intensity_negative;
    number_ms2_total = other.number_ms2_total;
    number_ms2_empty = other.number_ms2_empty;
    empty_channels.clear();
    empty_channels.insert(other.empty_channels.begin(), other.empty_channels.end());
  }

  IsobaricQuantifierStatistics& IsobaricQuantifierStatistics::operator=(const IsobaricQuantifierStatistics& rhs)
  {
    if (this == &rhs)
      return *this;

    channel_count = rhs.channel_count;
    iso_number_ms2_negative = rhs.iso_number_ms2_negative;
    iso_number_reporter_negative = rhs.iso_number_reporter_negative;
    iso_number_reporter_different = rhs.iso_number_reporter_different;
    iso_solution_different_intensity = rhs.iso_solution_different_intensity;
    iso_total_intensity_negative = rhs.iso_total_intensity_negative;
    number_ms2_total = rhs.number_ms2_total;
    number_ms2_empty = rhs.number_ms2_empty;
    empty_channels.clear();
    empty_channels.insert(rhs.empty_channels.begin(), rhs.empty_channels.end());

    return *this;
  }

} // namespace
