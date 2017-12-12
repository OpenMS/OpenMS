// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ISOBARICQUANTIFIERSTATISTICS_H
#define OPENMS_ANALYSIS_QUANTITATION_ISOBARICQUANTIFIERSTATISTICS_H

#include <OpenMS/CONCEPT/Types.h>

#include <map>

namespace OpenMS
{
  class String;

  /**
    @brief Statistics for quantitation performance and comparison of NNLS vs. naive method (aka matrix inversion)
   */
  struct OPENMS_DLLAPI IsobaricQuantifierStatistics
  {
    /**
     @brief Create stats object.
     */
    IsobaricQuantifierStatistics();

    /**
     @brief Reset statistics object.
     */
    void reset();

    Size channel_count; ///< 4plex, 6plex, or 8 plex?!
    Size iso_number_ms2_negative; ///< number of MS2 spectra where one or more channels had negative solution
    Size iso_number_reporter_negative; ///< number of channels where naive solution was negative
    Size iso_number_reporter_different; ///< number of channels >0 where naive solution was different; happens when naive solution is negative in other channels
    double iso_solution_different_intensity; ///< absolute intensity difference between both solutions (for channels > 0)
    double iso_total_intensity_negative; ///< only for spectra where naive solution is negative
    Size number_ms2_total; ///< total number of MS2 spectra
    Size number_ms2_empty; ///< number of empty MS2 (no reporters at all)
    std::map<String, Size> empty_channels; ///< Channel_ID -> Missing; indicating the number of empty channels from all MS2 scans, i.e., numbers are between number_ms2_empty and number_ms2_total

    /// Copy c'tor
    IsobaricQuantifierStatistics(const IsobaricQuantifierStatistics& other);

    /// Assignment operator
    IsobaricQuantifierStatistics& operator=(const IsobaricQuantifierStatistics& rhs);
  };
} // namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ISOBARICQUANTIFIERSTATISTICS_H
