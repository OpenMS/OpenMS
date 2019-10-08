// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{

class MSSpectrum;

class OPENMS_DLLAPI Deisotoper
{
  public:

  /** @brief Detect isotopic clusters in a mass spectrum.

    Deisotoping is based on C13 abundance and will try to identify a simple
    model based on the C12-C13 distance and charge state. This is often a good
    approximation for peptide fragment ion spectra but may not work well for
    other spectra. The algorithm will consider each peak (starting from the
    right of a spectrum) and, for each peak, attempt to add isotopic peaks to
    its envelope until either no peak is found, the maximum number of isotopic
    peaks is reached or (only when using @p use_decreasing_model) the intensity
    of the peak is higher than the previous peak.

    Deisotoping is done in-place and if @p annotate_charge is true,
    an additional IntegerDataArray "charge" will be appended. If 
    @p annotate_iso_peak_count is true, an additional IntegerDataArray 
    "iso_peak_count" containing the number of isotopic peaks will be 
    appended.
    Existing DataArrays are kept and shrunken to the peaks which
    remain in the spectrum.

   * @param [spectra] Input spectra (sorted by m/z)
   * @param [fragment_tolerance] The tolerance used to match isotopic peaks
   * @param [fragment_unit_ppm] Whether ppm or m/z is used as tolerance
   * @param [min_charge] The minimum charge considered
   * @param [max_charge] The maximum charge considered
   * @param [keep_only_deisotoped] Only monoisotopic peaks of fragments with isotopic pattern are retained
   * @param [min_isopeaks] The minimum number of isotopic peaks (at least 2) required for an isotopic cluster
   * @param [max_isopeaks] The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster
   * @param [make_single_charged] Convert deisotoped monoisotopic peak to single charge
   * @param [annotate_charge] Annotate the charge to the peaks in the IntegerDataArray: "charge" (0 for unknown charge)
   * @param [annotate_iso_peak_count] Annotate the number of isotopic peaks in a pattern for each monoisotopic peak in the IntegerDataArray: "iso_peak_count"
   * @param [use_decreasing_model] Use a simple averagine model that expects heavier isotopes to have less intensity. If false, no intensity checks are applied.
   * @param [add_up_intensity] Sum up the total intensity of each isotopic pattern into the intensity of the reported monoisotopic peak
   *
   * Note: If @p make_single_charged is selected, the original charge (>=1) gets annotated.
   */
  static void deisotopeAndSingleCharge(MSSpectrum& spectrum,
            double fragment_tolerance,
            bool fragment_unit_ppm,
            int min_charge = 1,
            int max_charge = 3,
            bool keep_only_deisotoped = false,
            unsigned int min_isopeaks = 3,
            unsigned int max_isopeaks = 10,
            bool make_single_charged = true,
            bool annotate_charge = false,
            bool annotate_iso_peak_count = false,
            bool use_decreasing_model = true,
            bool add_up_intensity = false);
};

}
