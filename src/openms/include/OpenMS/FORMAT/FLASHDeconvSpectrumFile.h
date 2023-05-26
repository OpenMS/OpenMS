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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

#include <iomanip>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO
**/

  class OPENMS_DLLAPI FLASHDeconvSpectrumFile
  {
  public:
    /**
            @brief write the header in the tsv output file (spectrum level)
            @param fs file stream to the output file
            @param ms_level ms level of the spectrum
            @param detail if set true, detailed information of the mass (e.g., peak list for the mass) is written
            @param dummy if set true, dummy and qvalue information will be written.
       */
    static void writeDeconvolvedMassesHeader(std::fstream& fs,
                                             uint ms_level,
                                             bool detail,
                                             bool dummy);

    /**
          @brief write the deconvolved masses in the output file (spectrum level)
          @param dspec deconvolved spectrum to write
          @param target_spec target spectrum only used for dummy spectrum output
          @param fs file stream to the output file
          @param file_name the output file name that the deconvolved masses will be written.
          @param avg averagine information to calculate monoisotopic and average mass difference within this function. In PeakGroup (peaks of DeconvolvedSpectrum) only monoisotopic mass is recorded. To write both monoisotopic and average masses, their mass difference should be calculated using this averagine information.
          @param tol mass tolerance
          @param write_detail if this is set, more detailed information on each mass will be written in the output file.
          @param dummy if set true, dummy and qvalue information will be written.
          Default MS1 headers are:
            FileName, ScanNum, TargetDummyType, RetentionTime, MassCountInSpec, AverageMass, MonoisotopicMass,
            SumIntensity, MinCharge, MaxCharge,
            PeakCount, IsotopeCosine, ChargeScore, MassSNR, ChargeSNR, RepresentativeCharge, RepresentativeMzStart, RepresentativeMzEnd, Qscore, PerChargeIntensity, PerIsotopeIntensity

          Default MS2 headers include MS1 headers plus:
            PrecursorScanNum, PrecursorMz, PrecursorIntensity, PrecursorCharge, PrecursorSNR, PrecursorMonoisotopicMass, PrecursorQscore

          Detailed MS1 and MS2 headers include all corresponding headers above plus:
            PeakMZs, PeakIntensities, PeakCharges, PeakMasses, PeakIsotopeIndices, PeakPPMErrors
        */
    static void writeDeconvolvedMasses(DeconvolvedSpectrum& dspec,
                                       DeconvolvedSpectrum& target_spec,
                                       std::fstream& fs,
                                       const String& file_name,
                                       const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                       double tol,
                                       bool write_detail,
                                       bool record_dummy);

    /**
      @brief write the deconvolved masses TopFD output (*.msalign)
      @param dspec deconvolved spectrum to write
      @param fs file stream to the output file
      @param snr_threshold SNR threshold to filter out low SNR precursors. Even if a PeakGroup has a high deconvolution quality, it should be still discarded for identification when its precursor SNR (SNR within the isolation window) is too low.
      @param min_ms_level min ms level of the dataset
      @param randomize_precursor_mass if set, a random number between -100 to 100 is added to precursor mass
      @param randomize_fragment_mass if set, a random number between -100 to 100 is added to fragment mass
    */
    //      @param avg averagine information to calculate monoisotopic and average mass difference
    static void writeTopFD(DeconvolvedSpectrum& dspec, std::fstream& fs,
                           double snr_threshold = 1.0,
                           const uint min_ms_level = 1,
                           bool randomize_precursor_mass = false,
                           bool randomize_fragment_mass = false);

  private:

    /// number of minimum peak count in topFD msalign file
    static const int topFD_min_peak_count_ = 3;
    /// number of maximum peak count in topFD msalign file
    static const int topFD_max_peak_count_ = 500;

  };
}// namespace OpenMS