//--------------------------------------------------------------------------
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <iomanip>

namespace OpenMS
{
  class PeakGroup;

  /**
       @brief A class representing a deconvolved spectrum.
       DeconvolvedSpectrum consists of PeakGroups representing masses.
       For MSn n>1, a PeakGroup representing the precursor mass is also added in this class. Properly assigning a precursor mass
       from the original precursor peak and its deconvolution result is very important in top down proteomics. This assignment is
       performed here for conventional datasets. But for FLASHIda acquired datasets, the assignment is already done by FLASHIda.
       So this class simply use the results from FLASHIda log file for assignment. The parsing of FLASHIda log file is done
       in FLASHDeconv tool class.
       It also contains functions to write in different formats and a function to export this class into MSSpectrum.
  */
  class OPENMS_DLLAPI DeconvolvedSpectrum :
      private std::vector<PeakGroup>
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    using std::vector<PeakGroup>::push_back;
    using std::vector<PeakGroup>::operator[];
    using std::vector<PeakGroup>::size;
    using std::vector<PeakGroup>::begin;
    using std::vector<PeakGroup>::end;
    using std::vector<PeakGroup>::swap;
    using std::vector<PeakGroup>::empty;
    using std::vector<PeakGroup>::reserve;
    using std::vector<PeakGroup>::clear;

    /// default constructor
    DeconvolvedSpectrum() = default;

    /**
       @brief Constructor for DeconvolvedSpectrum. Takes the spectrum and scan number calculated from outside
       @param spectrum spectrum for which the deconvolution will be performed
       @param scan_number scan number of the spectrum: this argument is put here for real time case where scan number should be input separately.
  */
    explicit DeconvolvedSpectrum(const MSSpectrum& spectrum, const int scan_number);

    /// default deconstructor
    ~DeconvolvedSpectrum() = default;

    /// copy constructor
    DeconvolvedSpectrum(const DeconvolvedSpectrum& ) = default;

    /// move constructor
    DeconvolvedSpectrum(DeconvolvedSpectrum&& other) = default;

    /// assignment operator
    DeconvolvedSpectrum& operator=(const DeconvolvedSpectrum& deconvolved_spectrum) = default;

    /**
        @brief write the header in the tsv output file (spectrum level)
        @param fs file stream to the output file
        @param ms_level ms level of the spectrum
        @param detail if set true, detailed information of the mass (e.g., peak list for the mass) is written
   */
    static void writeDeconvolvedMassesHeader(std::fstream& fs,
                                              const int ms_level,
                                              const bool detail);

    /**
      @brief write the deconvolved masses in the output file (spectrum level)
      @param fs file stream to the output file
      @param file_name the output file name that the deconvolved masses will be written.
      @param avg averagine information to calculate monoisotopic and average mass difference within this function. In PeakGroup (peaks of DeconvolvedSpectrum) only monoisotopic mass is recorded. To write both monoisotopic and average masses, their mass difference should be calculated using this averagine information.
      @param write_detail if this is set, more detailed information on each mass will be written in the output file.
      Default MS1 headers are:
        FileName, ScanNum, RetentionTime, MassCountInSpec, AverageMass, MonoisotopicMass,
        SumIntensity, MinCharge, MaxCharge,
        PeakCount, IsotopeCosine, ChargeScore, MassSNR, ChargeSNR, RepresentativeCharge, RepresentativeMzStart, RepresentativeMzEnd, QScore, PerChargeIntensity, PerIsotopeIntensity

      Default MS2 headers include MS1 headers plus:
        PrecursorScanNum, PrecursorMz, PrecursorIntensity, PrecursorCharge, PrecursorSNR, PrecursorMonoisotopicMass, PrecursorQScore

      Detailed MS1 and MS2 headers include all corresponding headers above plus:
        PeakMZs, PeakIntensities, PeakCharges, PeakMasses, PeakIsotopeIndices, PeakPPMErrors
    */
    void writeDeconvolvedMasses(std::fstream& fs,
                                 const String& file_name,
                                 const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                 const bool write_detail);

    /**
      @brief write the deconvolved masses TopFD output (*.msalign)
      @param fs file stream to the output file
      @param avg averagine information to calculate monoisotopic and average mass difference
      @param snr_threshold SNR threshold to filter out low SNR precursors. Even if a PeakGroup has a high deconvolution quality, it should be still discarded for identification when its precursor SNR (SNR within the isolation window) is too low.
      @param decoy_harmonic_factor this factor will be multiplied to precursor mass and charge. To generate decoy spectra
      @param decoy_precursor_offset this value will be added to precursor mass. To generate decoy spectra
    */
    void writeTopFD(std::fstream& fs,
                    const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                    const double snr_threshold = 1.0,
                    const double decoy_harmonic_factor = 1.0,
                    const double decoy_precursor_offset = .0);

    /// Convert DeconvolutedSpectrum to MSSpectrum (e.g., used to store in mzML format).
    /// @param to_charge the charge of each peak in mzml output.
    MSSpectrum toSpectrum(const int to_charge);

    /**
    @brief register the precursor peak as well as the precursor peak group (or mass) if possible for MSn (n>1) spectrum.
    Given a precursor peak (found in the original MS n-1 Spectrum) the masses containing the precursor peak are searched.
    If multiple masses are detected, the one with the best QScore is selected. For the selected mass, its corresponding peak group (along with precursor peak) is registered.
    If no such mass exists, only the precursor peak is registered.
    @param survey_scans the candidate precursor spectra - the user may allow search of previous N survey scans.
    @param is_positive if MS mode is positive
    @param isolation_window_size_ default isolation window size for precursor.
    @param precursor_map_for_real_time_acquisition this contains the deconvolved mass information from FLASHIda runs.
    */
    bool registerPrecursor(const std::vector<DeconvolvedSpectrum>& survey_scans,
                           const bool is_positive, double isolation_window_size_,
                           const std::map<int, std::vector<std::vector<double>>>& precursor_map_for_real_time_acquisition);

    /// original spectrum getter
    const MSSpectrum& getOriginalSpectrum() const;

    /// get precursor peak group for MSn (n>1) spectrum. It returns an empty peak group if no peak group is registered (by registerPrecursor)
    PeakGroup getPrecursorPeakGroup() const;

    /// precursor charge getter (set in registerPrecursor)
    int getPrecursorCharge() const;

    const Precursor getPrecursor() const;

    /// get possible max mass of the deconvolved masses - for MS1, max mass specified by user
    /// for MSn, min value between max mass specified by the user and precursor mass
    /// @param max_mass the max mass specified by the user
    double getCurrentMaxMass(const double max_mass) const;

    /// get possible min mass of the deconvolved masses - for MS1, min mass specified by user
    /// for MSn, 50.0
    /// @param min_mass the min mass specified by the user
    double getCurrentMinMass(const double min_mass) const;

    /// get possible max charge of the deconvolved masses - for MS1, max charge specified by user
    /// for MSn, min value between max charge specified by the user and precursor charge
    /// @param max_abs_charge the max absolute value of the charge specified by the user
    int getCurrentMaxAbsCharge(const int max_abs_charge) const;

    /// get scan number of the original spectrum
    int getScanNumber() const;

    /// get precursor scan number - only if it is registered. Otherwise return 0
    int getPrecursorScanNumber() const;

    /// get activation method
    //std::string getActivation_method();

  private:
    /// the original raw spectrum (not deconvolved)
    MSSpectrum spec_;
    /// precursor peakGroup (or mass)
    PeakGroup precursor_peak_group_;
    /// precursor raw peak (not deconvolved one)
    Precursor precursor_peak_;
    /// activation method for file output
    std::string activation_method_;
    /// scan number and precursor scan number
    int scan_number_ = 0, precursor_scan_number_ = 0;
    /// number of minimum peak count in topFD msalign file
    int topFD_min_peak_count_ = 3;
    /// number of maximum peak count in topFD msalign file
    int topFD_max_peak_count_ = 500;
  };
}
