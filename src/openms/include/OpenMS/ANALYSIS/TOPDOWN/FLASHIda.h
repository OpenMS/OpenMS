// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>

namespace OpenMS
{
  /**
   * @brief FLASHIda class for real time deconvolution
   * This class contains functions to perform deconvolution (by FLASHDeconvAlgorithm) for the spectrum received from Thermo iAPI.
   * Also precursor selection is done in this class.
   * The functions in this class are invoked in C# Thermo iAPI side through the functions in FLASHIdaBridgeFunctions class
   * @see FLASHIdaBridgeFunctions
   * @reference: https://stackoverflow.com/questions/31417688/passing-a-vector-array-from-unmanaged-c-to-c-sharp
   */
  class OPENMS_DLLAPI FLASHIda
  {
  public:
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// constructor that takes string input argument
    explicit FLASHIda(char *arg);

    /// destructor
    ~FLASHIda() = default;

    /// copy constructor
    FLASHIda(const FLASHIda& ) = default;

    /// move constructor
    FLASHIda(FLASHIda&& other) = default;

    /// assignment operator
    FLASHIda& operator=(const FLASHIda& fd) = default;

    /**
           @brief get peak groups (deconvolved masses) from input spectrum, specified by mzs and intensities (due to C# interface it is necessary)
           @param mzs mz values of the input spectrum
           @param intensities intensities of the input spectrum
           @param length length of mzs and ints
           @param rt Retention time in seconds
           @param ms_level ms level
           @param name spectrum name
           @return number of peak groups
      */
    int getPeakGroups(const double *mzs,
                      const double *intensities,
                      const int length,
                      const double rt,
                      const int ms_level,
                      const char *name);

    /**
           @brief get isolation windows using FLASHDeconv algorithm. Many parameters are in primitive types so they can be passed to C# FLASHIda side.
           All parameters are for isolation windows.
           @param window_start window start mzs
           @param window_end window end mzs
           @param qscores QScores of windows
           @param charges charges of windows
           @param min_charges minimum charges
           @param max_charges maximum charges
           @param mono_masses monoisotopic masses
           @param charge_cos charge cosine scores
           @param charge_snrs charge SNRs or precursor SNRs
           @param iso_cos mass cosine scores
           @param snrs mass SNRs
           @param charge_scores charge distribution scores
           @param ppm_errors average PPM errors
           @param precursor_intensities precursor peak intensities
           @param peakgroup_intensities precursor mass intensities
      */
    void getIsolationWindows(double *window_start,
                             double *window_end,
                             double *qscores,
                             int *charges,
                             int *min_charges,
                             int *max_charges,
                             double *mono_masses,
                             double *charge_cos,
                             double *charge_snrs,
                             double *iso_cos,
                             double *snrs, double *charge_scores,
                             double *ppm_errors,
                             double *precursor_intensities,
                             double *peakgroup_intensities);

    /**
           @brief parse FLASHIda log file
           @param in_log_file input log file
           @return parsed information : scan number - percursor information
    **/
    static std::map<int, std::vector<std::vector<double>>> parseFLASHIdaLog(const String& in_log_file);

  private:
    /// PeakGroup comparator for soring by QScore
    /*struct
    {
      bool operator()(const PeakGroup& a, const PeakGroup& b) const
      {
        return a.getQScore() > b.getQScore();
      }
    } QscoreComparator_;
*/
    /// Maps that are necessary for mass exclusion
    std::unordered_map<int, double> tqscore_exceeding_mz_rt_map_; /// integer mz value vs. retention time with tqscore exceeding total qscore threshold
    std::unordered_map<int, double> tqscore_exceeding_mass_rt_map_; /// integer mass value vs. retention time with tqscore exceeding total qscore threshold
    std::unordered_map<int, double> all_mass_rt_map_; /// mz value vs. retention time for all acquired precursors
    std::unordered_map<int, double> mass_qscore_map_; /// mass value vs. retention time for all acquired precursors

    /**
         @brief discard peak groups using mass exclusion
         @param ms_level MS level
         @param rt Retention time
    */
    void filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt);

    /**
         @brief generate MSSpectrum class using mzs and intensities. mzs and intensities and other information are
         provided by Thermo iAPI
         @param mzs m/z values
         @param ints intensities
         @param length number of peaks
         @param rt Retention time
         @param ms_level MS level
         @param name spectrum name
    */
    static MSSpectrum makeMSSpectrum_(const double *mzs,
                                      const double *ints,
                                      const int length,
                                      const double rt,
                                      const int ms_level,
                                      const char *name);

    /// deconvolved spectrum that contains the peak group
    DeconvolvedSpectrum deconvolved_spectrum_;
    /// peakGroup charges to be triggered
    std::vector<int> trigger_charges;
    /// peakGroup isolation window ranges
    std::vector<double> trigger_left_isolation_mzs_;
    std::vector<double> trigger_right_isolation_mzs_;

    /// FLASHDeconvAlgorithm class for deconvolution
    FLASHDeconvAlgorithm fd_;

    /// total QScore threshold
    double tqscore_threshold = .9;

    /// q score threshold - determined from C# side
    double qscore_threshold_;
    /// retention time window - determined from C# side
    double rt_window_;
    /// how many masses will be selected per ms level? - determined from C# side
    IntList mass_count_;

    int targeting_mode_ = 0; /// 0 no targeting 1 inclusive 2 exclusive
    /// maps for global inclusion/exclusion targeting
    std::map<double, std::vector<double>> target_mass_rt_map_;
    std::vector<double> target_masses_; /// current target masses

    /// precursor SNR threshold
    double snr_threshold_ = 1.0;

    /// mass tolerance
    DoubleList tol_;
  };
}
