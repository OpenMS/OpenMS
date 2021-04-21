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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>

namespace OpenMS {
    /**
     * @brief FLASHIda class for real time deconvolution
     *
     * @see FLASHIdaBridgeFunctions
     * @reference: FeatureFinderAlgorithmPickedHelperStructs
     * @reference: https://stackoverflow.com/questions/31417688/passing-a-vector-array-from-unmanaged-c-to-c-sharp
     */
    class OPENMS_DLLAPI FLASHIda {
    public:
        typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
        typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

        /// constructor that takes string input argument
        explicit FLASHIda(char *arg);

        /// destructor
        ~FLASHIda() = default;

        /// copy constructor
        FLASHIda(const FLASHIda &) = default;

        /// move constructor
        FLASHIda(FLASHIda &&other) = default;

        /// assignment operator
        FLASHIda &operator=(const FLASHIda &fd) = default;

        /**
               @brief get peak groups from input spectrum, specified by mzs and intensities (due to C# interface it is necessary)
               @param mzs mz values of the input spectrum
               @param intensities intensities of the input spectrum
               @param length length of mzs and ints
               @param rt Retention time in seconds
               @param ms_level ms level
               @param name spectrum name
          */
        int getPeakGroups(const double *mzs,
                          const double *intensities,
                          const int length,
                          const double rt,
                          const int ms_level,
                          const char *name);

        /**
               @brief get isolation windows
               @param window_start window start mzs
               @param window_end windo end mzs
               @param qscores QScores of windows
               @param charges charges of windows
               @param avg_masses average masses of windows
          */
        void getIsolationWindows(double *wstart,
                                 double *wend,
                                 double *qscores,
                                 int *charges,
                                 int *min_charges,
                                 int *max_charges,
                                 double *mono_masses,
                                 double *chare_cos,
                                 double *charge_snrs,
                                 double *iso_cos,
                                 double *snrs, double *charge_scores,
                                 double *ppm_errors,
                                 double *precursor_intensities,
                                 double *peakgroup_intensities);

    private:

        /// PeakGroup comparator for soring by QScore
        struct {
            bool operator()(const PeakGroup &a, const PeakGroup &b) const {
                return a.getQScore() > b.getQScore();
            }
        } QscoreComparator_;

        /// Selected integer masses - necessary for mass exclusion
        std::unordered_map<int, double> mz_rt_map_;
        /// Selected integer masses - necessary for mass exclusion
        std::unordered_map<int, double> mass_rt_map_;
        ///
        std::unordered_map<int, double> all_mass_rt_map_;
        std::unordered_map<int, double> mass_qscore_map_;

        /// discard peak groups using mass exclusion
        void filterPeakGroupsUsingMassExclusion_(const MSSpectrum &spec, const int ms_level, const double rt);

        /// generate MSSpectrum class using mzs and intensities
        static MSSpectrum
        makeMSSpectrum_(const double *mzs, const double *ints, const int length, const double rt, const int ms_level,
                        const char *name);

        /// deconvoluted spectrum that contains the peak group
        DeconvolutedSpectrum deconvoluted_spectrum_;
        /// peakGroup charges to be triggered
        std::vector<int> trigger_charges;

        /// FLASHDeconvAlgorithm class for deconvolution
        FLASHDeconvAlgorithm fd_;

        double error_threshold_ = .9;

        /// q score threshold - determined from C# side
        double qscore_threshold_;
        /// retention time window - determined from C# side
        double rt_window_;
        /// how many masses will be selected per ms level? - determined from C# side
        IntList mass_count_;
        /// minimum isolation window width divided by two
        const double min_isolation_window_half_ = .6;

        std::map<int, std::vector<double>> target_nominal_masses_;
        std::set<double> target_masses_;
        double charge_snr_threshold_ = 1.0;
        //const double snr_threshold = 0.0;
        //const double isotope_cosine_threshold = 0;
    };
}
