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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES//Matrix.h>
#include "boost/dynamic_bitset.hpp"
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>

namespace OpenMS
{
  /**
  @brief FLASHDeocnv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset
  From MSSPectrum, this class outputs DeconvolutedSpectrum.
  Deconvolution takes three steps:
   i) decharging and select candidate masses - speed up via binning
   ii) collecting isotopes from the candidate masses and deisotope - peak groups are defined here
   iii) scoring and filter out low scoring masses (i.e., peak groups)
  @ingroup Topdown
*/

  class OPENMS_DLLAPI FLASHDeconvAlgorithm :
      public DefaultParamHandler
  {
  public:
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    FLASHDeconvAlgorithm();

    /// default destructor
    ~FLASHDeconvAlgorithm() override;

    /// copy constructor
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm &) = default;

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm &&other) = default;

    /// assignment operator
    FLASHDeconvAlgorithm &operator=(const FLASHDeconvAlgorithm &fd);

    /**
      @brief main deconvolution function that generates the deconvoluted spectrum from the original spectrum.
      @param spec the original spectrum
      @param survey_scans the survey scans to assign precursor mass to the deconvoluted spectrum.
      @param scan_number scan number can be retrieved from the spectrum in most cases.
      But this parameter is put for real time deconvolution where scan number may be put separately.
      @param precursor_map_for_FLASHIda deconvoluted precursor information from FLASHIda
      @return the deconvoluted spectrum (as DeconvolutedSpectrum class)
 */
    DeconvolutedSpectrum &getDeconvolutedSpectrum(const MSSpectrum &spec,
                                                  const std::vector<DeconvolutedSpectrum> &survey_scans,
                                                  const int scan_number,
                                                  const std::map<int, std::vector<std::vector<double>>> &precursor_map_for_FLASHIda);

    /// get calculated averagine
    PrecalculatedAveragine getAveragine();


    void setTargetMasses(const std::vector<double> &masses, int ms_level);

    /** @brief precalculate averagine (for predifined mass bins) to speed up averagine generation
        @param use_RNA_averagine if set, averagine for RNA (nucleotides) is calcualted
     */
    void calculateAveragine(const bool use_RNA_averagine);

    /// convert double to nominal mass
    static int getNominalMass(const double mass);

    /** Examine charge intensity distribution of each peak group
        @per_charge_intensity per charge intensity - aggregated through isotope indices
        @charge_range max charge range (current_max_charge_ - minCharge)
        @return calculated charge fit score (0 - 1)
     */
    static double getChargeFitScore(const std::vector<double> &per_charge_intensity, const int charge_range);

    /** @brief Examine intensity distribution over iostope indices. Also determines the most plausible isotope index or, monoisotopic mono_mass
        @param mono_mass monoisotopic mass
        @param per_isotope_intensities per isotope intensity - aggregated through charges
        @param offset output offset between input monoisotopic mono_mass and determined monoisotopic mono_mass
        @param avg precalculated averagine
        @param use_shape_diff use shape diff function for more accurate monoisotopic mass determination
        @return calculated cosine similar score
     */
    static double getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass,
                                                           const std::vector<double> &per_isotope_intensities,
                                                           int &offset,
                                                           const PrecalculatedAveragine &avg,
                                                           bool use_shape_diff = true);


  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconv parameters

    /// range of RT subject to analysis (in seconds)
    double min_rt_, max_rt_;
    /// range of mz subject to analysis
    double min_mz_, max_mz_;
    /// min charge and max charge subject to analysis, set by users
    int min_abs_charge_, max_abs_charge_;
    /// is positive mode
    bool is_positive_;
    /// when a spectrum is deconvoluted, the deconvoluted masses within the this rt window are favorably considered.
    double rt_window_;
    /// mass ranges of deconvolution, set by users
    double min_mass_, max_mass_;
    /// min charge: 1 for MSn n>1; otherwise just min_abs_charge_
    int current_min_charge_;
    /// max charge: controlled by precursor charge for MSn n>1; otherwise just max_abs_charge_
    int current_max_charge_;
    /// max mass is controlled by precursor mass for MSn n>1; otherwise just max_mass
    double current_max_mass_;
    /// max mass is max_mass for MS1 and 50 for MS2
    double current_min_mass_;
    /// peak intensity threshold subject to analysis
    double intensity_threshold_;
    /// minimum number of peaks supporting a mass
    IntList min_support_peak_count_;
    /// tolerance in ppm for each MS level
    DoubleList tolerance_;
    /// bin size for first stage of mass selection - for fast convolution, binning is used
    DoubleList bin_width_;
    /// cosine threshold between observed and theoretical isotope patterns for each MS level
    DoubleList min_isotope_cosine_;
    /// max mass count per spectrum for each MS level
    IntList max_mass_count_;

    /// precalculated averagine distributions for fast averagine generation
    FLASHDeconvHelperStructs::PrecalculatedAveragine avg_;
    /// The data structures for spectra overlapping.
    std::vector<std::vector<Size>> prev_mass_bins_ms1_;
    std::vector<std::map<int, std::vector<Size>>> prev_mass_bins_ms2_;
    std::vector<double> prev_rts_ms1_;
    std::vector<double> prev_rts_ms2_;

    /// mass bins that are targeted for FLASHIda global targeting mode
    boost::dynamic_bitset<> target_mass_bins_;
    std::vector<double> target_masses_;

    /// harmonic charge factors that will be considered for harmonic mass reduction. For example, 2 is for 1/2 charge harmonic component reduction
    const std::vector<int> harmonic_charges_{2, 3, 5, 7};
    /// Stores log mz peaks
    std::vector<LogMzPeak> log_mz_peaks_;
    /// deconvoluted_spectrum_ stores the decovnoluted mass peak groups
    DeconvolutedSpectrum deconvoluted_spectrum_;

    /// mass_bins_ stores the selected bins for this spectrum + overlapped spectrum (previous a few spectra).
    boost::dynamic_bitset<> mass_bins_;
    /// mz_bins_ stores the binned log mz peaks
    boost::dynamic_bitset<> mz_bins_;
    /// mz_bins_for_edge_effect_ stores the binned log mz peaks, considering edge effect
    boost::dynamic_bitset<> mz_bins_for_edge_effect_;

    /// This stores the "universal pattern"
    std::vector<double> filter_;
    /// This stores the patterns for harmonic reduction
    Matrix<double> harmonic_filter_matrix_;

    /// This stores the "universal pattern" in binned dimension
    std::vector<int> bin_offsets_;
    /// This stores the patterns for harmonic reduction in binned dimension
    Matrix<int> harmonic_bin_offset_matrix_;

    /// minimum mass and mz values representing the first bin of massBin and mzBin, respectively: to save memory space
    double mass_bin_min_value_;
    double mz_bin_min_value_;

    /// current ms Level
    int ms_level_;

    /// high and low charges are differently deconvoluted. This vale determines the (includisve) threshold for low charge.
    const int low_charge_ = 10; // inclusive

    /// allowed maximum peak count per spectrum - intensity based.
    const int max_peak_count_ = 30000;//30000

    /** @brief static function that converts bin to value
        @param bin bin number
        @param min_value minimum value (corresponding to bin number = 0)
        @param bin_width bin width
        @return value corresponding to bin
     */
    static double getBinValue_(const Size bin, const double min_value, const double bin_width);

    /** @brief static function that converts value to bin
        @param value value
        @param min_value minimum value (corresponding to bin number = 0)
        @param bin_width bin width
        @return bin corresponding to value
     */
    static Size getBinNumber_(const double value, const double min_value, const double bin_width);

    ///static function that retruns average PPM error of the input peak group
    static float getAvgPPMError_(PeakGroup pg);

    ///generate log mz peaks from the input spectrum
    void updateLogMzPeaks_(const MSSpectrum *spec);

    /** @brief generate mz bins and intensity per mz bin from log mz peaks
        @param bin_number number of mz bins
        @param mz_bin_intensities intensity per mz bin
     */
    void updateMzBins_(const Size &bin_number, std::vector<float> &mz_bin_intensities);

    ///this function takes the previous deconvolution results (from ovelapped spectra) for sensitive deconvolution of the current spectrum
    void unionPrevMassBins_();

    ///Generate peak groups from the input spectrum
    void generatePeakGroupsFromSpectrum_();

    /** @brief Update mass_bins_. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses. This function does not perform deisotoping
        @param mz_intensities per mz bin intensity
        @return a matrix containing charge ranges for all found masses
     */
    Matrix<int> updateMassBins_(const std::vector<float> &mz_intensities);

    /** @brief Subfunction of updateMassBins_.
        @param mass_intensities per mass bin intensity
        @return a matrix containing charge ranges for all found masses
     */
    Matrix<int> filterMassBins_(const std::vector<float> &mass_intensities);

    /** @brief Subfunction of updateMassBins_. It select candidate masses and update mass_bins_ using the universal pattern, eliminate possible harmonic masses
        @param mass_intensitites mass bin intensities which are updated in te function
        @param mz_intensities mz bin intensities
     */
    void updateCandidateMassBins_(std::vector<float> &mass_intensitites, const std::vector<float> &mz_intensities);

    /** @brief For selected masses in mass_bins_, select the peaks from the original spectrum. Also isotopic peaks are clustered in this function.
        @param per_mass_abs_charge_ranges charge range per mass
     */
    void getCandidatePeakGroups_(const Matrix<int> &per_mass_abs_charge_ranges);

    /// Make the universal pattern.
    void setFilters_();

    /// function for peak group scoring and filtering
    void scoreAndFilterPeakGroups_();

    void removeHarmonicsPeakGroups_();

    /// filter out overlapping masses
    void removeOverlappingPeakGroups_(const double tol, const int iso_length = 1);

    /**@brief Calculate per charge and per isotope intensites from peak groups
     * @param per_isotope_intensity per isotope intensities being calculated
     * @param per_charge_intensity per charge intensities being calculated
     * @param max_isotope_count maximum isotope count
     * @param pg peak group
     */
    std::vector<int> calculatePerChargeIsotopeIntensity_(
        std::vector<double> &per_isotope_intensity,
        std::vector<double> &per_charge_intensity,
        const int max_isotope_count,
        PeakGroup &pg);

    ///Filter out masses with low isotope cosine scores, only retaining current_max_mass_count masses
    void filterPeakGroupsByIsotopeCosine_(const int current_max_mass_count);

    ///Check intensity ratio between charges.
    bool checkChargeDistribution_(const std::vector<double> &per_charge_intensity);

    /** calculate cosine between two vectors a and b with index offset off
     * @param a vector a
     * @param b vector b
     * @param off index offset
     */
    static double getCosine_(const std::vector<double> &a, const std::vector<double> &b, const int off = 0);


    /** calculate cosine between two vectors a and b with additional parameters for fast calculation
     * @param a vector a
     * @param a_start non zero start index of a
     * @param a_end non zero end index of a
     * @param b vector b
     * @param b_size size of b
     * @param offset element index offset between a and b
     */
    static double getCosine_(const std::vector<double> &a,
                             const int &a_start,
                             const int &a_end,
                             const IsotopeDistribution &b,
                             const int &b_size,
                             const int offset);

    /** calculate difference between two vectors a and b with additional parameters for fast calculation
     * @param a vector a
     * @param a_start non zero start index of a
     * @param a_end non zero end index of a
     * @param b vector b
     * @param b_size size of b
     * @param max_offset maximum absolute offset value
     * @param offset element index offset between a and b
     */
    static double getShapeDiff_(const std::vector<double> &a,
                                const int &a_start,
                                const int &a_end,
                                const IsotopeDistribution &b,
                                const int &b_size,
                                const int max_offset,
                                const int offset);

  };
}
