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

namespace OpenMS
{
  class DeconvolutedSpectrum;

  /**
  @brief FLASHDeocnv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset
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
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm& ) = default;

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& other) = default;

    /// assignment operator
    FLASHDeconvAlgorithm& operator= (const FLASHDeconvAlgorithm& fd);

    /**
      @brief main deconvolution function.
      @param dspec spectrum
      @param scanNumber scan number
      @param specIndex index for each spectrum, for file output
      @param massIndex index for each mass, for file output
 */
    void fillPeakGroupsInDeconvolutedSpectrum(DeconvolutedSpectrum& deconvoluted_spectrum,
                                              const int scan_number);

    /// get calculated averagine
    PrecalculatedAveragine getAveragine();

    /** calculate averagine
        @useRNAavg if set, averagine for RNA (nucleotides) is calcualted
     */
    void calculateAveragine(const bool use_RNA_averagine);

    /// convert double to nominal mass
    static int getNominalMass(const double mass);

    /** Examine charge intensity distribution of each peak group
        @perChargeIntensity per charge intensity - aggregated through isotope indices
        @chargeRange max charge range (current_max_charge_ - minCharge)
     */
    static double getChargeFitScore(const std::vector<double>& per_charge_intensity, const int charge_range);

    /** examine intensity distribution over iostope indices. Also determines the most plausible isotope index or, monoisotopic mass
        @mass monoisotopic mass
        @perIsotopeIntensities per isotope intensity - aggregated through charges
        @offset output offset between input monoisotopic mass and determined monoisotopic mass
        @avg precalculated averagine
     */
    static double getIsotopeCosineAndDetermineIsotopeIndex(const double mass,
                                                           const std::vector<double>& per_isotope_intensities,
                                                           int& offset,
                                                           const PrecalculatedAveragine& avg);


  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconv parameters

    /// range of RT subject to analysis (in seconds)
    double min_rt_, max_rt_;
    /// range of mz subject to analysis
    double min_mz_, max_mz_;
    /// min charge and max charge subject to analysis
    int min_charge_, max_charge_;
    /// when a spectrum is deconvoluted, the deconvoluted masses within the this rt window are favorably considered.
    double rt_window_;
    /// mass ranges of deconvolution
    double min_mass_, max_mass_;
    /// max charge is controlled by precursor charge for MSn
    int current_max_charge_;
    /// max mass is controlled by precursor mass for MSn
    double current_max_mass_;
    /// peak intensity threshold subject to analysis
    double intensity_threshold_;
    /// minimum number of peaks supporting a mass
    IntList min_support_peak_count_;
    /// tolerance in ppm for each MS level
    DoubleList tolerance_;
    /// bin size for first stage of mass selection
    DoubleList bin_width_;
    /// cosine threshold between observed and theoretical isotope patterns for each MS level
    DoubleList min_isotope_cosine_;
    /// max mass count per spectrum for each MS level
    IntList max_mass_count_;
    /// number of min mass per spec
    IntList min_mass_count_;

    /// precalculated averagine distributions
    FLASHDeconvHelperStructs::PrecalculatedAveragine avg_;
    /// The data structures for spectra overlapping.
    std::vector<std::vector<Size>> prev_mass_bin_vector_;
    std::vector<double> prev_minbin_logmass_vector_;
    std::vector<double> prev_rt_vector_;

    /// harmonic charge factors that will be considered. For example, 2 is for 1/2 charge harmonic component
    const std::vector<int> harmonic_charges_{2, 3, 5};
    /// Stores log mz peaks
    std::vector<LogMzPeak> log_mz_peaks_;
    /// dspec stores the decovnoluted mass peak groups
    DeconvolutedSpectrum *deconvoluted_spectrum_;

    /// mass_bins_ stores the selected bins for this spectrum + overlapped spectrum (previous a few spectra).
    boost::dynamic_bitset<> mass_bins_;
    /// mz_bins_ stores the binned log mz peaks
    boost::dynamic_bitset<> mz_bins_;
    /// mz_bins_for_edge_effect_ stores the bins to consider edge effect of log mz peak binning
    boost::dynamic_bitset<> mz_bins_for_edge_effect_;

    /// This stores the "universal pattern"
    std::vector<double> filter_;
    /// This stores the patterns for harmonic reduction
    Matrix<double> harmonic_filter_matrix_;

    /// This stores the "universal pattern" in binned dimenstion
    std::vector<int> bin_offsets_;
    /// This stores the patterns for harmonic reduction in binned dimenstion
    Matrix<int> harmonic_bin_offset_matrix_;

    /// minimum mass and mz values representing the first bin of massBin and mzBin, respectively
    double mass_bin_min_value_;
    double mz_bin_min_value_;

    /// current ms Level
    int ms_level_;

    ///static fucntion that converts bin to value
    static double getBinValue_(const Size bin, const double min_value, const double bin_width);

    ///static function that converts value to bin
    static Size getBinNumber_(const double value, const double min_value, const double bin_width);

    ///generate log mz peaks from the input spectrum
    void updateLogMzPeaks_(const MSSpectrum *spec);

    /** Generate mz bins and intensity per mz bin from log mz peaks
        @binNumber number of mz bins
        @mzBinIntensities intensity per mz bin
     */
    void updateMzBins_(const Size& bin_number, std::vector<float>& mz_bin_intensities);

    ///this function takes the previous deconvolution results (from ovelapped spectra) for sensitive deconvolution of the current spectrum
    void unionPrevMassBins_();

    ///Generate peak groups from the input spectrum
    void generatePeakGroupsFromSpectrum_();

    /** Update mass bins. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses
        @mzIntensities per mz bin intensity
        Returns charge range per mass
     */
    Matrix<int> updateMassBins_(const std::vector<float>& mz_intensities);

    /** Subfunction of updateMassBins_.
        @candidateMassBinsForThisSpectrum first candidate mass bins. They are updated using mass intensities in this function
        @massIntensities per mass bin intensity
        Returns charge range per mass
     */
    Matrix<int> filterMassBins_(const std::vector<float>& mass_intensities);

    /** Update mass bins and mass bin intensities. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses
        @massIntensitites mass bin intensities
        @mzIntensities mz bin intensities
     */
    void updateCandidateMassBins_(std::vector<float>& mass_intensitites, const std::vector<float>& mz_intensities);

    /** From selected candidate mass bins, select the peaks. It searches the original spectrum to select peaks. Also isotopic peaks are selected in this function.
        @chargeRanges charge range per mass
     */
    void getCandidatePeakGroups_(const Matrix<int>& charge_ranges);

    /// Make the universal pattern..
    void setFilters_();

    /// function for peak group scoring and filtering
    void scoreAndFilterPeakGroups_();

    /// filter out overlapping masses
    void removeOverlappingPeakGroups_(const double tol);

    void removeOverlappingPeakGroupsWithNominalMass_();

    /// filter out possible harmonics
    void removeHarmonicPeakGroups_(const double tol);

    /// test function to assign a peak to a single mass
    void reassignPeaksinPeakGroups_();

    /** calculate per charge and per isotope intensites from peak groups
     * @perIsotopeIntensity per isotope intensities
     * @perChargeIntensity per charge intensities
     * @maxIsotopeCount maximum isotope count
     * @pg peak groups
     */
    std::vector<int> calculatePerChargeIsotopeIntensity_(
        std::vector<double>& per_isotope_intensity,
        std::vector<double>& per_charge_intensity,
        const int max_isotope_count,
        PeakGroup& pg);

    ///Filter out masses with low isotope cosine scores
    void filterPeakGroupsByIsotopeCosine_(const int current_max_mass_count);

    ///Filter out masses with low QScores
    void filterPeakGroupsByQScore_(const int current_max_mass_count);

    ///For MS1, check intensity ratio between charges.
    bool checkChargeDistribution_(const std::vector<double>& per_charge_intensity);

    /** calculate cosine between two vectors a and b with index offset off
     * @a vector a
     * @b vector b
     * @off index offset
     */
    static double getCosine_(const std::vector<double>& a, const std::vector<double>& b, const int off = 0);


    /** calculate cosine between two vectors a and b with additional parameters for fast calculation
     * @a vector a
     * @aStart non zero start index of a
     * @aEnd non zero end index of a
     * @b vector b
     * @bSize size of b
     * @bNorm precalculated L2 norm of b
     * @index offset
     */
    static double getCosine_(const std::vector<double>& a,
                             const int& a_start,
                             const int& a_end,
                             const IsotopeDistribution& b,
                             const int& b_size,
                             const double& b_norm,
                             const int offset);
  };
}