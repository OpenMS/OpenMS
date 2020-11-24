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
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm &) = default;

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm &&other) = default;

    /// assignment operator
    FLASHDeconvAlgorithm &operator=(const FLASHDeconvAlgorithm &fd);

    /**
      @brief main deconvolution function.
      @param dspec spectrum
      @param scanNumber scan number
      @param specIndex index for each spectrum, for file output
      @param massIndex index for each mass, for file output
 */
    void fillPeakGroupsInDeconvolutedSpectrum(DeconvolutedSpectrum &dspec,
                                              int scanNumber,
                                              int &specIndex,
                                              int &massIndex);

    /// get calculated averagine
    PrecalculatedAveragine getAveragine();

    /** calculate averagine
        @useRNAavg if set, averagine for RNA (nucleotides) is calcualted
     */
    void calculateAveragine(bool useRNAavg);

    /// convert double to nominal mass
    static int getNominalMass(double m);

    /** Examine charge intensity distribution of each peak group
        @perChargeIntensity per charge intensity - aggregated through isotope indices
        @chargeRange max charge range (currentMaxCharge - minCharge)
     */
    static double getChargeFitScore(const std::vector<double> &perChargeIntensity, int chargeRange);

    /** examine intensity distribution over iostope indices. Also determines the most plausible isotope index or, monoisotopic mass
        @mass monoisotopic mass
        @perIsotopeIntensities per isotope intensity - aggregated through charges
        @offset output offset between input monoisotopic mass and determined monoisotopic mass
        @avg precalculated averagine
     */
    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           const std::vector<double> &perIsotopeIntensities,
                                                           int &offset,
                                                           PrecalculatedAveragine &avg);


  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconv parameters

    /// range of RT subject to analysis (in seconds)
    double minRT, maxRT;
    /// range of mz subject to analysis
    double minMz, maxMz;
    /// min charge and max charge subject to analysis
    int minCharge, maxCharge;
    /// when a spectrum is deconvoluted, the deconvoluted masses in the spectra within the overlapped scans are favorably considered.
    int numOverlappedScans;
    /// mass ranges of deconvolution
    double minMass, maxMass;
    /// max charge is controlled by precursor charge for MSn
    int currentMaxCharge;
    /// max mass is controlled by precursor mass for MSn
    double currentMaxMass;
    /// peak intensity threshold subject to analysis
    double intensityThreshold;
    /// minimum number of peaks supporting a mass
    IntList minSupportPeakCount;
    /// tolerance in ppm for each MS level
    DoubleList tolerance;
    /// bin size for first stage of mass selection
    DoubleList binWidth;
    /// cosine threshold between observed and theoretical isotope patterns for each MS level
    DoubleList minIsotopeCosine;
    /// max mass count per spectrum for each MS level
    IntList maxMassCount;

    /// precalculated averagine distributions
    FLASHDeconvHelperStructs::PrecalculatedAveragine avg;
    /// The data structures for spectra overlapping.
    std::vector<std::vector<Size>> prevMassBinVector;
    std::vector<double> prevMinBinLogMassVector;

    /// harmonic charge factors that will be considered. For example, 2 is for 1/2 charge harmonic component
    const std::vector<int> hCharges{2, 3, 5};
    /// Stores log mz peaks
    std::vector<LogMzPeak> logMzPeaks;
    /// dspec stores the decovnoluted mass peak groups
    DeconvolutedSpectrum *deconvolutedSpectrum;

    /// massBinsForThisSpectrum stores the selected bins only for this spectrum
    boost::dynamic_bitset<> massBinsForThisSpectrum;
    /// massBins stores the selected bins for this spectrum + overlapped spectrum (previous a few spectra).
    boost::dynamic_bitset<> massBins;
    /// mzBins stores the binned log mz peaks
    boost::dynamic_bitset<> mzBins;
    /// mzBinsForEdgeEffect stores the bins to consider edge effect of log mz peak binning
    boost::dynamic_bitset<> mzBinsForEdgeEffect;

    /// This stores the "universal pattern"
    std::vector<double> filter;
    /// This stores the patterns for harmonic reduction
    Matrix<double> harmonicFilter;

    /// This stores the "universal pattern" in binned dimenstion
    std::vector<int> binOffsets;
    /// This stores the patterns for harmonic reduction in binned dimenstion
    Matrix<int> hBinOffsets;

    /// minimum mass and mz values representing the first bin of massBin and mzBin, respectively
    double massBinMinValue;
    double mzBinMinValue;

    /// current ms Level
    int msLevel;

    ///static fucntion that converts bin to value
    static double getBinValue(Size bin, double minV, double binWidth);

    ///static function that converts value to bin
    static Size getBinNumber(double v, double minV, double binWidth);

    ///generate log mz peaks from the input spectrum
    void updateLogMzPeaks(const MSSpectrum *spec);

    /** Generate mz bins and intensity per mz bin from log mz peaks
        @binNumber number of mz bins
        @mzBinIntensities intensity per mz bin
     */
    void updateMzBins(Size &binNumber, std::vector<float> &mzBinIntensities);

    ///this function takes the previous deconvolution results (from ovelapped spectra) for sensitive deconvolution of the current spectrum
    void unionPrevMassBins();

    ///Generate peak groups from the input spectrum
    void generatePeakGroupsFromSpectrum();

    /** Update mass bins. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses
        @mzIntensities per mz bin intensity
        Returns charge range per mass
     */
    Matrix<int> updateMassBins(const std::vector<float> &mzIntensities);

    /** Subfunction of updateMassBins.
        @candidateMassBinsForThisSpectrum first candidate mass bins. They are updated using mass intensities in this function
        @massIntensities per mass bin intensity
        Returns charge range per mass
     */
    Matrix<int> updateMassBins_(const boost::dynamic_bitset<> &candidateMassBinsForThisSpectrum,
                                const std::vector<float> &massIntensities);

    /** Update mass bins and mass bin intensities. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses
        @massIntensitites mass bin intensities
        @mzIntensities mz bin intensities
     */
    boost::dynamic_bitset<> getCandidateMassBinsForThisSpectrum(std::vector<float> &massIntensitites,
                                                                const std::vector<float> &mzIntensities);

    /** From selected candidate mass bins, select the peaks. It searches the original spectrum to select peaks. Also isotopic peaks are selected in this function.
        @chargeRanges charge range per mass
     */
    void getCandidatePeakGroups(const Matrix<int> &chargeRanges);

    /// if the logmz peaks are empty..
    bool empty();

    /// Make the universal pattern..
    void setFilters();

    /// function for peak group scoring and filtering
    void scoreAndFilterPeakGroups();

    /// filter out overlapping masses
    void removeOverlappingPeakGroups(double tol);

    /// filter out possible harmonics
    void removeHarmonicPeakGroups(double tol);

    /// test function to assign a peak to a single mass
    void reassignPeaksinPeakGroups();

    /** calculate per charge and per isotope intensites from peak groups
     * @perIsotopeIntensity per isotope intensities
     * @perChargeIntensity per charge intensities
     * @maxIsotopeCount maximum isotope count
     * @pg peak groups
     */
    std::vector<int> calculatePerChargeIsotopeIntensity(
        std::vector<double> &perIsotopeIntensity,
        std::vector<double> &perChargeIntensity,
        int maxIsotopeCount,
        PeakGroup &pg);

    ///Filter out masses with low isotope cosine scores
    void filterPeakGroupsByIsotopeCosine(int currentMaxMassCount);

    ///Filter out masses with low QScores
    void filterPeakGroupsByQScore(int currentMaxMassCount);

    ///For MS1, check intensity ratio between charges.
    bool checkChargeDistribution(const std::vector<double> &perChargeIntensity);

    /** calculate cosine between two vectors a and b with index offset off
     * @a vector a
     * @b vector b
     * @off index offset
     */
    static double getCosine(const std::vector<double> &a, const std::vector<double> &b, int off = 0);


    /** calculate cosine between two vectors a and b with additional parameters for fast calculation
     * @a vector a
     * @aStart non zero start index of a
     * @aEnd non zero end index of a
     * @b vector b
     * @bSize size of b
     * @bNorm precalculated L2 norm of b
     * @index offset
     */
    static double getCosine(const std::vector<double> &a,
                            int &aStart,
                            int &aEnd,
                            const IsotopeDistribution &b,
                            int &bSize,
                            double &bNorm,
                            int offset);

  };
}