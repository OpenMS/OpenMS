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
      @brief main deconvolution function
      @param spec spectrum
      @param scanNumber scan number
      @param specIndex index for each spectrum, for file output
      @param massIndex index for each mass, for file output
 */
    void getPeakGroups(DeconvolutedSpectrum &spec,
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

    // examine intensity distribution over charges
    static double getChargeFitScore(double *perChargeIntensity, int chargeRange);

    // examine intensity distribution over iostope indices. Also determines the most plausible isotope index or, monoisotopic mass
    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           std::vector<double> &perIsotopeIntensities,
                                                           int &offset,
                                                           PrecalculatedAveragine &avg);


  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconv parameters
    double minRT, maxRT;
    double minMz, maxMz;
    // min charge and max charge of deconvolution
    int minCharge, maxCharge;
    // when a spectrum is deconvoluted, the deconvoluted masses in the spectra within the overlapped scans are favorably considered.
    int numOverlappedScans;
    // mass ranges of deconvolution
    double minMass, maxMass;
    int currentMaxCharge;
    double currentMaxMass;

    double intensityThreshold;
    // minimum number of peaks supporting a mass
    IntList minSupportPeakCount;
    // tolerance in ppm for each MS level
    DoubleList tolerance;
    // bin size for first stage of mass selection
    DoubleList binWidth;
    // cosine threshold between observed and theoretical isotope patterns for each MS level
    DoubleList minIsotopeCosine;
    // cosien thereshold between charge distribution and fit gaussian
    //double minChargeCosine;
    // max mass count per spectrum for each MS level
    IntList maxMassCount;

    /// precalculated averagine distributions
    FLASHDeconvHelperStructs::PrecalculatedAveragine avg;
    ///The data structures for spectra overlapping.
    std::vector<std::vector<Size>> prevMassBinVector;
    std::vector<double> prevMinBinLogMassVector;

    const std::vector<int> hCharges{2, 3, 5,};
    //Stores log mz peaks
    std::vector<LogMzPeak> logMzPeaks;
    //Bins for mass and mz. Bin size is deteremine by Parameter.tolerance
    //peakGroups stores the decovnoluted mass peak groups
    std::vector<PeakGroup> peakGroups;

    //massBinsForThisSpectrum stores the selected bins only for this spectrum
    boost::dynamic_bitset<> massBinsForThisSpectrum;
    //massBins stores the selected bins for this spectrum + overlapped spectrum (previous a few spectra).
    boost::dynamic_bitset<> massBins;
    //mzBins stores the binned log mz peaks
    boost::dynamic_bitset<> mzBins;

    //This stores the "universal pattern"
    std::vector<double> filter;
    //This stores the patterns for harmonic reduction
    Matrix<double> harmonicFilter;

    //This stores the "universal pattern" in binned dimenstion
    std::vector<int> binOffsets;
    //This stores the patterns for harmonic reduction in binned dimenstion
    Matrix<int> hBinOffsets;

    double massBinMinValue;
    double mzBinMinValue;
    int msLevel;

    //static fucntion that converts bin to value
    static double getBinValue(Size bin, double minV, double binWidth);

    //static function that converts value to bin
    static Size getBinNumber(double v, double minV, double binWidth);

    //generate log mz peaks from the input spectrum
    void updateLogMzPeaks(MSSpectrum *spec);

    //generate mz bins from log mz peaks
    void updateMzBins(Size &binNumber, std::vector<float> &mzBinIntensities);

    //this function takes the previous deconvolution results (from ovelapped spectra) for sensitive deconvolution of the current spectrum
    void unionPrevMassBins();

    //Main function
    void generatePeakGroupsFromSpectrum();

    //Update mass bins from mz bins and universal pattern. It select candidate mass bins using the pattern, eliminate possible harmonic masses
    Matrix<int> updateMassBins(//float *massIntensities,
        std::vector<float> &mzIntensities);

    //Subfunction of updateMassBins.
    Matrix<int> updateMassBins_(boost::dynamic_bitset<> &candidateMassBinsForThisSpectrum,
                                std::vector<float> &massIntensities,
                                long &binStart, long &binEnd);

    //Subfunction of updateMassBins.
    boost::dynamic_bitset<> getCandidateMassBinsForThisSpectrum(std::vector<float> &massIntensitites,
                                                                std::vector<float> &mzIntensities);

    //From selected candidate mass bins, reselect the peaks. It researches the original spectrum to select peaks.
    //Also isotopic peaks are selected in this function.
    void getCandidatePeakGroups(//float *sumLogIntensities,
        Matrix<int> &chargeRanges);

    bool empty();

    //Make the universal pattern here..
    void setFilters();

    /////////


    //the main function of this class
    void scoreAndFilterPeakGroups();

    //filter out overlapping masses
    void removeOverlappingPeakGroups(double tol);

    //filter out possible harmonics
    void removeHarmonicPeakGroups(double tol);

    void reassignPeaksinPeakGroups();

    //From peaks distributions over charge and isotope are calculated
    std::vector<int> updatePerChargeIsotopeIntensity(
        std::vector<double> &perIsotopeIntensity,
        std::vector<double> &perChargeIntensity,
        int maxIsotopeCount,
        PeakGroup &pg);

    //Filter out masses with low isotope cosine scores
    void filterPeakGroupsByIsotopeCosine(int currentMaxMassCount);

    //Filter out masses with low QScores
    void filterPeakGroupsByQScore(int currentMaxMassCount);

    //For MS1, check intensity ratio between charges.
    bool checkChargeDistribution(std::vector<double> &perChargeIntensity);

    //cosine function
    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    //cosine function
    static double getCosine(const double *a, double *b, Size size);

    //cosine function for fast calculatoin
    static double getCosine(std::vector<double> &a,
                            int &aStart,
                            int &aEnd,
                            IsotopeDistribution &b,
                            int &bSize,
                            double &bNorm,
                            int offset);
  };
}