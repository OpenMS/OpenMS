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
    ~FLASHDeconvAlgorithm() = default;

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
    FLASHDeconvHelperStructs::PrecalculatedAveragine getAveragine();

    /** calculate averagine
        @useRNAavg if set, averagine for RNA (nucleotides) is calcualted
     */
    void calculateAveragine(bool useRNAavg);

    /// convert double to nominal mass
    static int getNominalMass(double m);

  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconv parameters
    // min charge and max charge of deconvolution
    int minCharge, maxCharge;
    // when a spectrum is deconvoluted, the deconvoluted masses in the spectra within the overlapped scans are favorably considered.
    int numOverlappedScans;
    // mass ranges of deconvolution
    double minMass, maxMass;
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
    double minChargeCosine;
    // max mass count per spectrum for each MS level
    IntList maxMassCount;

    /// precalculated averagine distributions
    FLASHDeconvHelperStructs::PrecalculatedAveragine avg;
    ///The data structures for spectra overlapping.
    std::vector<std::vector<Size>> prevMassBinVector;
    std::vector<double> prevMinBinLogMassVector;
  };
}