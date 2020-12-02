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

namespace OpenMS
{
  /**
   * @brief FLASHIda class for real time deconvolution
   *
   * @see FLASHIdaBridgeFunctions
   * @reference: FeatureFinderAlgorithmPickedHelperStructs
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
    FLASHIda(const FLASHIda &) = default;

    /// move constructor
    FLASHIda(FLASHIda &&other) = default;

    /// assignment operator
    FLASHIda &operator=(const FLASHIda &fd) = default;

    /**
           @brief get peak groups from input spectrum, specified by mzs and intensities (due to C# interface it is necessary)
           @param mzs mz values of the input spectrum
           @param ints intensities of the input spectrum
           @param length length of mzs and ints
           @param rt Retention time
           @param msLevel ms level
           @param spectrum name
      */
    int getPeakGroups(double *mzs,
                      double *ints,
                      int length,
                      double rt,
                      int msLevel,
                      char *name);

    /**
           @brief get isolation windows
           @param wstart window start mzs
           @param wend windo end mzs
           @param qScores q scores of windows
           @param charges charges of windows
           @avgMasses average masses of windows
      */
    void getIsolationWindows(double *wstart, double *wend, double *qScores, int *charges, double *avgMasses);

  private:

    /// PeakGroup comparator for soring by QScore
    struct
    {
      bool operator()(const PeakGroup& a, const PeakGroup& b) const
      {
        return a.getQScore() > b.getQScore();
      }
    } QscoreComparator;

    /// Selected integer masses - necessary for mass exclusion
    std::map<int, std::vector<double>> selected; // int mass, rt, qscore
    /// precalculated averagine for fast selection
    PrecalculatedAveragine avg;
    /// discard peak groups using mass exclusion
    void filterPeakGroupsUsingMassExclusion(MSSpectrum &spec, int msLevel);
    /// generate MSSpectrum class using mzs and intensities
    static MSSpectrum makeMSSpectrum(double *mzs, double *ints, int length, double rt, int msLevel, char *name);
    /// deconvoluted spectrum that contains the peak groups
    DeconvolutedSpectrum deconvolutedSpectrum;
    /// FLASHDeconvAlgorithm class for deconvolution
    FLASHDeconvAlgorithm fd;
    /// q score threshold - determined from C# side
    double qScoreThreshold;
    /// retention time window - determined from C# side
    double RTwindow;
    /// how many masses will be selected per ms level? - determined from C# side
    IntList massCount;
  };
}