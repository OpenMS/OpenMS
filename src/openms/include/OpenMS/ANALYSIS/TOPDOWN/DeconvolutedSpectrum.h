//--------------------------------------------------------------------------
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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iomanip>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

namespace OpenMS
{
  class PeakGroup;

  class OPENMS_DLLAPI DeconvolutedSpectrum
  {

  public: // TODO protect and public devide..
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    //typedef FLASHDeconvHelperStructs::hash_LogMzPeak hash_LogMzPeak;

    DeconvolutedSpectrum();

    explicit DeconvolutedSpectrum(MSSpectrum &s, int n);

    ~DeconvolutedSpectrum();

    void writeDeconvolutedMasses(std::fstream &fs, Parameter &param);

    void writeAttCsv(std::fstream &fs, int msLevel, double qScoreThreshold, int numMaxMS2);

    void writeMassList(std::fstream &fs, double retDelta, double qScoreThreshold, int numMaxMS2);

    void writeTopFD(std::fstream &fs, int id);

    MSSpectrum toSpectrum();

    static void writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail);

    static void writeAttCsvHeader(std::fstream &fs);

    static void writeThermoInclusionHeader(std::fstream &fs);

    bool empty() const;

    void clearPeakGroupsChargeInfo();

    MSSpectrum *spec;
    std::vector<PeakGroup> peakGroups;
    //std::vector<LogMzPeak> peaks;
    //std::vector<double> mzs; // sorted peaks from the original spectrum

    PeakGroup *precursorPeakGroup = nullptr;
    Precursor precursorPeak;
    //double originalPrecursorIntensity = 0;
    std::string activationMethod;

    bool registerPrecursor(DeconvolutedSpectrum &precursorSpectrum);

    int scanNumber, precursorScanNumber;
  };
}

