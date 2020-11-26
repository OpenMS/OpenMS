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

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <iomanip>

namespace OpenMS
{
  class PeakGroup;

  /**
       @brief A class representing a deconvoluted spectrum. Also contains deconvoluted precursro information for MSn n>1.
  */
  class OPENMS_DLLAPI DeconvolutedSpectrum :
      public std::vector<PeakGroup>
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    DeconvolutedSpectrum() = default;

    /**
       @brief Constructor for DeconvolutedSpectrum
       @param s spectrum
       @param n scan number
  */
    explicit DeconvolutedSpectrum(const MSSpectrum &s, int n);

    /// default deconstructor
    ~DeconvolutedSpectrum() = default;

    /// copy constructor
    DeconvolutedSpectrum(const DeconvolutedSpectrum &) = default;

    /// move constructor
    DeconvolutedSpectrum(DeconvolutedSpectrum &&other) = default;

    /// assignment operator
    DeconvolutedSpectrum &operator=(const DeconvolutedSpectrum &fd) = default;

    /**
        @brief write the header in the output file (spectrum level)
        @param fs file stream to the output file
        @param n the index to the mass
        @param detail if set true, detailed information of the mass (e.g., peak list for the mass) is written
   */
    static void writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail);

    /**
      @brief write the deconvoluted masses in the output file (spectrum level)
      @param fs file stream to the output file
      @param param FLASHDeconv paramter
    */
    void writeDeconvolutedMasses(std::fstream &fs,
                                 const String &fileName,
                                 const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                 bool writeDetail);

    /**
      @brief write the deconvoluted masses TopFD format
      @param fs file stream to the output file
      @param id the index to the spectrum. updated outside.
 */
    void writeTopFD(std::fstream &fs, int id);

    /// cast DeconvolutedSpectrum into MSSpectrum object to write mzml format
    MSSpectrum toSpectrum();

    /// write the header for Thermo Inclusion List header format
    static void writeThermoInclusionHeader(std::fstream &fs);

    /// for memory save... clear unnecessary information in mass tracing
    void clearPeakGroupsChargeInfo();

    /**
     @brief register the precusor info in this from the precursor DeconvolutedSpectrum
     @param precursorSpectrum the precursor DeconvolutedSpectrum
     */
    bool registerPrecursor(DeconvolutedSpectrum &precursorSpectrum);

    /// original spectrum setter
    MSSpectrum &getOriginalSpectrum();

    /// peakGroup getter
    PeakGroup getPrecursorPeakGroup();

    /// precursor charge getter : set in registerPrecursor
    int getPrecursorCharge();

    /// get max mass - which is min mass between max mass specified by users or mass determined by precursor mass for MSn
    double getCurrentMaxMass(double maxMass);

    /// get max charge - which is min charge between max charge specified by users or charge determined by precursor charge for MSn
    int getCurrentMaxCharge(int maxCharge);

  private:
    /// the original spectrum from which this is generated
    MSSpectrum spec;
    /// precursor peakGroup (or mass)
    PeakGroup *precursorPeakGroup = nullptr;
    /// precursor peak (not deconvoluted one)
    Precursor precursorPeak;
    /// activation method for file output
    std::string activationMethod;
    /// scan number and precursor scan number
    int scanNumber, precursorScanNumber;
  };
}