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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>

namespace OpenMS
{
   

  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(FLASHDeconvHelperStructs::PrecalculatedAveragine &a, Parameter &p) :
      param(p), avg(a)
  {
    prevMassBinMap = std::vector<std::vector<Size>>();
    prevMinBinLogMassMap = std::vector<double>();
    //peakGroups = std::vector<PeakGroup>();
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
    //std::vector<PeakGroup>().swap(peakGroups);
  }

  FLASHDeconvAlgorithm &FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm &fd)
  {
    //ALWAYS CHECK FOR SELF ASSIGNEMT!
    if (this == &fd)
    {
      return *this;
    }
    //...
    return *this;
  }

  int FLASHDeconvAlgorithm::getNominalMass(double &m)
  {
    return (int) (m * 0.999497 + .5);
  }

  void FLASHDeconvAlgorithm::getPeakGroups(DeconvolutedSpectrum &dspec, int& specIndex, int& massIndex)
  {

    auto* spec = dspec.spec;
    int msLevel = spec->getMSLevel();
    if (msLevel == 1 || dspec.precursorPeakGroup == nullptr) {
      param.currentMaxMass = param.maxMass;
      param.currentChargeRange = param.chargeRange;
    }
    else{
     param.currentChargeRange = dspec.precursorPeak.getCharge();
     param.currentMaxMass = dspec.precursorPeakGroup->monoisotopicMass;
    }

    auto sd = SpectrumDeconvolution(*spec, param);

    dspec.peakGroups = sd.getPeakGroupsFromSpectrum(prevMassBinMap,
                                               prevMinBinLogMassMap,
                                               avg, msLevel);


    if (dspec.empty())
    {
      return;
    }

    for (auto &pg : (dspec.peakGroups))
    {
      sort(pg.peaks.begin(), pg.peaks.end());
      pg.spec = spec;
      pg.specIndex = specIndex;
      pg.scanNumber = dspec.scanNumber;
      pg.massIndex = massIndex++;
    }

    specIndex++;
  }
}


