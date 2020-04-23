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

  // constructor

  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(int &specIndex, int &massIndex, FLASHDeconvHelperStructs::PrecalcularedAveragine &a, Parameter &p) :
      param(p), avg(a), specIndex(specIndex), massIndex(massIndex)
  {
    prevMassBinMap = std::vector<std::vector<Size>>();
    prevMinBinLogMassMap = std::vector<double>();

    prevMassBinMap.reserve(param.numOverlappedScans[0] * 10);
    prevMinBinLogMassMap.reserve(param.numOverlappedScans[0] * 10);
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
    std::vector<std::vector<Size>>().swap(prevMassBinMap);
    std::vector<double>().swap(prevMinBinLogMassMap);
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


  void FLASHDeconvAlgorithm::Deconvolution(DeconvolutedSpectrum &deconvolutedSpectrum, int &scanNumber)
  {
    auto sd = SpectrumDeconvolution(*deconvolutedSpectrum.spec, param);
    //std::vector<SpectrumDeconvolution::PeakGroup> peakGroups;
    if (sd.empty())
    {
      return;
    }

    deconvolutedSpectrum.peakGroups = sd.getPeakGroupsFromSpectrum(prevMassBinMap,
                                              prevMinBinLogMassMap,
                                              avg,
                                              deconvolutedSpectrum.spec->getMSLevel());// FLASHDeconvAlgorithm::Deconvolution (specCntr, qspecCntr, massCntr);

    if (deconvolutedSpectrum.empty())
    {
      return;
    }

    deconvolutedSpectrum.peaks = sd.logMzPeaks;
    deconvolutedSpectrum.specIndex = specIndex;
    deconvolutedSpectrum.scanNumber = scanNumber;
    deconvolutedSpectrum.massCntr = (int) deconvolutedSpectrum.peakGroups.size();

    for (auto &pg : deconvolutedSpectrum.peakGroups)
    {
      sort(pg.peaks.begin(), pg.peaks.end());
      pg.deconvSpec = &deconvolutedSpectrum;
      pg.massIndex = massIndex++;
    }
    specIndex++;
  }




  /*
    bool FLASHDeconvAlgorithm::checkSpanDistribution(int *mins, int *maxs, int range, int threshold)
    {
      int nonZeroStart = -1, nonZeroEnd = 0;
      int maxSpan = 0;

      for (int i = 0; i < range; i++)
      {
        if (maxs[i] >= 0)
        {
          if (nonZeroStart < 0)
          {
            nonZeroStart = i;
          }
          nonZeroEnd = i;
          maxSpan = std::max(maxSpan, maxs[i] - mins[i]);
        }
      }
      if (maxSpan <= 0)
      {
        return false;
      }

      int prevCharge = nonZeroStart;
      int n_r = 0;

      double spanThreshold = maxSpan / 1.5;//

      for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++)
      {
        if (maxs[k] < 0)
        {
          continue;
        }


        if (k - prevCharge == 1)
        {
          int intersectSpan = std::min(maxs[prevCharge], maxs[k])
                              - std::max(mins[prevCharge], mins[k]);

          if (spanThreshold <= intersectSpan)
          { //
            n_r++;
          }
        }
        prevCharge = k;
      }

      if (n_r < threshold)
      {
        return -100.0;
      }

      for (int i = 2; i < std::min(12, range); i++)
      {
        for (int l = 0; l < i; l++)
        {
          int t = 0;
          prevCharge = nonZeroStart + l;
          for (int k = prevCharge + i; k <= nonZeroEnd; k += i)
          {
            if (maxs[k] < 0)
            {
              continue;
            }


            if (k - prevCharge == i)
            {
              int intersectSpan = std::min(maxs[prevCharge], maxs[k])
                                  - std::max(mins[prevCharge], mins[k]);

              if (spanThreshold <= intersectSpan)
              {
                t++;
              }
            }
            prevCharge = k;
          }
          if (n_r <= t)
          {
            return false;
          }
        }
      }

      return true;
    }


   */



}


