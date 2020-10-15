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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include "OpenMS/ANALYSIS/TOPDOWN/QScore.h"
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

namespace OpenMS
{

  double QScore::getQScore(PeakGroup *pg, double intensity, int charge)
  {
    double score;
    if (pg == nullptr)
    { // all zero
      score = -1.3923 * log10(intensity + 1) + 8.3952;
    }
    else
    {
      score =
          -0.367 * log10(pg->perChargeInfo[charge][0] + 1e-3)
          - 1.4339 * log10(intensity + 1)
          + 1.3297 * log10(pg->perChargeInfo[charge][2] + 1)
          - 2.7199 * pg->perChargeInfo[charge][1]
          - 0.8385 * log10(pg->totalSNR + 1e-3)
          + 4.29 * pg->isotopeCosineScore
          - 0.4913 * pg->chargeCosineScore
          - 0.4561 * log10(pg->intensity + 1)
          + 1.424;
      /*
auto score =
        -0.367 * log10(pg.perChargeSNR[charge] + 1e-3)
        - 1.4339 * log10(preChargeMaxIntensity[j] + 1)
        + 1.3297 * log10(pg.perChargeSumInt[charge] + 1)
        -2.7199 * pg.perChargeICos[charge]
        - 0.8385 * log10(pg.totalSNR + 1e-3)
        + 4.29 * pg.isotopeCosineScore
        - 0.4913 * pg.chargeCosineScore
        - 0.4561 * log10(pg.intensity + 1)
        + 1.424;

          -0.5495 * log10(pg->perChargeSNR[peak->charge] + 1e-3)
              - 0.3263 * log10(peak->intensity + 1)
              - 0.9337 * log10(pg->totalSNR + 1e-3)
              + 2.8442 * pg->isotopeCosineScore
              - 0.9721 * pg->chargeCosineScore
              - 0.1885 * log10(pg->intensity + 1)
              + 1.9703;*/
    }

    return -score;
  }
}
