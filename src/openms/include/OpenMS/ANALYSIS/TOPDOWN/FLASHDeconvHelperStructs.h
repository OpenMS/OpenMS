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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

namespace OpenMS
{
  /**
   * @brief Wrapper struct for all the structs needed by the FLASHDeconv
   *
   * @see FLASHDeconv
   * @reference: FeatureFinderAlgorithmPickedHelperStructs
   */

  struct OPENMS_DLLAPI FLASHDeconvHelperStructs
  {
    /* FLASHDeconv parameter
    struct OPENMS_DLLAPI Parameter
    {
      int minCharge = 1;
      double minMass = 50;
      double maxMass = 100000;
      double currentMaxMass = 100000;
      DoubleList tolerance;

      double intensityThreshold = 0;// advanced parameters
      DoubleList minIsotopeCosine = {.75, .75};
      double minChargeCosine = .8;

      IntList minContinuousChargePeakCount = {3, 2};
      IntList maxMassCount = {-1, -1};
      unsigned int maxMSLevel = 100;//maxMSL;
      unsigned int currentMaxMSLevel = 100;//maxMSL;

      double RTwindow = 60.0;
      double minRTSpan = 10.0;

      int chargeRange = 100;
      int currentChargeRange = 100;

      bool useRNAavg = false;

      /// automatically set parameters...
      std::vector<int> hCharges{2, 3, 5,};
      DoubleList binWidth;
      UInt minNumOverLappedScans = 15;
      int numOverlappedScans = 15;

      //void print();

      double chargeMass = Constants::PROTON_MASS_U;
    };*/

    struct OPENMS_DLLAPI PrecalculatedAveragine
    {
      std::vector<IsotopeDistribution> isotopes;
      std::vector<double> norms;
      std::vector<double> averageMassDelta;
      std::vector<Size> leftIndices;
      std::vector<Size> rightIndices;
      int maxIsotopeCount;

      double massInterval;
      double minMass;

      PrecalculatedAveragine();

      PrecalculatedAveragine(double m,
                             double M,
                             double delta,
                             CoarseIsotopePatternGenerator *generator,
                             bool useRNAavg);

      IsotopeDistribution get(double mass);
      double getNorm(double mass);
      Size getLeftIndex(double mass);
      Size getRightIndex(double mass);
      double getAverageMassDelta(double mass);

    };

    struct OPENMS_DLLAPI LogMzPeak
    {
      double mz = 0;
      double intensity = 0;
      double logMz = 0;
      double mass = .0;
      int charge = 0;
      int isotopeIndex = -1;

      LogMzPeak();

      explicit LogMzPeak(Peak1D &peak, bool positive);

      explicit LogMzPeak(double mz, bool positive);

      LogMzPeak(LogMzPeak &peak, int c, int i);

      ~LogMzPeak();

      double getUnchargedMass();

      bool operator<(const LogMzPeak &a) const;

      bool operator>(const LogMzPeak &a) const;

      bool operator==(const LogMzPeak &other) const;

    };

    static PrecalculatedAveragine calculateAveragines(double maxMass, bool useRNAavg);

    static double getLogMz(double mz, bool positive);

    static double getChargeMass(bool positive);
  };
}


