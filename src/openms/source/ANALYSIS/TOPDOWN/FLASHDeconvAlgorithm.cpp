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
  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(MSExperiment &m, Parameter &p): map(m), param(p)
  {
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
  }

  // FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm&)
  // {
  // }

  FLASHDeconvAlgorithm& FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm& fd)
  {
    //ALWAYS CHECK FOR SELF ASSIGNEMT!
    if (this == &fd) return *this;
    //...
    return *this;
  }

  int FLASHDeconvAlgorithm::getNominalMass(double &m)
  {
    return (int) (m * 0.999497 + .5);
  }

  std::vector<FLASHDeconvAlgorithm::PeakGroup> FLASHDeconvAlgorithm::Deconvolution(int* specCntr, int* qspecCntr,
                                                                                   int* massCntr, int& specIndex, int& massIndex,
                                                                                   FLASHDeconvHelperStructs::PrecalcularedAveragine &avg)
  {
    //calculateAveragines(param);
    float prevProgress = .0;
    std::vector<PeakGroup> allPeakGroups;
    allPeakGroups.reserve(200000);
    //to overlap previous mass bins.

    std::vector<std::vector<Size>> prevMassBinMap;
    std::vector<double> prevMinBinLogMassMap;

    std::unordered_map<UInt,std::map<double, int>> peakChargeMap; // mslevel, mz -> maxCharge should be ordered
    std::unordered_map<UInt,std::unordered_map<double, double>> peakIntMap; // mslevel, mz -> intensity
    std::unordered_map<UInt,std::unordered_map<double, double>> peakMassMap; // mslevel, mz -> mass
    std::unordered_map<UInt,std::unordered_map<double, float>> peakSNRMap; // mslevel, mz -> snr

    std::unordered_map<UInt, int> scanNumberMap;
    std::unordered_map<UInt, int> specIndexMap;

    param.currentMaxMSLevel = 0;

    //if(param.currentMaxMSLevel == 0){
    for (auto it = map.begin(); it != map.end(); ++it)
    {
      auto msLevel = it->getMSLevel();
      param.currentMaxMSLevel = param.currentMaxMSLevel < msLevel? msLevel : param.currentMaxMSLevel;
    }

    param.currentMaxMSLevel = param.currentMaxMSLevel > param.maxMSLevel ? param.maxMSLevel : param.currentMaxMSLevel;
    //}

    prevMassBinMap.reserve(param.numOverlappedScans[0] * 10 );
    //prevMinBinLogMassMap[j] = std::vector<double>();
    prevMinBinLogMassMap.reserve(param.numOverlappedScans[0] * 10 );

    for(UInt j=1;j<=param.currentMaxMSLevel;j++){
      //prevMassBinMap[j] = std::vector<std::vector<Size>>();

      peakChargeMap[j] = std::map<double, int>();
      peakIntMap[j] = std::unordered_map<double, double>();
      peakMassMap[j] = std::unordered_map<double, double>();
      peakSNRMap[j] = std::unordered_map<double, float>();
    }

    auto prevChargeRanges = new int[param.currentMaxMSLevel];
    auto prevMaxMasses = new double[param.currentMaxMSLevel];

    std::fill_n(prevChargeRanges, param.currentMaxMSLevel, param.chargeRange);
    std::fill_n(prevMaxMasses, param.currentMaxMSLevel, param.maxMass);

    int scanNumber = 0;
    //double massMargin = 100.0;
    for (auto it = map.begin(); it != map.end(); ++it)
    {
      ++scanNumber;
      auto msLevel =  it->getMSLevel();
      if (msLevel> param.currentMaxMSLevel)
      {
        continue;
      }

      float progress = (float) (it - map.begin()) / map.size();
      if (progress > prevProgress + .01)
      {
        printProgress(progress); //
        prevProgress = progress;
      }

      specCntr[msLevel-1]++;

      auto sd = SpectrumDeconvolution(*it, param);
      auto precursorMsLevel = msLevel - 1;

      // to find precursor peaks with assigned charges..
      double preMz = 0, preMass= 0, preInt= 0;
      float preSNR = 0;
      int preCharge = 0;

      if (msLevel == 1)
      {
        param.currentChargeRange = param.chargeRange;
        param.currentMaxMass = param.maxMass;
        param.currentMaxMassCount = param.maxMassCount;
      }else{
        auto &subPeakChargeMap = peakChargeMap[precursorMsLevel];
        auto &subPeakIntMap = peakIntMap[precursorMsLevel];
        auto &subPeakMassMap = peakMassMap[precursorMsLevel];
        auto &subPeakSNRMap = peakSNRMap[precursorMsLevel];

        int mc = -1;
        float pint = 0;
        double mm = -1;
        // auto tmz = 0.0;
        for(auto &pre : it->getPrecursors()){
          auto startMz = pre.getIsolationWindowLowerOffset() > 100.0 ? pre.getIsolationWindowLowerOffset() : -pre.getIsolationWindowLowerOffset() + pre.getMZ();
          auto endMz = pre.getIsolationWindowUpperOffset() > 100.0 ? pre.getIsolationWindowUpperOffset() : pre.getIsolationWindowUpperOffset() + pre.getMZ();

          for(auto iter = subPeakChargeMap.begin(); iter != subPeakChargeMap.end(); ++iter)
          {
            auto& mz = iter->first;
            if(mz < startMz){
              continue;
            }
            if(mz > endMz){
              break;
            }
            //tmz=mz;
            auto &cint = subPeakIntMap[mz];
            if(pint < cint){
              pint = cint;
              mc = subPeakChargeMap[mz];
              mm = subPeakMassMap[mz];//mc * mz;
              preSNR = subPeakSNRMap[mz];
              preMz = mz;
              preMass = mm;
              preInt = cint;
              preCharge = mc;
            }
          }
        }
        if(mc > 0){
          param.currentChargeRange = mc  - param.minCharge + 1; //
          param.currentMaxMass = mm + avg.getAverageMassDelta(mm); // isotopie margin

          prevChargeRanges[msLevel - 1] = param.currentChargeRange;
          prevMaxMasses[msLevel-1] =  param.currentMaxMass;
        }else{
          param.currentChargeRange =  prevChargeRanges[msLevel - 1];
          param.currentMaxMass = prevMaxMasses[msLevel - 1];
        }
        param.currentMaxMassCount = 500;
      }

      if(sd.empty()){
        continue;
      }

      auto & peakGroups = sd.getPeakGroupsFromSpectrum(prevMassBinMap, prevMinBinLogMassMap ,avg, msLevel);// FLASHDeconvAlgorithm::Deconvolution (specCntr, qspecCntr, massCntr);

      if (peakGroups.empty())
      {
        continue;
      }

      qspecCntr[msLevel-1]++;
      specIndex++;

      if (msLevel < param.currentMaxMSLevel)
      {
        auto subPeakChargeMap = std::map<double,int>();
        auto subPeakIntMap = std::unordered_map<double,double>();
        auto subPeakMassMap = std::unordered_map<double,double>();
        auto subPeakSNRMap = std::unordered_map<double,float>();

        std::map<double,int>().swap(peakChargeMap[msLevel]);
        std::unordered_map<double,double>().swap(peakIntMap[msLevel]);
        std::unordered_map<double,double>().swap(peakMassMap[msLevel]);
        std::unordered_map<double,float>().swap(peakSNRMap[msLevel]);

        for (auto &pg : peakGroups)
        {
          for (auto &p : pg.peaks)
          {
            int mc = p.charge;
            if (subPeakSNRMap.find(p.mz) != subPeakSNRMap.end())
            {
              float csnr = pg.perChargeSNR[mc];
              if (csnr < subPeakSNRMap[p.mz]){
                continue;
              }
            }
            subPeakChargeMap[p.mz] = mc;
            subPeakIntMap[p.mz] = p.intensity;
            subPeakMassMap[p.mz] = pg.monoisotopicMass;
            subPeakSNRMap[p.mz] = pg.perChargeSNR[mc];
          }
        }
        specIndexMap[msLevel] = specIndex;
        scanNumberMap[msLevel] = scanNumber;

        peakChargeMap[msLevel] = subPeakChargeMap;
        peakIntMap[msLevel] = subPeakIntMap;
        peakMassMap[msLevel] = subPeakMassMap;
        peakSNRMap[msLevel] = subPeakSNRMap;
      }

      for (auto &pg : peakGroups)
      {
        massCntr[msLevel-1]++;
        massIndex++;
        pg.spec = &(*it);
        pg.massIndex = massIndex;
        pg.specIndex = specIndex;
        pg.scanNumber = scanNumber;
        pg.massCntr = (int) peakGroups.size();

        if (msLevel > 1)
        {
          pg.precursorCharge = preCharge;
          pg.precursorMonoMass = preMass;
          pg.precursorIntensity = preInt;
          pg.precursorMz = preMz;
          pg.precursorSpecIndex = specIndexMap[msLevel-1];
          pg.precursorScanNumber = scanNumberMap[msLevel -1];
          pg.precursorSNR = preSNR;
        }
        pg.perChargeSNR.clear();
        allPeakGroups.push_back(pg);
      }
      //allPeakGroups.reserve(allPeakGroups.size() + peakGroups.size());

    }
    printProgress(1); //
    std::cout<<std::endl;
    //allPeakGroups.shrink_to_fit();
    delete[] prevChargeRanges;
    delete[] prevMaxMasses;

    return allPeakGroups; //
  }

  void FLASHDeconvAlgorithm::printProgress(float progress)
  {
    //return; //
    int barWidth = 70;
    std::cout << "[";
    int pos = (int) (barWidth * progress);
    for (int i = 0; i < barWidth; ++i)
    {
      if (i < pos)
      {
        std::cout << "=";
      }
      else if (i == pos)
      {
        std::cout << ">";
      }
      else
      {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
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


