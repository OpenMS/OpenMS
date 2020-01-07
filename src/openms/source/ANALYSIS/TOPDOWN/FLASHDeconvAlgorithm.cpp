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



  std::vector<FLASHDeconvAlgorithm::PeakGroup> FLASHDeconvAlgorithm::Deconvolution(int &specCntr, int &qspecCntr,
                                         int &massCntr, FLASHDeconvHelperStructs::PrecalcularedAveragine &avg)
  {
    //calculateAveragines(param);
    float prevProgress = .0;
    std::vector<PeakGroup> allPeakGroups;
    allPeakGroups.reserve(200000);
    //to overlap previous mass bins.

    std::map<UInt, std::vector<std::vector<Size>>> prevMassBinMap;
    std::map<UInt,std::vector<double>> prevMinBinLogMassMap;
    std::map<UInt,std::map<double, int>> peakChargeMap; // mslevel, mz -> maxCharge

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      auto msLevel =  it->getMSLevel();
      if (msLevel> param.maxMSLevel)
      {
        continue;
      }
      //if(qspecCntr > 50){
      //  break; //
      //}

      float progress = (float) (it - map.begin()) / map.size();
      if (progress > prevProgress + .01)
      {
        printProgress(progress); //
        prevProgress = progress;
      }

      specCntr++;

      auto sd = SpectrumDeconvolution(*it, param);

      if(prevMassBinMap.find(msLevel) == prevMassBinMap.end()){
        prevMassBinMap[msLevel] = std::vector<std::vector<Size>>();
        prevMinBinLogMassMap[msLevel] = std::vector<double>();
      }
      auto precursorMsLevel = msLevel - 1;

      // to find precursor peaks with assigned charges..
      if (peakChargeMap.find(precursorMsLevel) == peakChargeMap.end())
      {
        param.currentChargeRange = param.chargeRange;
        param.currentMaxMass = param.maxMass;
      }else{
        auto &subPeakChargeMap = peakChargeMap[precursorMsLevel];
        int mc = -1;
        double mm = -1;
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
            int pc = subPeakChargeMap[mz];
            mc = mc > pc? mc : pc;
            double pm = mc * mz;
            mm = mm > pm ? mm : pm;
          }
        }
        if(mc > 0){
          param.currentChargeRange = mc;
          param.currentMaxMass = mm;
          //std::cout<<mc<< " " << mm << std::endl;
        }else{
          param.currentChargeRange = param.chargeRange;
          param.currentMaxMass = param.maxMass;
        }
      }

      auto & peakGroups = sd.getPeakGroupsFromSpectrum(prevMassBinMap[msLevel], prevMinBinLogMassMap[msLevel],avg, msLevel);// FLASHDeconvAlgorithm::Deconvolution (specCntr, qspecCntr, massCntr);

      auto subPeakChargeMap = std::map<double, int>();
      for (auto& pg : peakGroups)
      {
        for(auto& p : pg.peaks){
          int mc = p.charge;
          if (subPeakChargeMap.find(p.mz) == subPeakChargeMap.end()){
            int pc = subPeakChargeMap[p.mz];
            mc = mc > pc? mc : pc;
          }
          subPeakChargeMap[p.mz] = mc;
        }
      }
      peakChargeMap[msLevel] = subPeakChargeMap;

      if (peakGroups.empty())
      {
        continue;
      }

      qspecCntr++;

      //allPeakGroups.reserve(allPeakGroups.size() + peakGroups.size());
      for (auto &pg : peakGroups)
      {
        massCntr++;
        pg.spec = &(*it);
        pg.massIndex = massCntr;
        pg.specIndex = qspecCntr;
        pg.massCntr = (int) peakGroups.size();
        allPeakGroups.push_back(pg);
      }
    }
    printProgress(1); //
    //allPeakGroups.shrink_to_fit();
    return allPeakGroups; //
  }

  void FLASHDeconvAlgorithm::printProgress(float progress)
  {
    //return;
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

/*
  void
  FLASHDeconvAlgorithm::filterPeakGroupsByIntensity(std::vector<PeakGroup> &peakGroups, std::vector<double> &intensities, const Parameter &param)
  {
    if (param.maxMassCount < 0 || intensities.size() <= (Size) param.maxMassCount)
    {
      return;
    }
    Size mc = (Size) param.maxMassCount;
    sort(intensities.begin(), intensities.end());
    auto threshold = intensities[intensities.size() - mc];
    for (auto pg = peakGroups.begin(); pg != peakGroups.end();)
    {
      if (peakGroups.size() <= mc)
      {
        break;
      }
      if (pg->intensity < threshold)
      {
        pg = peakGroups.erase(pg);
        continue;
      }
      ++pg;
    }
  }
*/
}

