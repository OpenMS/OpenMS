//
// Created by Kyowon Jeong on 2/14/20.
//

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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvRealTimeAlgorithm.h>

namespace OpenMS
{
/*
  FLASHDeconvRealTimeAlgorithm::mzRange::MzRange(double mz1, double mz2, double rt, int mass): mz1(mz1), mz2(mz2), rt(rt), nominalMass(mass)
  {
  }


  FLASHDeconvRealTimeAlgorithm::ShortMassTrace::ShortMassTrace(double rt): lastRt(rt), color('x')
  {
    intensities = std::queue<double>();
  }

  void FLASHDeconvRealTimeAlgorithm::ShortMassTrace::setColor(char c)
  {
    color = c;
  }

  void FLASHDeconvRealTimeAlgorithm::ShortMassTrace::updateColor(){
    if(intensities.size() < 3){
      return;
    }

  }

  void FLASHDeconvRealTimeAlgorithm::ShortMassTrace::addIntensity(double intensity)
  {
    intensities.push(intensity);
  }

  void FLASHDeconvRealTimeAlgorithm::ShortMassTrace::trim()
  {
    while(intensities.size()>=3){
      intensities.pop();
    }
  }


  FLASHDeconvRealTimeAlgorithm::FLASHDeconvRealTimeAlgorithm(MSSpectrum &s, std::vector<int> greenMasses, Parameter &p): spec(s), param(p)
  {
    for (auto iter = massMemory.begin(); iter != massMemory.cend();)
    {
      auto v = iter->second;
      if(v.lastRt > spec.getRT() - deltaRtForMassMemory){// no delete
        v.trim();
        ++iter;
        continue;
      }
      massMemory.erase(iter++);
    }

    for (auto& m : greenMasses)
    {
      if(massMemory.find(m) == massMemory.end()){ // no mass
        auto mt = ShortMassTrace(spec.getRT());
        mt.setColor('g');
        //mt.addIntensity(0); //
        massMemory[m] = mt;
      }else{
        auto& mt= massMemory[m];
        mt.setColor('g');
      }
    }
  }

  FLASHDeconvRealTimeAlgorithm::~FLASHDeconvRealTimeAlgorithm(){

  }

  bool FLASHDeconvRealTimeAlgorithm::compareIntensity(const FLASHDeconvHelperStructs::PeakGroup &pg1,
                                                      const FLASHDeconvHelperStructs::PeakGroup &pg2)
  {
    if(pg1.intensity == pg2.intensity)
      return pg1.monoisotopicMass < pg2.monoisotopicMass;
    return pg1.intensity  < pg2.intensity ;
  }



  std::vector<FLASHDeconvRealTimeAlgorithm::MzRange>& FLASHDeconvRealTimeAlgorithm::Deconvolution(FLASHDeconvHelperStructs::PrecalcularedAveragine &avg){
    auto sd = SpectrumDeconvolution(spec, param);

    std::vector<std::vector<Size>> prevMassBinMap;
    std::vector<double> prevMinBinLogMassMap;
    UInt msLevel = 1;
    auto & peakGroups = sd.getPeakGroupsFromSpectrum(prevMassBinMap, prevMinBinLogMassMap, avg, msLevel);// FLASHDeconvAlgorithm::Deconvolution (specCntr, qspecCntr, massCntr);
    //std::sort(peakGroups.begin(),peakGroups.end(), &FLASHDeconvRealTimeAlgorithm::compareIntensity);

    auto mzRanges = std::vector<FLASHDeconvRealTimeAlgorithm::MzRange>();
    for (auto &pg : peakGroups)
    {
      auto nominalMass = FLASHDeconvAlgorithm::getNominalMass(pg.monoisotopicMass);
      ShortMassTrace mt;
      if(massMemory.find(nominalMass) != massMemory.end()){ // mass found
        mt = massMemory[nominalMass];
        mt.addIntensity(pg.intensity);
      }else{
        mt = ShortMassTrace(spec.getRT());
      }
      if(mt.color == 'g'){
        continue;
      }






    }


  }
*/


}
