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

#include "OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h"
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <utility>

namespace OpenMS
{
  MassFeatureTrace::MassFeatureTrace() :
      DefaultParamHandler("MassFeatureTrace")
  {
    Param mtd_defaults = MassTraceDetection().getDefaults();

    //mtd_defaults.setValue("mass_error_da", 1.5);
    mtd_defaults.setValue("min_trace_length", 10.0);

    //mtd_defaults.remove("mass_error_da");
    mtd_defaults.remove("chrom_peak_snr");
    defaults_.insert("", mtd_defaults);

    //defaults_.setValue("min_charge_cosine",
    //                   .5,
    //                   "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)");
    defaults_.setValue("min_isotope_cosine", .75, "cosine threshold between avg. and observed isotope pattern for MS1");
    defaultsToParam_();
  }

  MassFeatureTrace::~MassFeatureTrace()
  {
    for (auto &item : peakGroupMap)
    {
      std::unordered_map<double, PeakGroup>().swap(item.second);
    }
    std::unordered_map<double, std::unordered_map<double, PeakGroup>>().swap(peakGroupMap);
  }


  void MassFeatureTrace::findFeatures(const String &fileName, bool promexOut, int &featureCntr,
                                      int &featureIndex,
                                      std::fstream &fsf,
                                      std::fstream &fsp,
                                      PrecalculatedAveragine averagines)
  {
    MSExperiment map;
    std::map<int, MSSpectrum> indexSpecMap;
    int minCharge = INT_MAX;
    int maxCharge = INT_MIN;

    for (auto &item : peakGroupMap)
    {
      auto rt = item.first;
      MSSpectrum deconvSpec;
      deconvSpec.setRT(rt);
      for (auto &pg : item.second)
      {
        maxCharge = maxCharge > pg.second.maxCharge ? maxCharge : pg.second.maxCharge;
        minCharge = minCharge < pg.second.minCharge ? minCharge : pg.second.minCharge;

        Peak1D tp(pg.first, (float) pg.second.intensity);
        deconvSpec.push_back(tp);
      }
      map.addSpectrum(deconvSpec);
    }

    if (map.size() < 3)
    {
      return;
    }

    map.sortSpectra();
    MassTraceDetection mtdet;
    Param mtd_param = getParameters().copy("");
    //mtd_param.remove("min_charge_cosine");
    mtd_param.remove("min_isotope_cosine");

    mtdet.setParameters(mtd_param);
    std::vector<MassTrace> m_traces;

    mtdet.run(map, m_traces);  // m_traces : output of this function

    int chargeRange = maxCharge - minCharge + 1;

    //std::cout<<chargeRange << " " << averagines.maxIsotopeIndex<<std::endl;
    int tmp[] = {0, 0, 0, 0};

    for (auto &mt : m_traces)
    {
      int minFCharge = INT_MAX; // min feature charge
      int maxFCharge = INT_MIN; // max feature charge

      int minFIso = INT_MAX; // min feature isotope index
      int maxFIso = INT_MIN; // max feature isotope index

      auto perChargeIntensity = std::vector<double>(chargeRange + 1, 0);
      auto perChargeMaxIntensity = std::vector<double>(chargeRange + 1, 0);
      auto perChargeMz = std::vector<double>(chargeRange + 1, 0);
      auto perIsotopeIntensity = std::vector<double>(averagines.maxIsotopeIndex, 0);

      int minScanNum = (int) map.size() + 1000;
      int maxScanNum = 0;

      int repScan = 0, repCharge = 0;
      double maxIntensity = 0;
      double maxMass = .0;
      double maxIso = 0;
      boost::dynamic_bitset<> charges(chargeRange + 1);

      for (auto &p2 : mt)
      {
        auto &pgMap = peakGroupMap[p2.getRT()];
        auto &pg = pgMap[p2.getMZ()];
        auto scanNumber = pg.scanNumber;

        minFCharge = minFCharge < pg.minCharge ? minFCharge : pg.minCharge;
        maxFCharge = maxFCharge > pg.maxCharge ? maxFCharge : pg.maxCharge;

        minScanNum = minScanNum < scanNumber ? minScanNum : scanNumber;
        maxScanNum = maxScanNum > scanNumber ? maxScanNum : scanNumber;

        if (pg.intensity > maxIntensity)
        {
          maxIntensity = pg.intensity;
          repScan = scanNumber;

        }

        if (pg.isotopeCosineScore > maxIso)
        {
          maxIso = pg.isotopeCosineScore;
          maxMass = pg.monoisotopicMass;
        }

        for (auto &p : pg)
        {
          if (p.isotopeIndex < 0 || p.isotopeIndex >= averagines.maxIsotopeIndex || p.charge < minCharge ||
              p.charge >= chargeRange + minCharge + 1)
          {
            continue;
          }

          charges[p.charge - minCharge] = true;
          perChargeIntensity[p.charge - minCharge] += p.intensity;
          perIsotopeIntensity[p.isotopeIndex] += p.intensity;
          if (perChargeMaxIntensity[p.charge - minCharge] > p.intensity)
          {
            continue;
          }
          perChargeMaxIntensity[p.charge - minCharge] = p.intensity;
          perChargeMz[p.charge - minCharge] = p.mz;
        }
      }

      //double chargeScore = FLASHDeconvAlgorithm::getChargeFitScore(perChargeIntensity,
      //                                                            chargeRange);
      //if (chargeScore < minChargeCosine) //
      //{
      //  continue;
      //}

      int offset = 0;

      double mass = mt.getCentroidMZ();
      double isoScore = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                                       perIsotopeIntensity,
                                                                                       offset, averagines);
      if (isoScore < minIsotopeCosine)
      {
        continue;
      }

      //perIsotopeIntensity, param.maxIsotopeIndex


      if (offset != 0)
      {
        mass += offset * Constants::ISOTOPE_MASSDIFF_55K_U;
      }

      auto sumInt = .0;

      for (auto &p : mt)
      {
        sumInt += p.getIntensity();
      }
      //mt.computeSmoothedPeakArea()
      //mt.estimateFWHM(true);
      auto avgMass = averagines.getAverageMassDelta(mass) + mass;
      ++featureCntr;
      fsf << featureIndex++ << "\t" << fileName << "\t" << std::to_string(mass) << "\t"
          << std::to_string(avgMass) << "\t" // massdiff
          << mt.getSize() << "\t"
          << mt.begin()->getRT() << "\t"
          << mt.rbegin()->getRT() << "\t"
          << mt.getTraceLength() << "\t"
          << mt[mt.findMaxByIntPeak()].getRT() << "\t"
          << sumInt << "\t"
          << mt.getMaxIntensity(false) << "\t"
          << mt.computePeakArea() << "\t"
          << minFCharge << "\t"
          << maxFCharge << "\t"
          << charges.count() << "\t"
          << isoScore << "\t";

      for (int i = minFCharge; i <= maxFCharge; i++)
      {
        fsf << perChargeIntensity[i - minCharge];
        if (i < maxFCharge)
        {
          fsf << ";";
        }
      }
      fsf << "\t";
      int isoEndIndex = 0;

      for (int i = 0; i < averagines.maxIsotopeIndex; i++)
      {
        if (perIsotopeIntensity[i] == 0)
        {
          continue;
        }
        isoEndIndex = i;
      }
      for (int i = 0; i <= isoEndIndex; i++)
      {
        fsf << perIsotopeIntensity[i];
        if (i < isoEndIndex)
        {
          fsf << ";";
        }
      }
      /*fsf<< "\t";
      double tmp2[maxFCharge + 1][isoEndIndex + 1];
      for (int i = 0; i <= maxFCharge; i++)
      {
        for (int j = 0; j <= isoEndIndex; j++)
        {
          tmp2[i][j] = .0;
        }
      }
      for (auto &p2 : mt)
      {
        auto &pgMap = peakGroupMap[p2.getRT()];
        auto &pg = pgMap[p2.getMZ()];
        for (auto &p: pg)
        {
          tmp2[p.charge][p.isotopeIndex] += p.intensity;
        }
      }
      for (int i = 0; i <= maxFCharge; i++)
      {
        for (int j = 0; j <= isoEndIndex; j++)
        {
          fsf<<tmp2[i][j]<<",";
        }
        fsf<<";";
      }*/

      fsf << "\n";

      if (abs(avgMass - 44086) < 3)
      {
        tmp[0] = 1;
      }
      if (abs(avgMass - 44166) < 3)
      {
        tmp[1] = 1;
      }
      if (abs(avgMass - 44328) < 3)
      {
        tmp[2] = 1;
      }
      if (abs(avgMass - 44369) < 3)
      {
        tmp[3] = 1;
      }

      if (promexOut)
      {
        double maxChargeIntensity = 0;
        for (int c = 0; c < chargeRange; c++) // c is charge range!!
        {
          if (perChargeIntensity[c] > maxChargeIntensity)
          {
            maxChargeIntensity = perChargeIntensity[c];
            repCharge = c + minCharge;
          }
        }
        auto apex = mt[mt.findMaxByIntPeak()];
        //int si = rtSpecMap[(float) apex.getRT()];
        auto &spgMap = peakGroupMap[apex.getRT()];
        auto &spg = spgMap[apex.getMZ()];

        fsp << featureIndex << "\t" << minScanNum << "\t" << maxScanNum << "\t" << minFCharge << "\t"
            << maxFCharge << "\t" << std::to_string(mass) << "\t" << std::fixed << std::setprecision(2)
            << repScan << "\t" << repCharge << "\t" << perChargeMz[repCharge] << "\t" << sumInt << "\t"
            << spg.scanNumber << "\t" << spg.intensity << "\t"
            << mt.begin()->getRT() / 60.0 << "\t"
            << mt.rbegin()->getRT() / 60.0 << "\t"
            << mt.getTraceLength() / 60.0 << "\t";


        for (int j = 0; j < averagines.maxIsotopeIndex; ++j)
        {
          if (perIsotopeIntensity[j] <= 0)
          {
            continue;
          }
          fsp << j << "," << perIsotopeIntensity[j] << ";";
        }
        fsp << "\t" << isoScore << "\n";
        fsp << std::setprecision(0);
      }
    }

    //std::cout << "**" << tmp[0] << tmp[1] << tmp[2] << tmp[3] << std::endl;
  }

  void MassFeatureTrace::addDeconvolutedSpectrum(DeconvolutedSpectrum &deconvolutedSpectrum)
  {
    if (deconvolutedSpectrum.getOriginalSpectrum().getMSLevel() != 1)
    {
      return;
    }
    double rt = deconvolutedSpectrum.getOriginalSpectrum().getRT();
    peakGroupMap[rt] = std::unordered_map<double, PeakGroup>();
    auto &subMap = peakGroupMap[rt];
    for (auto &pg : deconvolutedSpectrum)
    {
      subMap[pg.monoisotopicMass] = pg;
    }
  }

  void MassFeatureTrace::writeHeader(std::fstream &fs)
  {
    fs << "FeatureIndex\tFileName\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime"
          "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
          "\tSumIntensity\tMaxIntensity\tFeatureArea\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tPerChargeIntensity\tPerIsotopeIntensity"
          "\n";
  }


  void MassFeatureTrace::writePromexHeader(std::fstream &fs)
  {
    fs << "FeatureID\tMinScan\tMaxScan\tMinCharge\tMaxCharge\t"
          "MonoMass\tRepScan\tRepCharge\tRepMz\tAbundance\tApexScanNum\tApexIntensity\tMinElutionTime\tMaxElutionTime\t"
          "ElutionLength\tEnvelope\tLikelihoodRatio"
          "\n";
  }

  void MassFeatureTrace::updateMembers_()
  {
    tol = param_.getValue("mass_error_ppm");
    //minChargeCosine = param_.getValue("min_charge_cosine");
    minIsotopeCosine = param_.getValue("min_isotope_cosine");
  }
}