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

namespace OpenMS
{

  MassFeatureTrace::MassFeatureTrace(Parameter &p, Param &mp, PrecalculatedAveragine &avg) :
      param(p), mtd_param(mp), averagines(avg)
  {

  }

  MassFeatureTrace::~MassFeatureTrace()
  {
    for (auto &item : peakGroupMap)
    {
      std::unordered_map<double, PeakGroup>().swap(item.second);
    }
    std::unordered_map<double, std::unordered_map<double, PeakGroup>>().swap(peakGroupMap);
  }


  void MassFeatureTrace::findFeatures(int &featureCntr,
                                      int &featureIndex,
                                      std::fstream &fsf,
                                      std::fstream &fsp)
  {
    MSExperiment map;
    std::map<int, MSSpectrum> indexSpecMap;

    for (auto &item : peakGroupMap)
    {
      auto rt = item.first;
      MSSpectrum deconvSpec;
      deconvSpec.setRT(rt);
      for(auto &pg : item.second){
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

    mtd_param.setValue("mass_error_ppm", param.tolerance[0] * 1e6, "");
    mtd_param.setValue("trace_termination_criterion", "outlier", "");

    mtd_param.setValue("reestimate_mt_sd", "false", "");
    mtd_param.setValue("quant_method", "area", "");
    mtd_param.setValue("noise_threshold_int", .0, "");

    //double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;
    mtd_param.setValue("min_sample_rate", 0.01, "");
    mtd_param.setValue("trace_termination_outliers", param.numOverlappedScans, "");
    mtd_param.setValue("min_trace_length", param.minRTSpan, "");
    //mtd_param.setValue("max_trace_length", 1000.0, "");
    mtdet.setParameters(mtd_param);

    std::vector<MassTrace> m_traces;

    mtdet.run(map, m_traces);  // m_traces : output of this function

    auto *perChargeIntensity = new double[param.chargeRange + param.minCharge + 1];
    auto *perChargeMaxIntensity = new double[param.chargeRange + param.minCharge + 1];
    auto *perChargeMz = new double[param.chargeRange + param.minCharge + 1];
    auto *perIsotopeIntensity = new double[param.maxIsotopeCount];

    for (auto &mt : m_traces)
    {
      //if(mt.getSize() < 3){
      //  continue;
      //}
      int minCharge = param.chargeRange + param.minCharge + 1;
      int maxCharge = 0;

      int minScanNum = (int) map.size() + 1000;
      int maxScanNum = 0;

      int repScan = 0, repCharge = 0;
      double maxIntensity = 0;
      double maxMass = .0;
      double maxIso = 0;
      boost::dynamic_bitset<> charges(param.chargeRange + param.minCharge + 1);
      std::fill_n(perChargeIntensity, param.chargeRange + param.minCharge + 1, 0);
      std::fill_n(perChargeMaxIntensity, param.chargeRange + param.minCharge + 1, 0);
      std::fill_n(perChargeMz, param.chargeRange + param.minCharge + 1, 0);
      std::fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);

      for (auto &p2 : mt)
      {
        auto &pgMap = peakGroupMap[p2.getRT()];
        auto &pg = pgMap[p2.getMZ()];
        auto scanNumber = pg.scanNumber;

        minCharge = minCharge < pg.minCharge ? minCharge : pg.minCharge;
        maxCharge = maxCharge > pg.maxCharge ? maxCharge : pg.maxCharge;

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


        for (auto &p : pg.peaks)
        {
          if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount || p.charge < 0 ||
              p.charge >= param.chargeRange + param.minCharge + 1)
          {
            continue;
          }

          charges[p.charge] = true;
          perChargeIntensity[p.charge] += p.intensity;
          perIsotopeIntensity[p.isotopeIndex] += p.intensity;
          if (perChargeMaxIntensity[p.charge] > p.intensity)
          {
            continue;
          }
          perChargeMaxIntensity[p.charge] = p.intensity;
          perChargeMz[p.charge] = p.mz;
        }

      }

      double chargeScore = PeakGroupScoring::getChargeFitScore(perChargeIntensity,
                                                               param.minCharge + param.chargeRange + 1);
      if (chargeScore < param.minChargeCosine) //
      {
        continue;
      }

      int offset = 0;

      double mass = mt.getCentroidMZ();
      double isoScore = PeakGroupScoring::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                                   perIsotopeIntensity,
                                                                                   param.maxIsotopeCount,
                                                                                   offset, averagines);
      if (isoScore < param.minIsotopeCosine[0])
      {
        continue;
      }

      //perIsotopeIntensity, param.maxIsotopeCount


      if (offset != 0)
      {
        mass += offset * Constants::ISOTOPE_MASSDIFF_55K_U;
      }

      auto sumInt = .0;
      for (auto &p : mt)
      {
        sumInt += p.getIntensity();
      }

      auto massDelta = averagines.getAverageMassDelta(mass);
      ++featureCntr;
      fsf << featureIndex++ << "\t" << param.fileName << "\t" << std::to_string(mass) << "\t"
          << std::to_string(mass + massDelta) << "\t" // massdiff
          << mt.getSize() << "\t"
          << mt.begin()->getRT() << "\t"
          << mt.rbegin()->getRT() << "\t"
          << mt.getTraceLength() << "\t"
          << mt[mt.findMaxByIntPeak()].getRT() << "\t"
          << sumInt << "\t"
          << mt.getMaxIntensity(false) << "\t"
          << minCharge << "\t"
          << maxCharge << "\t"
          << charges.count() << "\t"
          << isoScore << "\t"
          << chargeScore << "\n";


      if(param.promexOut){
        double maxChargeIntensity = 0;
        for (int c = 0; c < param.chargeRange + param.minCharge + 1; c++)
        {
          if (perChargeIntensity[c] > maxChargeIntensity)
          {
            maxChargeIntensity = perChargeIntensity[c];
            repCharge = c;
          }
        }
        auto apex = mt[mt.findMaxByIntPeak()];
        //int si = rtSpecMap[(float) apex.getRT()];
        auto &spgMap = peakGroupMap[apex.getRT()];
        auto &spg = spgMap[apex.getMZ()];

        fsp << featureIndex << "\t" << minScanNum << "\t" << maxScanNum << "\t" << minCharge << "\t"
            << maxCharge << "\t" << std::to_string(mass) << "\t" << std::fixed << std::setprecision(2)
            << repScan << "\t" << repCharge << "\t" << perChargeMz[repCharge] << "\t" << sumInt << "\t"
            << spg.scanNumber << "\t" << spg.intensity << "\t"
            << mt.begin()->getRT()/60.0 << "\t"
            << mt.rbegin()->getRT()/60.0 << "\t"
            << mt.getTraceLength()/60.0 << "\t";


        for (int j = 0; j < param.maxIsotopeCount; ++j)
        {
          if (perIsotopeIntensity[j] <= 0)
          {
            continue;
          }
          fsp << j << "," << perIsotopeIntensity[j] << ";";
        }
        fsp << "\t"<<isoScore<<"\n";
        fsp << std::setprecision(0);
      }
    }

    delete[] perIsotopeIntensity;
    delete[] perChargeMz;
    delete[] perChargeMaxIntensity;
    delete[] perChargeIntensity;
//    delete[] peakGroupMap;
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
          "\tSumIntensity\tMaxIntensity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tChargeIntensityCosineScore"
          "\n";
  }


  void MassFeatureTrace::writePromexHeader(std::fstream &fs)
  {
    fs << "FeatureID\tMinScan\tMaxScan\tMinCharge\tMaxCharge\t"
          "MonoMass\tRepScan\tRepCharge\tRepMz\tAbundance\tApexScanNum\tApexIntensity\tMinElutionTime\tMaxElutionTime\t"
          "ElutionLength\tEnvelope\tLikelihoodRatio"
          "\n";
  }

}
