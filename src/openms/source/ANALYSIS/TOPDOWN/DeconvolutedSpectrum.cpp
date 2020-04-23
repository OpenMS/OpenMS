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
//
// Created by Kyowon Jeong on 4/22/20.
//

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h"

namespace OpenMS
{
  DeconvolutedSpectrum::DeconvolutedSpectrum()
  {
  }

  DeconvolutedSpectrum::DeconvolutedSpectrum(MSSpectrum &s) :
      spec(&s)
  {
  }

  DeconvolutedSpectrum::~DeconvolutedSpectrum()
  {
    std::vector<LogMzPeak>().swap(peaks);
    std::vector<PeakGroup>().swap(peakGroups);
    //std::unordered_map<int, PeakGroup>().swap(peakGroupMap);
  }

  bool DeconvolutedSpectrum::empty() const
  {
    return peakGroups.empty();
  }

  /*void DeconvolutedSpectrum::updatePeakGroupMap()
  {
    if (!peakGroupMap.empty())
    {
      return;
    }
    //std::cout<< peaks.size()<<std::endl;

    for (auto &pg : peakGroups)
    {
      for (auto &p : pg.peaks)
      {
        auto position = p.index;
        if (peakGroupMap.find(position) == peakGroupMap.end())
        {
          peakGroupMap[position] = pg;
        }
        else
        {
          if (peakGroupMap[position].isotopeCosineScore < pg.isotopeCosineScore)
          {
            peakGroupMap[position] = pg;
          }
          continue;
        }
      }
    }
    //std::cout<<"2"<<std::endl;
  }*/


  void DeconvolutedSpectrum::writeDeconvolutedMasses(std::fstream &fs,
                                                     FLASHDeconvHelperStructs::Parameter &param)//, fstream &fsm, fstream &fsp)
  {
    if (empty())
    {
      return;
    }

    for (auto &pg : peakGroups)
    {
      if (pg.peaks.empty())
      {
        continue;
      }
      double &m = pg.monoisotopicMass;
      double &am = pg.avgMass;
      double &intensity = pg.intensity;
      int minCharge = param.chargeRange + param.minCharge;
      int maxCharge = -1;
      for (auto &p : pg.peaks)
      {
        minCharge = minCharge < p.charge ? minCharge : p.charge;
        maxCharge = maxCharge > p.charge ? maxCharge : p.charge;
      }

      fs << pg.massIndex << "\t" << specIndex << "\t" << param.fileName << "\t" << spec->getNativeID() << "\t"
         << spec->getMSLevel() << "\t"
         << massCntr << "\t"
         << std::to_string(am) << "\t" << std::to_string(m) << "\t" << intensity << "\t"
         << (maxCharge - minCharge + 1) << "\t" << minCharge << "\t" << maxCharge << "\t"
         << std::to_string(spec->getRT())
         << "\t" << pg.peaks.size() << "\t";
      fs << std::fixed << std::setprecision(2);
      fs << pg.maxSNRcharge << "\t" << pg.maxSNR << "\t" << pg.maxSNRminMz << "\t" << pg.maxSNRmaxMz << "\t";
      fs << std::fixed << std::setprecision(-1);
      if (spec->getMSLevel() > 1)
      {
        if (precursorPeakGroup == nullptr){
          fs << "N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\t";
        }else{
          fs << precursorPeakGroup->deconvSpec->specIndex << "\t" << precursorPeak->mz << "\t" << precursorPeak->charge
              << "\t" << precursorPeakGroup->perChargeSNR[precursorPeak->charge]  << "\t" << precursorPeak->intensity
              << "\t" << precursorPeakGroup->monoisotopicMass << "\t" << precursorPeakGroup->totalSNR << "\t" << precursorPeakGroup->isotopeCosineScore
              << "\t" << precursorPeakGroup->chargeCosineScore << "\t" <<precursorPeakGroup->intensity << "\t";
        }
      }

      if (param.writeDetail)
      {
        fs << std::fixed << std::setprecision(2);
        for (auto &p : pg.peaks)
        {
          fs << p.mz << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.charge << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.getUnchargedMass() << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.isotopeIndex << ";";
        }
        fs << "\t";

        for (auto &p : pg.peaks)
        {
          auto tm = pg.monoisotopicMass + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          auto diff = (tm / p.charge + Constants::PROTON_MASS_U - p.mz) / p.mz;

          fs << 1e6 * diff << ";";
        }
        fs << "\t";

        fs << std::fixed << std::setprecision(1);
        for (auto &p : pg.peaks)
        {
          fs << p.intensity << ";";
        }
        fs << "\t";
      }

      fs << std::fixed << std::setprecision(3);
      fs << pg.isotopeCosineScore;

      if (spec->getMSLevel() == 1)
      {
        fs << "\t" << pg.chargeCosineScore;
      }

      fs << "\t" << pg.totalSNR << "\n" << std::setprecision(-1); // TODO
    }

  }


  void DeconvolutedSpectrum::writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail)
  {
    if (detail)
    {
      if (n == 1)
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRChargeMzStart\tMaxSNRChargeMzEnd\t"
               "PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tChargeIntensityCosine\tMassSNR\n";
      }
      else
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRChargeMzStart\tMaxSNRChargeMzEnd\t"
               "PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorChargeSNR\tPrecursorIntensity\t"
               "PrecursorMonoMass\tPrecursorMassSNR\tPrecursorIsotopeCosine\tPrecursorChargeIntensityCosine\tPrecursorMassIntensity\t"
               "PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tMassSNR\n";
      }
    }
    else
    {
      if (n == 1)
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRChargeMzStart\tMaxSNRChargeMzEnd\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tChargeIntensityCosine\tMassSNR\n";
      }
      else
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRChargeMzStart\tMaxSNRChargeMzEnd\t"
               "PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorChargeSNR\tPrecursorIntensity\t"
               "PrecursorMonoMass\tPrecursorMassSNR\tPrecursorIsotopeCosine\tPrecursorChargeIntensityCosine\tPrecursorMassIntensity\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tMassSNR\n";
      }

    }
    //pg.maxSNRcharge << "\t" << pg.maxSNR << "\t" << pg.maxSNRminMz << "\t" << pg.maxSNRmaxMz
    //MinScan	MaxScan	 RepScan	RepCharge	RepMz ApexScanNum    Envelope

  }

  void DeconvolutedSpectrum::writeAttCsvHeader(std::fstream &fs){
    fs<<"ScanNumber,Charge,ChargeSNR,PeakIntensity,"
        "MonoMass,MassSNR,IsotopeCosine,ChargeIntensityCosine,MassIntensity,Class\n";
  }

  void DeconvolutedSpectrum::writeAttCsv(std::fstream &fs, int msLevel){
    if (msLevel>1)
    {
      fs << scanNumber << "," << precursorPeak->charge
         << "," << precursorPeakGroup->perChargeSNR[precursorPeak->charge] << "," << precursorPeak->intensity
         << "," << precursorPeakGroup->monoisotopicMass << "," << precursorPeakGroup->totalSNR << ","
         << precursorPeakGroup->isotopeCosineScore
         << "," << precursorPeakGroup->chargeCosineScore << "," << precursorPeakGroup->intensity << "0\n";
    }else{
      for (auto &pg : peakGroups)
      {
        LogMzPeak *peak = nullptr;
        for (auto &p : pg.peaks)
        {
          if (p.charge != pg.maxSNRcharge){
            continue;
          }
          if (p.mz > pg.maxSNRminMz && p.mz < pg.maxSNRmaxMz){
            peak = &p;
            break;
          }
        }
        if (peak == nullptr){
          continue;
        }

        fs<<scanNumber << "," << pg.maxSNRcharge
          << "," << pg.maxSNR  << "," << peak->intensity
          << "," << pg.monoisotopicMass << "," << pg.totalSNR << "," << pg.isotopeCosineScore
          << "," << pg.chargeCosineScore << "," <<pg.intensity << "0\n";
      }
    }
  }


  void DeconvolutedSpectrum::writeTopFD(std::fstream &fs, int id)//, fstream &fsm, fstream &fsp)
  {
    if (spec->getMSLevel() == 1 || precursorPeakGroup== nullptr)
    {
      return;
    }
    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << id << "\n"
       << "SCANS=" << scanNumber << "\n"
       << "RETENTION_TIME=" << spec->getRT() << "\n";

      fs << "ACTIVATION=" << activationMethod << "\n";

    fs << "MS_ONE_ID=" << precursorPeakGroup->deconvSpec->specIndex << "\n"
       << "MS_ONE_SCAN=" << precursorPeakGroup->deconvSpec->scanNumber << "\n"
       << "PRECURSOR_MZ="
       << std::to_string(precursorPeak->mz) << "\n"
       << "PRECURSOR_CHARGE=" << precursorPeak->charge << "\n"
       << "PRECURSOR_MASS=" << std::to_string(precursorPeakGroup->monoisotopicMass) << "\n"
       << "PRECURSOR_INTENSITY=" << precursorPeak->intensity << "\n";
    fs << std::setprecision(-1);

    for (auto &pg : peakGroups)
    {
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.monoisotopicMass) << "\t" << pg.intensity << "\t" << pg.maxSNRcharge
         //  << "\t" << log10(pg.precursorSNR+1e-10) << "\t" << log10(pg.precursorTotalSNR+1e-10)
         //  << "\t" << log10(pg.maxSNR + 1e-10) << "\t" << log10(pg.totalSNR + 1e-10)
         << "\n";
      fs << std::setprecision(-1);
    }

    fs << "END IONS\n\n";
  }

  bool DeconvolutedSpectrum::registerPrecursor(DeconvolutedSpectrum &precursorSpectrum)
  {
    //precursorSpectrum.updatePeakGroupMap();
    for (auto &p: spec->getPrecursors())
    {
      for (auto &act :  p.getActivationMethods())
      {
        activationMethod = Precursor::NamesOfActivationMethodShort[act];
        break;
      }
      auto startMz = p.getIsolationWindowLowerOffset() > 100.0 ?
                     p.getIsolationWindowLowerOffset() :
                     -p.getIsolationWindowLowerOffset() + p.getMZ();
      auto endMz = p.getIsolationWindowUpperOffset() > 100.0 ?
                   p.getIsolationWindowUpperOffset() :
                   p.getIsolationWindowUpperOffset() + p.getMZ();

      auto position =
          std::lower_bound(precursorSpectrum.peaks.begin(), precursorSpectrum.peaks.end(), LogMzPeak(p.getMZ())) -
          precursorSpectrum.peaks.begin();
      //std::cout<<position<< " " << precursorSpectrum.mzs.size() <<std::endl;
      double minDistance = 10000.0;

      for (; position < (int) precursorSpectrum.peaks.size(); position++)
      {
        auto mz = precursorSpectrum.peaks[position].mz;
        if (mz < startMz)
        {
          continue;
        }
        if (mz > endMz)
        {
          break;
        }

        auto distance = abs(p.getMZ() - mz);
        if (distance > minDistance)
        {
          continue;
        }

        minDistance = distance;
        precursorPeak = &(precursorSpectrum.peaks[position]);
        //precursorPeakGroup = &(peakGroupMap[position]);

        if (mz > p.getMZ())
        {
          break;
        }
      }
    }
    if (precursorPeak != nullptr)
    {
      double pmz = precursorPeak->mz;
      for (auto &pg:precursorSpectrum.peakGroups)
      {
        if (pg.peaks[0].mz > pmz || pg.peaks[pg.peaks.size() - 1].mz < pmz)
        {
          continue;
        }
        for (auto &p:pg.peaks)
        {
          if (p.mz < pmz)
          {
            continue;
          }
          if (p.mz > pmz)
          {
            break;
          }
          if (precursorPeakGroup == nullptr)
          {
            precursorPeakGroup = &pg;
            precursorPeak->charge = p.charge;
          }
          else if (pg.isotopeCosineScore > precursorPeakGroup->isotopeCosineScore)
          {
            precursorPeakGroup = &pg;
            precursorPeak->charge = p.charge;
          }
        }
      }
    }
    return precursorPeakGroup != nullptr;

  }


}