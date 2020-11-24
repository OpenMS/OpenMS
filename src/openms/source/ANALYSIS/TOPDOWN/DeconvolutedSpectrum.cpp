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

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h"

namespace OpenMS
{
  DeconvolutedSpectrum::DeconvolutedSpectrum(const MSSpectrum &s, int n) :
      scanNumber(n)
  {
    spec = s;
  }

  MSSpectrum DeconvolutedSpectrum::toSpectrum()
  {
    auto outSpec = MSSpectrum(spec);
    outSpec.clear(false);
    for (auto &pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      outSpec.emplace_back(pg.getMonoMass(), pg.getIntensity());
    }
    if (precursorPeakGroup != nullptr && !spec.getPrecursors().empty())
    {
      Precursor precursor(spec.getPrecursors()[0]);
      precursor.setCharge(precursorPeakGroup->getRepCharge());
      precursor.setMZ(precursorPeakGroup->getMonoMass());
      precursor.setIntensity(precursorPeakGroup->getIntensity());
      outSpec.getPrecursors().clear();
      outSpec.getPrecursors().emplace_back(precursor);
    }
    return outSpec;
  }

  void DeconvolutedSpectrum::writeDeconvolutedMasses(std::fstream &fs,
                                                     const String &fileName,
                                                     const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                                     bool writeDetail)//, fstream &fsm, fstream &fsp)
  {
    if (empty())
    {
      return;
    }

    for (auto &pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      const double m = pg.getMonoMass();
      const double am = pg.getMonoMass() + avg.getAverageMassDelta(m);
      const double intensity = pg.getIntensity();

      auto crange = pg.getChargeRange();

      fs << "\t" << fileName << "\t" << pg.getScanNumber() << "\t"
         << std::to_string(spec.getRT()) << "\t"
         << size() << "\t"
         << std::to_string(am) << "\t" << std::to_string(m) << "\t" << intensity << "\t"
         << std::get<0>(crange) << "\t" << std::get<1>(crange) << "\t"
         << pg.size() << "\t";

      if (writeDetail)
      {
        fs << std::fixed << std::setprecision(2);
        for (auto &p : pg)
        {
          fs << p.mz << ";";
        }
        fs << "\t";

        fs << std::fixed << std::setprecision(1);
        for (auto &p : pg)
        {
          fs << p.intensity << ";";
        }
        fs << "\t";
        fs << std::setprecision(-1);

        for (auto &p : pg)
        {
          fs << p.charge << ";";
        }
        fs << "\t";
        for (auto &p : pg)
        {
          fs << p.getUnchargedMass() << ";";
        }
        fs << "\t";
        for (auto &p : pg)
        {
          fs << p.isotopeIndex << ";";
        }
        fs << "\t";

        for (auto &p : pg)
        {
          auto tm = pg.getMonoMass() + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          auto diff = (tm / abs(p.charge) + FLASHDeconvHelperStructs::getChargeMass(p.charge > 0) - p.mz) / p.mz;
          fs << 1e6 * diff << ";";
        }
        fs << "\t";
      }
      if (spec.getMSLevel() > 1)
      {
        //PrecursorScanNum	PrecursorMz	PrecursorIntensity PrecursorCharge	PrecursorMonoMass		PrecursorQScore
        fs << precursorScanNumber << "\t" << std::to_string(precursorPeak.getMZ()) << "\t"
           << precursorPeak.getIntensity() << "\t"
           << precursorPeak.getCharge()
           << "\t";

        if (precursorPeakGroup == nullptr)
        {
          fs << "N/A\tN/A\t";
        }
        else
        {
          fs << std::to_string(precursorPeakGroup->getMonoMass()) << "\t"
             << precursorPeakGroup->getQScore() << "\t";
        }
      }
      fs << pg.getIsotopeCosine() << "\t";

      auto qrange = pg.getMzxQScoreMzRange();
      fs << pg.getSNR() << "\t"
         << pg.getRepCharge() << "\t" << std::to_string(std::get<0>(qrange)) << "\t"
         << std::to_string(std::get<1>(qrange)) << "\t"
         << pg.getQScore() << "\t" << std::setprecision(-1); //

      for (int i = std::get<0>(crange); i <= std::get<1>(crange); i++)
      {

        fs << pg.getChargeIntensity(i);

        if (i < std::get<1>(crange))
        {
          fs << ";";
        }
      }
      fs << "\t";
      int isoEndIndex = 0;

      for (auto &p : pg)
      {
        isoEndIndex = isoEndIndex < p.isotopeIndex ? p.isotopeIndex : isoEndIndex;
      }
      auto perIsotopeIntensity = new double[isoEndIndex + 1];
      std::fill_n(perIsotopeIntensity, isoEndIndex + 1, .0);
      for (auto &p : pg)
      {
        perIsotopeIntensity[p.isotopeIndex] += p.intensity;
      }

      for (int i = 0; i <= isoEndIndex; i++)
      {
        fs << perIsotopeIntensity[i];
        if (i < isoEndIndex)
        {
          fs << ";";
        }
      }

      fs << "\n";
    }
  }


  void DeconvolutedSpectrum::writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail)
  {
    if (detail)
    {
      if (n == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "IsotopeCosine\tMassSNR\tRepresentativeCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tMassSNR\tRepresentativeCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (n == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tMassSNR\tRepresentativeCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";

      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tMassSNR\tRepresentativeCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
  }

  /*  void DeconvolutedSpectrum::writeAttCsvHeader(std::fstream &fs)
    {
      fs
          << "ScanNumber,RetentionTime,PrecursorScanNumber,Charge,ChargeSNR,PeakIntensity,EnvIntensity,EnvIsotopeCosine,PeakMz,"
             "MonoMass,MassSNR,IsotopeCosine,MassIntensity,QScore,Class\n";
    }*/

  void DeconvolutedSpectrum::writeThermoInclusionHeader(std::fstream &fs)
  {
    fs << "Compound,Formula,Adduct,m/z,z,t start (min),t stop (min),Isolation Window (m/z),Normalized AGC Target (%)\n";
  }


  void DeconvolutedSpectrum::clearPeakGroupsChargeInfo()
  {
    for (auto &pg: *this)
    {
      pg.clearChargeInfo();
    }
  }

  void DeconvolutedSpectrum::writeTopFD(std::fstream &fs, int id)//, fstream &fsm, fstream &fsp)
  {
    auto msLevel = spec.getMSLevel();

    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << id << "\n"
       << "SCANS=" << scanNumber << "\n"
       << "RETENTION_TIME=" << spec.getRT() << "\n";

    fs << "ACTIVATION=" << activationMethod << "\n";

    if (msLevel > 1)
    {
      if (precursorPeakGroup != nullptr)
      {
        fs << "MS_ONE_ID=" << precursorPeakGroup->getScanNumber() << "\n"
           << "MS_ONE_SCAN=" << precursorPeakGroup->getScanNumber() << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursorPeak.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << precursorPeak.getCharge() << "\n"
           << "PRECURSOR_MASS=" << std::to_string(precursorPeakGroup->getMonoMass()) << "\n"
           << "PRECURSOR_INTENSITY=" << precursorPeak.getIntensity() << "\n";
      }
      else
      {
        fs << "MS_ONE_ID=" << 0 << "\n"
            << "MS_ONE_SCAN=" << 0 << "\n"
            << "PRECURSOR_MZ="
            << std::to_string(precursorPeak.getMZ()) << "\n"
            << "PRECURSOR_CHARGE=" << precursorPeak.getCharge() << "\n"
            << "PRECURSOR_MASS=" << std::to_string(precursorPeak.getMZ() * precursorPeak.getCharge()) << "\n"
            << "PRECURSOR_INTENSITY=" << precursorPeak.getIntensity() << "\n";
      }
    }
    fs << std::setprecision(-1);

    double scoreThreshold = 0;
    std::vector<double> scores;

    if (size() > 500)// max peak count for TopPic
    {
      scores.reserve(size());
      for (auto &pg : *this)
      {
        scores.push_back(pg.getIsotopeCosine());
      }
      std::sort(scores.begin(), scores.end());
      scoreThreshold = scores[scores.size() - 500];
      std::vector<double>().swap(scores);
    }

    int size = 0;
    for (auto &pg : *this)
    {
      if (pg.getIsotopeCosine() < scoreThreshold)
      {
        continue;
      }
      if (size >= 500)
      {
        break;
      }

      size++;
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.getMonoMass()) << "\t" << pg.getIntensity() << "\t" << pg.getRepCharge()
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
    for (auto &p: spec.getPrecursors())
    {
      for (auto &act :  p.getActivationMethods())
      {
        activationMethod = Precursor::NamesOfActivationMethodShort[act];
        break;
      }
      precursorPeak = p;
      precursorScanNumber = precursorSpectrum.scanNumber;
      auto startMz = p.getIsolationWindowLowerOffset() > 100.0 ?
                     p.getIsolationWindowLowerOffset() :
                     -p.getIsolationWindowLowerOffset() + p.getMZ();
      auto endMz = p.getIsolationWindowUpperOffset() > 100.0 ?
                   p.getIsolationWindowUpperOffset() :
                   p.getIsolationWindowUpperOffset() + p.getMZ();

      double maxSumIntensity = 0.0;
      for (auto &pg: precursorSpectrum)
      {
        if (pg[0].mz > endMz || pg[pg.size() - 1].mz < startMz)
        {
          continue;
        }

        double sumIntensity = .0;
        double maxIntensity = .0;
        LogMzPeak *tmp = nullptr;
        for (auto &pt:pg)
        {
          if (pt.mz < startMz)
          {
            continue;
          }
          if (pt.mz > endMz)
          {
            break;
          }
          sumIntensity += pt.intensity;

          if (pt.intensity < maxIntensity)
          {
            continue;
          }
          maxIntensity = pt.intensity;
          precursorPeak.setCharge(pt.charge);
          tmp = &pt;
        }

        if (sumIntensity <= maxSumIntensity || tmp == nullptr)
        {
          continue;
        }

        maxSumIntensity = sumIntensity;
        precursorPeakGroup = &pg;
      }
    }

    return precursorPeakGroup != nullptr;
  }

  MSSpectrum &DeconvolutedSpectrum::getOriginalSpectrum()
  {
    return spec;
  }

  PeakGroup DeconvolutedSpectrum::getPrecursorPeakGroup()
  {
    if (precursorPeakGroup == nullptr)
    {
      return PeakGroup();
    }
    return *precursorPeakGroup;
  }

  int DeconvolutedSpectrum::getPrecursorCharge()
  {
    return precursorPeak.getCharge();
  }

  double DeconvolutedSpectrum::getCurrentMaxMass(double maxMass)
  {
    if (spec.getMSLevel() == 1 || precursorPeakGroup == nullptr || precursorPeakGroup->empty())
    {
      return maxMass;
    }
    return precursorPeakGroup->getMonoMass();
  }

  int DeconvolutedSpectrum::getCurrentMaxCharge(int maxCharge)
  {
    if (spec.getMSLevel() == 1 || precursorPeakGroup == nullptr || precursorPeakGroup->empty())
    {
      return maxCharge;
    }
    return precursorPeak.getCharge();
  }

}