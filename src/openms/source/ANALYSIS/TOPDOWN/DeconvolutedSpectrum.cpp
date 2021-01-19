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
  DeconvolutedSpectrum::DeconvolutedSpectrum(const MSSpectrum& spectrum, const int scan_number) :
      scan_number(scan_number)
  {
    spec = spectrum;
  }

  MSSpectrum DeconvolutedSpectrum::toSpectrum(const int mzml_charge)
  {
    auto outSpec = MSSpectrum(spec);
    outSpec.clear(false);
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      outSpec.emplace_back(pg.getMonoMass(), pg.getIntensity());
    }
    if (precursor_peak_group != nullptr && !spec.getPrecursors().empty())
    {
      Precursor precursor(spec.getPrecursors()[0]);
      precursor.setCharge(precursor_peak_group->getRepCharge());
      precursor.setMZ(precursor_peak_group->getMonoMass() + mzml_charge * (mzml_charge >= 0 ? Constants::PROTON_MASS_U : Constants::ELECTRON_MASS_U));
      precursor.setIntensity(precursor_peak_group->getIntensity());
      outSpec.getPrecursors().clear();
      outSpec.getPrecursors().emplace_back(precursor);
    }
    return outSpec;
  }

  void DeconvolutedSpectrum::writeDeconvolutedMasses(std::fstream& fs,
                                                     const String& file_name,
                                                     const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                                     const bool write_detail)//, fstream& fsm, fstream& fsp)
  {
    if (empty())
    {
      return;
    }

    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      const double m = pg.getMonoMass();
      const double am = pg.getMonoMass() + avg.getAverageMassDelta(m);
      const double intensity = pg.getIntensity();

      auto crange = pg.getChargeRange();

      fs << file_name << "\t" << pg.getScanNumber() << "\t"
         << std::to_string(spec.getRT()) << "\t"
         << size() << "\t"
         << std::to_string(am) << "\t" << std::to_string(m) << "\t" << intensity << "\t"
         << std::get<0>(crange) << "\t" << std::get<1>(crange) << "\t"
         << pg.size() << "\t";

      if (write_detail)
      {
        fs << std::fixed << std::setprecision(2);
        for (auto& p : pg)
        {
          fs << p.mz << ";";
        }
        fs << "\t";

        fs << std::fixed << std::setprecision(1);
        for (auto& p : pg)
        {
          fs << p.intensity << ";";
        }
        fs << "\t";
        fs << std::setprecision(-1);

        for (auto& p : pg)
        {
          fs << p.charge << ";";
        }
        fs << "\t";
        for (auto& p : pg)
        {
          fs << p.getUnchargedMass() << ";";
        }
        fs << "\t";
        for (auto& p : pg)
        {
          fs << p.isotopeIndex << ";";
        }
        fs << "\t";

        for (auto& p : pg)
        {
          double tm = pg.getMonoMass() + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          double diff = (tm / abs(p.charge) + FLASHDeconvHelperStructs::getChargeMass(p.charge > 0) - p.mz) / p.mz;
          fs << 1e6 * diff << ";";
        }
        fs << "\t";
      }
      if (spec.getMSLevel() > 1)
      {
        //PrecursorScanNum	PrecursorMz	PrecursorIntensity PrecursorCharge	PrecursorMonoMass		PrecursorQScore
        fs << precursor_scan_number << "\t" << std::to_string(precursor_peak.getMZ()) << "\t"
           << precursor_peak.getIntensity() << "\t"
           << precursor_peak.getCharge()
           << "\t";

        if (precursor_peak_group == nullptr)
        {
          fs << "nan\tnan\t";
        }
        else
        {
          fs << std::to_string(precursor_peak_group->getMonoMass()) << "\t"
             << precursor_peak_group->getQScore() << "\t";
        }
      }
      fs << pg.getIsotopeCosine() << "\t" << pg.getChargeScore() <<"\t";

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

      for (auto& p : pg)
      {
        isoEndIndex = isoEndIndex < p.isotopeIndex ? p.isotopeIndex : isoEndIndex;
      }
      auto perIsotopeIntensity = std::vector<double>(isoEndIndex + 1, .0);
      for (auto& p : pg)
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


  void DeconvolutedSpectrum::writeDeconvolutedMassesHeader(std::fstream& fs, const int ms_level, const bool detail)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (ms_level == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";

      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
  }

  /*  void DeconvolutedSpectrum::writeAttCsvHeader(std::fstream& fs)
    {
      fs
          << "ScanNumber,RetentionTime,PrecursorScanNumber,Charge,ChargeSNR,PeakIntensity,EnvIntensity,EnvIsotopeCosine,PeakMz,"
             "MonoMass,MassSNR,IsotopeCosine,MassIntensity,QScore,Class\n";
    }*/

  void DeconvolutedSpectrum::writeThermoInclusionHeader(std::fstream& fs)
  {
    fs << "Compound,Formula,Adduct,m/z,z,t start (min),t stop (min),Isolation Window (m/z),Normalized AGC Target (%)\n";
  }


  void DeconvolutedSpectrum::clearPeakGroupsChargeInfo()
  {
    for (auto& pg: *this)
    {
      pg.clearChargeInfo();
    }
  }

  void DeconvolutedSpectrum::writeTopFD(std::fstream& fs, const int id, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg)//, fstream& fsm, fstream& fsp)
  {
    UInt msLevel = spec.getMSLevel();

    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << id << "\n"
       << "SCANS=" << scan_number << "\n"
       << "RETENTION_TIME=" << spec.getRT() << "\n";

    fs << "ACTIVATION=" << activation_method << "\n";

    if (msLevel > 1)
    {
      if (precursor_peak_group != nullptr)
      {
        fs << "MS_ONE_ID=" << precursor_scan_number << "\n"
           << "MS_ONE_SCAN=" << precursor_scan_number << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursor_peak.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << precursor_peak.getCharge() << "\n"
           << "PRECURSOR_MASS=" << std::to_string(precursor_peak_group->getMonoMass()) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak.getIntensity() << "\n";
      }
      else
      {
        double avgMass = (precursor_peak.getMZ() - FLASHDeconvHelperStructs::getChargeMass(precursor_peak.getCharge() > 0)) * abs(precursor_peak.getCharge());
        double mMass = avgMass - avg.getAverageMassDelta(avgMass);
        fs << "MS_ONE_ID=" << 0 << "\n"
           << "MS_ONE_SCAN=" << precursor_scan_number << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursor_peak.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << precursor_peak.getCharge() << "\n"
           << "PRECURSOR_MASS=" << std::to_string(mMass) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak.getIntensity() << "\n";
      }
    }
    fs << std::setprecision(-1);

    double scoreThreshold = 0;
    std::vector<double> scores;

    if (size() > 500)// max peak count for TopPic
    {
      scores.reserve(size());
      for (auto& pg : *this)
      {
        scores.push_back(pg.getIsotopeCosine());
      }
      std::sort(scores.begin(), scores.end());
      scoreThreshold = scores[scores.size() - 500];
      std::vector<double>().swap(scores);
    }

    int size = 0;
    for (auto& pg : *this)
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
                                                                                //  << "\t" << log10(pg.maxSNR + 1e-10) << "\t" << log10(pg.total_snr + 1e-10)
                                                                                << "\n";
      fs << std::setprecision(-1);
    }

    fs << "END IONS\n\n";
  }

  bool DeconvolutedSpectrum::registerPrecursor(DeconvolutedSpectrum& precursor_spectrum)
  {
    //precursor_spectrum.updatePeakGroupMap();
    for (auto& p: spec.getPrecursors())
    {
      for (auto& act :  p.getActivationMethods())
      {
        activation_method = Precursor::NamesOfActivationMethodShort[act];
        break;
      }
      precursor_peak = p;
      precursor_scan_number = precursor_spectrum.scan_number;
      double startMz = p.getIsolationWindowLowerOffset() > 100.0 ?
                     p.getIsolationWindowLowerOffset() :
                     -p.getIsolationWindowLowerOffset() + p.getMZ();
      double endMz = p.getIsolationWindowUpperOffset() > 100.0 ?
                   p.getIsolationWindowUpperOffset() :
                   p.getIsolationWindowUpperOffset() + p.getMZ();

      double maxSumIntensity = 0.0;
      for (auto& pg: precursor_spectrum)
      {
        std::sort(pg.begin(), pg.end());
        if (pg[0].mz > endMz || pg[pg.size() - 1].mz < startMz)
        {
          continue;
        }

        double sumIntensity = .0;
        double maxIntensity = .0;
        LogMzPeak *tmp = nullptr;
        for (auto& pt:pg)
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

          tmp = &pt;
        }

        if (sumIntensity <= maxSumIntensity || tmp == nullptr)
        {
          continue;
        }

        precursor_peak.setMZ(tmp->mz);
        precursor_peak.setIntensity(tmp->intensity);
        precursor_peak.setCharge(tmp->charge);
        maxSumIntensity = sumIntensity;
        precursor_peak_group = &pg;
      }
      if(precursor_peak_group != nullptr){
        break;
      }
    }

    return precursor_peak_group != nullptr;
  }

  MSSpectrum& DeconvolutedSpectrum::getOriginalSpectrum()
  {
    return spec;
  }

  PeakGroup DeconvolutedSpectrum::getPrecursorPeakGroup()
  {
    if (precursor_peak_group == nullptr)
    {
      return PeakGroup();
    }
    return *precursor_peak_group;
  }

  int DeconvolutedSpectrum::getPrecursorCharge()
  {
    return precursor_peak.getCharge();
  }

  double DeconvolutedSpectrum::getCurrentMaxMass(const double max_mass)
  {
    if (spec.getMSLevel() == 1 || precursor_peak_group == nullptr || precursor_peak_group->empty())
    {
      return max_mass;
    }
    return precursor_peak_group->getMonoMass();
  }

  int DeconvolutedSpectrum::getCurrentMaxCharge(const int max_charge)
  {
    if (spec.getMSLevel() == 1 || precursor_peak_group == nullptr || precursor_peak_group->empty())
    {
      return max_charge;
    }
    return precursor_peak.getCharge();
  }

}