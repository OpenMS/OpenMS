// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <random>

namespace OpenMS
{
  /**
   * @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO

   */

  void FLASHDeconvSpectrumFile::writeDeconvolvedMasses(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& file_name,
                                                       const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double tol,
                                                      const bool write_detail, const bool decoy)
  {
    static std::vector<uint> indices{};

    if (dspec.empty())
    {
      return;
    }

    while(indices.size() <= dspec.getOriginalSpectrum().getMSLevel())
    {
      indices.push_back(1);
    }
    uint& index = indices[dspec.getOriginalSpectrum().getMSLevel()-1];

    for (auto& pg : dspec)
    {
      const double mono_mass = pg.getMonoMass();
      const double avg_mass = pg.getMonoMass() + avg.getAverageMassDelta(mono_mass);
      const double intensity = pg.getIntensity();

      auto charge_range = pg.getAbsChargeRange();
      int min_charge = pg.isPositive() ? std::get<0>(charge_range) : -std::get<1>(charge_range);
      int max_charge = pg.isPositive() ? std::get<1>(charge_range) : -std::get<0>(charge_range);

      pg.setIndex(index);
      fs << index++ << "\t" << file_name << "\t" << pg.getScanNumber() << "\t";
      if (decoy)
      {
        fs << pg.getDecoyFlag() << "\t";
      }
      fs << std::to_string(dspec.getOriginalSpectrum().getRT()) << "\t" << dspec.size() << "\t" << std::to_string(avg_mass) << "\t" << std::to_string(mono_mass) << "\t" << intensity << "\t"
         << min_charge << "\t" << max_charge << "\t" << pg.size() << "\t";

      if (write_detail)
      {
        std::unordered_set<int> excluded_integer_mzs;
        auto noisy_peaks = pg.recruitAllPeaksInSpectrum(dspec.getOriginalSpectrum(), tol * 1e-6, avg, pg.getMonoMass(), excluded_integer_mzs);

        fs << std::fixed << std::setprecision(2);
        for (auto& p : pg)
        {
          fs << p.mz << " ";
        }

        fs << "\t";
        fs << std::fixed << std::setprecision(1);
        for (auto& p : pg)
        {
          fs << p.intensity << " ";
        }

        fs << "\t";
        fs << std::setprecision(-1);
        for (auto& p : pg)
        {
          fs << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }

        fs << "\t";
        for (auto& p : pg)
        {
          fs << p.getUnchargedMass() << " ";
        }

        fs << "\t";
        for (auto& p : pg)
        {
          fs << p.isotopeIndex << " ";
        }

        fs << "\t";
        for (auto& p : pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) - p.mz) / p.mz;
          fs << 1e6 * mass_error << " ";
        }

        fs << "\t";
        fs << std::fixed << std::setprecision(2);
        for (auto& np: noisy_peaks)
        {
          fs << np.mz << " ";
        }

        fs << "\t";
        fs << std::fixed << std::setprecision(1);
        for (auto& np: noisy_peaks)
        {
          fs << np.intensity << " ";
        }

        fs << "\t";
        fs << std::setprecision(-1);
        for (auto& np: noisy_peaks)
        {
          fs << (np.is_positive ? np.abs_charge : -np.abs_charge) << " ";
        }

        fs << "\t";
        for (auto& np: noisy_peaks)
        {
          fs << np.getUnchargedMass() << " ";
        }

        fs << "\t";
        for (auto& np: noisy_peaks)
        {
          fs << np.isotopeIndex << " ";
        }

        fs << "\t";
        for (auto& np: noisy_peaks)
        {
          double average_mass = pg.getMonoMass() + np.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / np.abs_charge + FLASHDeconvHelperStructs::getChargeMass(np.is_positive) - np.mz) / np.mz;
          fs << 1e6 * mass_error << " ";
        }
        fs << "\t";
      }
      if (dspec.getOriginalSpectrum().getMSLevel() > 1)
      {
        // PrecursorScanNum	PrecursorMz	PrecursorIntensity PrecursorCharge	PrecursorMonoMass		PrecursorQScore
        fs << dspec.getPrecursorScanNumber() << "\t" << std::to_string(dspec.getPrecursor().getMZ()) << "\t" << dspec.getPrecursor().getIntensity() << "\t" << dspec.getPrecursor().getCharge() << "\t";

        if (dspec.getPrecursorPeakGroup().empty())
        {
          fs << "nan\tnan\tnan\t";
          if (decoy)
            fs << "nan\tnan\tnan\tnan\t";
        }
        else
        {
          fs << dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) << "\t" << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) << "\t"
             << dspec.getPrecursorPeakGroup().getQScore() << "\t";
          if (decoy)
            fs << dspec.getPrecursorPeakGroup().getQvalue() << "\t" << dspec.getPrecursorPeakGroup().getQvalue(PeakGroup::DecoyFlag::isotope_decoy) << "\t"
               << dspec.getPrecursorPeakGroup().getQvalue(PeakGroup::DecoyFlag::noise_decoy) << "\t" << dspec.getPrecursorPeakGroup().getQvalue(PeakGroup::DecoyFlag::charge_decoy) << "\t";
        }
      }
      fs << pg.getIsotopeCosine() << "\t" << pg.getChargeIsotopeCosine(pg.getRepAbsCharge()) << "\t" << pg.getChargeScore() << "\t";

      auto max_qscore_mz_range = pg.getRepMzRange();
      fs << pg.getSNR() << "\t" << pg.getChargeSNR(pg.getRepAbsCharge()) << "\t"<< pg.getAvgPPMError() << "\t" << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t" << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t" << pg.getQScore();

      if (decoy)
      {
        fs << "\t" << pg.getQvalue() << "\t" << pg.getQvalue(PeakGroup::DecoyFlag::isotope_decoy) << "\t"
           << pg.getQvalue(PeakGroup::DecoyFlag::noise_decoy) << "\t" << pg.getQvalue(PeakGroup::DecoyFlag::charge_decoy);
      }

      if (write_detail)
      {
        fs << "\t" << std::setprecision(-1);

        for (int i = std::get<0>(charge_range); i <= std::get<1>(charge_range); i++)
        {
          fs << pg.getChargeIntensity(i);

          if (i < std::get<1>(charge_range))
          {
            fs << ";";
          }
        }
        fs << "\t";

        auto iso_intensities = pg.getIsotopeIntensities();
        for (size_t i = 0; i < iso_intensities.size(); i++)
        {
          fs << iso_intensities[i];
          if (i < iso_intensities.size() - 1)
          {
            fs << ";";
          }
        }
      }
      fs << "\n";
    }
  }

  void FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(std::fstream& fs, const uint ms_level, const bool detail, const bool decoy)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\t";
        if (decoy)
        {
          fs << "Decoy\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\t";
        if (decoy)
        {
          fs << "Qvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly\t";
        }
        fs << "PerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\t";
        if (decoy)
        {
          fs << "Decoy\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQScore\t";
        if (decoy)
        {
          fs << "PrecursorQvalue\tPrecursorQvalueWithIsotopeDecoyOnly\tPrecursorQvalueWithNoiseDecoyOnly\tPrecursorQvalueWithChargeDecoyOnly\t";
        }
        fs << "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\t";
        if (decoy)
        {
          fs << "Qvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly\t";
        }
        fs << "PerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\t";
        if (decoy)
        {
          fs << "Decoy\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\t";
        if (decoy)
        {
          fs << "Qvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly";
        }
        fs << "\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\t";
        if (decoy)
        {
          fs << "Decoy\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQScore\t";
        if (decoy)
        {
          fs << "PrecursorQvalue\tPrecursorQvalueWithIsotopeDecoyOnly\tPrecursorQvalueWithNoiseDecoyOnly\tPrecursorQvalueWithChargeDecoyOnly\t";
        }
        fs << "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\t";
        if (decoy)
        {
          fs << "Qvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly";
        }
         fs << "\n";
      }
    }
  }

  void FLASHDeconvSpectrumFile::writeTopFD(const DeconvolvedSpectrum& dspec, std::fstream& fs,
                                           const double snr_threshold, const double decoy_harmonic_factor,
                                           const double decoy_precursor_offset) 
  {
    UInt ms_level = dspec.getOriginalSpectrum().getMSLevel();

    if (ms_level > 1)
    {
      if (dspec.getPrecursorPeakGroup().empty() || dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) < snr_threshold)
      {
        return;
      }
    }

    if (dspec.size() < topFD_min_peak_count_)
    {
      return;
    }

    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << dspec.getScanNumber() << "\n"
       << "FRACTION_ID=" << 0 << "\n"
       << "SCANS=" << dspec.getScanNumber() << "\n"
       << "RETENTION_TIME=" << dspec.getOriginalSpectrum().getRT() << "\n"
       << "LEVEL=" << dspec.getOriginalSpectrum().getMSLevel() << "\n";

    if (ms_level > 1)
    {
      fs << "ACTIVATION=" << Precursor::NamesOfActivationMethodShort[dspec.getActivationMethod()] << "\n";
      fs << "MS_ONE_ID=" << dspec.getPrecursorScanNumber() << "\n"
         << "MS_ONE_SCAN=" << dspec.getPrecursorScanNumber() << "\n"
         << "PRECURSOR_MZ=" << std::to_string(dspec.getPrecursor().getMZ()) << "\n"
         << "PRECURSOR_CHARGE=" << (int)(dspec.getPrecursor().getCharge() * decoy_harmonic_factor) << "\n"
         << "PRECURSOR_MASS=" << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass() * decoy_harmonic_factor + decoy_precursor_offset) << "\n"
         << "PRECURSOR_INTENSITY=" << dspec.getPrecursor().getIntensity() << "\n";
    }

    fs << std::setprecision(-1);

    double qscore_threshold = 0;
    std::vector<double> qscores;

    if (dspec.size() > topFD_max_peak_count_) // max peak count for TopPic = 500
    {
      qscores.reserve(dspec.size());
      for (auto& pg : dspec)
      {
        qscores.push_back(pg.getQScore());
      }
      std::sort(qscores.begin(), qscores.end());
      qscore_threshold = qscores[qscores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(qscores);
    }

    int size = 0;
    for (auto& pg : dspec)
    {
      if (pg.getQScore() < qscore_threshold)
      {
        continue;
      }

      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.getMonoMass()) << "\t" << pg.getIntensity() << "\t" << (pg.isPositive() ? std::get<1>(pg.getAbsChargeRange()) : -std::get<1>(pg.getAbsChargeRange())) << "\n";
      fs << std::setprecision(-1);
      if (++size >= topFD_max_peak_count_)
      {
        break;
      }
    }
    fs << "END IONS\n\n";
  }
} // namespace OpenMS