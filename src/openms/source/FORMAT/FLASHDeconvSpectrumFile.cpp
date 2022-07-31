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

namespace OpenMS
{
  /**
    @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO
**/


  void FLASHDeconvSpectrumFile::writeDeconvolvedMasses(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& file_name, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                                      const bool write_detail)
  {
    if (dspec.empty())
    {
      return;
    }
    int index = 1;

    for (auto& pg : dspec)
    {
      const double mono_mass = pg.getMonoMass();
      const double avg_mass = pg.getMonoMass() + avg.getAverageMassDelta(mono_mass);
      const double intensity = pg.getIntensity();

      auto charge_range = pg.getAbsChargeRange();
      int min_charge = pg.isPositive() ? std::get<0>(charge_range) : -std::get<1>(charge_range);
      int max_charge = pg.isPositive() ? std::get<1>(charge_range) : -std::get<0>(charge_range);

      fs << index++ << "\t" << file_name << "\t" << pg.getScanNumber() << "\t" << pg.getDecoyIndex() << "\t" << std::to_string(dspec.getOriginalSpectrum().getRT()) << "\t" << dspec.size() << "\t"
         << std::to_string(avg_mass) << "\t" << std::to_string(mono_mass) << "\t" << intensity << "\t" << min_charge << "\t" << max_charge << "\t" << pg.size() << "\t";

      if (write_detail)
      {
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
        for (auto& p : pg.noisy_peaks)
        {
          fs << p.mz << " ";
        }

        fs << "\t";

        fs << std::fixed << std::setprecision(1);
        for (auto& p : pg.noisy_peaks)
        {
          fs << p.intensity << " ";
        }


        fs << "\t";
        fs << std::setprecision(-1);

        for (auto& p : pg.noisy_peaks)
        {
          fs << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }

        fs << "\t";
        for (auto& p : pg.noisy_peaks)
        {
          fs << p.getUnchargedMass() << " ";
        }

        fs << "\t";
        for (auto& p : pg.noisy_peaks)
        {
          fs << p.isotopeIndex << " ";
        }

        fs << "\t";

        for (auto& p : pg.noisy_peaks)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) - p.mz) / p.mz;
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
          fs << "nan\tnan\tnan\tnan\tnan\tnan\tnan\t";
        }
        else
        {
          fs << dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) << "\t" << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) << "\t"
             << dspec.getPrecursorPeakGroup().getQScore() << "\t" << dspec.getPrecursorPeakGroup().getQvalue() << "\t" << dspec.getPrecursorPeakGroup().getQvalueWithIsotopeDecoyOnly() << "\t"  << dspec.getPrecursorPeakGroup().getQvalueWithNoiseDecoyOnly() << "\t"   << dspec.getPrecursorPeakGroup().getQvalueWithChargeDecoyOnly() << "\t" ;
        }
      }
      fs << pg.getIsotopeCosine() << "\t" << pg.getChargeScore() << "\t";

      auto max_qscore_mz_range = pg.getMaxQScoreMzRange();
      fs << pg.getSNR() << "\t" << pg.getChargeSNR(pg.getRepAbsCharge()) << "\t" << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t" << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t" << pg.getQScore() << "\t" << pg.getQvalue() << "\t" << pg.getQvalueWithIsotopeDecoyOnly() << "\t" << pg.getQvalueWithNoiseDecoyOnly() << "\t" << pg.getQvalueWithChargeDecoyOnly()<<"\t"
         << std::setprecision(-1); //

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
      for (int i=0;i<iso_intensities.size();i++)
      {
        fs << iso_intensities[i];
        if (i < iso_intensities.size() - 1)
        {
          fs << ";";
        }
      }
      fs << "\n";
    }
  }

  void FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(std::fstream& fs, const int ms_level, const bool detail)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\tDecoy\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tQvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\tDecoy\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQScore\tPrecursorQvalue\tPrecursorQvalueWithIsotopeDecoyOnly\tPrecursorQvalueWithNoiseDecoyOnly\tPrecursorQvalueWithChargeDecoyOnly\t"
              "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tQvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\tDecoy\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
              "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tQvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\tDecoy\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQScore\tPrecursorQvalue\tPrecursorQvalueWithIsotopeDecoyOnly\tPrecursorQvalueWithNoiseDecoyOnly\tPrecursorQvalueWithChargeDecoyOnly\t"
              "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tQvalue\tQvalueWithIsotopeDecoyOnly\tQvalueWithNoiseDecoyOnly\tQvalueWithChargeDecoyOnly\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
  }


  void FLASHDeconvSpectrumFile::writeTopFD(const DeconvolvedSpectrum& dspec, std::fstream& fs,
                                           //                                       const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                           const double snr_threshold, const double decoy_harmonic_factor,
                                           const double decoy_precursor_offset) //, fstream& fsm, fstream& fsp)
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
      fs << "ACTIVATION=" << dspec.getActivationMethod() << "\n";
      // if (!dspec.getPrecursorPeakGroup().empty())
      // {
      fs << "MS_ONE_ID=" << dspec.getPrecursorScanNumber() << "\n"
         << "MS_ONE_SCAN=" << dspec.getPrecursorScanNumber() << "\n"
         << "PRECURSOR_MZ=" << std::to_string(dspec.getPrecursor().getMZ()) << "\n"
         << "PRECURSOR_CHARGE=" << (int)(dspec.getPrecursor().getCharge() * decoy_harmonic_factor) << "\n"
         << "PRECURSOR_MASS=" << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass() * decoy_harmonic_factor + decoy_precursor_offset) << "\n"
         << "PRECURSOR_INTENSITY=" << dspec.getPrecursor().getIntensity() << "\n";
      //}
      /*else
      {
        double average_mass =
            (dspec.getPrecursor().getMZ() -
             FLASHDeconvHelperStructs::getChargeMass(dspec.getPrecursor().getCharge() > 0)) *
            abs(dspec.getPrecursor().getCharge() * decoy_harmonic_factor);
        double mono_mass = average_mass - avg.getAverageMassDelta(average_mass) + decoy_precursor_offset;
        fs << "MS_ONE_ID=" << dspec.getPrecursorScanNumber() << "\n"
           << "MS_ONE_SCAN=" << dspec.getPrecursorScanNumber() << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(dspec.getPrecursor().getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << (int) (dspec.getPrecursor().getCharge() * decoy_harmonic_factor) << "\n"
           << "PRECURSOR_MASS=" << std::to_string(mono_mass) << "\n"
           << "PRECURSOR_INTENSITY=" << dspec.getPrecursor().getIntensity() << "\n";
      }*/
    }
    fs << std::setprecision(-1);

    double isotope_score_threshold = 0;
    std::vector<double> isotope_scores;

    if (dspec.size() > topFD_max_peak_count_) // max peak count for TopPic = 500
    {
      isotope_scores.reserve(dspec.size());
      for (auto& pg : dspec)
      {
        isotope_scores.push_back(pg.getIsotopeCosine());
      }
      std::sort(isotope_scores.begin(), isotope_scores.end());
      isotope_score_threshold = isotope_scores[isotope_scores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(isotope_scores);
    }

    int size = 0;
    for (auto& pg : dspec)
    {
      if (pg.getIsotopeCosine() < isotope_score_threshold)
      {
        continue;
      }
      if (size >= topFD_max_peak_count_)
      {
        break;
      }
      size++;
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.getMonoMass()) << "\t" << pg.getIntensity() << "\t" << (pg.isPositive() ? std::get<1>(pg.getAbsChargeRange()) : -std::get<1>(pg.getAbsChargeRange())) << "\n";
      fs << std::setprecision(-1);
      if (size >= topFD_max_peak_count_)
      {
        break;
      }
    }

    fs << "END IONS\n\n";
  }


} // namespace OpenMS