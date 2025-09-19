// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <random>

namespace OpenMS
{
  /**
   * @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO

   */

  inline std::default_random_engine generator_;
  inline std::uniform_real_distribution<double> distribution_(0.0,1.0);

  void FLASHDeconvSpectrumFile::writeDeconvolvedMasses(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& file_name, const FLASHHelperClasses::PrecalculatedAveragine& avg, const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg, double tol,
                                                       const bool write_detail, const bool report_decoy, const double noise_decoy_weight)
  {
    if (!report_decoy && dspec.isDecoy()) return;
    static std::vector<uint> indices {};
    std::stringstream ss;
    if (dspec.empty())
    {
      return;
    }

    while (indices.size() <= dspec.getOriginalSpectrum().getMSLevel())
    {
      indices.push_back(1);
    }
    uint& index = indices[dspec.getOriginalSpectrum().getMSLevel() - 1];

    std::stringstream precursor_ss;
    if (dspec.getOriginalSpectrum().getMSLevel() > 1)
    {
      precursor_ss << dspec.getPrecursorScanNumber() << "\t"
                   << (dspec.getPrecursorPeakGroup().getFeatureIndex() == 0 ? "nan" : std::to_string(dspec.getPrecursorPeakGroup().getFeatureIndex()))
                   << "\t" << std::to_string(dspec.getPrecursor().getMZ()) << "\t" << dspec.getPrecursor().getIntensity() << "\t"
                   << dspec.getPrecursor().getCharge() << "\t";

      if (dspec.getPrecursorPeakGroup().empty())
      {
        precursor_ss << "nan\tnan\tnan\tnan\t";
        if (report_decoy) precursor_ss << "nan\t";
      }
      else
      {
        precursor_ss << dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) << "\t"
                     << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) << "\t"
                     << std::to_string(dspec.getPrecursorPeakGroup().getQscore()) << "\t" << dspec.getPrecursorPeakGroup().getQscore2D() << "\t";
        if (report_decoy) { precursor_ss << dspec.getPrecursorPeakGroup().getQvalue() << "\t"; }
      }
    }

    for (int i = 0; i < dspec.size(); i++)
    {
      auto& pg = dspec[i];
      if (!report_decoy && pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target) continue;

      if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy)
      {
        double number = distribution_(generator_);
        if (number > noise_decoy_weight)
        {
          continue;
        }
        if (number * noise_decoy_weight > 1.0)
        {
          i --;
        }
      }

      const auto& avg_ = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy? decoy_avg : avg;
      const double mono_mass = pg.getMonoMass();
      const double avg_mass = pg.getMonoMass() + avg_.getAverageMassDelta(mono_mass);
      const double intensity = pg.getIntensity();

      const auto& charge_range = pg.getAbsChargeRange();
      int min_charge = pg.isPositive() ? std::get<0>(charge_range) : -std::get<1>(charge_range);
      int max_charge = pg.isPositive() ? std::get<1>(charge_range) : -std::get<0>(charge_range);

      pg.setIndex(index);
      ss << index++ << "\t" << file_name << "\t" << pg.getScanNumber() << "\t" << (pg.getFeatureIndex() == 0 ? "nan" : std::to_string(pg.getFeatureIndex())) << "\t";

      if (report_decoy)
      {
        ss << pg.getTargetDecoyType() << "\t";
      }
      ss << std::to_string(dspec.getOriginalSpectrum().getRT()/60.0) << "\t" << dspec.size() << "\t" << std::to_string(avg_mass) << "\t" << std::to_string(mono_mass) << "\t" << intensity << "\t"
         << min_charge << "\t" << max_charge << "\t" << pg.size() << "\t";

      if (write_detail)
      {
        auto noisy_peaks = pg.recruitAllPeaksInSpectrum(dspec.getOriginalSpectrum(), tol * 1e-6, avg_, pg.getMonoMass(), false);

        std::sort(noisy_peaks.begin(), noisy_peaks.end());
        ss << std::fixed << std::setprecision(2);
        for (auto& p : pg)
        {
          ss << std::to_string(p.mz) << " ";
        }

        ss << "\t";
        ss << std::fixed << std::setprecision(1);
        for (auto& p : pg)
        {
          ss << p.intensity << " ";
        }

        ss << "\t";
        ss << std::setprecision(-1);
        for (auto& p : pg)
        {
          ss << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }

        ss << "\t";
        for (auto& p : pg)
        {
          ss << p.getUnchargedMass() << " ";
        }

        ss << "\t";
        for (auto& p : pg)
        {
          ss << p.isotopeIndex << " ";
        }

        ss << "\t";
        ss << std::setprecision(2);
        for (auto& p : pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / p.abs_charge + FLASHHelperClasses::getChargeMass(p.is_positive) - p.mz) / p.mz;
          ss << 1e6 * mass_error << " ";
        }
        ss << std::setprecision(-1);
        ss << "\t";
        ss << std::fixed << std::setprecision(2);
        for (auto& np : noisy_peaks)
        {
          ss << std::to_string(np.mz) << " ";
        }

        ss << "\t";
        ss << std::fixed << std::setprecision(1);
        for (auto& np : noisy_peaks)
        {
          ss << np.intensity << " ";
        }

        ss << "\t";
        ss << std::setprecision(-1);
        for (auto& np : noisy_peaks)
        {
          ss << (np.is_positive ? np.abs_charge : -np.abs_charge) << " ";
        }

        ss << "\t";
        for (auto& np : noisy_peaks)
        {
          ss << np.getUnchargedMass() << " ";
        }

        ss << "\t";
        for (auto& np : noisy_peaks)
        {
          ss << np.isotopeIndex << " ";
        }

        ss << "\t";
        ss << std::setprecision(2);
        for (auto& np : noisy_peaks)
        {
          double average_mass = pg.getMonoMass() + np.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / np.abs_charge + FLASHHelperClasses::getChargeMass(np.is_positive) - np.mz) / np.mz;
          ss << 1e6 * mass_error << " ";
        }
        ss << std::setprecision(-1);
        ss << "\t";
      }

      if (dspec.getOriginalSpectrum().getMSLevel() > 1)
      {
        ss << precursor_ss.str();
      }
      ss << pg.getIsotopeCosine() << "\t" << pg.getChargeIsotopeCosine(pg.getRepAbsCharge()) << "\t" << pg.getChargeScore() << "\t";

      const auto& max_qscore_mz_range = pg.getRepMzRange();
      ss << pg.getSNR() << "\t" << pg.getChargeSNR(pg.getRepAbsCharge()) << "\t" << pg.getAvgPPMError() << "\t" << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t" << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t" << std::to_string(pg.getQscore()) << "\t"
         << std::to_string(pg.getQscore2D());

      if (report_decoy)
      {
        ss << "\t" << pg.getQvalue();
      }

      if (write_detail)
      {
        ss << "\t" << std::setprecision(-1);

        for (int j = min_charge; j <= max_charge; j++)
        {
          ss << pg.getChargeIntensity(j);

          if (j < max_charge)
          {
            ss << ";";
          }
        }
        ss << "\t";

        const auto& iso_intensities = pg.getIsotopeIntensities();
        for (size_t j = 0; j < iso_intensities.size(); j++)
        {
          ss << iso_intensities[j];
          if (j < iso_intensities.size() - 1)
          {
            ss << ";";
          }
        }
      }
      ss << "\n";
    }
    fs << ss.str();
  }

  void FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(std::fstream& fs, const uint ms_level, const bool detail, const bool report_decoy)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue\t";
        }
        fs << "PerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "PrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQscore\tPrecursorQscore2D\t";
        if (report_decoy)
        {
          fs << "PrecursorQvalue\t";
        }
        fs << "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue\t";
        }
        fs << "PerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue";
        }
        fs << "\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "PrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQscore\tPrecursorQscore2D\t";
        if (report_decoy)
        {
          fs << "PrecursorQvalue\t";
        }
        fs << "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue\t";
        }
        fs << "\n";
      }
    }
  }

  void FLASHDeconvSpectrumFile::writeIsobaricQuantification(std::fstream& fs, std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
  {
    std::stringstream ss;
    ss << "Scan\tPrecursorScan\tPrecursorMZ\tRetentionTime\tMassCountInSpec\tQuantifiedChannelCount\tPrecursorQscore\tActivation\tActivationEnergy"
          "\tPrecursorCharge\tIsolationLeft\tIsolationRight\tPrecursorMass\tPrecursorSNR";

    Size channel_count = 0;
    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.isDecoy() || dspec.getOriginalSpectrum().getMSLevel() == 1)
        continue;
      auto quant = dspec.getQuantities();
      if (!quant.empty())
      {
        channel_count = quant.quantities.size();
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tNormalizedQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tMergedQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tNormalizedMergedQuantForCh" << (i + 1);
        }
        ss << "\n";
        break;
      }
    }
    if (channel_count == 0) return;

    for (auto& dspec : deconvolved_spectra)
    {
      bool precursor_found = ! dspec.getPrecursorPeakGroup().empty();
      if (dspec.isDecoy() || dspec.getOriginalSpectrum().getMSLevel() == 1) continue;
      int scan = dspec.getScanNumber();
      auto quant = dspec.getQuantities();
      // if (quant.empty())
      //   continue;
      int ch_count = 0;
      for (auto q : quant.quantities) if (q != 0) ch_count++;

      ss << scan << "\t" << dspec.getPrecursorScanNumber() << "\t" << dspec.getPrecursor().getMZ() << "\t" << dspec.getOriginalSpectrum().getRT()
         << "\t" << dspec.size() << "\t" << ch_count << "\t" << (precursor_found ? std::to_string(dspec.getPrecursorPeakGroup().getQscore2D()) : "NaN") << "\t"
         << (dspec.getActivationMethod() < Precursor::ActivationMethod::SIZE_OF_ACTIVATIONMETHOD
               ? Precursor::NamesOfActivationMethodShort[dspec.getActivationMethod()]
               : "N/A")
         << "\t" << dspec.getPrecursor().getActivationEnergy() << "\t" << dspec.getPrecursor().getCharge() << "\t"
         << dspec.getPrecursor().getMZ() - dspec.getPrecursor().getIsolationWindowLowerOffset() << "\t"
         << dspec.getPrecursor().getMZ() + dspec.getPrecursor().getIsolationWindowUpperOffset() << "\t"
         << (precursor_found ? std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) : "NaN") << "\t"
         << (precursor_found ? std::to_string(dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursorCharge())) : "NaN");
      double sum = 0;

      if (quant.empty())
      {
        for (Size i = 0; i < channel_count * 4; i++)
        {
          ss << "\t" << 0;
        }
        ss << "\n";
      }
      else
      {
        for (auto q : quant.quantities)
        {
          ss << "\t" << std::to_string(q);
          sum += q;
        }
        for (auto q : quant.quantities)
        {
          ss << "\t" << std::to_string(q / sum);
        }
        sum = 0;
        for (auto q : quant.merged_quantities)
        {
          ss << "\t" << std::to_string(q);
          sum += q;
        }
        for (auto q : quant.merged_quantities)
        {
          ss << "\t" << std::to_string(q / sum);
        }
        ss << "\n";
      }
    }
    fs << ss.str();
  }

  void FLASHDeconvSpectrumFile::writeMzML(const MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const String& deconvolved_mzML_file, const String& annotated_mzML_file,
                                          int mzml_charge, DoubleList tols)
  {
    if (deconvolved_mzML_file.empty() && annotated_mzML_file.empty())
      return;

    MSExperiment deconvolved_map;
    MSExperiment annotated_map;

    if (!deconvolved_mzML_file.empty())
    {
      deconvolved_map = MSExperiment(map);
      deconvolved_map.clear(false);
    }
    if (!annotated_mzML_file.empty())
    {
      annotated_map = MSExperiment(map);
      annotated_map.clear(false);
    }

    for (auto & deconvolved_spectrum : deconvolved_spectra)
    {
      if (deconvolved_spectrum.isDecoy()) continue;
      auto deconvolved_mzML = deconvolved_spectrum.toSpectrum(mzml_charge, tols[deconvolved_spectrum.getOriginalSpectrum().getMSLevel() - 1], false);
      if (!deconvolved_mzML_file.empty())
      {
        if (deconvolved_mzML.empty())
          continue;
        deconvolved_map.addSpectrum(deconvolved_mzML);
      }
      if (!annotated_mzML_file.empty())
      {
        auto anno_spec = MSSpectrum(deconvolved_spectrum.getOriginalSpectrum());
        anno_spec.sortByPosition();
        std::stringstream val {};

        for (auto& pg : deconvolved_spectrum)
        {
          if (pg.empty() || pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target)
          {
            continue;
          }
          val << std::to_string(pg.getMonoMass()) << ":";
          size_t pindex = 0;
          for (size_t k = 0; k < pg.size(); k++)
          {
            auto& p = pg[k];

            // Increment the pindex while it points to values less than p.mz
            while (pindex < anno_spec.size() - 1 && anno_spec[pindex].getMZ() < p.mz)
            {
              pindex++;
            }
            size_t nearest_index = pindex;
            if (anno_spec[pindex].getMZ() != p.mz)
            {
              if (pindex > 0 && std::abs(anno_spec[pindex - 1].getMZ() - p.mz) < std::abs(anno_spec[pindex].getMZ() - p.mz))
              {
                nearest_index = pindex - 1;
              }
            }

            val << nearest_index;

            if (k < pg.size() - 1)
            {
              val << ",";
            }
          }
          val << ";";
        }
        if (anno_spec.empty()) continue;
        anno_spec.setMetaValue("DeconvMassPeakIndices", val.str());
        annotated_map.addSpectrum(anno_spec);
      }
    }

    if (!deconvolved_mzML_file.empty())
    {
      MzMLFile mzml_file;
      mzml_file.store(deconvolved_mzML_file, deconvolved_map);
    }

    if (!annotated_mzML_file.empty())
    {
      MzMLFile mzml_file;
      mzml_file.store(annotated_mzML_file, annotated_map);
    }
  }

  void FLASHDeconvSpectrumFile::writeTopFDHeader(std::fstream& fs, const Param& param)
  {
    fs << "#FLASHDeconv generated msalign file\n";
    fs << "####################### Parameters ######################\n";
    for (const auto& p : param)
    {
      fs << "#" << p.name << ": " << p.value << "\n";
    }
    fs << "####################### Parameters ######################\n";
  }

  void FLASHDeconvSpectrumFile::writeTopFD(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& filename, const double qval_threshold, const uint min_ms_level,
                                           const bool randomize_precursor_mass, const bool randomize_fragment_mass)
  {
    std::stringstream ss;
    UInt ms_level = dspec.getOriginalSpectrum().getMSLevel();
    if (ms_level > min_ms_level)
    {
      if (dspec.getPrecursorPeakGroup().empty() ||
          dspec.getPrecursorPeakGroup().getQvalue() > qval_threshold)
      {
        return;
      }
    }

    if (dspec.size() < topFD_min_peak_count_)
    {
      return;
    }
    ss << std::fixed << std::setprecision(2);
    ss << "BEGIN IONS\n"
       << "FILE_NAME=" << filename << "\n"
       << "SPECTRUM_ID=" << dspec.getScanNumber() - 1 << "\n"
       << "TITLE=" << (dspec.isDecoy()? "DScan_" : "Scan_") << dspec.getScanNumber() << "\n"
       //<< "NATIVE_ID=" << dspec.getOriginalSpectrum().getNativeID() << "\n"
       //<< "FRACTION_ID=" << 0 << "\n"
       << "SCANS=" << dspec.getScanNumber() << "\n"
       << "RETENTION_TIME=" << dspec.getOriginalSpectrum().getRT() << "\n"
       << "LEVEL=" << dspec.getOriginalSpectrum().getMSLevel() << "\n";


    if (ms_level > 1)
    {
      double precursor_mass = dspec.getPrecursorPeakGroup().getMonoMass();

      ss << "MS_ONE_ID=" << dspec.getPrecursorScanNumber() - 1 << "\n"
         << "MS_ONE_SCAN=" << dspec.getPrecursorScanNumber() << "\n"
        << "PRECURSOR_WINDOW_BEGIN=" << -dspec.getPrecursor().getIsolationWindowLowerOffset() + dspec.getPrecursor().getMZ() << "\n"
         << "PRECURSOR_WINDOW_END=" << dspec.getPrecursor().getIsolationWindowUpperOffset() + dspec.getPrecursor().getMZ() << "\n";

        if (dspec.getActivationMethod() < Precursor::ActivationMethod::SIZE_OF_ACTIVATIONMETHOD)
      {
        ss << "ACTIVATION=" << Precursor::NamesOfActivationMethodShort[dspec.getActivationMethod()] << "\n";
      }
       ss << "PRECURSOR_MZ=" << std::to_string(dspec.getPrecursor().getMZ()) << "\n"
         << "PRECURSOR_CHARGE=" << (int)(dspec.getPrecursorCharge()) << "\n"
         << "PRECURSOR_MASS=" << std::to_string(precursor_mass + (randomize_precursor_mass ? (((double)rand() / (RAND_MAX)) * 200.0 - 100.0) : .0)) << "\n" // random number between 0 and 100.
         << "PRECURSOR_INTENSITY=" << dspec.getPrecursor().getIntensity() << "\n"
         << "PRECURSOR_FEATURE_ID=" << dspec.getPrecursorPeakGroup().getFeatureIndex() << "\n"
         << "PRECURSOR_QSCORE=" << dspec.getPrecursorPeakGroup().getQscore2D() << "\n";
    }

    ss << std::setprecision(-1);

    double qscore_threshold = 0;
    std::vector<double> qscores;

    if (dspec.size() > topFD_max_peak_count_) // max peak count for TopPic = 500
    {
      qscores.reserve(dspec.size());
      for (auto& pg : dspec)
      {
        qscores.push_back(pg.getQscore2D());
      }
      std::sort(qscores.begin(), qscores.end());
      qscore_threshold = qscores[qscores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(qscores);
    }

    int size = 0;
    for (auto& pg : dspec)
    {
      if (pg.getQscore2D() < qscore_threshold || pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target)
      {
        continue;
      }

      ss << std::fixed << std::setprecision(2);
      ss << std::to_string(pg.getMonoMass() + (randomize_fragment_mass ? (((double)rand() / (RAND_MAX)) * 200.0 - 100.0) : .0)) << "\t" << pg.getIntensity() << "\t"
         << (pg.isPositive() ? std::get<1>(pg.getAbsChargeRange()) : -std::get<1>(pg.getAbsChargeRange())) << "\n";
      ss << std::setprecision(-1);
      if (++size >= topFD_max_peak_count_)
      {
        break;
      }
    }
    ss << "END IONS\n\n";
    fs << ss.str();
  }
} // namespace OpenMS
