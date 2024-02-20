// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
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

  void FLASHDeconvSpectrumFile::writeDeconvolvedMasses(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& file_name, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double tol,
                                                       const bool write_detail, const bool report_decoy, const double noise_decoy_weight)
  {
    static std::vector<uint> indices {};

    if (dspec.empty())
    {
      return;
    }

    while (indices.size() <= dspec.getOriginalSpectrum().getMSLevel())
    {
      indices.push_back(1);
    }
    uint& index = indices[dspec.getOriginalSpectrum().getMSLevel() - 1];

    for (int i = 0; i < dspec.size(); i++)
    {
      auto& pg = dspec[i];

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

      const double mono_mass = pg.getMonoMass();
      const double avg_mass = pg.getMonoMass() + avg.getAverageMassDelta(mono_mass);
      const double intensity = pg.getIntensity();

      auto charge_range = pg.getAbsChargeRange();
      int min_charge = pg.isPositive() ? std::get<0>(charge_range) : -std::get<1>(charge_range);
      int max_charge = pg.isPositive() ? std::get<1>(charge_range) : -std::get<0>(charge_range);

      pg.setIndex(index);
      fs << index++ << "\t" << file_name << "\t" << pg.getScanNumber() << "\t" << (pg.getFeatureIndex() == 0 ? "nan" : std::to_string(pg.getFeatureIndex())) << "\t";

      if (report_decoy)
      {
        fs << pg.getTargetDecoyType() << "\t";
      }
      fs << std::to_string(dspec.getOriginalSpectrum().getRT()) << "\t" << dspec.size() << "\t" << std::to_string(avg_mass) << "\t" << std::to_string(mono_mass) << "\t" << intensity << "\t"
         << min_charge << "\t" << max_charge << "\t" << pg.size() << "\t";

      if (write_detail)
      {
        auto noisy_peaks = pg.recruitAllPeaksInSpectrum(dspec.getOriginalSpectrum(), tol * 1e-6, avg, pg.getMonoMass());

        std::sort(noisy_peaks.begin(), noisy_peaks.end());
        fs << std::fixed << std::setprecision(2);
        for (auto& p : pg)
        {
          fs << std::to_string(p.mz) << " ";
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
        fs << std::setprecision(2);
        for (auto& p : pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) - p.mz) / p.mz;
          fs << 1e6 * mass_error << " ";
        }
        fs << std::setprecision(-1);
        fs << "\t";
        fs << std::fixed << std::setprecision(2);
        for (auto& np : noisy_peaks)
        {
          fs << std::to_string(np.mz) << " ";
        }

        fs << "\t";
        fs << std::fixed << std::setprecision(1);
        for (auto& np : noisy_peaks)
        {
          fs << np.intensity << " ";
        }

        fs << "\t";
        fs << std::setprecision(-1);
        for (auto& np : noisy_peaks)
        {
          fs << (np.is_positive ? np.abs_charge : -np.abs_charge) << " ";
        }

        fs << "\t";
        for (auto& np : noisy_peaks)
        {
          fs << np.getUnchargedMass() << " ";
        }

        fs << "\t";
        for (auto& np : noisy_peaks)
        {
          fs << np.isotopeIndex << " ";
        }

        fs << "\t";
        fs << std::setprecision(2);
        for (auto& np : noisy_peaks)
        {
          double average_mass = pg.getMonoMass() + np.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / np.abs_charge + FLASHDeconvHelperStructs::getChargeMass(np.is_positive) - np.mz) / np.mz;
          fs << 1e6 * mass_error << " ";
        }
        fs << std::setprecision(-1);
        fs << "\t";
      }
      if (dspec.getOriginalSpectrum().getMSLevel() > 1)
      {
        // PrecursorScanNum	PrecursorMz	PrecursorIntensity PrecursorCharge	PrecursorMonoMass		PrecursorQscore
        fs << dspec.getPrecursorScanNumber() << "\t" << (dspec.getPrecursorPeakGroup().getFeatureIndex() == 0 ? "nan" : std::to_string(dspec.getPrecursorPeakGroup().getFeatureIndex())) << "\t"
           << std::to_string(dspec.getPrecursor().getMZ()) << "\t" << dspec.getPrecursor().getIntensity() << "\t" << dspec.getPrecursor().getCharge() << "\t";

        if (dspec.getPrecursorPeakGroup().empty())
        {
          fs << "nan\tnan\tnan\tnan\t";
          if (report_decoy)
            fs << "nan\t";
        }
        else
        {
          fs << dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) << "\t" << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) << "\t"
             << std::to_string(dspec.getPrecursorPeakGroup().getQscore()) << "\t" << dspec.getPrecursorPeakGroup().getQscore2D() << "\t";
          if (report_decoy)
          {
            fs << dspec.getPrecursorPeakGroup().getQvalue() << "\t";
          }
        }
      }
      fs << pg.getIsotopeCosine() << "\t" << pg.getChargeIsotopeCosine(pg.getRepAbsCharge()) << "\t" << pg.getChargeScore() << "\t";

      auto max_qscore_mz_range = pg.getRepMzRange();
      fs << pg.getSNR() << "\t" << pg.getChargeSNR(pg.getRepAbsCharge()) << "\t" << pg.getAvgPPMError() << "\t" << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t" << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t" << std::to_string(pg.getQscore()) << "\t"
         << std::to_string(pg.getQscore2D());


      if (report_decoy)
      {
        fs << "\t" << pg.getQvalue();
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
    fs << "Scan\tPrecursorScan\tPrecursorMass\tPrecursorSNR";
    bool begin = true;
    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.getOriginalSpectrum().getMSLevel() == 1)
        continue;
      if (dspec.getPrecursorPeakGroup().empty())
        continue;
      int scan = dspec.getScanNumber();
      auto quant = dspec.getQuantities();
      if (quant.empty())
        continue;
      if (begin)
      {
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          fs << "\tQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          fs << "\tNormalizedQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          fs << "\tMergedQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          fs << "\tNormalizedMergedQuantForCh" << (i + 1);
        }
        fs << "\n";
        begin = false;
      }
      fs << scan << "\t" << dspec.getPrecursorScanNumber() << "\t" << dspec.getPrecursorPeakGroup().getMonoMass() << "\t" << dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursorCharge());
      double sum = 0;
      for (auto q : quant.quantities)
      {
        fs << "\t" << std::to_string(q);
        sum += q;
      }
      for (auto q : quant.quantities)
      {
        fs << "\t" << std::to_string(q / sum);
      }
      sum = 0;
      for (auto q : quant.merged_quantities)
      {
        fs << "\t" << std::to_string(q);
        sum += q;
      }
      for (auto q : quant.merged_quantities)
      {
        fs << "\t" << std::to_string(q / sum);
      }
      fs << "\n";
    }
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

    uint current_min_ms_level = 0;
    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      if (current_min_ms_level == 0 || current_min_ms_level > ms_level)
        current_min_ms_level = ms_level;
    }

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      auto deconvolved_mzML = deconvolved_spectrum.toSpectrum(mzml_charge, current_min_ms_level, tols[deconvolved_spectrum.getOriginalSpectrum().getMSLevel() - 1], false);
      if (!deconvolved_mzML_file.empty())
      {
        if (deconvolved_mzML.empty())
          continue;
        deconvolved_map.addSpectrum(deconvolved_mzML);
      }
      if (!annotated_mzML_file.empty())
      {
        auto anno_spec = MSSpectrum(deconvolved_spectrum.getOriginalSpectrum());

        std::stringstream val {};

        for (auto& pg : deconvolved_spectrum)
        {
          val << std::to_string(pg.getMonoMass()) << ":";
          for (size_t k = 0; k < pg.size(); k++)
          {
            auto& p = pg[k];
            auto pindex = anno_spec.findNearest(p.mz);
            val << pindex;
            if (k < pg.size() - 1)
            {
              val << ",";
            }
          }
          val << ";";
        }
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

  void FLASHDeconvSpectrumFile::writeDLMatrixHeader(std::fstream& fs)
  {
    for (int i = 0; i < 3; i++)
    {
      String prefix = "Set" + std::to_string(i + 1) + "_";
      PeakGroup peak_group;
      int dlm = peak_group.getIsotopeRangeForDL() * peak_group.getChargeRangeForDL() / peak_group.getBinWidthDL();

      for (int j = 0; j < dlm; j++)
      {
        fs << prefix << j << ",";
      }
    }
    fs << "Class\n";
  }

  template<class BidiIter>
  BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random)
  {
    size_t left = std::distance(begin, end);
    while (num_random--)
    {
      BidiIter r = begin;
      std::advance(r, rand() % left);
      std::swap(*begin, *r);
      ++begin;
      --left;
    }
    return begin;
  }

  void FLASHDeconvSpectrumFile::writeDLMatrix(std::vector<DeconvolvedSpectrum>& dspecs, double tol, std::fstream& fs, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg)
  {
    String cns[] = {"T", "D1", "D2", "D3"};
    // int count = 2000000;
    // class,ID,group,data
    std::vector<std::vector<PeakGroup>> grouped(4);

    for (auto& dspec : dspecs)
    {
      for (auto& pg : dspec)
      {
        if (pg.size() == 0)
        {
          continue;
        }
        int cl = pg.getTargetDecoyType();
        if (cl < 0 || cl >= (int)grouped.size())
        {
          continue;
        }

        pg.calculateDLMatrices(dspec.getOriginalSpectrum(), tol, avg);

        auto dlmatrix = pg.getDLMatrix(0).asVector();
        grouped[cl].push_back(pg);
      }
    }

    /*
    for(auto& g : grouped)
    {
      if((int)g.size() < count)
      {
         continue;
      }
      random_unique(g.begin(), g.end(), count);
    }
  */
    for (auto& g : grouped)
    {
      // int cntr = 0;
      for (auto& pg : g)
      {
        int cl = pg.getTargetDecoyType();

        for (int i = 0; i < 3; i++)
        {
          auto dlm = pg.getDLMatrix(i).asVector();

          for (double v : dlm)
          {
            fs << v << ",";
          }
        }
        fs << cns[cl] << "\n";
        // if(++cntr >= count)
        //{
        //  break;
        // }
      }
    }
  }

  void FLASHDeconvSpectrumFile::writeTopFDHeader(std::fstream& fs, const Param& param)
  {
    fs << "#FLASHDeconv generated msalign file\n";
    fs << "####################### Parameters ######################\n";
    fs << "#Maximum charge:                              \t" << param.getValue("max_charge") << "\n";
    fs << "#Maximum monoisotopic mass:                   \t" << param.getValue("max_mass") << " Dalton\n";
    fs << "#Peak error tolerance:                        \t" << param.getValue("tol") << " ppm\n";
    fs << "####################### Parameters ######################\n";
  }

  void FLASHDeconvSpectrumFile::writeTopFD(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& filename, const double snr_threshold, const double qval_threshold, const uint min_ms_level,
                                           const bool randomize_precursor_mass, const bool randomize_fragment_mass)
  {
    UInt ms_level = dspec.getOriginalSpectrum().getMSLevel();
    if (ms_level > min_ms_level)
    {
      if (dspec.getPrecursorPeakGroup().empty() || dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) < snr_threshold ||
          dspec.getPrecursorPeakGroup().getQvalue() > qval_threshold)
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
       << "FILE_NAME=" << filename << "\n"
       << "NATIVE_ID=" << dspec.getOriginalSpectrum().getNativeID() << "\n"
       << "FRACTION_ID=" << 0 << "\n"
       << "SCANS=" << dspec.getScanNumber() << "\n"
       << "RETENTION_TIME=" << dspec.getOriginalSpectrum().getRT() << "\n"
       << "LEVEL=" << dspec.getOriginalSpectrum().getMSLevel() << "\n";


    if (ms_level > 1)
    {
      double precursor_mass = dspec.getPrecursorPeakGroup().getMonoMass();
      if (dspec.getActivationMethod() < Precursor::ActivationMethod::SIZE_OF_ACTIVATIONMETHOD)
      {
        fs << "ACTIVATION=" << Precursor::NamesOfActivationMethodShort[dspec.getActivationMethod()] << "\n";
      }
      fs << "MS_ONE_ID=" << dspec.getPrecursorScanNumber() << "\n"
         << "MS_ONE_SCAN=" << dspec.getPrecursorScanNumber() << "\n"
         << "PRECURSOR_MZ=" << std::to_string(dspec.getPrecursor().getMZ()) << "\n"
         << "PRECURSOR_CHARGE=" << (int)(dspec.getPrecursorCharge()) << "\n"
         << "PRECURSOR_MASS=" << std::to_string(precursor_mass + (randomize_precursor_mass ? (((double)rand() / (RAND_MAX)) * 200.0 - 100.0) : .0)) << "\n" // random number between 0 and 100.
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
        qscores.push_back(pg.getQscore2D());
      }
      std::sort(qscores.begin(), qscores.end());
      qscore_threshold = qscores[qscores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(qscores);
    }

    int size = 0;
    for (auto& pg : dspec)
    {
      if (pg.getQscore2D() < qscore_threshold)
      {
        continue;
      }

      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.getMonoMass() + (randomize_fragment_mass ? (((double)rand() / (RAND_MAX)) * 200.0 - 100.0) : .0)) << "\t" << pg.getIntensity() << "\t"
         << (pg.isPositive() ? std::get<1>(pg.getAbsChargeRange()) : -std::get<1>(pg.getAbsChargeRange())) << "\n";
      fs << std::setprecision(-1);
      if (++size >= topFD_max_peak_count_)
      {
        break;
      }
    }
    fs << "END IONS\n\n";
  }
} // namespace OpenMS
