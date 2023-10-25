// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownIsobaricQuantifier.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  inline const Size max_peak_count_ = 3e4;

  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm() : DefaultParamHandler("FLASHDeconvAlgorithm"), ProgressLogger()
  {
    //   registerIntOption_(
    //      "report_FDR", "<0: Do not report 1: report>", 0,
    //      "Report qvalues (roughly, point-wise FDR) for deconvolved masses in the tsv files for deconvolved spectra. Dummy masses to calculate qvalues and FDR are also reported. Beta version.", false, false);
    //    setMinInt_("report_FDR", 0);
    //    setMaxInt_("report_FDR", 1);

    //    registerIntOption_("use_RNA_averagine", "", 0, "If set to 1, RNA averagine model is used", false, true);
    //    setMinInt_("use_RNA_averagine", 0);
    //    setMaxInt_("use_RNA_averagine", 1);

    //    registerIntOption_("preceding_MS1_count", "<number>", 3,
    //                       "Specifies the number of preceding MS1 spectra for MS2 precursor determination. "
    //                       "In TDP, the precursor peak of a MS2 spectrum may not belong to any "
    //                       "deconvolved masses in the MS1 spectrum immediately preceding the MS2 spectrum. "
    //                       "Increasing this parameter to N allows for the search for the deconvolved masses in the N preceding MS1 spectra from the MS2 spectrum"
    //                       ", increasing the chance that its precursor is deconvolved.",
    //                       false, true);
    //
    //    setMinInt_("preceding_MS1_count", 1);

//    registerIntOption_("max_MS_level", "<number>", 3, "Maximum MS level (inclusive) for deconvolution.", false, true);
//    setMinInt_("max_MS_level", 1);
//
//    registerIntOption_("forced_MS_level", "", 0,
//                       "If set to an integer N, MS level of all spectra will be set to N regardless of original MS level. Useful when deconvolving datasets containing only MS2 spectra.", false, true);
//    setMinInt_("forced_MS_level", 0);
//
//    registerIntOption_("merging_method", "<0: None 1: gaussian averaging 2: block method>", 0,
//                       "Method for spectra merging before deconvolution. 0: No merging "
//                       "1: Average gaussian method to perform moving gaussian averaging of spectra per MS level. Effective to increase proteoform ID sensitivity (in particular for Q-TOF datasets). "
//                       "2: Block method to perform merging of all spectra into a single one per MS level (e.g., for NativeMS datasets)",
//                       false);
//    setMinInt_("merging_method", 0);
//    setMaxInt_("merging_method", 2);
//    registerIntOption_("target_precursor_charge", "<Target precursor charge>", 0,
    //                       "Charge state of the target precursor. All precursor charge is fixed to this value. "
    //                       "This parameter is useful for targeted studies where MS2 spectra are generated from a fixed precursor (e.g., Native-MS). "
    //                       "This option also gives the maximum charge and masses (together with precursor m/z) of fragment ions, which overrides -Algorithm:max_charge and -Algorithm:max_mass.",
    //                       false, false);
    //
    //    registerDoubleOption_("target_precursor_mz", "<m/z value>", 0.0,
    //                          "Target precursor m/z value. This option must be used with -target_precursor_charge option. Otherwise it will be ignored. "
    //                          "If -target_precursor_charge option is used but this option is not used, the precursor m/z value written in MS2 spectra will be used by default. "
    //                          "Together with -target_precursor_charge, this option overrides -Algorithm:max_mass.",
    //                          false, false);

    // Param mf_param = getParam_().copy("FeatureTracing:", true);
    //
    //    if (((double)mf_param.getValue("mass_error_ppm")) < 0)
    //    {
    //      mf_param.setValue("mass_error_ppm", tols[0]);
    //    }
    //    mf_param.setValue("noise_threshold_int", .0);
    //    mf_param.setValue("reestimate_mt_sd", "false");
    //    mf_param.setValue("trace_termination_criterion", "outlier");
    //    mf_param.setValue("trace_termination_outliers", 20);
    //    mf_param.setValue("chrom_peak_snr", .0);
    //
    //    if (((double)mf_param.getValue("min_qscore")) < 0)
    //    {
    //      mf_param.setValue("min_qscore", qscores[0]);
    //    }

    defaults_.insert("Deconvolution", SpectralDeconvolution().getDefaults());
    //subsections_.emplace_back("FeatureTracing");
    defaults_.insert("FeatureTracing", MassFeatureTrace().getDefaults());
    defaults_.insert("IsobaricQuantification", TopDownIsobaricQuantifier().getDefaults());

    defaultsToParam_();
  }

  void FLASHDeconvAlgorithm::updateMembers_()
  {
    forced_ms_level_ = 0;
    max_ms_level_ = 1;

    current_min_ms_level_ = max_ms_level_;
    current_max_ms_level_ = 0;
    target_precursor_charge_ = 0;
    max_charge_ = 0;
    tols_ = DoubleList();
    preceding_MS1_count_ = param_.getValue("preceding_MS1_count");
    use_RNA_averagine_;
    report_decoy_;
  }

  int FLASHDeconvAlgorithm::scan_map_(MSExperiment& map)
  {
    // read input dataset once to count spectra
    double gradient_rt = .0;
    for (auto& it : map)
    {
      gradient_rt = std::max(gradient_rt, it.getRT());
      // if forced_ms_level > 0, force MS level of all spectra to 1.
      if (forced_ms_level_ > 0)
      {
        it.setMSLevel(forced_ms_level_);
      }

      if (it.empty())
      {
        continue;
      }

      if (it.getMSLevel() > max_ms_level_)
      {
        continue;
      }

      uint ms_level = it.getMSLevel();
      current_max_ms_level_ = current_max_ms_level_ < ms_level ? ms_level : current_max_ms_level_;
      current_min_ms_level_ = current_min_ms_level_ > ms_level ? ms_level : current_min_ms_level_;

      if (ms_level > 1 && target_precursor_charge_ != 0)
      {
        if (it.getPrecursors().empty())
        {
          if (target_precursor_mz_ == 0)
          {
            OPENMS_LOG_INFO << "Target precursor charge is set but no precursor is found in MS2 spectra. Specify target precursor m/z with -target_precursor_mz option" << std::endl;
            return -1;
          }
          else
          {
            Precursor precursor;
            precursor.setCharge(target_precursor_charge_);
            precursor.setMZ(target_precursor_mz_);
            it.setPrecursors({precursor});
          }
        }
        else
        {
          it.getPrecursors()[0].setCharge(target_precursor_charge_);
          if (target_precursor_mz_ != 0)
          {
            it.getPrecursors()[0].setMZ(target_precursor_mz_);
          }
        }
      }
    }
    // Max MS Level is adjusted according to the input dataset
    current_max_ms_level_ = current_max_ms_level_ > max_ms_level_ ? max_ms_level_ : current_max_ms_level_;
    return 0;
  }

  void FLASHDeconvAlgorithm::filterLowPeaks_(MSExperiment& map, Size count)
  {
    for (auto& it : map)
    {
      double threshold;
      if (it.getType(false) == SpectrumSettings::CENTROID)
      {
        if (it.size() <= count)
        {
          return;
        }
        it.sortByIntensity(true);
        threshold = it[count].getIntensity();
      }
      else
      {
        if (it.size() <= count)
        {
          continue;
        }

        it.sortByIntensity(true);
        double max_intensity = log10(it[0].getIntensity());
        double min_intensity = 0;
        for (auto& p : it)
        {
          if (p.getIntensity() <= 0)
          {
            break;
          }
          min_intensity = log10(p.getIntensity());
        }
        Size bin_size = 500;
        std::vector<int> freq(bin_size + 1, 0);
        for (auto& p : it)
        {
          if (p.getIntensity() <= 0)
          {
            break;
          }
          Size bin = round((log10(p.getIntensity()) - min_intensity) / (max_intensity - min_intensity) * bin_size);
          freq[bin]++;
        }

        auto mod_bin = std::distance(freq.begin(), std::max_element(freq.begin(), freq.end())); // most frequent intensity is the threshold to distinguish between signal and noise

        threshold =
          3.0 * (pow(10.0, (double)mod_bin / bin_size * (max_intensity - min_intensity) +
                             min_intensity)); // multiply by 3 to the most frequent intensity to make sure more signal component remains. Later this could be determined to use signal-to-noise ratio.
      }
      // pop back the low intensity peaks using threshold
      while (!it.empty() && it[it.size() - 1].getIntensity() <= threshold)
      {
        it.pop_back();
      }

      it.sortByPosition();
    }
  }

  void FLASHDeconvAlgorithm::mergeSpectra_(MSExperiment& map)
  {
    filterLowPeaks_(map, max_peak_count_);

    // if a merged spectrum is analyzed, replace the input dataset with the merged one
    if (merge_spec_ == 1)
    {
      OPENMS_LOG_INFO << "Merging spectra using gaussian averaging... " << std::endl;
      SpectraMerger merger;
      merger.setLogType(CMD);
      Param sm_param = merger.getDefaults();
      sm_param.setValue("average_gaussian:precursor_mass_tol", tols_[0]);
      sm_param.setValue("average_gaussian:precursor_max_charge", abs(max_charge_));
      sm_param.setValue("mz_binning_width", 1.0);

      merger.setParameters(sm_param);
      map.sortSpectra();

      for (uint tmp_ms_level = current_min_ms_level_; tmp_ms_level <= current_max_ms_level_; tmp_ms_level++)
      {
        merger.average(map, "gaussian", (int)tmp_ms_level);
      }
    }
    else if (merge_spec_ == 2)
    {
      OPENMS_LOG_INFO << "Merging spectra into a single spectrum per MS level... " << std::endl;
      SpectraMerger merger;
      merger.setLogType(CMD);
      Param sm_param = merger.getDefaults();
      sm_param.setValue("block_method:rt_max_length", .0);
      // rt_max_length = 0 TODO what is the meaning of this??
      map.sortSpectra();
      sm_param.setValue("mz_binning_width", 1.0);

      for (int ml = 1; ml <= (int)current_max_ms_level_; ml++)
      {
        // sm_param.setValue("mz_binning_width", tols[ml - 1] / 2.0);
        sm_param.setValue("block_method:ms_levels", IntList {ml});
        merger.setParameters(sm_param);
        merger.mergeSpectraBlockWise(map);
      }
    }
  }

  int FLASHDeconvAlgorithm::runFD_(const MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
  {
    auto last_deconvolved_spectra = std::unordered_map<UInt, std::vector<DeconvolvedSpectrum>>();

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      int scan_number = map.getSourceFiles().empty() ? -1 : SpectrumLookup::extractScanNumber(it->getNativeID(), map.getSourceFiles()[0].getNativeIDTypeAccession());

      if (scan_number < 0)
      {
        scan_number = (int)std::distance(map.begin(), it) + 1;
      }

      nextProgress();

      if (it->empty())
      {
        continue;
      }

      uint ms_level = it->getMSLevel();
      if (ms_level > current_max_ms_level_)
      {
        continue;
      }

      if (ms_level > 1 && target_precursor_charge_ != 0 && it->getPrecursors().empty())
      {
        OPENMS_LOG_INFO << "Target precursor charge is set but no precursor m/z is found in MS2 spectra. Specify target precursor m/z with -target_precursor_mz option" << std::endl;
        return -1;
      }

      // for MS>1 spectrum, register precursor
      std::vector<DeconvolvedSpectrum> precursor_specs;

      if (ms_level > 1 && last_deconvolved_spectra.find(ms_level - 1) != last_deconvolved_spectra.end())
      {
        precursor_specs = (last_deconvolved_spectra[ms_level - 1]);
      }

      sd_.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_ida_);
      auto& deconvolved_spectrum = sd_.getDeconvolvedSpectrum();
      if (deconvolved_spectrum.empty())
      {
        continue;
      }

      if (ms_level > 1 && target_precursor_charge_ != 0)
        setTargetPrecursorCharge_(deconvolved_spectrum, *it);

      if (it->getMSLevel() > 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
      {
        ms2scan_to_precursor_peak_group_map_[scan_number] = deconvolved_spectrum.getPrecursorPeakGroup();
      }

      if (ms_level < current_max_ms_level_)
      {
        if ((int)last_deconvolved_spectra[ms_level].size() >= preceding_MS1_count_)
        {
          last_deconvolved_spectra.erase(last_deconvolved_spectra.begin());
        }
        last_deconvolved_spectra[ms_level].push_back(deconvolved_spectrum);
      }

      if (report_decoy_)
      {
#pragma omp parallel sections default(none) shared(sd_charge_decoy_, sd_noise_decoy_, sd_isotope_decoy_, it, precursor_specs, scan_number, precursor_map_for_ida_)
        {
#pragma omp section
          sd_charge_decoy_.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_ida_);
#pragma omp section
          sd_noise_decoy_.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_ida_);
#pragma omp section
          sd_isotope_decoy_.performSpectrumDeconvolution(*it, precursor_specs, scan_number, precursor_map_for_ida_);
        }
        DeconvolvedSpectrum dummy_deconvolved_spectrum(scan_number);
        deconvolved_spectrum.sortByQscore();
        float qscore_threshold_for_dummy = deconvolved_spectrum[deconvolved_spectrum.size() - 1].getQscore();
        dummy_deconvolved_spectrum.setOriginalSpectrum(*it);
        dummy_deconvolved_spectrum.reserve(sd_isotope_decoy_.getDeconvolvedSpectrum().size() + sd_charge_decoy_.getDeconvolvedSpectrum().size() + sd_noise_decoy_.getDeconvolvedSpectrum().size());

        for (auto& pg : sd_charge_decoy_.getDeconvolvedSpectrum())
        {
          if (pg.getQscore() < qscore_threshold_for_dummy)
          {
            continue;
          }
          dummy_deconvolved_spectrum.push_back(pg);
        }

        for (auto& pg : sd_isotope_decoy_.getDeconvolvedSpectrum())
        {
          if (pg.getQscore() < qscore_threshold_for_dummy)
          {
            continue;
          }
          dummy_deconvolved_spectrum.push_back(pg);
        }

        for (auto& pg : sd_noise_decoy_.getDeconvolvedSpectrum())
        {
          if (pg.getQscore() < qscore_threshold_for_dummy)
          {
            continue;
          }
          dummy_deconvolved_spectrum.push_back(pg);
        }

        deconvolved_spectrum.sort();
        dummy_deconvolved_spectrum.sort();

        deconvolved_spectra.push_back(dummy_deconvolved_spectrum);
      }
      deconvolved_spectra.push_back(deconvolved_spectrum);
    }
    endProgress();
    return 0;
  }

  void FLASHDeconvAlgorithm::setTargetPrecursorCharge_(DeconvolvedSpectrum& deconvolved_spectrum, const MSSpectrum& it)
  {
    auto precursor = it.getPrecursors()[0];
    double target_precursor_mass = (precursor.getMZ() - FLASHDeconvHelperStructs::getChargeMass(target_precursor_charge_ > 0)) * std::abs(target_precursor_charge_);
    // precursor.setCharge(target_precursor_charge);
    PeakGroup precursorPeakGroup(1, std::abs(target_precursor_charge_), target_precursor_charge_ > 0);
    precursorPeakGroup.push_back(FLASHDeconvHelperStructs::LogMzPeak());
    precursorPeakGroup.setMonoisotopicMass(target_precursor_mass);
    precursorPeakGroup.setSNR(1.0);

    precursorPeakGroup.setChargeSNR(std::abs(target_precursor_charge_), 1.0);
    precursorPeakGroup.setQscore(1.0);
    deconvolved_spectrum.setPrecursor(precursor);
    deconvolved_spectrum.setPrecursorPeakGroup(precursorPeakGroup);
  }

  const FLASHDeconvHelperStructs::PrecalculatedAveragine& FLASHDeconvAlgorithm::getAveragine()
  {
    return sd_.getAveragine();
  }

  void FLASHDeconvAlgorithm::run(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHDeconvHelperStructs::MassFeature>& deconvolved_features)
  {
    // initialize
    precursor_map_for_ida_ = FLASHIda::parseFLASHIdaLog(ida_log_file_); // ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color

    if (!precursor_map_for_ida_.empty())
    {
      precursor_SNR_threshold_ = .0;
    }

    uint current_max_ms_level = 0;
    uint current_min_ms_level = 1000;

    // scan the dataset and extract necessary basic information.
    if (scan_map_(map) != 0)
      return;

    // FLASHIda information reading and Parameter parsing
    if (!ida_log_file_.empty())
    {
      preceding_MS1_count_ = 50; // if FLASHIda log file exists, keep up to 50 survey scans.
    }

    sd_ = SpectralDeconvolution();

    Param sd_param = param_.copy("Deconvolution", true);
    tols_ = sd_param.getValue("tol");

    if (target_precursor_charge_ != 0)
    {
      sd_param.setValue("max_charge", target_precursor_charge_);
      if (target_precursor_mz_ != 0)
      {
        sd_param.setValue("max_mass", std::abs(target_precursor_charge_) * (target_precursor_mz_ - FLASHDeconvHelperStructs::getChargeMass(target_precursor_charge_ > 0)));
      }
    }

    // merge spectra if the merging option is turned on (> 0)
    mergeSpectra_(map);

    sd_.setParameters(sd_param);
    sd_.calculateAveragine(use_RNA_averagine_);
    auto avg = sd_.getAveragine();

    if (report_decoy_)
    {
      sd_charge_decoy_.setParameters(sd_param);
      sd_charge_decoy_.setAveragine(avg);
      sd_charge_decoy_.setTargetDummyType(PeakGroup::TargetDecoyType::charge_decoy, sd_.getDeconvolvedSpectrum()); // charge

      sd_noise_decoy_.setParameters(sd_param);
      sd_noise_decoy_.setAveragine(avg);
      sd_noise_decoy_.setTargetDummyType(PeakGroup::TargetDecoyType::noise_decoy, sd_.getDeconvolvedSpectrum()); // noise

      sd_isotope_decoy_.setParameters(sd_param);
      sd_isotope_decoy_.setAveragine(avg);
      sd_isotope_decoy_.setTargetDummyType(PeakGroup::TargetDecoyType::isotope_decoy, sd_.getDeconvolvedSpectrum()); // isotope
    }

    setLogType(CMD);
    startProgress(0, (SignedSize)map.size(), "running FLASHDeconv");
    deconvolved_spectra.reserve(map.size() * 4);

    // run FLASHDeconv here and get deconvolved spectra
    if (runFD_(map, deconvolved_spectra) != 0)
      return;
    // feature tracing here and update FeatureQScores
    runFeatureFinding_(deconvolved_spectra, deconvolved_features);

    updatePrecursorQScores_(deconvolved_spectra);
    Qvalue::updatePeakGroupQvalues(deconvolved_spectra);

    TopDownIsobaricQuantifier quantifier;
    Param quant_param = param_.copy("IsobaricQuantification", true);
    quantifier.setParameters(quant_param);
    // Isobaric quant run
    quantifier.quantify(map, deconvolved_spectra, deconvolved_features);
  }

  void FLASHDeconvAlgorithm::updatePrecursorQScores_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
  {
    // update precursor feature QScores and qvalues
    std::map<int, DeconvolvedSpectrum> scan_fullscan;

    for (auto& dspec : deconvolved_spectra)
    {
      int scan = dspec.getScanNumber();
      scan_fullscan[scan] = dspec;
    }

    for (auto& dspec : deconvolved_spectra)
    {
      auto& precursor_pg = dspec.getPrecursorPeakGroup();
      if (precursor_pg.empty())
        continue;

      int pscan = precursor_pg.getScanNumber();
      if (scan_fullscan.find(pscan) == scan_fullscan.end())
        continue;

      auto fullscan = scan_fullscan[pscan];

      auto iter = std::lower_bound(fullscan.begin(), fullscan.end(), precursor_pg);
      if (precursor_pg.getMonoMass() == iter->getMonoMass())
      {
        precursor_pg.setFeatureIndex(iter->getFeatureIndex());
        precursor_pg.setQscore(iter->getQscore());
        if (iter->getFeatureIndex() > 0)
          precursor_pg.setFeatureQscore(iter->getFeatureQscore());
      }
      else
      {
        // std::cout<<iter->getMonoMass() << " " <<  precursor_pg.getMonoMass() << std::endl;
        precursor_pg.setFeatureIndex(0);
      }
    }
  }

  void FLASHDeconvAlgorithm::runFeatureFinding_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHDeconvHelperStructs::MassFeature>& deconvolved_features)
  {
    if (merge_spec_ == 2)
      return;

    auto mass_tracer = MassFeatureTrace();
    auto decoy_mass_tracer = MassFeatureTrace();

    Param mf_param = param_.copy("FeatureTracing", true);

    mass_tracer.setParameters(mf_param); // maybe go to set param
    decoy_mass_tracer.setParameters(mf_param);

    // deconvolved_spectra
    deconvolved_features = mass_tracer.findFeatures(sd_.getAveragine(), deconvolved_spectra, 1, false);

    if (report_decoy_)
    {
      // remove the decoy peak groups overlapping with target features.
      for (auto& dspec : deconvolved_spectra)
      {
        if (!dspec.isDecoy())
          continue;

        double rt = dspec.getOriginalSpectrum().getRT();
        auto filtered_dspec = dspec;
        filtered_dspec.clear();

        boost::dynamic_bitset<> remove_index;
        remove_index.resize(dspec.size());

        for (auto& feature : deconvolved_features)
        {
          double frt_begin = feature.mt.begin()->getRT();
          double frt_end = feature.mt.rbegin()->getRT();
          if (rt < frt_begin || rt > frt_end)
            continue;

          double f_mass = feature.mt.getCentroidMZ() + feature.iso_offset * Constants::ISOTOPE_MASSDIFF_55K_U;
          for (uint i = 0; i < dspec.size(); i++)
          {
            if (remove_index[i])
              continue;
            auto pg = dspec[i];
            double pg_mass = pg.getMonoMass();
            if (abs(f_mass - pg_mass) < 1e-6 * tols_[0] * pg_mass)
            {
              remove_index[i] = true;
            }
          }
        }

        for (uint i = 0; i < dspec.size(); i++)
        {
          if (remove_index[i])
            continue;
          filtered_dspec.push_back(dspec[i]);
        }
        dspec = filtered_dspec;
      }
      auto decoy_deconvolved_features = decoy_mass_tracer.findFeatures(sd_.getAveragine(), deconvolved_spectra, 1, true);
      deconvolved_features.insert(deconvolved_features.end(), decoy_deconvolved_features.begin(), decoy_deconvolved_features.end());
    }
  }
} // namespace OpenMS