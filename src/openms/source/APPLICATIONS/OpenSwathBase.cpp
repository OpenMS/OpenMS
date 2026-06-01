// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest, Justin Sing$
// $Authors: Hannes Roest, Justin Sing$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/OpenSwathBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/TargetedDataFileLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

#include <algorithm>
#include <iostream>
#include <utility>

using namespace std;

namespace OpenMS
{
  // Helper function to check if ion mobility data is in CCS format
  // CCS data is not directly compatible with OpenSwath which expects 1/K0 (inverse reduced ion mobility)
  static void warnIfCCSData(const std::vector<OpenSwath::SwathMap>& swath_maps)
  {
    for (const auto& swath_map : swath_maps)
    {
      if (swath_map.sptr && swath_map.sptr->getNrSpectra() > 0)
      {
        // Get the first spectrum and check its data arrays for CCS indicators
        auto spectrum = swath_map.sptr->getSpectrumById(0);
        if (spectrum)
        {
          for (const auto& arr : spectrum->getDataArrays())
          {
            // Check for CCS CV term (MS:1002954) or square angstrom unit (UO:0000324)
            if (arr->description.find("MS:1002954") != std::string::npos ||
                arr->description.find("UO:0000324") != std::string::npos ||
                arr->description.find("collision cross section") != std::string::npos)
            {
              OPENMS_LOG_WARN << "Warning: Ion mobility data appears to be in CCS (Collisional Cross Section) format. "
                              << "OpenSwath expects ion mobility in 1/K0 (inverse reduced ion mobility) units. "
                              << "Results may be incorrect if the spectral library also uses 1/K0 units. "
                              << "Consider converting the data or using a CCS-based spectral library." << std::endl;
              return; // Only warn once
            }
          }
        }
        return; // Only check first map with spectra
      }
    }
  }

  TOPPOpenSwathBase::TOPPOpenSwathBase(String name, String description, bool official, const std::vector<Citation>& citations) :
    TOPPBase(name, description, official, citations)
  {
  }

  TOPPOpenSwathBase::~TOPPOpenSwathBase() = default;

  // Private

  void TOPPOpenSwathBase::loadSwathFiles_(const StringList& file_list,
                       const bool split_file,
                       const String& tmp,
                       const String& readoptions,
                       std::shared_ptr<ExperimentalSettings > & exp_meta,
                       std::vector< OpenSwath::SwathMap > & swath_maps,
                       std::vector<String> & swath_map_sources,
                       Interfaces::IMSDataConsumer* plugin_consumer)
  {
    SwathFile swath_file;
    swath_file.setLogType(log_type_);

    if (split_file)
    {
      // Treat the list as parts of a single split mzML
      // (one file per SWATH window). Keep for backward compatibility.
      swath_maps = swath_file.loadSplit(file_list, tmp, exp_meta, readoptions);
      // assign the same source filename (first entry) for these maps
      swath_map_sources.clear();
      for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
    }
    else if (file_list.size() > 1)
    {
      // Treat the list as multiple independent experiment files
      // (e.g. multiple SRM/MRM / experiment mzMLs). Load each file individually
      // and append their swath maps so they can be processed in a single run.
      for (const auto & f : file_list)
      {
        FileTypes::Type in_file_type = FileHandler::getTypeByFileName(f);
        if (in_file_type == FileTypes::MZML)
        {
          // Dispatch to targeted loader which will choose SRM/MRM or SWATH path
          auto maps = TargetedDataFileLoader::loadFile(f, tmp, exp_meta, readoptions, plugin_consumer);
          for (auto & m : maps)
          {
            swath_maps.push_back(m);
            swath_map_sources.push_back(f);
          }
        }
        else if (in_file_type == FileTypes::MZXML)
        {
          auto maps = swath_file.loadMzXML(f, tmp, exp_meta, readoptions);
          for (auto & m : maps)
          {
            swath_maps.push_back(m);
            swath_map_sources.push_back(f);
          }
        }
        else if (in_file_type == FileTypes::SQMASS)
        {
          auto maps = swath_file.loadSqMass(f, exp_meta);
          for (auto & m : maps)
          {
            swath_maps.push_back(m);
            swath_map_sources.push_back(f);
          }
        }
#ifdef WITH_OPENTIMS
        else if (in_file_type == FileTypes::BRUKER_TDF)
        {
          auto maps = swath_file.loadBrukerTdf(f, tmp, exp_meta, readoptions);
          for (auto & m : maps)
          {
            swath_maps.push_back(m);
            swath_map_sources.push_back(f);
          }
        }
#endif
#ifdef WITH_THERMO_RAW
        else if (in_file_type == FileTypes::RAW)
        {
          // Load .raw fully via ThermoRawFile (in-process), then build SwathMaps
          // directly from the in-memory MSExperiment — no temp mzML round-trip.
          auto exp_raw = std::make_shared<PeakMap>();
          FileHandler().loadExperiment(f, *exp_raw, {FileTypes::RAW});
          auto maps = swath_file.loadFromMSExperiment(exp_raw, tmp, exp_meta, readoptions);
          for (auto & m : maps)
          {
            swath_maps.push_back(m);
            swath_map_sources.push_back(f);
          }
        }
#endif
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Input file needs to have a supported extension "
                                           "(mzML, mzXML, sqMass"
#ifdef WITH_OPENTIMS
                                           ", .d (Bruker TDF)"
#endif
#ifdef WITH_THERMO_RAW
                                           ", raw (Thermo)"
#endif
                                           ")");
        }
      }
    }
    else
    {
      FileTypes::Type in_file_type = FileHandler::getTypeByFileName(file_list[0]);
      if (in_file_type == FileTypes::MZML)
      {
        // Dispatch to targeted loader which will choose SRM/MRM or SWATH path
        swath_maps = TargetedDataFileLoader::loadFile(file_list[0], tmp, exp_meta, readoptions, plugin_consumer);
        // single-file: all maps originate from file_list[0]
        swath_map_sources.clear();
        for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
      }
      else if (in_file_type == FileTypes::MZXML)
      {
        swath_maps = swath_file.loadMzXML(file_list[0], tmp, exp_meta, readoptions);
        swath_map_sources.clear();
        for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
      }
      else if (in_file_type == FileTypes::SQMASS)
      {
        swath_maps = swath_file.loadSqMass(file_list[0], exp_meta);
        swath_map_sources.clear();
        for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
      }
#ifdef WITH_OPENTIMS
      else if (in_file_type == FileTypes::BRUKER_TDF)
      {
        swath_maps = swath_file.loadBrukerTdf(file_list[0], tmp, exp_meta, readoptions);
        swath_map_sources.clear();
        for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
      }
#endif
#ifdef WITH_THERMO_RAW
      else if (in_file_type == FileTypes::RAW)
      {
        // Load .raw fully via ThermoRawFile (in-process), then build SwathMaps
        // directly from the in-memory MSExperiment — no temp mzML round-trip.
        auto exp_raw = std::make_shared<PeakMap>();
        FileHandler().loadExperiment(file_list[0], *exp_raw, {FileTypes::RAW});
        swath_maps = swath_file.loadFromMSExperiment(exp_raw, tmp, exp_meta, readoptions);
        swath_map_sources.clear();
        for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
      }
#endif
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Input file needs to have a supported extension "
                                         "(mzML, mzXML, sqMass"
#ifdef WITH_OPENTIMS
                                         ", .d (Bruker TDF)"
#endif
#ifdef WITH_THERMO_RAW
                                         ", raw (Thermo)"
#endif
                                         ")");
      }
    }
  }

  // Protected

  bool TOPPOpenSwathBase::loadSwathFiles(const StringList& file_list,
                      std::shared_ptr<ExperimentalSettings >& exp_meta,
                      std::vector< OpenSwath::SwathMap >& swath_maps,
                      std::vector<String> & swath_map_sources,
                      const bool split_file,
                      const String& tmp,
                      const String& readoptions,
                      const String& swath_windows_file,
                      const double min_upper_edge_dist,
                      const bool force,
                      const bool sort_swath_maps,
                      const bool prm,
                      Interfaces::IMSDataConsumer* plugin_consumer)
  {
    // (i) Load files
    loadSwathFiles_(file_list, split_file, tmp, readoptions, exp_meta, swath_maps, swath_map_sources, plugin_consumer);

    // (ii) Check consistency: all input files must be of the same type (SRM/MRM or DIA/PRM)
    if (!swath_maps.empty())
    {
      bool first_has_spectra = (swath_maps[0].sptr->getNrSpectra() > 0);
      String first_file_type = first_has_spectra ? "DIA/PRM (spectra-based)" : "SRM/MRM (chromatogram-only)";
      for (Size i = 1; i < swath_maps.size(); ++i)
      {
        bool current_has_spectra = (swath_maps[i].sptr->getNrSpectra() > 0);
        String current_file_type = current_has_spectra ? "DIA/PRM (spectra-based)" : "SRM/MRM (chromatogram-only)";
        if (current_has_spectra != first_has_spectra)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Mixed input file types detected. All input files must be of the same type. "
                                           "First file '" + swath_map_sources[0] + "' is " + first_file_type + ", "
                                           "but file '" + swath_map_sources[i] + "' is " + current_file_type + ". "
                                           "Please ensure all input files are either chromatogram-only (SRM/MRM) or spectra-based (DIA/PRM).");
        }
      }
    }

    // (iii) Allow the user to specify the SWATH windows
    if (!swath_windows_file.empty())
    {
      SwathWindowLoader::annotateSwathMapsFromFile(swath_windows_file, swath_maps, sort_swath_maps, force);
    }

    for (Size i = 0; i < swath_maps.size(); i++)
    {
      OPENMS_LOG_DEBUG << "Found swath map " << i
                       << " with lower " << swath_maps[i].lower
                       << " and upper " << swath_maps[i].upper
                       << " and im Lower bounds of " << swath_maps[i].imLower
                       << " and im Upper bounds of " << swath_maps[i].imUpper
                       << " and " << swath_maps[i].sptr->getNrSpectra()
                       << " spectra." << std::endl;
    }

    // (iii) Sanity check: there should be no overlap between the windows:
    std::vector<std::pair<double, double>> sw_windows;
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      if (!swath_maps[i].ms1)
      {
        sw_windows.push_back(std::make_pair(swath_maps[i].lower, swath_maps[i].upper));
      }
    }
    // sort by lower bound (first entry in pair)
    std::sort(sw_windows.begin(), sw_windows.end());

    // Auto-detect PASEF/IM data: check if any non-MS1 swath map has ion mobility bounds
    bool has_im_windows = std::any_of(swath_maps.begin(), swath_maps.end(),
      [](const OpenSwath::SwathMap& m) { return !m.ms1 && m.imLower >= 0 && m.imUpper >= 0; });

    for (Size i = 1; i < sw_windows.size(); i++)
    {
      double lower_map_end = sw_windows[i-1].second - min_upper_edge_dist;
      double upper_map_start = sw_windows[i].first;
      OPENMS_LOG_DEBUG << "Extraction will go up to " << lower_map_end << " and continue at " << upper_map_start << std::endl;

      if (prm) {continue;} // skip next step as expect them to overlap and have gaps...

      if (upper_map_start - lower_map_end > 0.01)
      {
        OPENMS_LOG_WARN << "Extraction will have a gap between " << lower_map_end << " and " << upper_map_start << std::endl;
        if (!force)
        {
          OPENMS_LOG_ERROR << "Extraction windows have a gap. Will abort (override with -force)" << std::endl;
          return false;
        }
      }

      if (has_im_windows) {continue;} // IM data present: m/z overlap expected across IM dimension

      if (lower_map_end - upper_map_start > 0.01)
      {
        OPENMS_LOG_WARN << "Extraction will overlap between " << lower_map_end << " and " << upper_map_start << "!\n"
                        << "This will lead to multiple extraction of the transitions in the overlapping region "
                        << "which will lead to duplicated output. It is very unlikely that you want this." << "\n"
                        << "Please fix this by providing an appropriate extraction file with -swath_windows_file" << std::endl;
        if (!force)
        {
          OPENMS_LOG_ERROR << "Extraction windows overlap. Will abort (override with -force)" << std::endl;
          return false;
        }
      }
    }

    // Check if ion mobility data is in CCS format and warn user
    warnIfCCSData(swath_maps);

    return true;
  }

  void TOPPOpenSwathBase::prepareChromOutput(Interfaces::IMSDataConsumer ** chromatogramConsumer,
                          const std::shared_ptr<ExperimentalSettings>& exp_meta,
                          const OpenSwath::LightTargetedExperiment& transition_exp,
                          const String& out_chrom,
                          const UInt64 run_id,
                          const String& source_file)
  {
    if (!out_chrom.empty())
    {
      const FileTypes::Type out_chrom_type = FileHandler::getType(out_chrom);
      if (out_chrom_type == FileTypes::SQMASS)
      {
        bool full_meta = false; // can lead to very large files in memory
        bool lossy_compression = true;
        *chromatogramConsumer = new MSDataSqlConsumer(out_chrom, run_id, 500, full_meta, lossy_compression);
      }
      else if (out_chrom_type == FileTypes::CHROMPARQUET)
      {
        auto * chromConsumer = new MSChromatogramParquetConsumer(out_chrom, run_id, source_file, transition_exp);
        Size expected_chromatograms = transition_exp.transitions.size();
        chromConsumer->setExpectedSize(0, expected_chromatograms);
        *chromatogramConsumer = chromConsumer;
      }
      else
      {
        PlainMSDataWritingConsumer * chromConsumer = new PlainMSDataWritingConsumer(out_chrom);
        int expected_chromatograms = transition_exp.transitions.size();
        chromConsumer->setExpectedSize(0, expected_chromatograms);
        chromConsumer->setExperimentalSettings(*exp_meta);
        chromConsumer->getOptions().setWriteIndex(true);  // ensure that we write the index
        chromConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));

        // prepare data structures for lossy compression
        MSNumpressCoder::NumpressConfig npconfig_mz;
        MSNumpressCoder::NumpressConfig npconfig_int;
        npconfig_mz.estimate_fixed_point = true; // critical
        npconfig_int.estimate_fixed_point = true; // critical
        npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
        npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
        npconfig_mz.setCompression("linear");
        npconfig_int.setCompression("slof");
        npconfig_mz.linear_fp_mass_acc = 0.05; // set the desired RT accuracy in seconds

        chromConsumer->getOptions().setNumpressConfigurationMassTime(npconfig_mz);
        chromConsumer->getOptions().setNumpressConfigurationIntensity(npconfig_int);
        chromConsumer->getOptions().setCompression(true);

        *chromatogramConsumer = chromConsumer;
      }
    }
    else
    {
      *chromatogramConsumer = new NoopMSDataWritingConsumer("");
    }
  }


  void TOPPOpenSwathBase::prepareMobilogramOutput(std::unique_ptr<MobilogramParquetConsumer>& mobilogramConsumer,
                                                  const std::shared_ptr<ExperimentalSettings>& /*exp_meta*/,
                                                  const OpenSwath::LightTargetedExperiment& transition_exp,
                                                  const String& out_mobilogram,
                                                  const UInt64 run_id,
                                                  const String& source_file)
  {
    if (!out_mobilogram.empty())
    {
      const FileTypes::Type out_mob_type = FileHandler::getType(out_mobilogram);
      if (out_mob_type == FileTypes::MOBILPARQUET)
      {
        mobilogramConsumer = std::make_unique<MobilogramParquetConsumer>(out_mobilogram, run_id, source_file, transition_exp);
        // estimate expected mobilograms from transitions
        Size expected = transition_exp.transitions.size();
        mobilogramConsumer->setExpectedSize(expected);
      }
      else
      {
        // Unknown mobilogram format: currently only MOBILPARQUET supported
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Unsupported mobilogram output type: " + FileTypes::typeToName(out_mob_type));
      }
    }
    else
    {
      mobilogramConsumer.reset();
    }
  }

  OpenSwath::LightTargetedExperiment TOPPOpenSwathBase::loadTransitionList(const FileTypes::Type& tr_type,
                                                        const String& tr_file,
                                                        const Param& tsv_reader_param)
  {
    OpenSwath::LightTargetedExperiment transition_exp;
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    if (tr_type == FileTypes::TRAML)
    {
      progresslogger.startProgress(0, 1, "Load TraML file");
      TargetedExperiment targeted_exp;
      FileHandler().loadTransitions(tr_file, targeted_exp, {FileTypes::TRAML});
      OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, transition_exp);
      progresslogger.endProgress();
    }
    else if (tr_type == FileTypes::PQP)
    {
      progresslogger.startProgress(0, 1, "Load PQP file");
      TransitionPQPFile().convertPQPToTargetedExperiment(tr_file.c_str(), transition_exp);
      progresslogger.endProgress();
    }
    else if (tr_type == FileTypes::TSV)
    {
      progresslogger.startProgress(0, 1, "Load TSV file");
      TransitionTSVFile tsv_reader;
      tsv_reader.setParameters(tsv_reader_param);
      tsv_reader.convertTSVToTargetedExperiment(tr_file.c_str(), tr_type, transition_exp);
      progresslogger.endProgress();
    }
    else if (tr_type == FileTypes::OSWPQ)
    {
      progresslogger.startProgress(0, 1, "Load Parquet library");
      TransitionParquetFile().convertParquetToTargetedExperiment(tr_file, transition_exp);
      progresslogger.endProgress();
    }
    else
    {
      OPENMS_LOG_ERROR << "Provide valid TraML, TSV, PQP or OSWPQ transition file." << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Need to provide valid input file.");
    }
    return transition_exp;
  }
} //  end NS OpenMS