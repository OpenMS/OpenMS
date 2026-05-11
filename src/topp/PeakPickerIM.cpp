// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeakPickerIM PeakPickerIM

@brief A tool for peak detection in the ion mobility dimension for mzML and Bruker .d files.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; PeakPickerIM &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
</tr>
</table>
</center>

This tool applies peak picking in the ion mobility dimension to raw LC-IMS-MS data.
The input file can be an mzML file containing ion mobility data in concatenated format
(where each spectrum contains an ion mobility float data array) or a Bruker TimsTOF .d
directory (requires OpenMS built with WITH_OPENTIMS).

Three peak picking methods are available:
- @b mobilogram: Picks peaks along the ion mobility dimension using a peak picker.
- @b cluster: Clusters peaks in the ion mobility dimension.
- @b traces: Picks peaks using ion mobility elution profiles.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PeakPickerIM.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PeakPickerIM.html

For the parameters of the algorithm section see the algorithm documentation: @ref OpenMS::PeakPickerIM "PeakPickerIM"

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPickerIM : public TOPPBase
{
public:
  TOPPPeakPickerIM() :
      TOPPBase("PeakPickerIM", "Applies PeakPickerIM to an mzML or Bruker .d file", false)
  {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (mzML or Bruker .d)");
    setValidFormats_("in", { "mzML",
#ifdef WITH_OPENTIMS
      "d",
#endif
    });

    registerOutputFile_("out", "<file>", "", "Output mzML file");
    setValidFormats_("out", { "mzML" });

    registerStringOption_("processOption", "<name>", "inmemory",
                          "Whether to load all data and process them in-memory or process on-the-fly (lowmemory) without loading the whole file into memory first",
                          false, true);
    setValidStrings_("processOption", { "inmemory", "lowmemory" } );

    registerStringOption_("method", "<name>", "mobilogram",
                          "Method to pick peaks in IM dimension", false, true);
    setValidStrings_("method", { "mobilogram", "cluster", "traces" } );

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters for PeakPickerIM (organized into pickIMTraces, pickIMCluster, pickIMElutionProfiles).");

#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerStringOption_("bruker:export_mode", "<mode>", "frame", "Export mode: 'auto' detects DDA/DIA acquisition type, "
      "'frame' returns raw 4D frames without signal processing.", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "frame"});
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0, "m/z recalibration tolerance (0 = library default)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false", "Enable m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
    registerDoubleOption_("bruker:ms1_centroid_mz_ppm", "<float>", 10.0,
      "MS1 m/z linking tolerance in ppm. HillBased default 10 ppm tuned for "
      "detector-centroided TIMS-PASEF MS1; Greedy2D additionally requires "
      "ms1_centroid_im_pct > 0. Set to 0 to disable MS1 centroiding.",
      false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 0.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. Both this and ms1_centroid_mz_ppm "
      "must be > 0 to enable. Suggested value: 3.0.", false, true);
    setMinFloat_("bruker:ms1_centroid_im_pct", 0.0);
    registerIntOption_("bruker:dia_ms2_n_neighbors", "<int>", 0,
      "DIA MS2 frame aggregation: 0 = raw per-frame export, 1 = 3-frame sum, 2 = 5-frame sum. "
      "Switches the entire DIA-MS2 export pipeline regardless of ms2_centroid_algo.",
      false, true);
    setMinInt_("bruker:dia_ms2_n_neighbors", 0);
    registerIntOption_("bruker:dia_ms2_min_support", "<int>", 1,
      "DIA MS2 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point "
      "(center cell excluded from count). Applied after frame aggregation. Only effective when "
      "dia_ms2_n_neighbors > 0. Set to 0 to disable denoising (useful for pure centroiding "
      "without noise filtering).", false, true);
    setMinInt_("bruker:dia_ms2_min_support", 0);
    registerStringOption_("bruker:dia_ms2_centroid", "<toggle>", "false",
      "Apply 2D Gaussian smoothing + local maxima peak picking to the denoised DIA MS2 grid. "
      "Produces IM_CENTROIDED spectra with sub-bin (m/z, IM) precision. Only effective when dia_ms2_n_neighbors > 0.", false, true);
    setValidStrings_("bruker:dia_ms2_centroid", {"true", "false"});

    // Hill-based centroiding (IM-axis trace linking + valley splitting).
    registerStringOption_("bruker:ms1_centroid_algo", "<algo>", "off",
      "MS1 centroiding algorithm. 'off' = no IM-axis centroiding. 'greedy2d' = legacy 2D box "
      "clustering using ms1_centroid_mz_ppm/pct. 'hillbased' = IM-axis hill detection.",
      false, true);
    setValidStrings_("bruker:ms1_centroid_algo", {"off", "greedy2d", "hillbased"});
    registerStringOption_("bruker:ms2_centroid_algo", "<algo>", "off",
      "MS2 centroiding algorithm. Takes precedence over dia_ms2_centroid.",
      false, true);
    setValidStrings_("bruker:ms2_centroid_algo", {"off", "greedy2d", "hillbased"});
    registerDoubleOption_("bruker:ms2_centroid_mz_ppm", "<float>", 20.0,
      "HillBased MS2 m/z linking tolerance in ppm. Default 20.0 is DIA-PASEF-tuned. "
      "Set to 0 to refuse to run HillBased MS2 (the algo helper falls back to Off).",
      false, true);
    setMinFloat_("bruker:ms2_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:centroid_valley_factor", "<float>", 1.3,
      "HillBased: hill valley factor (hvf). Smaller = more aggressive splitting.",
      false, true);
    setMinFloat_("bruker:centroid_valley_factor", 1.0);
    registerIntOption_("bruker:ms1_centroid_min_hill_length", "<int>", 1,
      "HillBased MS1: minimum number of IM scans a hill must span. Default 1 keeps single-IM-scan ions.",
      false, true);
    setMinInt_("bruker:ms1_centroid_min_hill_length", 1);
    registerIntOption_("bruker:ms2_centroid_min_hill_length", "<int>", 2,
      "HillBased MS2: minimum number of IM scans a hill must span. Default 2 is "
      "DIA-PASEF-tuned. DDA-PASEF users should override to 1 (most DDA fragments "
      "are seen in only one IM scan; min=2 drops ~93% of DDA peaks).",
      false, true);
    setMinInt_("bruker:ms2_centroid_min_hill_length", 1);
    registerIntOption_("bruker:centroid_max_scan_gap", "<int>", 0,
      "HillBased: max consecutive empty IM scans a hill may bridge (0 = strict).",
      false, true);
    setMinInt_("bruker:centroid_max_scan_gap", 0);
    registerStringOption_("bruker:expose_hill_bounds", "<toggle>", "false",
      "HillBased: attach hill bounding-box arrays per centroided spectrum for visual QC.",
      false, true);
    setValidStrings_("bruker:expose_hill_bounds", {"true", "false"});
#endif
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      OpenMS::PeakPickerIM picker_defaults;
      Param p = picker_defaults.getDefaults();
      Param combined;
      combined.insert("pickIMTraces:",         p.copy("pickIMTraces:", true));
      combined.insert("pickIMCluster:",        p.copy("pickIMCluster:", true));
      combined.insert("pickIMElutionProfiles:",p.copy("pickIMElutionProfiles:", true));
      return combined;
    }
    return Param();
  }

#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate = (getStringOption_("bruker:calibrate") == "true");
    String mode = getStringOption_("bruker:export_mode");
    if (mode == "frame") c.export_mode = BrukerTimsFile::Config::FRAME;
    else c.export_mode = BrukerTimsFile::Config::AUTO;
    c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_im_pct"));
    c.dia_ms2_n_neighbors = getIntOption_("bruker:dia_ms2_n_neighbors");
    c.dia_ms2_min_support = getIntOption_("bruker:dia_ms2_min_support");
    c.dia_ms2_centroid = (getStringOption_("bruker:dia_ms2_centroid") == "true");

    using CA = BrukerTimsFile::Config::CentroidAlgo;
    auto parse_algo = [](const String& s) {
      if (s == "greedy2d")  return CA::GREEDY2D;
      if (s == "hillbased") return CA::HILL_BASED;
      return CA::OFF;
    };
    c.ms1_centroid_algo            = parse_algo(getStringOption_("bruker:ms1_centroid_algo"));
    c.ms2_centroid_algo            = parse_algo(getStringOption_("bruker:ms2_centroid_algo"));
    c.ms2_centroid_mz_ppm          = static_cast<float>(getDoubleOption_("bruker:ms2_centroid_mz_ppm"));
    c.centroid_valley_factor       = getDoubleOption_("bruker:centroid_valley_factor");
    c.ms1_centroid_min_hill_length = static_cast<Size>(getIntOption_("bruker:ms1_centroid_min_hill_length"));
    c.ms2_centroid_min_hill_length = static_cast<Size>(getIntOption_("bruker:ms2_centroid_min_hill_length"));
    c.centroid_max_scan_gap        = static_cast<Size>(getIntOption_("bruker:centroid_max_scan_gap"));
    c.expose_hill_bounds           = (getStringOption_("bruker:expose_hill_bounds") == "true");
    return c;
  }
#endif

  // -------------------- Low-memory consumer --------------------
  class Consumer : public MSDataWritingConsumer
  {
  public:
    Consumer(String filename, const String& method, const PeakPickerIM& pp) :
        MSDataWritingConsumer(std::move(filename)), pp_(pp), method_(method) {}

    void processSpectrum_(MapType::SpectrumType& spectrum) override
    {
      if (method_ == "mobilogram")
      {
        pp_.pickIMTraces(spectrum);
      }
      else if (method_ == "cluster")
      {
        pp_.pickIMCluster(spectrum);
      }
      else if (method_ == "traces")
      {
        pp_.pickIMElutionProfiles(spectrum);
      }
    }

    void processChromatogram_(MapType::ChromatogramType&) override {}

  private:
    PeakPickerIM pp_;
    String method_;
  };

  // -------------------- Format detection consumer (reads first MS1 spectrum only) --------------------
  class FormatDetector : public Interfaces::IMSDataConsumer
  {
  public:
    IMFormat detected_format = IMFormat::NONE;

    // Exception to abort after first spectrum (efficient early exit)
    struct FirstSpectrumRead : std::exception {};

    void consumeSpectrum(SpectrumType& s) override
    {
      if (s.getMSLevel() != 1) return; // Only check MS1 spectra (consistent with in-memory path)
      detected_format = IMTypes::determineIMFormat(s);
      throw FirstSpectrumRead(); // Abort after first MS1 spectrum
    }
    void consumeChromatogram(ChromatogramType&) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
    void setExpectedSize(size_t, size_t) override {}
  };

  // -------------------- Passthrough consumer (copies without processing) --------------------
  class PassthroughConsumer : public MSDataWritingConsumer
  {
  public:
    PassthroughConsumer(const String& filename) : MSDataWritingConsumer(filename) {}
    void processSpectrum_(MapType::SpectrumType&) override {} // No processing
    void processChromatogram_(MapType::ChromatogramType&) override {}
  };

  // -------------------- Helper for low-memory path --------------------
  ExitCodes doLowMemAlgorithm(const String& method, const PeakPickerIM& pp,
                              const String& input_file, const String& output_file)
  {
    MzMLFile mzml;
    mzml.setLogType(log_type_);

    // Step 1: Detect IMFormat by reading only the first spectrum (minimal I/O)
    IMFormat im_format = IMFormat::NONE;
    {
      FormatDetector detector;
      try
      {
        mzml.transform(input_file, &detector);
        // If we reach here, file has no MS1 spectra - format stays NONE
      }
      catch (const FormatDetector::FirstSpectrumRead&)
      {
        im_format = detector.detected_format;
      }
    }

    // Step 2: Validate format
    if (im_format == IMFormat::IM_SPECTRUM)
    {
      OPENMS_LOG_ERROR << "Error: Input data has single drift time per spectrum (IM_SPECTRUM format). "
                       << "PeakPickerIM requires per-peak IM arrays (IM_PEAK format)." << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    if (im_format == IMFormat::NONE)
    {
      OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                      << "No peak picking will be performed." << std::endl;
      // Pass through unchanged
      PassthroughConsumer passthrough(output_file);
      mzml.transform(input_file, &passthrough);
      return EXECUTION_OK;
    }

    // Step 3: Proceed with streaming processing
    Consumer pp_consumer(output_file, method, pp);
    pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));
    mzml.transform(input_file, &pp_consumer);
    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    const String input_file  = getStringOption_("in");
    const String output_file = getStringOption_("out");
    const String process_opt = getStringOption_("processOption");
    const String method      = getStringOption_("method");

    // Collect algorithm parameters from 'algorithm:' We strip and pass the remaining keys directly to PeakPickerIM.
    Param algo = getParam_().copy("algorithm:",true);

    PeakPickerIM picker;
    picker.setParameters(algo);

    // Detect input file type
    FileTypes::Type in_type = FileHandler::getType(input_file);

#ifdef WITH_OPENTIMS
    if (in_type == FileTypes::BRUKER_TDF)
    {
      if (process_opt == "lowmemory")
      {
        OPENMS_LOG_WARN << "Warning: 'lowmemory' processing is not yet supported for Bruker .d files. "
                        << "Data will be loaded fully into memory." << std::endl;
      }

      auto bruker_config = getBrukerConfig_();
      BrukerTimsFile tims_file;
      tims_file.setLogType(log_type_);

      PeakMap exp;
      tims_file.load(input_file, exp, bruker_config);

      // If built-in IM centroiding was enabled, BrukerTimsFile already produced
      // IM_CENTROIDED spectra — skip PeakPickerIM and write directly.
      bool builtin_centroiding = (bruker_config.ms1_centroid_mz_ppm > 0.0f
                                  && bruker_config.ms1_centroid_im_pct > 0.0f);
      if (builtin_centroiding)
      {
        OPENMS_LOG_INFO << "Built-in Bruker IM centroiding was applied during .d loading "
                        << "(ms1_centroid_mz_ppm=" << bruker_config.ms1_centroid_mz_ppm
                        << ", ms1_centroid_im_pct=" << bruker_config.ms1_centroid_im_pct
                        << "). Skipping PeakPickerIM algorithm." << std::endl;
        addDataProcessing_(exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
        MzMLFile().store(output_file, exp);
        return EXECUTION_OK;
      }

      // Check MS1 spectra for IM format
      IMFormat im_format = IMTypes::determineIMFormat(exp, 1);
      if (im_format == IMFormat::NONE)
      {
        OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                        << "No peak picking will be performed." << std::endl;
        MzMLFile().store(output_file, exp);
        return EXECUTION_OK;
      }
      if (im_format == IMFormat::IM_SPECTRUM)
      {
        OPENMS_LOG_ERROR << "Error: Input data has single drift time per spectrum (IM_SPECTRUM format). "
                         << "PeakPickerIM requires per-peak IM arrays (IM_PEAK format). "
                         << "Try using bruker:export_mode=frame." << std::endl;
        return ILLEGAL_PARAMETERS;
      }

      std::exception_ptr first_error = nullptr;
#pragma omp parallel for
      for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
      {
        try
        {
          MSSpectrum& spectrum = exp[static_cast<Size>(i)];
          // Skip already-centroided spectra (e.g., DIA MS2 with bruker:dia_ms2_centroid=true)
          if (spectrum.getIMPeakType() == IMPeakType::IM_CENTROIDED) continue;
          if (method == "mobilogram")       picker.pickIMTraces(spectrum);
          else if (method == "cluster")     picker.pickIMCluster(spectrum);
          else if (method == "traces")      picker.pickIMElutionProfiles(spectrum);
        }
        catch (...)
        {
#pragma omp critical
          { if (!first_error) first_error = std::current_exception(); }
        }
      }
      if (first_error) std::rethrow_exception(first_error);

      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
      MzMLFile().store(output_file, exp);
      return EXECUTION_OK;
    }
#endif

    if (process_opt == "lowmemory")
    {
      return doLowMemAlgorithm(method, picker, input_file, output_file);
    }
    else
    {
      PeakMap exp;
      MzMLFile mzml;
      mzml.load(input_file, exp);

      // Check MS1 spectra for IM format (PeakPickerIM works on per-peak IM data in MS1 frames)
      IMFormat im_format = IMTypes::determineIMFormat(exp, 1);
      if (im_format == IMFormat::NONE)
      {
        OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                        << "No peak picking will be performed." << std::endl;
        mzml.store(output_file, exp);
        return EXECUTION_OK;
      }
      if (im_format == IMFormat::IM_SPECTRUM)
      {
        OPENMS_LOG_ERROR << "Error: Input data has single drift time per spectrum (IM_SPECTRUM format). "
                         << "PeakPickerIM requires per-peak IM arrays (IM_PEAK format)." << std::endl;
        return ILLEGAL_PARAMETERS;
      }

      std::exception_ptr first_error = nullptr;
#pragma omp parallel for
      for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
      {
        try
        {
          MSSpectrum& spectrum = exp[static_cast<Size>(i)];
          if (method == "mobilogram")       picker.pickIMTraces(spectrum);
          else if (method == "cluster")     picker.pickIMCluster(spectrum);
          else if (method == "traces")      picker.pickIMElutionProfiles(spectrum);
        }
        catch (...)
        {
#pragma omp critical
          { if (!first_error) first_error = std::current_exception(); }
        }
      }
      if (first_error) std::rethrow_exception(first_error);

      // Annotate processing info (same as low-memory path)
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));

      mzml.store(output_file, exp);
      return EXECUTION_OK;
    }
  }
};

int main(int argc, const char** argv)
{
  TOPPPeakPickerIM tool;
  return tool.main(argc, argv);
}

/// @endcond

