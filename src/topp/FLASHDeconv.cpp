// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
@page TOPP_FLASHDeconv FLASHDeconv

  @brief FLASHDeconv performs ultrafast deconvolution of top-down proteomics MS datasets.

  FLASHDeconv takes an mzML file as input and outputs deconvolved feature list (.tsv) and
  deconvolved spectra files (.tsv, .mzML, .msalign, .feature).
  FLASHDeconv uses SpectralDeconvolution for spectral level deconvolution and MassFeatureTrace to detect mass features.
  For MSn spectra, the precursor masses (not peak m/zs) are determined by tracking MSn-1 spectra deconvolution information.

  See https://openms.de/FLASHDeconv for more information.


<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FLASHDeconv.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FLASHDeconv.html
*/
class TOPPFLASHDeconv : public TOPPBase
{
public:
  TOPPFLASHDeconv():
      TOPPBase("FLASHDeconv",
               "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
               true,
               {Citation {"Jeong K, Kim J, Gaikwad M et al.", "FLASHDeconv: Ultrafast, High-Quality Feature Deconvolution for Top-Down Proteomics",
                          "Cell Syst 2020 Feb 26;10(2):213-218.e6", "10.1016/j.cels.2020.01.003"}})
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file in mzML format. ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Default output tsv file containing deconvolved features");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFile_(
      "out_spec1", "<file>", "",
      "Output tsv file for deconvolved MS1 spectra. Use -out_spec2, ..., -out_spec4 for MS2, ..., MS4 spectra.", false);
    setValidFormats_("out_spec1", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_spec2", "<file>", "", "Output TSV files for deconvolved MS2 spectra.", false, true);
    setValidFormats_("out_spec2", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_spec3", "<file>", "", "Output TSV files for deconvolved MS3 spectra.", false, true);
    setValidFormats_("out_spec3", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_spec4", "<file>", "", "Output TSV files for deconvolved MS4 spectra.", false, true);
    setValidFormats_("out_spec4", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_mzml", "<file>", "", "Output mzML file containing deconvolved spectra (for all MS levels).", false);
    setValidFormats_("out_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_quant", "<file>", "", "Output tsv file with isobaric quantification results for MS2 spectra.", false);
    setValidFormats_("out_quant", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_annotated_mzml", "<file>", "",
                        "Output annotated mzML file with monoisotopic mass, charge, and isotope index metadata for peaks. Unannotated peaks are also retained without metadata.",
                        false);
    setValidFormats_("out_annotated_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_(
      "out_msalign1", "<file>", "",
      "Output msalign (TopFD and ProMex compatible) file for MS1 deconvolved spectra. Ensure filename ends with ms1.msalign for TopPIC GUI compatibility (e.g., result_ms1.msalign; refer to TopPIC input formats).",
      false);
    setValidFormats_("out_msalign1", ListUtils::create<String>("msalign"), false);

    registerOutputFile_("out_msalign2", "<file>", "",
                        "Output msalign (TopFD and ProMex compatible) file for MS2 deconvolved spectra. Ensure filename ends with ms2.msalign for TopPIC GUI compatibility (e.g., result_ms2.msalign; refer to TopPIC input formats).",
                        false, true);
    setValidFormats_("out_msalign2", ListUtils::create<String>("msalign"), false);

    registerOutputFile_("out_feature1", "<file>", "",
                        "Output feature file (TopFD compatible) for MS1 spectra. It is needed for TopPIC feature intensity output (refer to TopPIC input formats).",
                        false);

    setValidFormats_("out_feature1", ListUtils::create<String>("feature"), false);

    registerOutputFile_("out_feature2", "<file>", "",
                        "Output feature file (TopFD compatible) for MS2 spectra. It is needed for TopPIC feature intensity output (refer to TopPIC input formats).",
                        false, true);

    setValidFormats_("out_feature2", ListUtils::create<String>("feature"), false);

    registerFlag_("keep_empty_out", "Retain empty output files (e.g., *.tsv files with no features).");

    registerIntOption_("mzml_mass_charge", "<0:uncharged 1: +1 charged -1: -1 charged>", 0,
                       "Charge state of deconvolved masses in mzML output specified by -out_mzml.", false, true);
    setMinInt_("mzml_mass_charge", -1);
    setMaxInt_("mzml_mass_charge", 1);

    registerFlag_("write_detail",
                  "Include detailed peak information (m/z, intensity, charge, isotope index) for each deconvolved mass in the output spectrum tsv files specified by out_spec* options.",
                  false);

    registerDoubleOption_("min_mz", "<m/z value>", -1.0, "Specify the minimum m/z values for peaks considered during deconvolution. Negative values disable the threshold.", false, true);
    registerDoubleOption_("max_mz", "<m/z value>", -1.0, "Specify the maximum m/z values for peaks considered during deconvolution. Negative values disable the threshold.", false, true);
    registerDoubleOption_("min_rt", "<RT value>", -1.0, "Specify the minimum retention time (in minutes) for spectra considered during deconvolution. Negative values disable the threshold.", false, true);
    registerDoubleOption_("max_rt", "<RT value>", -1.0, "Specify the maximum retention time (in minutes) for spectra considered during deconvolution. Negative values disable the threshold.", false, true);

    registerIntOption_("max_ms_level", "<MS level>", -1, "Set the maximum MS level (inclusive) for deconvolution. Negative values disable the threshold.", false, true);

    registerSubsection_("FD", "FLASHDeconv algorithm parameters");
    registerSubsection_("SD", "Spectral deconvolution parameters");
    registerSubsection_("ft", "Feature tracing parameters");
    registerSubsection_("iq", "Isobaric quantification parameters");
  }

  /// Returns default parameters for each subsection.
  /// FLASHDeconvAlgorithm stores all parameters in a single Param object with prefixes (SD:, ft:, iq:).
  /// This function splits them into separate subsections for a cleaner CLI interface:
  ///   -FD:*  -> top-level algorithm parameters (after removing SD:, ft:, iq: prefixed params)
  ///   -SD:*  -> spectral deconvolution parameters
  ///   -ft:*  -> feature tracing parameters
  ///   -iq:*  -> isobaric quantification parameters
  Param getSubsectionDefaults_(const String& prefix) const override
  {
    auto fd_param = FLASHDeconvAlgorithm().getDefaults();

    if (prefix == "FD")
    {
      // Remove nested subsection params to get only top-level algorithm parameters
      fd_param.removeAll("SD:");
      fd_param.removeAll("ft:");
      fd_param.removeAll("iq:");
      return fd_param;
    }
    else if (prefix == "SD")
    {
      return fd_param.copy("SD:", true);  // Extract SD:* params, strip prefix
    }
    else if (prefix == "ft")
    {
      return fd_param.copy("ft:", true);  // Extract ft:* params, strip prefix
    }
    else if (prefix == "iq")
    {
      return fd_param.copy("iq:", true);  // Extract iq:* params, strip prefix
    }
    else { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix); }
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**) override
  {
    OPENMS_LOG_INFO << "Initializing ... " << endl;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String in_file = getStringOption_("in");
    String out_file = getStringOption_("out");
    bool keep_empty_out = getFlag_("keep_empty_out");
    auto out_spec_file
      = StringList {getStringOption_("out_spec1"), getStringOption_("out_spec2"), getStringOption_("out_spec3"), getStringOption_("out_spec4")};

    auto out_topfd_file = StringList {getStringOption_("out_msalign1"), getStringOption_("out_msalign2")};
    auto out_topfd_feature_file = StringList {getStringOption_("out_feature1"), getStringOption_("out_feature2")};

    String out_mzml_file = getStringOption_("out_mzml");
    String out_anno_mzml_file = getStringOption_("out_annotated_mzml");
    String out_quant_file = getStringOption_("out_quant");

    bool write_detail = getFlag_("write_detail");
    int mzml_charge = getIntOption_("mzml_mass_charge");
    double min_mz = getDoubleOption_("min_mz");
    double max_mz = getDoubleOption_("max_mz");
    double min_rt = getDoubleOption_("min_rt") * 60.0;
    double max_rt = getDoubleOption_("max_rt") * 60.0;
    int max_ms_level = getIntOption_("max_ms_level");

    std::map<uint, int> per_ms_level_spec_count;
    std::map<uint, int> per_ms_level_deconv_spec_count;
    std::map<uint, int> per_ms_level_mass_count;
    FLASHDeconvAlgorithm fd;

    // Reassemble parameters from CLI subsections back into format expected by FLASHDeconvAlgorithm.
    // getSubsectionDefaults_() split the params for CLI display; here we merge them back:
    //   FD:* params -> inserted without prefix (top-level algorithm params)
    //   SD:*, ft:*, iq:* params -> inserted with prefix preserved (nested params)
    Param fd_param;
    Param tmp_fd_param = getParam_().copy("FD:", true);  // copy FD:* params, strip "FD:" prefix
    fd_param.insert("", tmp_fd_param);
    bool report_decoy = tmp_fd_param.getValue("report_FDR") != "false";

    tmp_fd_param = getParam_().copy("SD:", false);  // copy SD:* params, keep "SD:" prefix
    fd_param.insert("", tmp_fd_param);
    DoubleList tols = tmp_fd_param.getValue("SD:tol");

    tmp_fd_param = getParam_().copy("ft:", false);  // copy ft:* params, keep "ft:" prefix
    fd_param.insert("", tmp_fd_param);

    tmp_fd_param = getParam_().copy("iq:", false);  // copy iq:* params, keep "iq:" prefix
    fd_param.insert("", tmp_fd_param);
    fd.setParameters(fd_param);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    constexpr double MAX_RANGE_VALUE = 1e7; // effectively unlimited upper bound for RT/m/z ranges

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.
    PeakFileOptions opt = mzml.getOptions();
    if (min_rt > 0 || max_rt > 0)
    {
      if (min_rt > 0 && max_rt < 0) max_rt = MAX_RANGE_VALUE;
      opt.setRTRange(DRange<1> {min_rt, max_rt});
    }
    if (min_mz > 0 || max_mz > 0)
    {
      if (min_mz > 0 && max_mz < 0) max_mz = MAX_RANGE_VALUE;
      opt.setMZRange(DRange<1> {min_mz, max_mz});
    }
    if (max_ms_level > 0)
    {
      IntList ms_levels;
      for (int msl = 1; msl <= max_ms_level; msl++)
        ms_levels.push_back(msl);
      opt.setMSLevels(ms_levels);
    }

    mzml.setLogType(log_type_);
    mzml.setOptions(opt);
    mzml.load(in_file, map);

    std::vector<DeconvolvedSpectrum> deconvolved_spectra;
    std::vector<FLASHHelperClasses::MassFeature> deconvolved_features;
    std::map<int, double> scan_rt_map;

    // Run FLASHDeconvAlgorithm here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    fd.run(map, deconvolved_spectra, deconvolved_features);
    tols = fd.getTolerances();
    // collect statistics for information
    for (const auto& it : map)
    {
      uint ms_level = it.getMSLevel();
      if (per_ms_level_spec_count.find(ms_level) == per_ms_level_spec_count.end()) per_ms_level_spec_count[ms_level] = 0;
      per_ms_level_spec_count[ms_level]++;
    }

    for (const auto& deconvolved_spectrum : deconvolved_spectra)
    {
      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      scan_rt_map[deconvolved_spectrum.getScanNumber()] = deconvolved_spectrum.getOriginalSpectrum().getRT();

      if (deconvolved_spectrum.empty()) continue;
      if (per_ms_level_deconv_spec_count.find(ms_level) == per_ms_level_deconv_spec_count.end()) per_ms_level_deconv_spec_count[ms_level] = 0;
      if (per_ms_level_mass_count.find(ms_level) == per_ms_level_mass_count.end()) per_ms_level_mass_count[ms_level] = 0;

      per_ms_level_deconv_spec_count[ms_level]++;
      per_ms_level_mass_count[ms_level] += (int)deconvolved_spectrum.size();

    }
    for (const auto& val : per_ms_level_deconv_spec_count)
    {
      OPENMS_LOG_INFO << "So far, FLASHDeconv found " << per_ms_level_mass_count[val.first] << " masses in " << val.second << " MS" << val.first
                      << " spectra out of " << per_ms_level_spec_count[val.first] << endl;
    }
    if (! deconvolved_features.empty()) { OPENMS_LOG_INFO << "Mass tracer found " << deconvolved_features.size() << " features" << endl; }

    OPENMS_LOG_INFO << "FLASHDeconv run complete. Now writing the results in output files ..." << endl;
    // Write output files
    // default feature deconvolution tsv output

    if (keep_empty_out || ! deconvolved_features.empty())
    {
      OPENMS_LOG_INFO << "writing feature tsv ..." << endl;
      ofstream out_stream;
      out_stream.open(out_file);
      if (!out_stream)
      {
        OPENMS_LOG_FATAL_ERROR << "Error: Could not open output file '" << out_file << "'" << endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      FLASHDeconvFeatureFile::writeHeader(out_stream, report_decoy);
      FLASHDeconvFeatureFile::writeFeatures(deconvolved_features, in_file, out_stream, report_decoy);
      out_stream.close();
    }
    // Per ms level spectrum deconvolution tsv output
    // Check if any spectrum output file is specified
    auto has_any_output = [](const StringList& files) {
      for (const auto& f : files) if (!f.empty()) return true;
      return false;
    };

    if (has_any_output(out_spec_file))
    {
      std::vector<ofstream> out_spec_streams = std::vector<ofstream>(out_spec_file.size());
      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing spectrum tsv for MS level " << (i + 1) << " ..." << endl;

        out_spec_streams[i].open(out_spec_file[i]);
        if (!out_spec_streams[i])
        {
          OPENMS_LOG_FATAL_ERROR << "Error: Could not open output file '" << out_spec_file[i] << "'" << endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(out_spec_streams[i], i + 1, write_detail, report_decoy);
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (ms_level > out_spec_file.size() || out_spec_file[ms_level - 1].empty()) continue;
        FLASHDeconvSpectrumFile::writeDeconvolvedMasses(deconvolved_spectrum, out_spec_streams[ms_level - 1], in_file, fd.getAveragine(), fd.getDecoyAveragine(),
                                                        tols[ms_level - 1], write_detail, report_decoy, fd.getNoiseDecoyWeight());
      }

      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        out_spec_streams[i].close();
      }
    }

    // mzML output
    if (! out_anno_mzml_file.empty() || ! out_mzml_file.empty())
    {
      FLASHDeconvSpectrumFile::writeMzML(map, deconvolved_spectra, out_mzml_file, out_anno_mzml_file, mzml_charge, tols);
    }

    // isobaric quantification output
    if (! out_quant_file.empty())
    {
      OPENMS_LOG_INFO << "writing quantification tsv ..." << endl;
      ofstream out_quant_stream;
      out_quant_stream.open(out_quant_file);
      if (!out_quant_stream)
      {
        OPENMS_LOG_FATAL_ERROR << "Error: Could not open output file '" << out_quant_file << "'" << endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      FLASHDeconvSpectrumFile::writeIsobaricQuantification(out_quant_stream, deconvolved_spectra);
      out_quant_stream.close();
    }

    // topFD feature output - it should be at the end since zero feature IDs are redefined for TopPIC feature indices.
    if (has_any_output(out_topfd_feature_file))
    {
      std::vector<ofstream> out_topfd_feature_streams;
      out_topfd_feature_streams = std::vector<ofstream>(out_topfd_feature_file.size());
      for (Size i = 0; i < out_topfd_feature_file.size(); i++)
      {
        if (out_topfd_feature_file[i].empty()
            || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing topfd *.feature for MS level " << (i + 1) << " ..." << endl;

        out_topfd_feature_streams[i].open(out_topfd_feature_file[i]);
        if (!out_topfd_feature_streams[i])
        {
          OPENMS_LOG_FATAL_ERROR << "Error: Could not open output file '" << out_topfd_feature_file[i] << "'" << endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }
        FLASHDeconvFeatureFile::writeTopFDFeatureHeader(out_topfd_feature_streams[i], i + 1);
        FLASHDeconvFeatureFile::writeTopFDFeatures(deconvolved_spectra, deconvolved_features, scan_rt_map, in_file, out_topfd_feature_streams[i],
                                                   i + 1);
        out_topfd_feature_streams[i].close();
      }
    }
    // topFD msalign output
    if (has_any_output(out_topfd_file))
    {
      auto out_topfd_streams = std::vector<ofstream>(out_topfd_file.size());

      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing topfd *.msalign for MS level " << (i + 1) << " ..." << endl;

        out_topfd_streams[i].open(out_topfd_file[i]);
        if (!out_topfd_streams[i])
        {
          OPENMS_LOG_FATAL_ERROR << "Error: Could not open output file '" << out_topfd_file[i] << "'" << endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }
        FLASHDeconvSpectrumFile::writeTopFDHeader(out_topfd_streams[i], getParam_().copy("SD:", true));
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (ms_level > out_topfd_file.size() || out_topfd_file[ms_level - 1].empty()) continue;

        FLASHDeconvSpectrumFile::writeTopFD(deconvolved_spectrum, out_topfd_streams[ms_level - 1], in_file, 1,
                                            per_ms_level_deconv_spec_count.begin()->first, false, false);
      }

      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        out_topfd_streams[i].close();
      }
    }

    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHDeconv tool;
  return tool.main(argc, argv);
}