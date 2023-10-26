// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownTagger.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>


#include <QFileInfo>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHDeconv TOPP_FLASHDeconv

  @brief FLASHDeconv performs ultrafast deconvolution of top down proteomics MS datasets.
  FLASHDeconv takes mzML file as input and outputs deconvolved feature list (.tsv) and
  deconvolved spectra files (.tsv, .mzML, .msalign, .ms1ft).
  FLASHDeconv uses SpectralDeconvolution for spectral level deconvolution and MassFeatureTracer to detect mass features.
  Also for MSn spectra, the precursor masses (not peak m/zs) should be determined and assigned in most cases. This assignment
  can be done by tracking MSn-1 spectra deconvolution information. Thus FLASHDeconv class keeps MSn-1 spectra deconvolution information
  for a certain period for precursor mass assignment in DeconvolvedSpectrum class.
  In case of FLASHIda runs, this precursor mass assignment is done by FLASHIda. Thus FLASHDeconv class simply parses the log file
  from FLASHIda runs and pass the parsed information to DeconvolvedSpectrum class.
*/
class TOPPFLASHDeconv : public TOPPBase
{
public:
  TOPPFLASHDeconv() : TOPPBase("FLASHDeconv", "Ultra-fast high-quality deconvolution enables online processing of top-down MS data")
  {
  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (mzML)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Default output tsv file containing deconvolved features");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFileList_("out_spec", "<file for MS1, file for MS2, ...>", {""}, "Output tsv files containing deconvolved spectra (per MS level)", false);
    setValidFormats_("out_spec", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_mzml", "<file>", "", "Output mzml file containing deconvolved spectra (of all MS levels)", false);
    setValidFormats_("out_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_quant", "<file>", "", "Output tsv file containing isobaric quantification results", false);
    setValidFormats_("out_quant", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_annotated_mzml", "<file>", "",
                        "Output mzml file containing annotated spectra. For each annotated peak, monoisotopic mass, charge, and isotope index are stored as meta data. Unannotated peaks are also "
                        "copied as well without meta data.",
                        false);
    setValidFormats_("out_annotated_mzml", ListUtils::create<String>("mzML"));

    registerOutputFileList_("out_msalign", "<file for MS1, file for MS2, ...>", {""},
                            "Output msalign (topFD and ProMex compatible) files containing deconvolved spectra (per MS level)."
                            " The file name for MSn should end with msn.msalign to be able to be recognized by TopPIC GUI. "
                            "For example, -out_msalign [name]_ms1.msalign [name]_ms2.msalign",
                            false);

    setValidFormats_("out_msalign", ListUtils::create<String>("msalign"), false);

    registerOutputFileList_("out_feature", "<file for MS1, file for MS2, ...>", {""},
                            "Output feature (topFD compatible) files containing deconvolved features (per MS level). "
                            "The feature files are necessary for TopPIC feature intensity output.",
                            false);
    setValidFormats_("out_feature", ListUtils::create<String>("feature"), false);

    registerDoubleOption_("min_precursor_snr", "<SNR value>", 1.0,
                          "Minimum precursor SNR (SNR within the precursor envelope range) for identification. Similar to precursor interference level, but far more stringent as it also considers "
                          "the isotope distribution shape of signal."
                          "When FLASHIda log file is used, this parameter is ignored. Applied only for topFD msalign outputs.",
                          false, false);

    setMinFloat_("min_precursor_snr", .0);

    registerDoubleOption_("min_precursor_qvalue", "<q value>", 1.0,
                          "Minimum precursor q value for identification. To use this threshold, set -report_FDR 1."
                          " Specify to, for instance, 0.01 to make sure your precursor deconvolution FDR is less than 0.01."
                          " Applied only for topFD msalign outputs, regardless of using FLASHIda.",
                          false, false);

    setMinFloat_("min_precursor_qvalue", .0);
    setMaxFloat_("min_precursor_qvalue", 1.0);


    registerIntOption_("mzml_mass_charge", "<0:uncharged 1: +1 charged -1: -1 charged>", 0, "Charge state of deconvolved masses in mzml output (specified by out_mzml)", false);

    setMinInt_("mzml_mass_charge", -1);
    setMaxInt_("mzml_mass_charge", 1);


    registerIntOption_("write_detail", "<1:true 0:false>", 0,
                       "To write peak information per deconvolved mass in detail or not in tsv files for deconvolved spectra. "
                       "If set to 1, all peak information (m/z, intensity, charge and isotope index) per mass is reported.",
                       false, true);

    setMinInt_("write_detail", 0);
    setMaxInt_("write_detail", 1);

    registerDoubleOption_("min_mz", "<m/z value>", -1.0, "If set to positive value, minimum m/z to deconvolve.", false, true);

    registerDoubleOption_("max_mz", "<m/z value>", -1.0, "If set to positive value, maximum m/z to deconvolve.", false, true);

    registerDoubleOption_("min_rt", "<RT value>", -1.0, "If set to positive value, minimum RT (in second) to deconvolve.", false, true);

    registerDoubleOption_("max_rt", "<RT value>", -1.0, "If set to positive value, maximum RT (in second) to deconvolve.", false, true);

    Param combined;
    auto fd_param = FLASHDeconvAlgorithm().getDefaults();

    fd_param.removeAll("sd:");
    fd_param.removeAll("ft:");
    fd_param.removeAll("iq:");
    combined.insert("fd:", fd_param);

    fd_param = FLASHDeconvAlgorithm().getDefaults();
    combined.insert("", fd_param.copy("sd:", false));
    combined.insert("", fd_param.copy("ft:", false));
    combined.insert("", fd_param.copy("iq:", false));

    auto tagger_param = TopDownTagger().getDefaults();
    tagger_param.setValue("tol", DoubleList {}, "ppm tolerances for tag generation. If not set, sd:tol will be used.");
    tagger_param.addTag("tol", "advanced");

    combined.insert("tagger:", tagger_param);

    combined.setSectionDescription("fd", "FLASHDeconv algorithm parameters (prefix fd:)");
    combined.setSectionDescription("sd", "Spectral deconvolution parameters (prefix sd:)");
    combined.setSectionDescription("ft", "Feature tracing parameters (prefix ft:)");
    combined.setSectionDescription("iq", "Isobaric quantification parameters (prefix iq:)");
    combined.setSectionDescription("tagger", "Tagger parameters (prefix tagger:)");

    registerFullParam_(combined);
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
    //String in_train_file {};

    auto out_spec_file = getStringList_("out_spec");
    String out_mzml_file = getStringOption_("out_mzml");
    String out_anno_mzml_file = getStringOption_("out_annotated_mzml");
    String out_quant_file = getStringOption_("out_quant");
    auto out_topfd_file = getStringList_("out_msalign");
    auto out_topfd_feature_file = getStringList_("out_feature");
    double topFD_SNR_threshold = getDoubleOption_("min_precursor_snr");
    double topFD_qval_threshold = getDoubleOption_("min_precursor_qvalue");

    bool write_detail = getIntOption_("write_detail") > 0;
    int mzml_charge = getIntOption_("mzml_mass_charge");
    double min_mz = getDoubleOption_("min_mz");
    double max_mz = getDoubleOption_("max_mz");
    double min_rt = getDoubleOption_("min_rt");
    double max_rt = getDoubleOption_("max_rt");
    std::map<uint, int> per_ms_level_spec_count;
    std::map<uint, int> per_ms_level_deconv_spec_count;
    std::map<uint, int> per_ms_level_mass_count;
    FLASHDeconvAlgorithm fd;
    Param tmp_fd_param = getParam_().copy("fd:", true);
    Param fd_param;
    fd_param.insert("", tmp_fd_param);
    bool report_decoy = (int)tmp_fd_param.getValue("report_FDR") > 0;
    topFD_SNR_threshold = tmp_fd_param.getValue("ida_log").toString().empty()? topFD_SNR_threshold : 0;

    tmp_fd_param = getParam_().copy("sd:", false);
    fd_param.insert("", tmp_fd_param);
    DoubleList tols = tmp_fd_param.getValue("sd:tol");

    tmp_fd_param = getParam_().copy("ft:", false);
    fd_param.insert("", tmp_fd_param);

    tmp_fd_param = getParam_().copy("iq:", false);
    fd_param.insert("", tmp_fd_param);
    fd.setParameters(fd_param);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.
    PeakFileOptions opt = mzml.getOptions();
    if (min_rt > 0 || max_rt > 0)
    {
      opt.setRTRange(DRange<1> {min_rt, max_rt});
    }
    if (min_mz > 0 || max_mz > 0)
    {
      opt.setMZRange(DRange<1> {min_mz, max_mz});
    }
    mzml.setLogType(log_type_);
    mzml.setOptions(opt);
    mzml.load(in_file, map);

    std::vector<DeconvolvedSpectrum> deconvolved_spectra;
    std::vector<FLASHDeconvHelperStructs::MassFeature> deconvolved_features;
    std::map<int, double> scan_rt_map;
    std::map<int, PeakGroup> msNscan_to_precursor_pg;

    // Run FLASHDeconvAlgorithm here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;
    fd.run(map, deconvolved_spectra, deconvolved_features);

    // collect statistics for information
    for (auto& it : map)
    {
      uint ms_level = it.getMSLevel();
      if (per_ms_level_spec_count.find(ms_level) == per_ms_level_spec_count.end())
        per_ms_level_spec_count[ms_level] = 0;
      per_ms_level_spec_count[ms_level]++;
    }

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      if (per_ms_level_deconv_spec_count.find(ms_level) == per_ms_level_deconv_spec_count.end())
        per_ms_level_deconv_spec_count[ms_level] = 0;
      if (per_ms_level_mass_count.find(ms_level) == per_ms_level_mass_count.end())
        per_ms_level_mass_count[ms_level] = 0;

      per_ms_level_deconv_spec_count[ms_level]++;
      per_ms_level_mass_count[ms_level] += (int)deconvolved_spectrum.size();
      scan_rt_map[deconvolved_spectrum.getScanNumber()] = deconvolved_spectrum.getOriginalSpectrum().getRT();
      if (ms_level > 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
      {
        msNscan_to_precursor_pg[deconvolved_spectrum.getScanNumber()] = deconvolved_spectrum.getPrecursorPeakGroup();
      }
    }

    for (auto& val : per_ms_level_deconv_spec_count)
    {
      OPENMS_LOG_INFO << "So far, FLASHDeconv found " << per_ms_level_mass_count[val.first] << " masses in " << val.second << " MS" << val.first << " spectra out of "
                      << per_ms_level_spec_count[val.first] << endl;
    }
    if (!deconvolved_features.empty())
    {
      OPENMS_LOG_INFO << "Mass tracer found " << deconvolved_features.size() << " features" << endl;
    }

    OPENMS_LOG_INFO << "FLASHDeconv run complete. Now writing the results in output files ..." << endl;

    // Run tagger
    TopDownTagger tagger;
    auto tagger_param = getParam_().copy("tagger:", true);
    if (((DoubleList)tagger_param.getValue("tol")).empty())
    {
      tagger_param.setValue("tol", tols);
    }
    tagger.setParameters(tagger_param);

    std::vector<std::string> tags;
    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      tagger.run(deconvolved_spectrum, tags);
      for (auto& tag : tags)
      {
        std::cout<<tag<<std::endl;
        std::reverse(tag.begin(), tag.end());
        std::cout<< tag <<std::endl;
      }
      tags.clear();
    }
    // Write output files
    // default feature deconvolution tsv output
    if (!deconvolved_features.empty())
    {
      fstream out_stream;
      out_stream.open(out_file, fstream::out);
      FLASHDeconvFeatureFile::writeHeader(out_stream, report_decoy);
      FLASHDeconvFeatureFile::writeFeatures(deconvolved_features, in_file, out_stream, report_decoy);
      out_stream.close();
    }
    // Per ms level spectrum deconvolution tsv output
    if (!out_spec_file.empty())
    {
      std::vector<fstream> out_spec_streams = std::vector<fstream>(out_spec_file.size());
      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end())
          continue;
        out_spec_streams[i].open(out_spec_file[i], fstream::out);
        FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(out_spec_streams[i], i + 1, write_detail, report_decoy);
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (out_spec_file[ms_level - 1].empty()) continue;
        FLASHDeconvSpectrumFile::writeDeconvolvedMasses(deconvolved_spectrum, deconvolved_spectrum, out_spec_streams[ms_level - 1], in_file, fd.getAveragine(), tols[ms_level - 1], write_detail,
                                                        report_decoy);
      }

      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end())
          continue;
        out_spec_streams[i].close();
      }
    }
    // topFD feature output
    if (!out_topfd_feature_file.empty())
    {
      std::vector<fstream> out_topfd_feature_streams;
      out_topfd_feature_streams = std::vector<fstream>(out_topfd_feature_file.size());
      for (Size i = 0; i < out_topfd_feature_file.size(); i++)
      {
        if (out_topfd_feature_file[i].empty() || per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end())
          continue;
        out_topfd_feature_streams[i].open(out_topfd_feature_file[i], fstream::out);
        FLASHDeconvFeatureFile::writeTopFDFeatureHeader(out_topfd_feature_streams[i], i + 1);
        FLASHDeconvFeatureFile::writeTopFDFeatures(deconvolved_features, msNscan_to_precursor_pg, scan_rt_map, in_file, out_topfd_feature_streams[i], i + 1);
        out_topfd_feature_streams[i].close();
      }
    }

    // topFD msalign output
    if (!out_topfd_file.empty())
    {
      auto out_topfd_streams = std::vector<fstream>(out_topfd_file.size());
      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end())
          continue;
        out_topfd_streams[i].open(out_topfd_file[i], fstream::out);
      }
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (out_topfd_file[ms_level - 1].empty()) continue;
        FLASHDeconvSpectrumFile::writeTopFD(deconvolved_spectrum, out_topfd_streams[ms_level - 1], topFD_SNR_threshold, topFD_qval_threshold, per_ms_level_deconv_spec_count.begin()->first, false,
                                            false);
      }

      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end())
          continue;
        out_topfd_streams[i].close();
      }
    }

    // isobaric quantification output
    if (!out_quant_file.empty())
    {
      fstream out_quant_stream;
      out_quant_stream.open(out_quant_file, fstream::out);
      FLASHDeconvSpectrumFile::writeIsobaricQuantification(out_quant_stream, deconvolved_spectra);
      out_quant_stream.close();
    }
    // mzML output
    if (!out_anno_mzml_file.empty() || !out_mzml_file.empty())
    {
      FLASHDeconvSpectrumFile::writeMzML(map, deconvolved_spectra, out_mzml_file, out_anno_mzml_file, mzml_charge, tols);
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