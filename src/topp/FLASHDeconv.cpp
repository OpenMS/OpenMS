// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------
#define USE_TAGGER //
// #define TRAIN_OUT

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#ifdef USE_TAGGER
  #include <OpenMS/ANALYSIS/TOPDOWN/TopDownTagger.h>
  #include <OpenMS/CHEMISTRY/AASequence.h>
  #include <OpenMS/FORMAT/FASTAFile.h>
#endif
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFileInfo>
#ifdef TRAIN_OUT
  #include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#endif
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
@page TOPP_FLASHDeconv FLASHDeconv

  @brief FLASHDeconv performs ultrafast deconvolution of top down proteomics MS datasets.
  FLASHDeconv takes mzML file as input and outputs deconvolved feature list (.tsv) and
  deconvolved spectra files (.tsv, .mzML, .msalign, .ms1ft).
  FLASHDeconv uses SpectralDeconvolution for spectral level deconvolution and MassFeatureTracer to detect mass features.
  Also for MSn spectra, the precursor masses (not peak m/zs) should be determined and assigned in most cases. This assignment
  can be done by tracking MSn-1 spectra deconvolution information. Thus FLASHDeconv class keeps MSn-1 spectra deconvolution information
  for a certain period for precursor mass assignment in DeconvolvedSpectrum class.
  In case of FLASHIda runs, this precursor mass assignment is done by FLASHIda. Thus FLASHDeconv class simply parses the log file
  from FLASHIda runs and pass the parsed information to DeconvolvedSpectrum class.

See https://openms.de/FLASHDeconv for more information.


<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FLASHDeconv.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FLASHDeconv.html
*/
class TOPPFLASHDeconv : public TOPPBase
{
public:
  TOPPFLASHDeconv() :
      TOPPBase("FLASHDeconv", "Ultra-fast high-quality deconvolution enables online processing of top-down MS data", true,
               {Citation {"Jeong K, Kim J, Gaikwad M et al.", "FLASHDeconv: Ultrafast, High-Quality Feature Deconvolution for Top-Down Proteomics", "Cell Syst 2020 Feb 26;10(2):213-218.e6",
                          "10.1016/j.cels.2020.01.003"}})
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

    registerOutputFile_("out_spec1", "<file>", "", "Output tsv file containing deconvolved MS1 spectra. Likewise, use -out_spec2, ..., -out_spec4 to specify tsv files for MS2, ..., MS4.", false);
    setValidFormats_("out_spec1", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_spec2", "<file>", "", "Output tsv file containing deconvolved MS2 spectra.", false, true);
    setValidFormats_("out_spec2", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_spec3", "<file>", "", "Output tsv file containing deconvolved MS3 spectra.", false, true);
    setValidFormats_("out_spec3", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_spec4", "<file>", "", "Output tsv file containing deconvolved MS4 spectra.", false, true);
    setValidFormats_("out_spec4", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_mzml", "<file>", "", "Output mzml file containing deconvolved spectra (of all MS levels)", false);
    setValidFormats_("out_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_quant", "<file>", "", "Output tsv file containing isobaric quantification results for MS2 only", false);
    setValidFormats_("out_quant", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_annotated_mzml", "<file>", "",
                        "Output mzml file containing annotated spectra. For each annotated peak, monoisotopic mass, charge, and isotope index are stored as meta data. Unannotated peaks are also "
                        "copied as well without meta data.",
                        false);
    setValidFormats_("out_annotated_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_msalign1", "<file>", "",
                        "Output msalign (topFD and ProMex compatible) file containing MS1 deconvolved spectra. Likewise, use -out_msalign2 for MS2 spectra."
                        " The file names for MS1 and MS2 should end with ms1.msalign and ms2.msalgin respectively to be able to be recognized by TopPIC GUI. ",
                        false);
    setValidFormats_("out_msalign1", ListUtils::create<String>("msalign"), false);

    registerOutputFile_("out_msalign2", "<file>", "",
                        "Output msalign (topFD and ProMex compatible) file containing MS2 deconvolved spectra."
                        " The file name should end with ms2.msalign to be able to be recognized by TopPIC GUI. ",
                        false, true);
    setValidFormats_("out_msalign2", ListUtils::create<String>("msalign"), false);
    //
    //    registerOutputFile_("out_msalign3", "<file>", "",
    //                        "Output msalign (topFD and ProMex compatible) file containing MS3 deconvolved spectra."
    //                        " The file name should end with ms3.msalign to be able to be recognized by TopPIC GUI. ",
    //                        false);
    //    setValidFormats_("out_msalign3", ListUtils::create<String>("msalign"), false);
    //
    //    registerOutputFile_("out_msalign4", "<file>", "",
    //                        "Output msalign (topFD and ProMex compatible) file containing MS4 deconvolved spectra."
    //                        " The file name should end with ms4.msalign to be able to be recognized by TopPIC GUI. ",
    //                        false);
    //    setValidFormats_("out_msalign4", ListUtils::create<String>("msalign"), false);

    registerOutputFile_("out_feature1", "<file>", "",
                        "Output feature (topFD compatible) file containing MS1 deconvolved features. Likewise, use -out_feature2 for MS2 features. "
                        "The MS1 and MS2 feature files are necessary for TopPIC feature intensity output.",
                        false);

    setValidFormats_("out_feature1", ListUtils::create<String>("feature"), false);

    registerOutputFile_("out_feature2", "<file>", "",
                        "Output feature (topFD compatible) file containing MS2 deconvolved features. "
                        "The MS1 and MS2 feature files are necessary for TopPIC feature intensity output.",
                        false, true);

    setValidFormats_("out_feature2", ListUtils::create<String>("feature"), false);

    registerFlag_("keep_empty_out", "If set, empty output files (e.g., *.tsv file when no feature was generated) are kept.");

    registerIntOption_("mzml_mass_charge", "<0:uncharged 1: +1 charged -1: -1 charged>", 0, "Charge state of deconvolved masses in mzml output (specified by out_mzml)", false, true);
    setMinInt_("mzml_mass_charge", -1);
    setMaxInt_("mzml_mass_charge", 1);

    registerFlag_("write_detail",
                  "To write peak information per deconvolved mass in detail or not in tsv files for deconvolved spectra. "
                  "If set to 1, all peak information (m/z, intensity, charge and isotope index) per mass is reported.",
                  false);


    registerDoubleOption_("precursor_snr", "<snr value>", 1.0, "Precursor SNR threshold for TopFD MS2 msalign tsv files.", false, true);

    registerDoubleOption_("min_mz", "<m/z value>", -1.0, "If set to positive value, minimum m/z to deconvolve.", false, true);
    registerDoubleOption_("max_mz", "<m/z value>", -1.0, "If set to positive value, maximum m/z to deconvolve.", false, true);
    registerDoubleOption_("min_rt", "<RT value>", -1.0, "If set to positive value, minimum RT (in second) to deconvolve.", false, true);
    registerDoubleOption_("max_rt", "<RT value>", -1.0, "If set to positive value, maximum RT (in second) to deconvolve.", false, true);

    registerIntOption_("max_ms_level", "<MS level>", -1.0, "If set to positive value, maximum MS level (inclusive) to deconvolve.", false, true);

    registerSubsection_("FD", "FLASHDeconv algorithm parameters");
    registerSubsection_("SD", "Spectral deconvolution parameters");
    registerSubsection_("ft", "Feature tracing parameters");
    registerSubsection_("iq", "Isobaric quantification parameters");
#ifdef USE_TAGGER
    registerSubsection_("tagger", "Tagger parameters");
#endif
  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "FD")
    {
      auto fd_param = FLASHDeconvAlgorithm().getDefaults();
      fd_param.removeAll("SD:");
      fd_param.removeAll("ft:");
      fd_param.removeAll("iq:");
      return fd_param;
    }
    else if (prefix == "SD")
    {
      auto fd_param = FLASHDeconvAlgorithm().getDefaults();
      return fd_param.copy("SD:", true);
    }
    else if (prefix == "ft")
    {
      auto fd_param = FLASHDeconvAlgorithm().getDefaults();
      return fd_param.copy("ft:", true);
    }
    else if (prefix == "iq")
    {
      auto fd_param = FLASHDeconvAlgorithm().getDefaults();
      return fd_param.copy("iq:", true);
    }
#ifdef USE_TAGGER
    else if (prefix == "tagger")
    {
      auto tagger_param = TopDownTagger().getDefaults();
      tagger_param.remove("min_charge");
      tagger_param.remove("max_charge");

      tagger_param.setValue("fasta", "", "Target protein sequence database against which tags will be matched.");
      tagger_param.addTag("fasta", "input file");

      tagger_param.setValue("out_tag", "", "Output file containing tags.");
      tagger_param.addTag("out_tag", "output file");

      tagger_param.setValue("out_protein", "", "Output file containing matched proteins.");
      tagger_param.addTag("out_protein", "output file");

      return tagger_param;
    }
#endif
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix);
    }
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
    // String in_train_file {};
    bool keep_empty_out = getFlag_("keep_empty_out");
    auto out_spec_file = StringList {getStringOption_("out_spec1"), getStringOption_("out_spec2"), getStringOption_("out_spec3"), getStringOption_("out_spec4")};

    auto out_topfd_file = StringList {getStringOption_("out_msalign1"), getStringOption_("out_msalign2")};
    auto out_topfd_feature_file = StringList {getStringOption_("out_feature1"), getStringOption_("out_feature2")};

    String out_mzml_file = getStringOption_("out_mzml");
    String out_anno_mzml_file = getStringOption_("out_annotated_mzml");
    String out_quant_file = getStringOption_("out_quant");

    bool write_detail = getFlag_("write_detail");
    int mzml_charge = getIntOption_("mzml_mass_charge");
    double min_mz = getDoubleOption_("min_mz");
    double max_mz = getDoubleOption_("max_mz");
    double min_rt = getDoubleOption_("min_rt");
    double max_rt = getDoubleOption_("max_rt");
    int max_ms_level = getIntOption_("max_ms_level");

    std::map<uint, int> per_ms_level_spec_count;
    std::map<uint, int> per_ms_level_deconv_spec_count;
    std::map<uint, int> per_ms_level_mass_count;
    FLASHDeconvAlgorithm fd;
    Param tmp_fd_param = getParam_().copy("FD:", true);
    Param fd_param;
    fd_param.insert("", tmp_fd_param);
    bool report_decoy = tmp_fd_param.getValue("report_FDR") != "false";
    double topfd_snr_threshold = tmp_fd_param.getValue("ida_log").toString().empty() ? getDoubleOption_("precursor_snr") : .0;

    tmp_fd_param = getParam_().copy("SD:", false);
    fd_param.insert("", tmp_fd_param);
    DoubleList tols = tmp_fd_param.getValue("SD:tol");

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

#ifdef USE_TAGGER
    // Run tagger
    TopDownTagger tagger;

    auto tagger_param = getParam_().copy("tagger:", true);
    if ((int)tagger_param.getValue("max_tag_count") > 0 && !deconvolved_spectra.empty())
    {
      OPENMS_LOG_INFO << "Finding sequence tags from deconvolved spectra ..." << endl;

      String fastaname = tagger_param.getValue("fasta").toString();
      String out_tag = tagger_param.getValue("out_tag").toString();
      String out_protein_tag = tagger_param.getValue("out_protein").toString();

      tagger_param.remove("out_tag");
      tagger_param.remove("out_protein");
      tagger_param.remove("fasta");
      tagger.setParameters(tagger_param);

      DeconvolvedSpectrum dspec_for_tagging;
      for (const auto& dspec : deconvolved_spectra)
      {
        if (dspec.isDecoy())
          continue;
        for (const auto& pg : dspec)
          dspec_for_tagging.push_back(pg);
      }

      if (deconvolved_spectra.size() > 1)
      {
        dspec_for_tagging.sort();
        SpectralDeconvolution::removeOverlappingPeakGroups(dspec_for_tagging,
                                                           1e-6 * tols[deconvolved_spectra[0].getOriginalSpectrum().getMSLevel() - 1]); // merged peak groups have scan number information!
      }

      tagger.run(dspec_for_tagging, tols[dspec_for_tagging.getOriginalSpectrum().getMSLevel() - 1]);
      tagger.runMatching(fastaname);

      if (!out_protein_tag.empty())
      {
        fstream out_tagger_stream = fstream(out_protein_tag, fstream::out);
        out_tagger_stream << "ProteinIndex\tProteinAccession\ttProteinDescription\tMatchedAminoAcidCount\tCoverage(%)\tProteinScore\tProteinQvalue\tTagIndices\n";

        for (const auto& hit : tagger.getProteinHits())
        {
          String tagindices = "";
          for (const auto& tag : tagger.getTags(hit))
          {
            if (!tagindices.empty())
              tagindices += ";";
            tagindices += std::to_string(tagger.getTagIndex(tag));
          }

          out_tagger_stream << tagger.getProteinIndex(hit) << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("MatchedAA") << "\t" << 100.0 * hit.getCoverage()
                            << "\t" << hit.getScore() << "\t" << std::to_string((double)hit.getMetaValue("qvalue")) << "\t" << tagindices << "\n";
        }

        out_tagger_stream.close();
      }

      if (!out_tag.empty())
      {
        fstream out_tagger_stream = fstream(out_tag, fstream::out);
        out_tagger_stream << "Scan\tTagIndex\tProteinIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tLength\tDeNovoScore\tmzs\n";

        for (int n = 0; n <= tagger.getProteinHits().size(); n++)
        {
          for (const auto& tag : tagger.getTags())
          {
            auto hits = tagger.getProteinHits(tag);
            if (n < tagger.getProteinHits().size())
            {
              bool found = false;
              for (const auto& hit : hits)
              {
                if (n == tagger.getProteinIndex(hit))
                  found = true;
              }
              if (!found)
                continue;
            }

            if (n == tagger.getProteinHits().size() && !hits.empty())
              continue;

            String acc = "";
            String description = "";
            String hitindices = "";
            for (const auto& hit : hits)
            {
              if (!acc.empty())
                acc += ";";
              if (!description.empty())
                description += ";";
              if (!hitindices.empty())
                hitindices += ";";
              acc += hit.getAccession();
              description += hit.getDescription();
              hitindices += std::to_string(tagger.getProteinIndex(hit));
            }

            out_tagger_stream << dspec_for_tagging.getScanNumber() << "\t" << tagger.getTagIndex(tag) << "\t" << hitindices << "\t" << acc << "\t" << description << "\t" << tag.getSequence() << "\t"
                              << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

            for (const auto& mz : tag.getMzs())
            {
              out_tagger_stream << std::to_string(mz) << ",";
            }
            out_tagger_stream << "\n";
          }
        }

        out_tagger_stream.close();
      }
    }
#endif
    OPENMS_LOG_INFO << "FLASHDeconv run complete. Now writing the results in output files ..." << endl;

    // Write output files
    // default feature deconvolution tsv output
    if (keep_empty_out || !deconvolved_features.empty())
    {
      OPENMS_LOG_INFO << "writing feature tsv ..." << endl;
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
#ifdef TRAIN_OUT
      std::vector<fstream> out_train_streams = std::vector<fstream>(out_spec_file.size());
#endif

      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || (!keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing spectrum tsv for MS level " << (i + 1) << " ..." << endl;

        out_spec_streams[i].open(out_spec_file[i], fstream::out);
        FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(out_spec_streams[i], i + 1, write_detail, report_decoy);
#ifdef TRAIN_OUT
        out_train_streams[i].open(out_spec_file[i] + "_train.csv", fstream::out);
        Qscore::writeAttCsvForQscoreTrainingHeader(out_train_streams[i]);
#endif
      }

      std::map<int, DeconvolvedSpectrum> target_spec_map;
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (out_spec_file[ms_level - 1].empty())
          continue;
        if (deconvolved_spectrum.isDecoy())
          continue;
        if (report_decoy)
          target_spec_map[deconvolved_spectrum.getScanNumber()] = deconvolved_spectrum;
        FLASHDeconvSpectrumFile::writeDeconvolvedMasses(deconvolved_spectrum, out_spec_streams[ms_level - 1], in_file, fd.getAveragine(), tols[ms_level - 1], write_detail, report_decoy, 1.0);
#ifdef TRAIN_OUT
        Qscore::writeAttCsvForQscoreTraining(deconvolved_spectrum, out_train_streams[ms_level - 1]);
#endif
      }

      if (report_decoy)
      {
        double noise_decoy_weight = fd.getNoiseDecoyWeight();
        for (auto& deconvolved_spectrum : deconvolved_spectra)
        {
          uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
          if (out_spec_file[ms_level - 1].empty())
            continue;
          if (!deconvolved_spectrum.isDecoy())
            continue;
          FLASHDeconvSpectrumFile::writeDeconvolvedMasses(deconvolved_spectrum, out_spec_streams[ms_level - 1], in_file, fd.getAveragine(), tols[ms_level - 1], write_detail, report_decoy,
                                                          noise_decoy_weight);
#ifdef TRAIN_OUT
          Qscore::writeAttCsvForQscoreTraining(deconvolved_spectrum, out_train_streams[ms_level - 1]);
#endif
        }
      }

      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || (!keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        out_spec_streams[i].close();
#ifdef TRAIN_OUT
        out_train_streams[i].close();
#endif
      }
    }
    // topFD feature output
    if (!out_topfd_feature_file.empty())
    {
      std::vector<fstream> out_topfd_feature_streams;
      out_topfd_feature_streams = std::vector<fstream>(out_topfd_feature_file.size());
      for (Size i = 0; i < out_topfd_feature_file.size(); i++)
      {
        if (out_topfd_feature_file[i].empty() || (!keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing topfd *.feature for MS level " << (i + 1) << " ..." << endl;

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
        if (out_topfd_file[i].empty() || (!keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing topfd *.msalign for MS level " << (i + 1) << " ..." << endl;

        out_topfd_streams[i].open(out_topfd_file[i], fstream::out);
        FLASHDeconvSpectrumFile::writeTopFDHeader(out_topfd_streams[i], getParam_().copy("SD:", true));
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (out_topfd_file[ms_level - 1].empty())
          continue;

        FLASHDeconvSpectrumFile::writeTopFD(deconvolved_spectrum, out_topfd_streams[ms_level - 1], in_file, topfd_snr_threshold, 1, per_ms_level_deconv_spec_count.begin()->first, false, false);
      }

      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || (!keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        out_topfd_streams[i].close();
      }
    }

    // isobaric quantification output
    if (!out_quant_file.empty())
    {
      OPENMS_LOG_INFO << "writing quantification tsv ..." << endl;
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
