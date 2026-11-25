// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------
//#define TRAIN_OUT  //
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHugePages.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <limits>
#if defined(__GLIBC__)
  #include <malloc.h>
#endif
#ifdef TRAIN_OUT
  #include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>
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
  deconvolved spectra files (.tsv, .mzML, .msalign, .feature).
  FLASHDeconv uses SpectralDeconvolution for spectral level deconvolution and MassFeatureTracer to detect mass features.
  Also for MSn spectra, the precursor masses (not peak m/zs) should be determined and assigned in most cases. This assignment
  can be done by tracking MSn-1 spectra deconvolution information. Thus FLASHDeconv class keeps MSn-1 spectra deconvolution information
  for a certain period for precursor mass assignment in DeconvolvedSpectrum class.
  In case of FLASHIda runs, FLASHIda decided which precursor each MSn scan was acquired for. Given the run's scan_commands.tsv
  (-FD:scan_commands), FLASHDeconv joins each spectrum to its command by the tracking id in the spectrum's scan description, and
  takes its own deconvolved mass for the commanded precursor from the MS1 spectrum that command was decided from.

See https://openms.de/FLASHDeconv for more information.

On Linux with glibc, FLASHDeconv (like FLASHTnT) re-executes itself once at start-up to put its heap on transparent huge pages
(GLIBC_TUNABLES=glibc.malloc.hugetlb=1; ~20 % faster at high thread counts, identical output). Set FLASH_NO_HUGETLB to prevent this. With
OpenMP thread binding (OMP_PROC_BIND, OMP_PLACES, GOMP_CPU_AFFINITY, KMP_AFFINITY) it is skipped; set GLIBC_TUNABLES yourself there.


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
               {Citation {"Jeong K, Kim J, Gaikwad M et al.", "FLASHDeconv: Ultrafast, High-Quality Feature Deconvolution for Top-Down Proteomics",
                          "Cell Syst 2020 Feb 26;10(2):213-218.e6", "10.1016/j.cels.2020.01.003"}})
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "",
                       "Input file: mzML (profile or centroid), or a Thermo .raw file (recognised by its extension; read with Thermo's peak "
                       "picking, so centroided, as ThermoRawFileParser converts it by default; needs OpenMS built with WITH_THERMO_RAW and the "
                       ".NET 8 runtime). "
                       "A Bruker timsTOF .d directory, or a zipped one (.d.zip), is recognised by its name and read with OpenMS's timsTOF reader (DDA-PASEF or DIA-PASEF): the spectra are summed over the ion mobility scans and centroided, the spectra are numbered by retention time (scan=1, 2, ...), and the ion mobility values are not kept (needs OpenMS built with WITH_OPENTIMS, the default).");
    setValidFormats_("in", ListUtils::create<std::string>("mzML,raw,d"));

    registerOutputFile_("out", "<file>", "", "Default output tsv file containing deconvolved features");
    setValidFormats_("out", ListUtils::create<std::string>("tsv"));

    registerOutputFile_(
      "out_spec1", "<file>", "",
      "Output tsv file for deconvolved MS1 spectra. Use -out_spec2, ..., -out_spec4 for MS2, ..., MS4 spectra.", false);
    setValidFormats_("out_spec1", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_spec2", "<file>", "", "Output TSV files for deconvolved MS2 spectra.", false, true);
    setValidFormats_("out_spec2", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_spec3", "<file>", "", "Output TSV files for deconvolved MS3 spectra.", false, true);
    setValidFormats_("out_spec3", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_spec4", "<file>", "", "Output TSV files for deconvolved MS4 spectra.", false, true);
    setValidFormats_("out_spec4", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_mzml", "<file>", "", "Output mzML file containing deconvolved spectra (for all MS levels).", false);
    setValidFormats_("out_mzml", ListUtils::create<std::string>("mzML"));

    registerOutputFile_("out_quant", "<file>", "", "Output tsv file with isobaric quantification results for MS2 spectra.", false);
    setValidFormats_("out_quant", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_annotated_mzml", "<file>", "",
                        "Output annotated mzML file: the input spectra that have deconvolved masses (all their peaks), each with the meta value "
                        "DeconvMassPeakIndices listing its masses as monoisotopic mass:indices of its peaks (in the m/z-sorted spectrum), separated by ';'.",
                        false);
    setValidFormats_("out_annotated_mzml", ListUtils::create<std::string>("mzML"));

    registerOutputFile_(
      "out_msalign1", "<file>", "",
      "Output msalign (TopFD and ProMex compatible) file for MS1 deconvolved spectra. Ensure filename ends with ms1.msalign for TopPIC GUI compatibility (e.g., result_ms1.msalign; refer to TopPIC input formats).",
      false);
    setValidFormats_("out_msalign1", ListUtils::create<std::string>("msalign"), false);

    registerOutputFile_("out_msalign2", "<file>", "",
                        "Output msalign (TopFD and ProMex compatible) file for MS2 deconvolved spectra. Ensure filename ends with ms2.msalign for TopPIC GUI compatibility (e.g., result_ms2.msalign; refer to TopPIC input formats). "
                        "TopPIC (1.8 to 2.1) reads it as TopFD's: run TopPIC with -D (--disable-post-match), since FLASHDeconv writes no TopFD SQLite database, "
                        "and give -out_feature2 as <name>_ms2.feature next to <name>_ms2.msalign (or run TopPIC with -x).",
                        false, true);
    setValidFormats_("out_msalign2", ListUtils::create<std::string>("msalign"), false);

    registerOutputFile_("out_feature1", "<file>", "",
                        "Output feature file (TopFD compatible) for MS1 spectra. It is needed for TopPIC feature intensity output (refer to TopPIC input formats).",
                        false);

    setValidFormats_("out_feature1", ListUtils::create<std::string>("feature"), false);

    registerOutputFile_("out_feature2", "<file>", "",
                        "Output feature file (TopFD compatible) for MS2 spectra. It is needed for TopPIC feature intensity output (refer to TopPIC input formats).",
                        false, true);

    setValidFormats_("out_feature2", ListUtils::create<std::string>("feature"), false);

    registerFlag_("keep_empty_out", "Retain empty output files (e.g., *.tsv files with no features).");

    registerIntOption_("mzml_mass_charge", "<0:uncharged 1: +1 charged -1: -1 charged>", 0,
                       "Charge state of deconvolved masses in mzML output specified by -out_mzml.", false, true);
    setMinInt_("mzml_mass_charge", -1);
    setMaxInt_("mzml_mass_charge", 1);

    registerFlag_("quant",
                  "Quantify every MS1 mass from its raw signal (FLASHQuant). Input: the MS1 spectra of -in as loaded, centroided (as for "
                  "FLASHDeconv itself). The deconvolved masses are the starting point: the masses of one feature, or a mass without a feature, "
                  "take the raw mass traces their peaks lie on (FLASHQuant's mass trace and elution peak detection), and the traces no other mass "
                  "claims at their other charges and isotopes at their apex; traces that several masses share are resolved as in FLASHQuant (elution model and isotope pattern of each mass, non-negative least squares). Output: the "
                  "columns Abundance and AbundanceTraceCount at the end of the -out feature tsv and of the -out_spec1 tsv (one value per feature, "
                  "and each MS1 mass the value of its feature): FLASHQuant's quantity (area of each trace down to 10% of its apex, summed; "
                  "intensity x seconds; AbundanceTraceCount: the traces summed). A feature without a raw trace of at least 3 peaks, or whose traces "
                  "the resolution gave to other masses, gets its FeatureQuantity on the scale of the run's abundances (trace count 0; NaN when "
                  "no feature of the run has a raw abundance to set the scale); a mass "
                  "without a feature and without a raw trace (a peak seen in one scan) gets NaN. Either ion mode (that of the masses). With "
                  "fewer than 3 MS1 spectra with peaks there are no raw traces: every Abundance is NaN (with a warning); MS1 spectra without peaks "
                  "(e.g. every other scan) are left out of the mass traces. Without the flag, nothing changes: "
                  "the same masses, features and output files as before.");

    registerFlag_("write_detail",
                  "Include detailed peak information (m/z, intensity, charge, isotope index) for each deconvolved mass in the output spectrum tsv files specified by out_spec* options.",
                  false);

    //registerDoubleOption_("precursor_snr", "<snr value>", 1.0, "Precursor SNR threshold for TopFD MS2 msalign TSV files.", false, true);

    registerDoubleOption_("min_mz", "<m/z value>", -1.0, "Specify the minimum m/z values for peaks considered during deconvolution. Negative values disable the threshold.", false, true);
    registerDoubleOption_("max_mz", "<m/z value>", -1.0, "Specify the maximum m/z values for peaks considered during deconvolution. Negative values disable the threshold.", false, true);
    registerDoubleOption_("min_rt", "<RT value>", -1.0, "Specify the minimum retention time (in minutes) for spectra considered during deconvolution. Negative values disable the threshold.", false, true);
    registerDoubleOption_("max_rt", "<RT value>", -1.0, "Specify the maximum retention time (in minutes) for spectra considered during deconvolution. Negative values disable the threshold.", false, true);

    registerIntOption_("max_ms_level", "<MS level>", -1.0, "Set the maximum MS level (inclusive) for deconvolution. Negative values disable the threshold.", false, true);

    registerSubsection_("FD", "FLASHDeconv algorithm parameters");
    registerSubsection_("SD", "Spectral deconvolution parameters");
    registerSubsection_("ft", "Feature tracing parameters");
    registerSubsection_("iq", "Isobaric quantification parameters");
  }

  Param getSubsectionDefaults_(const std::string& prefix) const override
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
    else { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix); }
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**) override
  {
    OPENMS_LOG_INFO << "Initializing ... " << endl;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    std::string in_file = getStringOption_("in");
    std::string out_file = getStringOption_("out");
    // String in_train_file {};
    bool keep_empty_out = getFlag_("keep_empty_out");
    auto out_spec_file
      = StringList {getStringOption_("out_spec1"), getStringOption_("out_spec2"), getStringOption_("out_spec3"), getStringOption_("out_spec4")};

    auto out_topfd_file = StringList {getStringOption_("out_msalign1"), getStringOption_("out_msalign2")};
    auto out_topfd_feature_file = StringList {getStringOption_("out_feature1"), getStringOption_("out_feature2")};

    std::string out_mzml_file = getStringOption_("out_mzml");
    std::string out_anno_mzml_file = getStringOption_("out_annotated_mzml");
    std::string out_quant_file = getStringOption_("out_quant");

    bool write_detail = getFlag_("write_detail");
    const bool quant = getFlag_("quant");
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
    Param tmp_fd_param = getParam_().copy("FD:", true);
    Param fd_param;
    fd_param.insert("", tmp_fd_param);
    // precursor_FDR makes the decoy masses (its decoy precursors) and turns report_FDR on
    bool report_decoy = tmp_fd_param.getValue("report_FDR") != "false" || tmp_fd_param.getValue("precursor_FDR") == "true";

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
      if (min_rt > 0 && max_rt < 0) max_rt = 1e7;
      opt.setRTRange(DRange<1> {min_rt, max_rt});
    }
    if (min_mz > 0 || max_mz > 0)
    {
      if (min_mz > 0 && max_mz < 0) max_mz = 1e7;
      opt.setMZRange(DRange<1> {min_mz, max_mz});
    }
    if (max_ms_level > 0)
    {
      IntList ms_levels;
      for (int msl = 1; msl <= max_ms_level; msl++)
        ms_levels.push_back(msl);
      opt.setMSLevels(ms_levels);
    }

    FLASHHelperClasses::loadExperiment(in_file, map, opt, log_type_);

    std::vector<DeconvolvedSpectrum> deconvolved_spectra;
    std::vector<FLASHHelperClasses::MassFeature> deconvolved_features;
    std::map<int, double> scan_rt_map;

    // -quant: the MS1 spectra as loaded (FLASHDeconv removes low peaks from map in place)
    PeakMap ms1_raw;
    if (quant)
    {
      for (const auto& spec : map)
      {
        if (spec.getMSLevel() == 1) ms1_raw.addSpectrum(spec);
      }
      ms1_raw.sortSpectra(true);
    }

    // Run FLASHDeconvAlgorithm here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    fd.run(map, deconvolved_spectra, deconvolved_features);
    tols = fd.getTolerances();

    std::vector<std::vector<FLASHQuantAlgorithm::MassAbundance>> mass_abundances;
    std::map<uint, FLASHQuantAlgorithm::MassAbundance> feature_abundances;
    if (quant)
    {
      OPENMS_LOG_INFO << "Quantifying the MS1 masses from their raw signal (FLASHQuant) ..." << endl;
      FLASHQuantAlgorithm fq;
      fq.setLogType(log_type_);
      fq.runFromDeconvolvedMasses(ms1_raw, deconvolved_spectra, fd.getAveragine(), mass_abundances, feature_abundances, &deconvolved_features);
      ms1_raw.clear(true);
    }
    // collect statistics for information
    for (auto& it : map)
    {
      uint ms_level = it.getMSLevel();
      if (per_ms_level_spec_count.find(ms_level) == per_ms_level_spec_count.end()) per_ms_level_spec_count[ms_level] = 0;
      per_ms_level_spec_count[ms_level]++;
    }

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      scan_rt_map[deconvolved_spectrum.getScanNumber()] = deconvolved_spectrum.getOriginalSpectrum().getRT();

      if (deconvolved_spectrum.empty()) continue;
      if (per_ms_level_deconv_spec_count.find(ms_level) == per_ms_level_deconv_spec_count.end()) per_ms_level_deconv_spec_count[ms_level] = 0;
      if (per_ms_level_mass_count.find(ms_level) == per_ms_level_mass_count.end()) per_ms_level_mass_count[ms_level] = 0;

      per_ms_level_deconv_spec_count[ms_level]++;
      per_ms_level_mass_count[ms_level] += (int)deconvolved_spectrum.size();

    }
    for (auto& val : per_ms_level_deconv_spec_count)
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
      fstream out_stream;
      out_stream.open(out_file, fstream::out);

      FLASHDeconvFeatureFile::writeHeader(out_stream, report_decoy, quant);
      FLASHDeconvFeatureFile::writeFeatures(deconvolved_features, in_file, out_stream, report_decoy, quant ? &feature_abundances : nullptr);
      out_stream.close();
    }
    // Per ms level spectrum deconvolution tsv output
    if (! out_spec_file.empty())
    {
      std::vector<fstream> out_spec_streams = std::vector<fstream>(out_spec_file.size());
#ifdef TRAIN_OUT
      std::vector<fstream> out_train_streams = std::vector<fstream>(out_spec_file.size());
#endif

      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing spectrum tsv for MS level " << (i + 1) << " ..." << endl;

        out_spec_streams[i].open(out_spec_file[i], fstream::out);

        FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(out_spec_streams[i], i + 1, write_detail, report_decoy, quant && i == 0);
#ifdef TRAIN_OUT
        out_train_streams[i].open(out_spec_file[i] + "_train.csv", fstream::out);
        PeakGroupScoring::writeAttCsvForQscoreTrainingHeader(out_train_streams[i]);
#endif
      }


#ifdef TRAIN_OUT
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (out_spec_file[ms_level - 1].empty()) continue;
        FLASHDeconvSpectrumFile::writeDeconvolvedMasses(deconvolved_spectrum, out_spec_streams[ms_level - 1], in_file, fd.getAveragine(), fd.getDecoyAveragine(),
                                                        tols[ms_level - 1], write_detail, report_decoy, fd.getNoiseDecoyWeight(ms_level));
        PeakGroupScoring::writeAttCsvForQscoreTraining(deconvolved_spectrum, out_train_streams[ms_level - 1]);
      }
#else
      // the same rows in the same order, formatted in parallel (with report_FDR, the noise-decoy carries are replayed serially up front)
      FLASHDeconvSpectrumFile::writeDeconvolvedMasses(deconvolved_spectra, out_spec_streams, out_spec_file, in_file, fd.getAveragine(), fd.getDecoyAveragine(), tols,
                                                      write_detail, report_decoy, fd.getNoiseDecoyWeights((uint)tols.size()), quant ? &mass_abundances : nullptr);
#endif

      for (Size i = 0; i < out_spec_file.size(); i++)
      {
        if (out_spec_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        out_spec_streams[i].close();
#ifdef TRAIN_OUT
        out_train_streams[i].close();
#endif
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
      fstream out_quant_stream;
      out_quant_stream.open(out_quant_file, fstream::out);
      FLASHDeconvSpectrumFile::writeIsobaricQuantification(out_quant_stream, deconvolved_spectra);
      out_quant_stream.close();
    }

    // topFD feature output - it should be at the end since zero feature IDs are redefined for TopPIC feature indices.
    if (! out_topfd_feature_file.empty())
    {
      std::vector<fstream> out_topfd_feature_streams;
      out_topfd_feature_streams = std::vector<fstream>(out_topfd_feature_file.size());
      for (Size i = 0; i < out_topfd_feature_file.size(); i++)
      {
        if (out_topfd_feature_file[i].empty()
            || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing topfd *.feature for MS level " << (i + 1) << " ..." << endl;

        out_topfd_feature_streams[i].open(out_topfd_feature_file[i], fstream::out);
        FLASHDeconvFeatureFile::writeTopFDFeatureHeader(out_topfd_feature_streams[i], i + 1);
        FLASHDeconvFeatureFile::writeTopFDFeatures(deconvolved_spectra, deconvolved_features, scan_rt_map, in_file, out_topfd_feature_streams[i],
                                                   i + 1);
        out_topfd_feature_streams[i].close();
      }
    }
    // topFD msalign output
    if (! out_topfd_file.empty())
    {
      auto out_topfd_streams = std::vector<fstream>(out_topfd_file.size());

      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        OPENMS_LOG_INFO << "writing topfd *.msalign for MS level " << (i + 1) << " ..." << endl;

        out_topfd_streams[i].open(out_topfd_file[i], fstream::out);
        FLASHDeconvSpectrumFile::writeTopFDHeader(out_topfd_streams[i], getParam_().copy("SD:", true));
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
        if (out_topfd_file[ms_level - 1].empty()) continue;

        FLASHDeconvSpectrumFile::writeTopFD(deconvolved_spectrum, out_topfd_streams[ms_level - 1], in_file, 1,
                                            per_ms_level_deconv_spec_count.begin()->first);
      }

      for (Size i = 0; i < out_topfd_file.size(); i++)
      {
        if (out_topfd_file[i].empty() || (! keep_empty_out && per_ms_level_deconv_spec_count.find(i + 1) == per_ms_level_deconv_spec_count.end()))
          continue;
        out_topfd_streams[i].close();
      }
    }

    // All outputs are written: the deconvolved spectra (their peak groups) and the input spectra are released in parallel instead of one
    // after the other by the destructors at the end of this function (each element is moved into a local that is destroyed).
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < (SignedSize)deconvolved_spectra.size(); i++)
    {
      DeconvolvedSpectrum released(std::move(deconvolved_spectra[i]));
    }
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < (SignedSize)map.size(); i++)
    {
      MSSpectrum released(std::move(map[i]));
    }

    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  useTransparentHugePages(argv);
#if defined(__GLIBC__)
  // Keep freed heap memory in the per-thread arenas instead of trimming it back to the kernel after every
  // spectrum: trim, refault and heap regrowth serialise on the process-wide mmap lock at high thread counts.
  mallopt(M_TRIM_THRESHOLD, std::numeric_limits<int>::max());
  mallopt(M_TOP_PAD, 256 << 20);
#endif
  TOPPFLASHDeconv tool;
  return tool.main(argc, argv);
}