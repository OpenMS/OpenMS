// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHTaggerFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFileInfo>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
@page

<B>The command line parameters of this tool are:</B>
@verbinclude TOPPFLASHTagger.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPPFLASHTagger.html
*/
class TOPPFLASHTagger : public TOPPBase
{
public:
  TOPPFLASHTagger() :
      TOPPBase("FLASHTagger", "FLASHTagger to generate de novo sequence tags from deconvolved spectrum.", false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (mzML)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("fasta", "<file>", "", "Input proteome database file (fasta)");
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_protein", "<file>", "", "Default output protein level tsv file containing matched proteins");
    setValidFormats_("out_protein", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_tag", "<file>", "", "Default output tag level tsv file containing matched tags");
    setValidFormats_("out_tag", ListUtils::create<String>("tsv"));

    registerDoubleOption_("min_mz", "<m/z value>", -1.0, "If set to positive value, minimum m/z to deconvolve.", false, true);
    registerDoubleOption_("max_mz", "<m/z value>", -1.0, "If set to positive value, maximum m/z to deconvolve.", false, true);
    registerDoubleOption_("min_rt", "<RT value>", -1.0, "If set to positive value, minimum RT (in second) to deconvolve.", false, true);
    registerDoubleOption_("max_rt", "<RT value>", -1.0, "If set to positive value, maximum RT (in second) to deconvolve.", false, true);

    registerIntOption_("max_ms_level", "<MS level>", -1.0, "If set to positive value, maximum MS level (inclusive) to deconvolve.", false, true);

    registerSubsection_("Tagger", "FLASHTagger algorithm parameters");
    registerSubsection_("FD", "FLASHDeconv algorithm parameters");
    registerSubsection_("SD", "Spectral deconvolution parameters");
  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "FD")
    {
      auto fd_param = FLASHDeconvAlgorithm().getDefaults();
      fd_param.remove("report_FDR");
      fd_param.remove("allowed_isotope_error");
      fd_param.remove("preceding_MS1_count");
      fd_param.remove("isolation_window");
      fd_param.remove("forced_MS_level");
      fd_param.remove("merging_method");
      fd_param.remove("ida_log");
      fd_param.removeAll("SD:");
      fd_param.removeAll("ft:");
      fd_param.removeAll("iq:");
      return fd_param;
    }
    else if (prefix == "SD")
    {
      auto fd_param = FLASHDeconvAlgorithm().getDefaults();
      auto sd_param = fd_param.copy("SD:", true);
      sd_param.remove("max_qvalue");
      return sd_param;
    }
    else if (prefix == "Tagger")
    {
      auto tagger_param = FLASHTaggerAlgorithm().getDefaults();
      return tagger_param;
    }
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
    String in_fasta = getStringOption_("fasta");

    String out_tag_file = getStringOption_("out_tag");
    String out_protein_file = getStringOption_("out_protein");

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

    tmp_fd_param = getParam_().copy("SD:", false);
    fd_param.insert("", tmp_fd_param);
    DoubleList tols = tmp_fd_param.getValue("SD:tol");

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
    // Run tagger
    FLASHTaggerAlgorithm tagger;

    auto tagger_param = getParam_().copy("Tagger:", true);
    if ((int)tagger_param.getValue("max_tag_count") > 0 && !deconvolved_spectra.empty() && tols.size() > 1)
    {
      OPENMS_LOG_INFO << "Finding sequence tags from deconvolved MS2 spectra ..." << endl;
      tagger.setParameters(tagger_param);

      tagger.run(deconvolved_spectra, tols[1]);
      tagger.runMatching(in_fasta);

      if (!out_protein_file.empty())
      {
        fstream out_tagger_stream = fstream(out_protein_file, fstream::out);
        FLASHTaggerFile::writeProteinHeader(out_tagger_stream);
        FLASHTaggerFile::writeProteins(tagger, out_tagger_stream);
        out_tagger_stream.close();
      }

      if (!out_tag_file.empty())
      {
        fstream out_tagger_stream = fstream(out_tag_file, fstream::out);
        FLASHTaggerFile::writeTagHeader(out_tagger_stream);
        FLASHTaggerFile::writeTags(tagger, out_tagger_stream);

        out_tagger_stream.close();
      }
    }

    OPENMS_LOG_INFO << "FLASHTagger run complete. Now writing the results in output files ..." << endl;

    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHTagger tool;
  return tool.main(argc, argv);
}
