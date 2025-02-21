// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFileInfo>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
// We do not want this class to show up in the docu:

class TOPPFLASHTnT : public TOPPBase
{
public:
  TOPPFLASHTnT():
      TOPPBase("FLASHTnT",
               "FLASHTnT to generate de novo sequence tags from TDP spectrum and match them against proteome DB for proteoform identification.",
               false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (deconv.mzML from FLASHDeconv mzML output).", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("fasta", "<file>", "", "Input proteome database file (fasta)", true);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_pro", "<file>", "", "Output proteoform-level tsv file containing proteoform IDs. ", true);
    setValidFormats_("out_pro", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_prsm", "<file>", "", "Output PrSM-level tsv file containing PrSMs.");
    setValidFormats_("out_prsm", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_tag", "<file>", "", "Output tag-level tsv file containing matched tags.");
    setValidFormats_("out_tag", ListUtils::create<String>("tsv"));

    // Register PrSM-level FDR parameter
    registerDoubleOption_("prsm_fdr", "Specifies the PrSM-level FDR.", 1.0, "Specifies the PrSM-level FDR.", false);
    setMinFloat_("prsm_fdr", 0.0);

    // Register proteoform-level FDR parameter
    registerDoubleOption_("pro_fdr", "Specifies the proteoform-level FDR.", 1.0, "Specifies the proteoform-level FDR.", false);
    setMinFloat_("pro_fdr", 0.0);

    // Register single-hit-only option
    registerFlag_("only_single_hit", "Allows only a single hit per spectrum.", "Allows only a single hit per spectrum.");
    //setValidStrings_("only_single_hit", {"true", "false"});

    // Register underdetermined proteoform discard option
    registerFlag_("discard_underdetermined",
                          "Discards underdetermined proteoform IDs (e.g., those without exact precursor masses or start/end positions).");
    //setValidStrings_("discard_underdetermined", {"true", "false"});

    // Register decoy retention option
    registerFlag_("keep_decoy", "Retains decoy hits in the results.");
    //setValidStrings_("keep_decoy", {"true", "false"});

    // Register ion type parameter
    registerStringList_("ion_type", "Specifies ion types to consider.", ListUtils::toStringList(std::vector<std::string>{"b", "y"}), "Specifies ion types to consider.", false);
    setValidStrings_("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

    registerSubsection_("tag", "Tag algorithm parameters");
    registerSubsection_("ex", "Extension algorithm parameters");

  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "tag")
    {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("tag:", true);
    }
    else if (prefix == "ex") {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("ex:", true);
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
    String in_fasta = getStringOption_("fasta");

    String out_tag_file = getStringOption_("out_tag");
    String out_prsm_file = getStringOption_("out_prsm");
    String out_pro_file = getStringOption_("out_pro");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.

    mzml.setLogType(log_type_);
    mzml.load(in_file, map);

    auto tnt_param = getParam_();
    tnt_param.remove("in");
    tnt_param.remove("fasta");
    tnt_param.remove("out_tag");
    tnt_param.remove("out_prsm");
    tnt_param.remove("out_pro");
    tnt_param.remove("log");
    tnt_param.remove("debug");
    tnt_param.remove("threads");
    tnt_param.remove("no_progress");
    tnt_param.remove("force");
    tnt_param.remove("test");
    double max_mod_mass = tnt_param.getValue("ex:max_mod_mass");
    int max_mod_count = tnt_param.getValue("ex:max_mod_count");
    double pro_fdr = tnt_param.getValue("pro_fdr");
    OPENMS_LOG_INFO << "Finding sequence tags from deconvolved spectra ..." << endl;

    FASTAFile fasta_file;
    std::vector<FASTAFile::FASTAEntry> fasta_entry;
    fasta_file.load(in_fasta, fasta_entry);
    std::vector<ProteinHit> proteoform_hits;
    // Run here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    fstream out_tagger_stream;
    fstream out_prsm_stream;
    fstream out_pro_stream;

    if (! out_tag_file.empty())
    {
      out_tagger_stream = fstream(out_tag_file, fstream::out);
      FLASHTnTFile::writeTagHeader(out_tagger_stream);
    }

    if (! out_prsm_file.empty())
    {
      out_prsm_stream = fstream(out_prsm_file, fstream::out);
      FLASHTnTFile::writePrSMHeader(out_prsm_stream);
    }

    if (! out_pro_file.empty())
    {
      out_pro_stream = fstream(out_pro_file, fstream::out);
      FLASHTnTFile::writeProHeader(out_pro_stream);
    }

    FLASHTnTAlgorithm tnt;
    tnt.setParameters(tnt_param);
    tnt.run(map, fasta_entry);
    tnt.getProteoforms(proteoform_hits);

    OPENMS_LOG_INFO << "FLASHTnT run complete. Now writing the results in output files ..." << endl;
    if (! out_tag_file.empty())
    {
      FLASHTnTFile::writeTags(tnt, max_mod_count * max_mod_mass + 1, out_tagger_stream);
      out_tagger_stream.close();
    }
    if (! out_prsm_file.empty())
    {
      FLASHTnTFile::writePrSMs(proteoform_hits, out_prsm_stream);
      out_prsm_stream.close();
    }
    if (! out_pro_file.empty())
    {
      FLASHTnTFile::writeProteoforms(proteoform_hits, out_pro_stream, pro_fdr);
      out_pro_stream.close();
    }

    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHTnT tool;
  return tool.main(argc, argv);
}
