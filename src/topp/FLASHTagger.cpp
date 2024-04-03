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
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHTaggerFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFileInfo>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
// We do not want this class to show up in the docu:

class TOPPFLASHTagger : public TOPPBase
{
public:
  TOPPFLASHTagger(): TOPPBase("FLASHTagger", "FLASHTagger to generate de novo sequence tags from TDP spectrum.", false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (deconv.mzML from FLASHDeconv mzML output)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("fasta", "<file>", "", "Input proteome database file (fasta)");
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_protein", "<file>", "", "Default output protein level tsv file containing matched proteins");
    setValidFormats_("out_protein", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_tag", "<file>", "", "Default output tag level tsv file containing matched tags");
    setValidFormats_("out_tag", ListUtils::create<String>("tsv"));

    registerSubsection_("Tagger", "FLASHTagger algorithm parameters");
  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "Tagger")
    {
      auto tagger_param = FLASHTaggerAlgorithm().getDefaults();
      return tagger_param;
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
    String out_protein_file = getStringOption_("out_protein");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.

    mzml.setLogType(log_type_);
    mzml.load(in_file, map);

    std::vector<DeconvolvedSpectrum> deconvolved_spectra;
    double tol(10);
    // Run FLASHDeconvAlgorithm here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    // collect statistics for information
    for (int index = 0; index < map.size(); index++)
    {
      auto spec = map[index];
      if (spec.getMSLevel() != 2) continue;
      int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);
      DeconvolvedSpectrum dspec(scan);
      dspec.setOriginalSpectrum(spec);
      String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();
      int tol_loc_s = deconv_meta_str.find("tol=") + 4;
      int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);
      tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));

      int q_loc_s = deconv_meta_str.find("qscore=") + 7;
      int q_loc_e = deconv_meta_str.find(";", q_loc_s);
      auto q_str = deconv_meta_str.substr(q_loc_s, q_loc_e - q_loc_s);
      Size pos = 0;
      std::vector<double> qscores;
      while (true)
      {
        Size pos_t = q_str.find(",", pos);
        if (pos_t == String::npos) break;
        auto token = q_str.substr(pos, pos_t - pos);
        qscores.push_back(stod(token));
        pos = pos_t + 1;
      }

      int s_loc_s = deconv_meta_str.find("snr=") + 4;
      int s_loc_e = deconv_meta_str.find(";", s_loc_s);
      auto s_str = deconv_meta_str.substr(s_loc_s, s_loc_e - s_loc_s);
      pos = 0;
      std::vector<float> snrs;
      while (true)
      {
        Size pos_t = s_str.find(",", pos);
        if (pos_t == String::npos) break;
        auto token = s_str.substr(pos, pos_t - pos);
        snrs.push_back(stof(token));
        pos = pos_t + 1;
      }

      for (int i = 0; i < spec.size(); i++)
      {
        PeakGroup peak;
        peak.setQscore(qscores[i]);
        peak.setSNR(snrs[i]);
        peak.setMonoisotopicMass(spec[i].getMZ());
        peak.setScanNumber(scan);
        dspec.push_back(peak);
      }
      dspec.sort();
      deconvolved_spectra.push_back(dspec);
    }
    // Run tagger
    FLASHTaggerAlgorithm tagger;

    auto tagger_param = getParam_().copy("Tagger:", true);

    if ((int)tagger_param.getValue("max_tag_count") > 0)
    {
      OPENMS_LOG_INFO << "Finding sequence tags from deconvolved MS2 spectra ..." << endl;
      tagger.setParameters(tagger_param);
      tagger.run(deconvolved_spectra, tol); //
      OPENMS_LOG_INFO << "Matching sequence tags against database ..." << endl;
      tagger.runMatching(in_fasta);
      OPENMS_LOG_INFO << "FLASHTagger run complete. Now writing the results in output files ..." << endl;

      if (! out_protein_file.empty())
      {
        fstream out_tagger_stream = fstream(out_protein_file, fstream::out);
        FLASHTaggerFile::writeProteinHeader(out_tagger_stream);
        FLASHTaggerFile::writeProteins(tagger, out_tagger_stream);
        out_tagger_stream.close();
      }

      if (! out_tag_file.empty())
      {
        fstream out_tagger_stream = fstream(out_tag_file, fstream::out);
        FLASHTaggerFile::writeTagHeader(out_tagger_stream);
        FLASHTaggerFile::writeTags(tagger, out_tagger_stream);

        out_tagger_stream.close();
      }
    }
    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHTagger tool;
  return tool.main(argc, argv);
}
