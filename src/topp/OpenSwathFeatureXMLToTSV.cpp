// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

#include <fstream>
#include <clocale>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathFeatureXMLToTSV OpenSwathFeatureXMLToTSV

@brief Converts a featureXML to a mProphet tsv

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; OpenSwathFeatureXMLToTSV &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> Downstream data analysis </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathConfidenceScoring </td>
        </tr>
    </table>
</CENTER>

Creates a tsv that is compatible as input to mProphet.
Furthermore it creates the columns "decoy" and
"transition_group_id" which are required by mProphet.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathFeatureXMLToTSV.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathFeatureXMLToTSV.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

std::map<String, std::vector<const ReactionMonitoringTransition *> > peptide_transition_map;

void write_out_header(std::ostream &os, FeatureMap &feature_map, /* String main_var_name,  */ std::vector<String> &meta_value_names, bool short_format)
{
  std::vector<String> meta_value_names_tmp;

  os <<  "transition_group_id"  << "\t"
     <<  "run_id"  << "\t"
     <<  "filename"  << "\t"
     <<  "RT"  << "\t"
     <<  "id"  << "\t"
     <<  "Sequence"  << "\t"
     <<  "FullPeptideName"  << "\t"
     <<  "Charge"  << "\t"
     <<  "m/z"  << "\t"
     <<  "Intensity"  << "\t"
     <<  "ProteinName"  << "\t"
     <<  "decoy"  << "\t";

  // get all meta values from the first feature
  feature_map[0].getKeys(meta_value_names_tmp);
  for (Size i = 0; i < meta_value_names_tmp.size(); i++)
  {
    if (meta_value_names_tmp[i] != "PeptideRef" && meta_value_names_tmp[i] != "PrecursorMZ")
    {
      meta_value_names.push_back(meta_value_names_tmp[i]);
    }
  }
  std::sort(meta_value_names.begin(), meta_value_names.end());
  for (Size i = 0; i < meta_value_names.size(); i++)
  {
    os << meta_value_names[i] << "\t";
  }

  if (!short_format)
  {
    os << "Peak_Area" << "\t";
    os << "Peak_Apex" << "\t";
    os << "Fragment_Annotation" << "\t";
    os << "ProductMZ";
  }
  else
  {
    os << "aggr_Peak_Area" << "\t";
    os << "aggr_Peak_Apex" << "\t";
    os << "aggr_Fragment_Annotation";
  }
  os << std::endl;
}

void write_out_body_(std::ostream &os, Feature *feature_it, TargetedExperiment &transition_exp,
                     std::vector<String> &meta_value_names, int run_id, bool short_format, String identifier, String filename)
{

  String peptide_ref = feature_it->getMetaValue("PeptideRef");
  String precursor_mz = feature_it->getMetaValue("PrecursorMZ");

  String sequence;
  String full_peptide_name = "NA";
  String protein_name = "NA";
  String decoy = "NA";
  String charge = "NA";

  if (!transition_exp.hasPeptide(peptide_ref))
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Did not find the peptide " + peptide_ref + " in the targeted experiment.");
  }

  const OpenMS::TargetedExperiment::Peptide &pep = transition_exp.getPeptideByRef(peptide_ref);

  sequence = pep.sequence;
  if (!pep.protein_refs.empty())
  {
    // For now just take the first one, assuming the protein name is the id
    protein_name = pep.protein_refs[0];
  }

  // handle charge
  if (pep.hasCVTerm("MS:1000041"))
  {
    charge = pep.getCVTerms().at("MS:1000041")[0].getValue().toString();
  }
  else if (pep.hasCharge())
  {
    charge = (String)pep.getChargeState();
  }
  if (charge == "NA" && !full_peptide_name.empty())
  {
    // deal with FullPeptideNames like PEPTIDE/2
    std::vector<String> substrings;
    full_peptide_name.split("/", substrings);
    if (substrings.size() == 2)
    {
      charge = substrings[1];
    }
  }

  // handle decoy tag
  if (peptide_transition_map.find(peptide_ref) != peptide_transition_map.end() && !peptide_transition_map[peptide_ref].empty())
  {
    const ReactionMonitoringTransition *transition = peptide_transition_map[peptide_ref][0];
#if 1
    const auto& terms = transition->getCVTerms();
    if (terms.find("decoy") != terms.end())
    {
      decoy = transition->getCVTerms().at("decoy")[0].getValue().toString();
    }
    else if (terms.find("MS:1002007") != terms.end() )    // target SRM transition
    {
      decoy = "0";
    }
    else if (terms.find("MS:1002008") != terms.end() )    // decoy SRM transition
    {
      decoy = "1";
    }
    else if (terms.find("MS:1002007") != terms.end() && terms.find("MS:1002008") != terms.end() )    // both == illegal
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Peptide " + peptide_ref + " cannot be target and decoy at the same time.");
    }
    else
#endif
    if (transition->getDecoyTransitionType() == ReactionMonitoringTransition::UNKNOWN)
    {
      // assume its target
      decoy = "0";
    }
    else if (transition->getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
    {
      decoy = "0";
    }
    else if (transition->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
    {
      decoy = "1";
    }
    else
    {
      // assume its target
      decoy = "0";
    }
  }
  else
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Did not find the peptide " + peptide_ref + " in the targeted experiment.");
  }

  if (pep.metaValueExists("full_peptide_name"))
  {
    full_peptide_name = pep.getMetaValue("full_peptide_name");
  }

  // adjust peptide ref with current file identifier
  peptide_ref += "_" + identifier;

  String line = "";
  // Start writing out
  line += peptide_ref + "\t" + (String)run_id + "\t" + (String)filename + "\t" + feature_it->getRT() + "\tf_" + feature_it->getUniqueId() + "\t";
  line += sequence + "\t" + full_peptide_name + "\t";
  line += (String)charge + "\t";
  line += precursor_mz + "\t";
  line += (String)feature_it->getIntensity() + "\t";
  line += protein_name + "\t";
  line += decoy + "\t";

  String meta_values = "";
  for (Size i = 0; i < meta_value_names.size(); i++)
  {
    meta_values += feature_it->getMetaValue(meta_value_names[i]).toString() + "\t";
  }

  // Write out the individual transition
  if (short_format)
  {
    String aggr_Peak_Area = "";
    String aggr_Peak_Apex = "";
    String aggr_Fragment_Annotation = "";
    for (Feature& sub : feature_it->getSubordinates())
    {
      aggr_Peak_Area += String(sub.getIntensity()) + ";";

      if (sub.metaValueExists("peak_apex_int"))
      {
        aggr_Peak_Apex += String((double)sub.getMetaValue("peak_apex_int")) + ";";
      }
      else
      {
        aggr_Peak_Apex += "NA;";
      }

      aggr_Fragment_Annotation += (String)sub.getMetaValue("native_id") + ";";
    }

    // remove the last semicolon
    if (!feature_it->getSubordinates().empty())
    {
      aggr_Peak_Area = aggr_Peak_Area.substr(0, aggr_Peak_Area.size() - 1);
      aggr_Peak_Apex = aggr_Peak_Apex.substr(0, aggr_Peak_Apex.size() - 1);
      aggr_Fragment_Annotation = aggr_Fragment_Annotation.substr(0, aggr_Fragment_Annotation.size() - 1);
    }
    os << line << meta_values << aggr_Peak_Area << "\t" << aggr_Peak_Apex << "\t" << aggr_Fragment_Annotation << std::endl;
  }
  else
  {
    for (Feature& sub : feature_it->getSubordinates())
    {
      os.precision(writtenDigits(double()));
      String apex = "NA";
      if (sub.metaValueExists("peak_apex_int"))
      {
        apex = String((double)sub.getMetaValue("peak_apex_int"));
      }
      os << line << meta_values << String(sub.getIntensity()) << "\t" << apex << "\t" << (String)sub.getMetaValue("native_id") << "\t" << String(sub.getMZ()) << std::endl;
    }
  }
}

Feature *find_best_feature(const std::vector<Feature *> &features, String score_)
{
  double best_score = -std::numeric_limits<double>::max();
  Feature *best_feature = nullptr;

  for (Size i = 0; i < features.size(); i++)
  {
    double  score = features[i]->getMetaValue(score_).toString().toDouble();
    if (score > best_score)
    {
      best_feature = features[i];
      best_score = score;
    }
  }
  return best_feature;
}

void write_out_body_best_score(std::ostream &os, FeatureMap &feature_map,
                               TargetedExperiment &transition_exp, std::vector<String> &meta_value_names,
                               int run_id, bool short_format, String best_score, String filename)
{

  // for each peptide reference search for the best feature
  typedef std::map<String, std::vector<Feature *> > PeptideFeatureMapType;
  PeptideFeatureMapType peptide_feature_map;
  for (Feature& feature : feature_map)
  {
    String peptide_ref = feature.getMetaValue("PeptideRef");
    peptide_feature_map[peptide_ref].push_back(&feature);
  }

  for (PeptideFeatureMapType::iterator peptide_it = peptide_feature_map.begin(); peptide_it != peptide_feature_map.end(); ++peptide_it)
  {
    if (peptide_it->second.size() > 1)
    {
      //std::cout << "Warning " << peptide_it->first << " has " << peptide_it->second.size() << " features!" << std::endl;
      // for (Size j =0; j<peptide_it->second.size(); j++)
      //     std::cout <<  *peptide_it->second.at(j) << std::endl;
    }
  }

  for (PeptideFeatureMapType::iterator peptide_it = peptide_feature_map.begin(); peptide_it != peptide_feature_map.end(); ++peptide_it)
  {
    Feature *bestfeature = find_best_feature(peptide_it->second, best_score);
    if (bestfeature == nullptr)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Did not find best feature for peptide " + peptide_it->first);
    }
    write_out_body_(os, bestfeature, transition_exp, meta_value_names, run_id, short_format, feature_map.getIdentifier(), filename);
  }
}


class TOPPOpenSwathFeatureXMLToTSV
: public TOPPBase, public ProgressLogger
{
public:

  TOPPOpenSwathFeatureXMLToTSV() :
    TOPPBase("OpenSwathFeatureXMLToTSV", "Converts a featureXML to a mProphet tsv.", true)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("featureXML"));

    registerInputFile_("tr", "<file>", "", "TraML transition file");
    setValidFormats_("tr", ListUtils::create<String>("traML"));
    //registerStringOption_("main_var_name","<varname>","xx_lda_prelim_score","Name of the main variable", false);

    registerOutputFile_("out", "<file>", "", "tsv output file (mProphet compatible)");
    setValidFormats_("out", ListUtils::create<String>("csv"));

    registerFlag_("short_format", "Whether to write short (one peptide per line) or long format (one transition per line).");

    registerStringOption_("best_scoring_peptide", "<varname>", "", "If only the best scoring feature per peptide should be printed, give the variable name", false);
  }

  void write_out_body(std::ostream &os, FeatureMap &feature_map,
                      TargetedExperiment &transition_exp, std::vector<String> &meta_value_names,
                      int run_id, bool short_format, String filename)
  {

    Size progress = 0;
    startProgress(0, feature_map.size(), "writing out features");
    for (Feature& feature : feature_map)
    {
      setProgress(progress++);
      write_out_body_(os, &feature, transition_exp, meta_value_names, run_id, short_format, feature_map.getIdentifier(), filename);
    }
    endProgress();
  }

  ExitCodes main_(int, const char **) override
  {
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");
    String out = getStringOption_("out");
    //String main_var_name = getStringOption_("main_var_name");
    String best_scoring = getStringOption_("best_scoring_peptide");
    bool short_format = getFlag_("short_format");

    setLogType(log_type_);

    TargetedExperiment transition_exp;
    FileHandler().loadTransitions(tr_file, transition_exp, {FileTypes::TRAML});

    startProgress(0, transition_exp.getTransitions().size(), "indexing transitions peaks");
    for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
    {
      setProgress(i);
      const ReactionMonitoringTransition *transition = &transition_exp.getTransitions()[i];

      {
        peptide_transition_map[transition->getPeptideRef()].push_back(&transition_exp.getTransitions()[i]);
      }
    }
    endProgress();

    std::ofstream os(out.c_str());
    //set high precision for writing of floating point numbers
    os.precision(writtenDigits(double()));
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out);
    }

    // write the csv header (we need to know which parameters are in the map to do that)
    if (file_list.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No input files given ");
    }
    FeatureMap feature_map;
    // feature_file.load() resets the locale to the user's (Don't know where, maybe QT or Xerces)
    // Somehow even our variable OpenMS::Internal::OpenMS_locale is overwritten
    // Create copy here and reset it later. TODO this needs to be fixed more thouroughly.
    String locale_before = String(OpenMS::Internal::OpenMS_locale);
    FileHandler().loadFeatures(file_list[0], feature_map, {FileTypes::FEATUREXML}, log_type_);
    setlocale(LC_ALL, locale_before.c_str());
    if (feature_map.getIdentifier().empty())
    {
      feature_map.setIdentifier("run0");
    }
    std::vector<String> meta_value_names;

    if (feature_map.empty() && file_list.size() > 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature map " + file_list[0] + " is empty.");
    }
    else if (feature_map.empty())
    {
      std::cout << "Warning: Feature map " + file_list[0] + " is empty." << std::endl;
      return EXECUTION_OK;
    }

    write_out_header(os, feature_map, /* main_var_name, */ meta_value_names, short_format);
    String filename;
    filename = file_list[0];
    if (getFlag_("test"))
    {
      filename = "testfile.file";
    }
    // write out the one we just loaded
    if (best_scoring.empty())
    {
      write_out_body(os, feature_map, transition_exp, meta_value_names, 0, short_format, filename);
    }
    else
    {
      write_out_body_best_score(os, feature_map, transition_exp, meta_value_names, 0, short_format, best_scoring, filename);
    }

    // start with the second in the list (we just wrote out the first one)
    for (Size i = 1; i < file_list.size(); ++i)
    {
      FileHandler().loadFeatures(file_list[i], feature_map, {FileTypes::FEATUREXML}, log_type_);
      if (feature_map.getIdentifier().empty())
      {
        feature_map.setIdentifier("run" + (String)i);
      }

      if (feature_map.empty())
      {
        continue;
      }

      filename = file_list[i];
      if (getFlag_("test"))
      {
        filename = "testfile.file";
      }

      if (best_scoring.empty())
      {
        write_out_body(os, feature_map, transition_exp, meta_value_names, boost::numeric_cast<int>(i), short_format, filename);
      }
      else
      {
        write_out_body_best_score(os, feature_map, transition_exp, meta_value_names, boost::numeric_cast<int>(i), short_format, best_scoring, filename);
      }
    }

    os.close();
    return EXECUTION_OK;

  }

};

int main(int argc, const char **argv)
{
  TOPPOpenSwathFeatureXMLToTSV tool;
  int code = tool.main(argc, argv);
  return code;

}

/// @endcond
