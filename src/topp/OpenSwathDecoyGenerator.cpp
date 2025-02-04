// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathDecoyGenerator OpenSwathDecoyGenerator

@brief Generates decoys according to different models for a specific TraML

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; OpenSwathDecoyGenerator &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_OpenSwathAssayGenerator </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathWorkflow </td>
        </tr>
    </table>
</CENTER>

This module generates "decoy" transitions from a set of real or "target"
transitions. The idea is to use the decoy transitions in a statistical scoring
process to estimate the false hits in an SRM / SWATH experiment.  The tool
operates on @ref OpenMS::TraMLFile "TraML" files, which can come from @ref
TOPP_TargetedFileConverter or any other tool.

There are multiple methods to create the decoy transitions, the simplest ones
are reverse and pseudo-reverse which reverse the sequence either completely or
leaving the last (tryptic) AA untouched respectively.

Another decoy generation method is "shuffle" which uses an algorithm similar
to the one described in Lam, Henry, et al. (2010). "Artificial decoy spectral
libraries for false discovery rate estimation in spectral library searching in
proteomics".  Journal of Proteome Research 9, 605-610. It shuffles the amino
acid sequence (excluding N-/C-terminal and K/R/P) and shuffles the fragment 
ion intensities accordingly. If the new sequence does not reach a threshold of
identity within a set number of trials, a random amino acid (not N-/C-terminal
or modified) is "mutated" to a random other amino acid.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathDecoyGenerator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathDecoyGenerator.html
*/

// TODO: could theoretical also produce an annotation in the TraML of what it thinks the ion is?

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathDecoyGenerator
: public TOPPBase
{
public:

  TOPPOpenSwathDecoyGenerator() :
    TOPPBase("OpenSwathDecoyGenerator", "Generates decoys according to different models for a specific TraML", true)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content\n", false);
    String formats("tsv,mrm,pqp,TraML");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));

    formats = "tsv,pqp,TraML";
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\n", false);
    setValidStrings_("out_type", ListUtils::create<String>(formats));

    registerStringOption_("method", "<type>", "shuffle", "Decoy generation method", false);
    setValidStrings_("method", ListUtils::create<String>(String("shuffle,pseudo-reverse,reverse,shift")));

    registerStringOption_("decoy_tag", "<type>", "DECOY_", "decoy tag", false);

    registerDoubleOption_("min_decoy_fraction", "<double>", 0.8, "Minimum fraction of decoy / target peptides and proteins", false, true);
    registerDoubleOption_("aim_decoy_fraction", "<double>", 1.0, "Number of decoys the algorithm should generate (if unequal to 1, the algorithm will randomly select N peptides for decoy generation)", false, true);

    registerIntOption_("shuffle_max_attempts", "<int>", 30, "shuffle: maximum attempts to lower the amino acid sequence identity between target and decoy for the shuffle algorithm", false, true);
    registerDoubleOption_("shuffle_sequence_identity_threshold", "<double>", 0.5, "shuffle: target-decoy amino acid sequence identity threshold for the shuffle algorithm", false, true);

    registerDoubleOption_("shift_precursor_mz_shift", "<double>", 0.0, "shift: precursor ion MZ shift in Thomson for shift decoy method", false, true);
    registerDoubleOption_("shift_product_mz_shift", "<double>", 20, "shift: fragment ion MZ shift in Thomson for shift decoy method", false, true);

    registerDoubleOption_("product_mz_threshold", "<double>", 0.025, "MZ threshold in Thomson for fragment ion annotation", false, true);
    registerStringOption_("allowed_fragment_types", "<type>", "b,y", "allowed fragment types", false, true);
    registerStringOption_("allowed_fragment_charges", "<type>", "1,2,3,4", "allowed fragment charge states", false, true);
    registerFlag_("enable_detection_specific_losses", "set this flag if specific neutral losses for detection fragment ions should be allowed", true);
    registerFlag_("enable_detection_unspecific_losses", "set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for detection fragment ions should be allowed", true);
    registerStringOption_("switchKR", "<true/false>", "true", "Whether to switch terminal K and R (to achieve different precursor mass)", false);
    setValidStrings_("switchKR", ListUtils::create<String>(String("true,false")));

    registerFlag_("separate", "set this flag if decoys should not be appended to targets.", true);
  }

  ExitCodes main_(int, const char **) override
  {
    FileHandler fh;

    //input file type
    String in = getStringOption_("in");
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    //output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    String method = getStringOption_("method");
    String decoy_tag = getStringOption_("decoy_tag");

    double min_decoy_fraction = getDoubleOption_("min_decoy_fraction");
    double aim_decoy_fraction = getDoubleOption_("aim_decoy_fraction");

    Int max_attempts = getIntOption_("shuffle_max_attempts");
    double identity_threshold = getDoubleOption_("shuffle_sequence_identity_threshold");

    double precursor_mz_shift = getDoubleOption_("shift_precursor_mz_shift");
    double product_mz_shift = getDoubleOption_("shift_product_mz_shift");

    double product_mz_threshold = getDoubleOption_("product_mz_threshold");
    String allowed_fragment_types_string = getStringOption_("allowed_fragment_types");
    String allowed_fragment_charges_string = getStringOption_("allowed_fragment_charges");
    bool enable_detection_specific_losses = getFlag_("enable_detection_specific_losses");
    bool enable_detection_unspecific_losses = getFlag_("enable_detection_unspecific_losses");
    bool switchKR = getStringOption_("switchKR") == "true";

    bool separate = getFlag_("separate");

    std::vector<String> allowed_fragment_types;
    allowed_fragment_types_string.split(",", allowed_fragment_types);

    std::vector<String> allowed_fragment_charges_string_vector;
    std::vector<size_t> allowed_fragment_charges;
    allowed_fragment_charges_string.split(",", allowed_fragment_charges_string_vector);
    for (size_t i = 0; i < allowed_fragment_charges_string_vector.size(); i++)
    {
      size_t charge = std::atoi(allowed_fragment_charges_string_vector.at(i).c_str());
      allowed_fragment_charges.push_back(charge);
    }

    TargetedExperiment targeted_merged;
    {
      TargetedExperiment targeted_exp;
      TargetedExperiment targeted_decoy;

      // Load data
      OPENMS_LOG_INFO << "Loading targets from file: " << in << std::endl;
      if (in_type == FileTypes::TSV || in_type == FileTypes::MRM)
      {
        const char* tr_file = in.c_str();
        Param reader_parameters = getParam_().copy("algorithm:", true);
        TransitionTSVFile tsv_reader = TransitionTSVFile();
        tsv_reader.setLogType(log_type_);
        tsv_reader.setParameters(reader_parameters);
        tsv_reader.convertTSVToTargetedExperiment(tr_file, in_type, targeted_exp);
        tsv_reader.validateTargetedExperiment(targeted_exp);
      }
      else if (in_type == FileTypes::PQP)
      {
        const char* tr_file = in.c_str();
        TransitionPQPFile pqp_reader = TransitionPQPFile();
        Param reader_parameters = getParam_().copy("algorithm:", true);
        pqp_reader.setLogType(log_type_);
        pqp_reader.setParameters(reader_parameters);
        pqp_reader.convertPQPToTargetedExperiment(tr_file, targeted_exp);
        pqp_reader.validateTargetedExperiment(targeted_exp);
      }
      else if (in_type == FileTypes::TRAML)
      {
        FileHandler().loadTransitions(in, targeted_exp, {FileTypes::TRAML});
      }

      MRMDecoy decoys = MRMDecoy();
      decoys.setLogType(ProgressLogger::CMD);

      OPENMS_LOG_INFO << "Generate decoys" << std::endl;
      decoys.generateDecoys(targeted_exp, targeted_decoy, method,
                            aim_decoy_fraction, switchKR, decoy_tag, max_attempts,
                            identity_threshold, precursor_mz_shift,
                            product_mz_shift, product_mz_threshold,
                            allowed_fragment_types, allowed_fragment_charges,
                            enable_detection_specific_losses,
                            enable_detection_unspecific_losses);


      // Check if we have enough peptides left
      OPENMS_LOG_INFO << "Number of target peptides: " << targeted_exp.getPeptides().size() << std::endl;
      OPENMS_LOG_INFO << "Number of decoy peptides: " << targeted_decoy.getPeptides().size() << std::endl;
      OPENMS_LOG_INFO << "Number of target proteins: " << targeted_exp.getProteins().size() << std::endl;
      OPENMS_LOG_INFO << "Number of decoy proteins: " << targeted_decoy.getProteins().size() << std::endl;

      if ((float)targeted_decoy.getPeptides().size() / (float)targeted_exp.getPeptides().size() < min_decoy_fraction || (float)targeted_decoy.getProteins().size() / (float)targeted_exp.getProteins().size() < min_decoy_fraction)
      {
         throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The number of decoys for peptides or proteins is below the threshold of " + String(min_decoy_fraction * 100) + "% of the number of targets.");
      }

      if (separate)
      {
        OPENMS_LOG_INFO << "Writing only decoys to file: " << out << std::endl;
        targeted_merged = std::move(targeted_decoy);
      }
      else
      {
        OPENMS_LOG_INFO << "Writing targets and decoys to file: " << out << std::endl;
        targeted_merged = std::move(targeted_exp);
        targeted_merged += std::move(targeted_decoy);
      }

    }

    if (out_type == FileTypes::TSV)
    {
      const char* tr_file = out.c_str();
      TransitionTSVFile tsv_reader = TransitionTSVFile();
      tsv_reader.setLogType(log_type_);
      tsv_reader.convertTargetedExperimentToTSV(tr_file, targeted_merged);
    }
    if (out_type == FileTypes::PQP)
    {
      const char * tr_file = out.c_str();
      TransitionPQPFile pqp_reader = TransitionPQPFile();
      pqp_reader.setLogType(log_type_);
      pqp_reader.convertTargetedExperimentToPQP(tr_file, targeted_merged);
    }
    else if (out_type == FileTypes::TRAML)
    {
      FileHandler().storeTransitions(out, targeted_merged, {FileTypes::TRAML});
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char **argv)
{
  TOPPOpenSwathDecoyGenerator gen;
  return gen.main(argc, argv);
}

/// @endcond
