// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
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
    StringList formats = {"tsv", "mrm", "pqp", "TraML"};
    formats.push_back("oswpq");
    setValidFormats_("in", formats);
    setValidStrings_("in_type", formats);

    formats = {"tsv", "pqp", "TraML"};
    formats.push_back("oswpq");
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", formats);
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\n", false);
    setValidStrings_("out_type", formats);

    registerStringOption_("method", "<type>", "shuffle", "Decoy generation method", false);
    setValidStrings_("method", ListUtils::create<std::string>(std::string("shuffle,pseudo-reverse,reverse,shift")));

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
    setValidStrings_("switchKR", ListUtils::create<std::string>(std::string("true,false")));

    registerFlag_("separate", "set this flag if decoys should not be appended to targets.", true);
  }

  ExitCodes main_(int, const char **) override
  {
    FileHandler fh;

    //input file type
    std::string in = getStringOption_("in");
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(std::string("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    //output file names and types
    std::string out = getStringOption_("out");
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

    std::string method = getStringOption_("method");
    std::string decoy_tag = getStringOption_("decoy_tag");

    double min_decoy_fraction = getDoubleOption_("min_decoy_fraction");
    double aim_decoy_fraction = getDoubleOption_("aim_decoy_fraction");

    Int max_attempts = getIntOption_("shuffle_max_attempts");
    double identity_threshold = getDoubleOption_("shuffle_sequence_identity_threshold");

    double precursor_mz_shift = getDoubleOption_("shift_precursor_mz_shift");
    double product_mz_shift = getDoubleOption_("shift_product_mz_shift");

    double product_mz_threshold = getDoubleOption_("product_mz_threshold");
    std::string allowed_fragment_types_string = getStringOption_("allowed_fragment_types");
    std::string allowed_fragment_charges_string = getStringOption_("allowed_fragment_charges");
    bool enable_detection_specific_losses = getFlag_("enable_detection_specific_losses");
    bool enable_detection_unspecific_losses = getFlag_("enable_detection_unspecific_losses");
    bool switchKR = getStringOption_("switchKR") == "true";

    bool separate = getFlag_("separate");

    std::vector<std::string> allowed_fragment_types;
    StringUtils::split(allowed_fragment_types_string, ",", allowed_fragment_types);

    std::vector<std::string> allowed_fragment_charges_string_vector;
    std::vector<size_t> allowed_fragment_charges;
    StringUtils::split(allowed_fragment_charges_string, ",", allowed_fragment_charges_string_vector);
    for (size_t i = 0; i < allowed_fragment_charges_string_vector.size(); i++)
    {
      size_t charge = std::atoi(allowed_fragment_charges_string_vector.at(i).c_str());
      allowed_fragment_charges.push_back(charge);
    }

    OpenSwathLibraryPreparation::DecoyGeneratorParameters decoy_parameters;
    decoy_parameters.method = method;
    decoy_parameters.decoy_tag = decoy_tag;
    decoy_parameters.min_decoy_fraction = min_decoy_fraction;
    decoy_parameters.aim_decoy_fraction = aim_decoy_fraction;
    decoy_parameters.shuffle_max_attempts = max_attempts;
    decoy_parameters.shuffle_sequence_identity_threshold = identity_threshold;
    decoy_parameters.shift_precursor_mz_shift = precursor_mz_shift;
    decoy_parameters.shift_product_mz_shift = product_mz_shift;
    decoy_parameters.product_mz_threshold = product_mz_threshold;
    decoy_parameters.allowed_fragment_types = allowed_fragment_types;
    decoy_parameters.allowed_fragment_charges = allowed_fragment_charges;
    decoy_parameters.enable_detection_specific_losses = enable_detection_specific_losses;
    decoy_parameters.enable_detection_unspecific_losses = enable_detection_unspecific_losses;
    decoy_parameters.switch_kr = switchKR;
    decoy_parameters.separate = separate;

    OpenSwathLibraryPreparation preparation;
    preparation.setLogType(log_type_);
    preparation.generateDecoys(in, in_type, out, out_type, decoy_parameters, getParam_().copy("algorithm:", true));

    return EXECUTION_OK;
  }

};

int main(int argc, const char **argv)
{
  TOPPOpenSwathDecoyGenerator gen;
  return gen.main(argc, argv);
}

/// @endcond
