// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <iostream>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathAssayGenerator OpenSwathAssayGenerator

@brief Generates filtered and optimized assays using TraML files.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; OpenSwathAssayGenerator &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1>  </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathDecoyGenerator  </td>
        </tr>
    </table>
</CENTER>

This module generates assays for targeted proteomics using a set of rules
that was found to improve the sensitivity and selectivity for detection
of typical peptides (Schubert et al., 2015). The tool operates on @ref
OpenMS::TraMLFile "TraML" files, which can come from @ref
TOPP_TargetedFileConverter or any other tool. In a first step, the tool will
annotate all transitions according to the predefined criteria. In a second
step, the transitions will be filtered to improve sensitivity for detection
of peptides.

Optionally, theoretical identification transitions can be generated when the
TraML will be used for IPF scoring in OpenSWATH, see @ref OpenMS::MRMAssay
"MRMAssay" for more information on the algorithm. This is recommended if
post-translational modifications are scored with OpenSWATH.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathAssayGenerator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathAssayGenerator.html


*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOpenSwathAssayGenerator :
  public TOPPBase
{
public:

  TOPPOpenSwathAssayGenerator() :
    TOPPBase("OpenSwathAssayGenerator", "Generates assays according to different models for a specific TraML", true)
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

    registerIntOption_("min_transitions", "<int>", 6, "minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "maximal number of transitions", false);
    registerStringOption_("allowed_fragment_types", "<type>", "b,y", "allowed fragment types", false);
    registerStringOption_("allowed_fragment_charges", "<type>", "1,2,3,4", "allowed fragment charge states", false);
    registerFlag_("enable_detection_specific_losses", "set this flag if specific neutral losses for detection fragment ions should be allowed");
    registerFlag_("enable_detection_unspecific_losses", "set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for detection fragment ions should be allowed");

    registerDoubleOption_("precursor_mz_threshold", "<double>", 0.025, "MZ threshold in Thomson for precursor ion selection", false);
    registerDoubleOption_("precursor_lower_mz_limit", "<double>", 400, "lower MZ limit for precursor ions", false);
    registerDoubleOption_("precursor_upper_mz_limit", "<double>", 1200, "upper MZ limit for precursor ions", false);
    registerDoubleOption_("product_mz_threshold", "<double>", 0.025, "MZ threshold in Thomson for fragment ion annotation", false);
    registerDoubleOption_("product_lower_mz_limit", "<double>", 350, "lower MZ limit for fragment ions", false);
    registerDoubleOption_("product_upper_mz_limit", "<double>", 2000, "upper MZ limit for fragment ions", false);

    registerInputFile_("swath_windows_file", "<file>", "", "Tab separated file containing the SWATH windows for exclusion of fragment ions falling into the precursor isolation window: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped.", false, false);
    setValidFormats_("swath_windows_file", ListUtils::create<String>("txt"));

    registerInputFile_("unimod_file", "<file>", "", "(Modified) Unimod XML file (http://www.unimod.org/xml/unimod.xml) describing residue modifiability", false, false);
    setValidFormats_("unimod_file", ListUtils::create<String>("xml"));

    registerFlag_("enable_ipf", "IPF: set this flag if identification transitions should be generated for IPF. Note: Requires setting 'unimod_file'.");
    registerIntOption_("max_num_alternative_localizations", "<int>", 10000, "IPF: maximum number of site-localization permutations", false, true);
    registerFlag_("disable_identification_ms2_precursors", "IPF: set this flag if MS2-level precursor ions for identification should not be allowed for extraction of the precursor signal from the fragment ion data (MS2-level).", true);
    registerFlag_("disable_identification_specific_losses", "IPF: set this flag if specific neutral losses for identification fragment ions should not be allowed", true);
    registerFlag_("enable_identification_unspecific_losses", "IPF: set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for identification fragment ions should be allowed", true);
    registerFlag_("enable_swath_specifity", "IPF: set this flag if identification transitions without precursor specificity (i.e. across whole precursor isolation window instead of precursor MZ) should be generated.", true);

  }

  ExitCodes main_(int, const char**) override
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

    Int min_transitions = getIntOption_("min_transitions");
    Int max_transitions = getIntOption_("max_transitions");
    String allowed_fragment_types_string = getStringOption_("allowed_fragment_types");
    String allowed_fragment_charges_string = getStringOption_("allowed_fragment_charges");
    bool enable_detection_specific_losses = getFlag_("enable_detection_specific_losses");
    bool enable_detection_unspecific_losses = getFlag_("enable_detection_unspecific_losses");
    bool enable_identification_specific_losses = !getFlag_("disable_identification_specific_losses");
    bool enable_identification_unspecific_losses = getFlag_("enable_identification_unspecific_losses");
    bool enable_identification_ms2_precursors = !getFlag_("disable_identification_ms2_precursors");
    bool enable_ipf = getFlag_("enable_ipf");
    bool enable_swath_specifity = getFlag_("enable_swath_specifity");
    size_t max_num_alternative_localizations = getIntOption_("max_num_alternative_localizations");
    double precursor_mz_threshold = getDoubleOption_("precursor_mz_threshold");
    double precursor_lower_mz_limit = getDoubleOption_("precursor_lower_mz_limit");
    double precursor_upper_mz_limit = getDoubleOption_("precursor_upper_mz_limit");
    double product_mz_threshold = getDoubleOption_("product_mz_threshold");
    double product_lower_mz_limit = getDoubleOption_("product_lower_mz_limit");
    double product_upper_mz_limit = getDoubleOption_("product_upper_mz_limit");
    String swath_windows_file = getStringOption_("swath_windows_file");

    String unimod_file = getStringOption_("unimod_file");
    bool is_test = getFlag_("test");

    // Set specific seed for test mode
    int uis_seed = -1;
    bool disable_decoy_transitions = false;
    if (is_test)
    {
      uis_seed = 42;
      disable_decoy_transitions = true;
    }

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

    // Require Unimod XML file when running IPF to prevent accidential mistakes
    if (enable_ipf && unimod_file.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please provide a valid Unimod XML file for IPF.");
    }

    // Load Unimod file
    if (!unimod_file.empty())
    {
      if (!ModificationsDB::isInstantiated()) // We need to ensure that ModificationsDB was not instantiated before!
      {
        ModificationsDB* ptr = ModificationsDB::initializeModificationsDB(unimod_file, String(""), String(""));
        OPENMS_LOG_INFO << "Unimod XML: " << ptr->getNumberOfModifications() << " modification types and residue specificities imported from file: " << unimod_file << std::endl;
      }
      else
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ModificationsDB has been instantiated before and can not be generated from the provided Unimod XML file.");
      }
    }

    std::vector<std::pair<double, double> > swathes;
    // Check swath window input
    if (!swath_windows_file.empty())
    {
      OPENMS_LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      OPENMS_LOG_INFO << "Read Swath maps file with " << swath_prec_lower.size() << " windows." << std::endl;
      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        swathes.push_back(std::make_pair(swath_prec_lower[i], swath_prec_upper[i]));
        OPENMS_LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    TargetedExperiment targeted_exp;

    // Load data
    OPENMS_LOG_INFO << "Loading " << in << std::endl;
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

    MRMAssay assays = MRMAssay();
    assays.setLogType(ProgressLogger::CMD);

    OPENMS_LOG_INFO << "Annotating transitions" << std::endl;
    assays.reannotateTransitions(targeted_exp, precursor_mz_threshold, product_mz_threshold, allowed_fragment_types, allowed_fragment_charges, enable_detection_specific_losses, enable_detection_unspecific_losses);

    OPENMS_LOG_INFO << "Annotating detecting transitions" << std::endl;
    assays.restrictTransitions(targeted_exp, product_lower_mz_limit, product_upper_mz_limit, swathes);
    assays.detectingTransitions(targeted_exp, min_transitions, max_transitions);

    if (enable_ipf)
    {
      std::vector<std::pair<double, double> > uis_swathes;

      if (!enable_swath_specifity)
      {
        int num_precursor_windows = static_cast<int>(Math::round((precursor_upper_mz_limit - precursor_lower_mz_limit) / precursor_mz_threshold));
        for (int i = 0; i < num_precursor_windows; i++)
        {
          uis_swathes.push_back(std::make_pair((precursor_lower_mz_limit+(i*precursor_mz_threshold)),(precursor_lower_mz_limit+((i+1)*precursor_mz_threshold))));
        }
      }
      else
      {
        uis_swathes = swathes;
      }
      
      OPENMS_LOG_INFO << "Generating identifying transitions for IPF" << std::endl;
      assays.uisTransitions(targeted_exp, allowed_fragment_types, allowed_fragment_charges, enable_identification_specific_losses, enable_identification_unspecific_losses, enable_identification_ms2_precursors, product_mz_threshold, uis_swathes, -4, max_num_alternative_localizations, uis_seed, disable_decoy_transitions);
      std::vector<std::pair<double, double> > empty_swathes;
      assays.restrictTransitions(targeted_exp, product_lower_mz_limit, product_upper_mz_limit, empty_swathes);
    }

    OPENMS_LOG_INFO << "Writing assays " << out << std::endl;
    if (out_type == FileTypes::TSV)
    {
      const char* tr_file = out.c_str();
      TransitionTSVFile tsv_reader = TransitionTSVFile();
      tsv_reader.setLogType(log_type_);
      tsv_reader.convertTargetedExperimentToTSV(tr_file, targeted_exp);
    }
    if (out_type == FileTypes::PQP)
    {
      const char * tr_file = out.c_str();
      TransitionPQPFile pqp_reader = TransitionPQPFile();
      pqp_reader.setLogType(log_type_);
      pqp_reader.convertTargetedExperimentToPQP(tr_file, targeted_exp);
    }
    else if (out_type == FileTypes::TRAML)
    {
      FileHandler().storeTransitions(out, targeted_exp, {FileTypes::TRAML});
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPOpenSwathAssayGenerator gen;
  return gen.main(argc, argv);
}

/// @endcond
