// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_TargetedFileConverter TargetedFileConverter

@brief Converts different spectral libraries / transition files for targeted proteomics and metabolomics analysis.

Can convert multiple formats to and from TraML (standardized transition format). The following formats are supported:

      <ul>
        <li> @ref OpenMS::TraMLFile "TraML" </li>
        <li> @ref OpenMS::TransitionTSVFile "OpenSWATH TSV transition lists" </li>
        <li> @ref OpenMS::TransitionPQPFile "OpenSWATH PQP SQLite files" </li>
        <li> SpectraST MRM transition lists </li>
        <li> Skyline transition lists </li>
        <li> Spectronaut transition lists </li>
      </ul>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_TargetedFileConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_TargetedFileConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPTargetedFileConverter :
  public TOPPBase
{
public:

  TOPPTargetedFileConverter() :
    TOPPBase("TargetedFileConverter", "Converts different transition files for targeted proteomics / metabolomics analysis.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.\n "
                                           "See http://www.openms.de/current_doxygen/html/TOPP_TargetedFileConverter.html for format of OpenSWATH transition TSV file or SpectraST MRM file.");
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    StringList formats{"tsv", "mrm" ,"pqp", "TraML"};
    setValidFormats_("in", formats);
    setValidStrings_("in_type", formats);

    formats = { "tsv", "pqp", "TraML" };
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", formats);
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: not all conversion paths work or make sense.", false);
    setValidStrings_("out_type", formats);

    registerSubsection_("algorithm", "Algorithm parameters section");
    registerFlag_("legacy_traml_id", "PQP to TraML: Should legacy TraML IDs be used?", true);

  }

  Param getSubsectionDefaults_(const String&) const override
  {
    return TransitionTSVFile().getDefaults();
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

    bool legacy_traml_id = getFlag_("legacy_traml_id");

    //--------------------------------------------------------------------------- 
    // Start Conversion
    //--------------------------------------------------------------------------- 
    TargetedExperiment targeted_exp;
    if (in_type == FileTypes::TSV || in_type == FileTypes::MRM)
    {
      Param reader_parameters = getParam_().copy("algorithm:", true);
      TransitionTSVFile tsv_reader;
      tsv_reader.setLogType(log_type_);
      tsv_reader.setParameters(reader_parameters);
      tsv_reader.convertTSVToTargetedExperiment(in.c_str(), in_type, targeted_exp);
      tsv_reader.validateTargetedExperiment(targeted_exp);
    }
    else if (in_type == FileTypes::PQP)
    {
      TransitionPQPFile pqp_reader;
      Param reader_parameters = getParam_().copy("algorithm:", true);
      pqp_reader.setLogType(log_type_);
      pqp_reader.setParameters(reader_parameters);
      pqp_reader.convertPQPToTargetedExperiment(in.c_str(), targeted_exp, legacy_traml_id);
      pqp_reader.validateTargetedExperiment(targeted_exp);
    }
    else if (in_type == FileTypes::TRAML)
    {
      FileHandler().loadTransitions(in, targeted_exp, {FileTypes::TRAML});
    }

    if (out_type == FileTypes::TSV)
    {
      TransitionTSVFile tsv_reader;
      tsv_reader.setLogType(log_type_);
      tsv_reader.convertTargetedExperimentToTSV(out.c_str(), targeted_exp);
    }
    if (out_type == FileTypes::PQP)
    {
      TransitionPQPFile pqp_reader;
      pqp_reader.setLogType(log_type_);
      pqp_reader.convertTargetedExperimentToPQP(out.c_str(), targeted_exp);
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

  TOPPTargetedFileConverter tool;
  return tool.main(argc, argv);
}

/// @endcond
