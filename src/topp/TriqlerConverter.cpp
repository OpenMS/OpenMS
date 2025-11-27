// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/regex.hpp>
#include <OpenMS/FORMAT/TriqlerFile.h>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_TriqlerConverter TriqlerConverter

@brief Converter to input for Triqler

This util consumes an ID-mapped consensusXML file and OpenMS experimental design in TSV format to create a CSV file which can subsequently be used as input for the python tool Triqler [1].

[1] The, M. & Käll, L. (2019). Integrated identification and quantification error probabilities for shotgun proteomics. Molecular & Cellular Proteomics, 18 (3), 561-570.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_TriqlerConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_TriqlerConverter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPTriqlerConverter final :
  public TOPPBase
{
public:
  TOPPTriqlerConverter() :
          TOPPBase("TriqlerConverter", "Converter to input for Triqler")
  {
  }

protected:
    // this function will be used to register the tool parameters
    // it gets automatically called on tool execution
    void registerOptionsAndFlags_() override
    {
      // Input consensusXML
      registerInputFile_(param_in, "<in>", "", "Input consensusXML with peptide intensities",
                         true, false);
      setValidFormats_(param_in, ListUtils::create<String>("consensusXML"), true);

      registerInputFile_(param_in_design, "<in_design>", "", "Experimental Design file", true,
                         false);
      setValidFormats_(param_in_design, ListUtils::create<String>("tsv"), true);

      registerStringOption_(param_Triqler_condition, "<Triqler_condition>", "Triqler_Condition",
                            "Which column in the condition table should be used for Triqler 'Condition'", false, false);

      // advanced option to overwrite MS file annotations in consensusXML
      registerInputFileList_(param_reannotate_filenames, "<file(s)>", StringList(),
                             "Overwrite MS file names in consensusXML", false, true);
      setValidFormats_(param_reannotate_filenames, ListUtils::create<String>("mzML"), true);                             

      // Output CSV file
      registerOutputFile_(param_out, "<out>", "", "Input CSV file for Triqler.", true, false);
      setValidFormats_(param_out, ListUtils::create<String>("csv"));
    }

    // the main_ function is called after all parameters are read
    ExitCodes main_(int, const char **) final 
    {
      try
      {
        // Input file, must be consensusXML
        const String arg_in(getStringOption_(param_in));
        const FileTypes::Type in_type(FileHandler::getType(arg_in));

        fatalErrorIf_(
                in_type != FileTypes::CONSENSUSXML,
                "Input type is not consensusXML!",
                ILLEGAL_PARAMETERS);
        // Tool arguments
        const String arg_method = getStringOption_(param_method);
        const String arg_out = getStringOption_(param_out);

        // Experimental Design file
        const String arg_in_design = getStringOption_(param_in_design);
        const ExperimentalDesign design = ExperimentalDesignFile::load(arg_in_design, false);
        ExperimentalDesign::SampleSection sampleSection = design.getSampleSection();

        ConsensusMap consensus_map;
        FileHandler().loadConsensusFeatures(arg_in, consensus_map, {FileTypes::CONSENSUSXML});

        StringList reannotate_filenames = getStringList_(param_reannotate_filenames);
        String condition = getStringOption_(param_Triqler_condition);
        String retention_time_summarization_method = getStringOption_(param_retention_time_summarization_method);

        TriqlerFile TriqlerFile;

        TriqlerFile.storeLFQ(arg_out, consensus_map, 
                              design,
                              reannotate_filenames,
                              condition);

        return EXECUTION_OK;
      }
      catch (const ExitCodes &exit_code)
      {
        return exit_code;
      }

    }

    static const String param_in;
    static const String param_in_design;
    static const String param_method;
    static const String param_Triqler_condition;
    static const String param_out;    
    static const String param_retention_time_summarization_method;
    static const String param_reannotate_filenames;

private:
    static void fatalErrorIf_(const bool error_condition, const String &message, const int exit_code)
    {
      if (error_condition)
      {
        OPENMS_LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
        throw exit_code;
      }
    }
};

const String TOPPTriqlerConverter::param_in = "in";
const String TOPPTriqlerConverter::param_in_design = "in_design";
const String TOPPTriqlerConverter::param_method = "method";
const String TOPPTriqlerConverter::param_Triqler_condition = "Triqler_condition";
const String TOPPTriqlerConverter::param_out = "out";
const String TOPPTriqlerConverter::param_retention_time_summarization_method = "retention_time_summarization_method";
const String TOPPTriqlerConverter::param_reannotate_filenames = "reannotate_filenames";

// the actual main function needed to create an executable
int main(int argc, const char **argv)
{
  TOPPTriqlerConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
