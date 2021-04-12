// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_TargetedFileConverter TargetedFileConverter

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
  @verbinclude UTILS_TargetedFileConverter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_TargetedFileConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPTargetedFileConverter :
  public TOPPBase
{
public:

  TOPPTargetedFileConverter() :
    TOPPBase("TargetedFileConverter", "Converts different transition files for targeted proteomics / metabolomics analysis.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.\n "
                                           "See http://www.openms.de/current_doxygen/html/UTILS_TargetedFileConverter.html for format of OpenSWATH transition TSV file or SpectraST MRM file.");
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
      writeLog_("Error: Could not determine input file type!");
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
      writeLog_("Error: Could not determine output file type!");
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
      TraMLFile traml;
      traml.load(in, targeted_exp);
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
      TraMLFile traml;
      traml.store(out, targeted_exp);
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
