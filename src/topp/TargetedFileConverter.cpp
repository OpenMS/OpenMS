// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>

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
        <li> @ref OpenMS::TransitionParquetFile "OpenSWATH Parquet library (.oswpq)" </li>
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
    formats.push_back("oswpq");
    setValidFormats_("in", formats);
    setValidStrings_("in_type", formats);

    formats = { "tsv", "pqp", "TraML" };
    formats.push_back("oswpq");
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", formats);
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: not all conversion paths work or make sense.", false);
    setValidStrings_("out_type", formats);

    registerSubsection_("algorithm", "Algorithm parameters section");
    registerFlag_("legacy_traml_id", "PQP to TraML: Should legacy TraML IDs be used?", true);

  }

  Param getSubsectionDefaults_(const std::string&) const override
  {
    return TransitionTSVFile().getDefaults();
  }

  ExitCodes main_(int, const char**) override
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

    bool legacy_traml_id = getFlag_("legacy_traml_id");

    //---------------------------------------------------------------------------
    // Start Conversion
    //---------------------------------------------------------------------------

    // Use memory-efficient Light path for TSV/PQP → TSV/PQP conversions
    bool use_light_path = (in_type == FileTypes::TSV || in_type == FileTypes::MRM || in_type == FileTypes::PQP
      || in_type == FileTypes::OSWPQ
      )
                       && (out_type == FileTypes::TSV || out_type == FileTypes::PQP
                       || out_type == FileTypes::OSWPQ
                       );

    if (use_light_path)
    {
      // Memory-efficient Light path for TSV/PQP workflows
      OpenSwath::LightTargetedExperiment light_exp;
      OpenSwathLibraryIDNormalizer::SourceIDMapping source_ids;
      const bool persistent_output = out_type == FileTypes::PQP || out_type == FileTypes::OSWPQ;

      if (in_type == FileTypes::TSV || in_type == FileTypes::MRM)
      {
        Param reader_parameters = getParam_().copy("algorithm:", true);
        TransitionTSVFile tsv_reader;
        tsv_reader.setLogType(log_type_);
        tsv_reader.setParameters(reader_parameters);
        tsv_reader.convertTSVToTargetedExperiment(in.c_str(), in_type, light_exp);
        if (persistent_output)
        {
          source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(light_exp);
        }
      }
      else if (in_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_reader;
        Param reader_parameters = getParam_().copy("algorithm:", true);
        pqp_reader.setLogType(log_type_);
        pqp_reader.setParameters(reader_parameters);
        // TSV conversion may explicitly request source/TRAML identifiers. Persistent
        // outputs, however, stay in the existing canonical PQP ID domain.
        const bool restore_source_ids = !persistent_output;
        pqp_reader.convertPQPToTargetedExperiment(in.c_str(), light_exp, restore_source_ids);
        if (persistent_output)
        {
          source_ids.precursor_canonical_to_source =
            pqp_reader.getPQPCurrentIDToTraMLIDMap(in.c_str(), "PRECURSOR");
          source_ids.precursor_source_to_canonical.reserve(source_ids.precursor_canonical_to_source.size());
          for (const auto& [canonical_id, source_id] : source_ids.precursor_canonical_to_source)
          {
            if (source_id.empty())
            {
              continue;
            }
            const auto [it, inserted] = source_ids.precursor_source_to_canonical.emplace(source_id, canonical_id);
            if (!inserted && it->second != canonical_id)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "PQP precursor TRAML_ID maps to multiple PRECURSOR.ID values",
                                            source_id);
            }
          }
          source_ids.transition_canonical_to_source =
            pqp_reader.getPQPCurrentIDToTraMLIDMap(in.c_str(), "TRANSITION");
          for (auto it = source_ids.transition_canonical_to_source.begin();
               it != source_ids.transition_canonical_to_source.end();)
          {
            if (it->second.empty())
            {
              it = source_ids.transition_canonical_to_source.erase(it);
            }
            else
            {
              ++it;
            }
          }
        }
      }
      else if (in_type == FileTypes::OSWPQ)
      {
        TransitionParquetFile parquet_reader;
        parquet_reader.convertParquetToTargetedExperiment(in, light_exp, persistent_output ? &source_ids : nullptr);
      }

      if (out_type == FileTypes::TSV)
      {
        TransitionTSVFile tsv_writer;
        tsv_writer.setLogType(log_type_);
        tsv_writer.convertLightTargetedExperimentToTSV(out.c_str(), light_exp);
      }
      else if (out_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_writer;
        pqp_writer.setLogType(log_type_);
        pqp_writer.convertLightTargetedExperimentToPQP(out.c_str(), light_exp, &source_ids);
      }
      else if (out_type == FileTypes::OSWPQ)
      {
        TransitionParquetFile parquet_writer;
        parquet_writer.convertLightTargetedExperimentToParquet(out, light_exp, &source_ids);
      }
    }
    else
    {
      // Heavy path for TraML conversions (maintains full metadata)
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
      else if (in_type == FileTypes::OSWPQ)
      {
        writeLogError_("Error: Parquet input is only supported for light-weight conversions.");
        return PARSE_ERROR;
      }
      else if (in_type == FileTypes::TRAML)
      {
        FileHandler().loadTransitions(in, targeted_exp, {FileTypes::TRAML});
      }

      if (out_type == FileTypes::TSV)
      {
        TransitionTSVFile tsv_writer;
        tsv_writer.setLogType(log_type_);
        tsv_writer.convertTargetedExperimentToTSV(out.c_str(), targeted_exp);
      }
      else if (out_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_writer;
        pqp_writer.setLogType(log_type_);
        pqp_writer.convertTargetedExperimentToPQP(out.c_str(), targeted_exp);
      }
      else if (out_type == FileTypes::TRAML)
      {
        FileHandler().storeTransitions(out, targeted_exp, {FileTypes::TRAML});
      }
      else if (out_type == FileTypes::OSWPQ)
      {
        OpenSwath::LightTargetedExperiment light_exp;
        OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);
        auto source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(light_exp);
        TransitionParquetFile parquet_writer;
        parquet_writer.convertLightTargetedExperimentToParquet(out, light_exp, &source_ids);
      }
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
