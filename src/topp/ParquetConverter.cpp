// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <filesystem>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ParquetConverter ParquetConverter

@brief Converts between XML and Parquet formats for feature and consensus maps.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; ParquetConverter &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Any tool producing featureXML or consensusXML </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Any tool consuming featureXML or consensusXML </td>
        </tr>
    </table>
</CENTER>

ParquetConverter converts between XML and Parquet directory formats for
FeatureMap and ConsensusMap data. The conversion direction is detected
automatically from the input:

- If the input is a featureXML or consensusXML file, the output is a
  Parquet directory containing multiple .parquet files (features/consensus
  features, PSMs, proteins, protein groups, search parameters).
- If the input is a Parquet directory (previously exported), the output
  is a featureXML or consensusXML file.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ParquetConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ParquetConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPParquetConverter :
  public TOPPBase
{
public:
  TOPPParquetConverter() :
    TOPPBase("ParquetConverter",
             "Converts between XML and Parquet formats for feature and consensus maps.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<path>", "", "Input file (featureXML, consensusXML) or Parquet directory");
    registerOutputFile_("out", "<path>", "", "Output directory (for Parquet export) or file (for XML import)", true);
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");

    if (in.empty() || out.empty())
    {
      OPENMS_LOG_ERROR << "Both -in and -out must be specified.\n";
      return ILLEGAL_PARAMETERS;
    }

    // Detect direction based on input
    namespace fs = std::filesystem;

    if (fs::is_directory(fs::path(std::string(in))))
    {
      // Verify directory contains parquet files before treating as import source
      bool has_parquet = false;
      for (const auto& entry : fs::directory_iterator(fs::path(std::string(in))))
      {
        if (entry.path().extension() == ".parquet")
        {
          has_parquet = true;
          break;
        }
      }
      if (has_parquet)
      {
        return importFromParquet_(in, out);
      }
      // Fall through to error if directory has no parquet files
    }

    if (in.hasSuffix(".featureXML"))
    {
      return exportFeatureMap_(in, out);
    }

    if (in.hasSuffix(".consensusXML"))
    {
      return exportConsensusMap_(in, out);
    }

    OPENMS_LOG_ERROR << "Cannot determine conversion direction. "
                     << "Input must be a .featureXML, .consensusXML file, or a Parquet directory.\n";
    return ILLEGAL_PARAMETERS;
  }

private:
  ExitCodes exportFeatureMap_(const String& in, const String& out)
  {
    OPENMS_LOG_INFO << "Loading featureXML..." << endl;
    FeatureMap fm;
    FeatureXMLFile().load(in, fm);

    OPENMS_LOG_INFO << "Exporting to Parquet directory: " << out << endl;
    if (!FeatureMapArrowIO::exportToParquet(fm, out))
    {
      OPENMS_LOG_ERROR << "Failed to export FeatureMap to Parquet.\n";
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    OPENMS_LOG_INFO << "Export completed successfully." << endl;
    return EXECUTION_OK;
  }

  ExitCodes exportConsensusMap_(const String& in, const String& out)
  {
    OPENMS_LOG_INFO << "Loading consensusXML..." << endl;
    ConsensusMap cm;
    ConsensusXMLFile().load(in, cm);

    OPENMS_LOG_INFO << "Exporting to Parquet directory: " << out << endl;
    if (!ConsensusMapArrowIO::exportToParquet(cm, out))
    {
      OPENMS_LOG_ERROR << "Failed to export ConsensusMap to Parquet.\n";
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    OPENMS_LOG_INFO << "Export completed successfully." << endl;
    return EXECUTION_OK;
  }

  ExitCodes importFromParquet_(const String& in, const String& out)
  {
    if (out.hasSuffix(".featureXML"))
    {
      OPENMS_LOG_INFO << "Importing FeatureMap from Parquet directory: " << in << endl;
      FeatureMap fm;
      if (!FeatureMapArrowIO::importFromParquet(in, fm))
      {
        OPENMS_LOG_ERROR << "Failed to import FeatureMap from Parquet.\n";
        return INPUT_FILE_NOT_READABLE;
      }

      OPENMS_LOG_INFO << "Writing featureXML..." << endl;
      FeatureXMLFile().store(out, fm);
      OPENMS_LOG_INFO << "Import completed successfully." << endl;
      return EXECUTION_OK;
    }

    if (out.hasSuffix(".consensusXML"))
    {
      OPENMS_LOG_INFO << "Importing ConsensusMap from Parquet directory: " << in << endl;
      ConsensusMap cm;
      if (!ConsensusMapArrowIO::importFromParquet(in, cm))
      {
        OPENMS_LOG_ERROR << "Failed to import ConsensusMap from Parquet.\n";
        return INPUT_FILE_NOT_READABLE;
      }

      OPENMS_LOG_INFO << "Writing consensusXML..." << endl;
      ConsensusXMLFile().store(out, cm);
      OPENMS_LOG_INFO << "Import completed successfully." << endl;
      return EXECUTION_OK;
    }

    OPENMS_LOG_ERROR << "Output file must have .featureXML or .consensusXML extension "
                     << "when importing from a Parquet directory.\n";
    return ILLEGAL_PARAMETERS;
  }
};

int main(int argc, const char** argv)
{
  TOPPParquetConverter tool;
  return tool.main(argc, argv);
}

/// @endcond
