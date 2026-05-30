// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QPXConverter QPXConverter

@brief Converts IdXML files to parquet format following QPX PSM specification.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; QPXConverter &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Any identification tool producing idXML </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> QPX analysis tools </td>
        </tr>
    </table>
</CENTER>

QPXConverter reads peptide and protein identifications from idXML files
and converts them to parquet format following the QPX PSM (Peptide
Spectrum Match) specification.

The output parquet file contains PSM data with structured columns following
the QPX PSM specification including ProForma peptidoform notation,
structured modifications with UniMod accessions, structured additional scores,
protein accessions as list, and metavalue columns.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QPXConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QPXConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPQPXConverter :
  public TOPPBase
{
public:
  TOPPQPXConverter() :
    TOPPBase("QPXConverter", "Converts IdXML files to parquet format following QPX PSM specification.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input idXML file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));

    registerOutputFile_("out", "<file>", "", "Output parquet file", true);
    setValidFormats_("out", ListUtils::create<String>("parquet"));

    registerFlag_("export_all_psms", "Export all PSMs per spectrum (not just the best hit)");
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const bool export_all = getFlag_("export_all_psms");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    vector<ProteinIdentification> protein_identifications;
    PeptideIdentificationList peptide_identifications;

    OPENMS_LOG_INFO << "Loading idXML file..." << endl;
    IdXMLFile idxml_file;
    idxml_file.load(in, protein_identifications, peptide_identifications);

    if (peptide_identifications.empty())
    {
      OPENMS_LOG_WARN << "No peptide identifications found in input file." << endl;
      return INPUT_FILE_EMPTY;
    }

    OPENMS_LOG_INFO << "Found " << peptide_identifications.size()
                    << " peptide identifications in "
                    << protein_identifications.size()
                    << " protein identification runs." << endl;

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Converting to parquet format..." << endl;
    if (!QPXFile::exportToParquet(protein_identifications, peptide_identifications, out, export_all))
    {
      OPENMS_LOG_ERROR << "Failed to write Parquet file: " << out << endl;
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    OPENMS_LOG_INFO << "Conversion completed successfully." << endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPQPXConverter tool;
  return tool.main(argc, argv);
}

/// @endcond
