// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/QuantmsIO.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QuantmsIOConverter QuantmsIOConverter

@brief Converts IdXML files to parquet format following quantms.io PSM specification.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; QuantmsIOConverter &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Any identification tool producing idXML </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> quantms.io analysis tools </td>
        </tr>
    </table>
</CENTER>

QuantmsIOConverter reads peptide and protein identifications from idXML files 
and converts them to parquet format following the quantms.io PSM (Peptide 
Spectrum Match) specification.

The output parquet file contains PSM data with columns including:
- sequence: peptide sequence
- spectrum_reference: spectrum identifier
- charge: precursor charge
- retention_time: retention time in seconds
- mass_to_charge: precursor m/z value
- score: identification score
- rank: peptide hit rank
- protein_accessions: associated protein accessions
- modifications: peptide modifications
- is_decoy: decoy database flag
- search_engine: search engine name
- search_engine_score_name: score type name

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QuantmsIOConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QuantmsIOConverter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPQuantmsIOConverter :
  public TOPPBase
{
public:
  TOPPQuantmsIOConverter() :
    TOPPBase("QuantmsIOConverter", "Converts IdXML files to parquet format following quantms.io PSM specification.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input idXML file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));

    registerOutputFile_("out", "<file>", "", "Output parquet file", true);
    setValidFormats_("out", ListUtils::create<String>("parquet"));
  }

  ExitCodes main_(int, const char**) override
  {
#ifdef WITH_PARQUET
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");

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
    QuantmsIO quantms_io;
    quantms_io.store(out, protein_identifications, peptide_identifications);

    OPENMS_LOG_INFO << "Conversion completed successfully." << endl;
    return EXECUTION_OK;

#else
    OPENMS_LOG_ERROR << "This tool requires OpenMS to be compiled with parquet support (WITH_PARQUET=ON)." << endl;
    return EXTERNAL_PROGRAM_ERROR;
#endif
  }
};

int main(int argc, const char** argv)
{
  TOPPQuantmsIOConverter tool;
  return tool.main(argc, argv);
}

/// @endcond