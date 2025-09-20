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

@brief Converts IdXML files to parquet or TSV format following quantms.io PSM specification.

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
and converts them to parquet or TSV format following the quantms.io PSM (Peptide 
Spectrum Match) specification.

The output format is determined by the file extension:
- .parquet: Apache Parquet format with embedded metadata
- .tsv/.csv/.txt: Tab-separated values format (metadata not preserved)

The output file contains PSM data with columns following the quantms.io PSM specification:
- sequence: unmodified peptide sequence
- peptidoform: peptide sequence with modifications
- modifications: peptide modifications (null for now)
- precursor_charge: precursor charge
- posterior_error_probability: PEP score from metavalues (nullable)
- is_decoy: decoy flag (0=target, 1=decoy) based on target_decoy metavalue
- calculated_mz: theoretical m/z from sequence
- observed_mz: experimental precursor m/z
- additional_scores: additional scores (null for now)
- mp_accessions: protein accessions (null for now)
- predicted_rt: predicted retention time (null for now)
- reference_file_name: reference file name
- cv_params: CV parameters (null for now)
- scan: scan identifier
- rt: retention time in seconds (nullable)
- ion_mobility: ion mobility value (nullable, null for now)
- num_peaks: number of peaks (nullable, null for now)
- mz_array: m/z values array (null for now)
- intensity_array: intensity values array (null for now)

Only the first peptide hit per peptide identification is processed (no rank field).
PEP scores are automatically detected from metavalues using known PEP score names.

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
    TOPPBase("QuantmsIOConverter", "Converts IdXML files to parquet or TSV format following quantms.io PSM specification.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input idXML file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));

    registerOutputFile_("out", "<file>", "", "Output file (parquet, TSV, or CSV format)", true);
    setValidFormats_("out", ListUtils::create<String>("parquet,tsv,csv,txt"));
  }

  ExitCodes main_(int, const char**) override
  {
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
    OPENMS_LOG_INFO << "Converting to output format..." << endl;
    QuantmsIO quantms_io;
    quantms_io.store(out, protein_identifications, peptide_identifications);

    OPENMS_LOG_INFO << "Conversion completed successfully." << endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPQuantmsIOConverter tool;
  return tool.main(argc, argv);
}

/// @endcond