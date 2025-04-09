// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn, Leon Kuchenbecker$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>

using std::map;
using std::pair;
using std::vector;

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDRipper IDRipper

@brief IDRipper splits the protein/peptide identifications of an idXML file into several idXML files according their annotated file origin.

<center>
<table>
<tr>
<th ALIGN = "center"> potential predecessor tools </td>
<td VALIGN="middle" ROWSPAN=3> &rarr; IDRipper&rarr;</td>
<th ALIGN = "center"> potential successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDFilter</td>
<td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDMapper</td>
</tr>
</table>
</center>

<B>Example</B>

<p>Assuming each peptide identification in a given idXML file is annotated with its file origin (e.g. IDRipper_test.idXML) :</p>

@p <tt>&lt;userParam type="string" name="file_origin" value="IDMerger1_test.idXML"/&gt;</tt> or <br />
@p <tt>&lt;userParam type="string" name="file_origin" value="IDMerger2_test.idXML"/&gt;</tt>

<p>Obviously the file contains protein/peptide identifications of IDMerger1_test.idXML and IDMerger2_test.idXML.</p>

<p>Calling IDRipper with an input file (here: @p -in IDRipper_test.idXML) and an output directory (via @p out) will
result in two idXML files stored in the specified directory and named according to their file origin.</p>

<p>In theory, merging files with @p IDMerger and splitting the resulting file with @p IDRipper will result in the original input files.

<B>NOTE: The meta value "file_origin" is removed by the <tt>IDSplitter</tt>!</B>

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDRipper.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDRipper.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDRipper : public TOPPBase
{
public:
  TOPPIDRipper() : TOPPBase("IDRipper", "Split protein/peptide identification file into several files according to identification run and annotated file origin.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file, in which the protein/peptide identifications must be tagged with 'file_origin'");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputPrefix_("out", "<directory>", "", "Path to the output directory to write the ripped files to.", true, false);
    registerFlag_("numeric_filenames", "Do not infer output filenames from spectra_data or file_origin but use the input filename with numeric suffixes.");
    registerFlag_("split_ident_runs", "Split different identification runs into separate files.");
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String file_name = getStringOption_("in");
    String out_dir = getStringOption_("out");
    bool numeric_filenames = getFlag_("numeric_filenames");
    bool split_ident_runs = getFlag_("split_ident_runs");

    String output_directory = QFileInfo(out_dir.toQString()).absoluteFilePath().toStdString();

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    FileHandler().loadIdentifications(file_name, proteins, peptides, {FileTypes::IDXML});

    // ensure protein and peptide identifications are presented, otherwise we don't have to rip anything anyhow
    if (proteins.empty() || peptides.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "idXML file has to store protein and peptide identifications!");
    }

    IDRipper::RipFileMap ripped;

    // rip the idXML-file into several idXML according to the annotated file origin
    IDRipper ripper;
    ripper.rip(ripped, proteins, peptides, numeric_filenames, split_ident_runs);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    for (IDRipper::RipFileMap::iterator it = ripped.begin(); it != ripped.end(); ++it)
    {
      const IDRipper::RipFileIdentifier& rfi = it->first;
      const IDRipper::RipFileContent& rfc = it->second;

      QString output = output_directory.toQString();

      String out_fname;
      if (numeric_filenames)
      {
        String s_ident_run_idx = split_ident_runs ? '_' + String(rfi.ident_run_idx) : "";
        String s_file_origin_idx = '_' + String(rfi.file_origin_idx);
        out_fname = QFileInfo(file_name.toQString()).completeBaseName().toStdString() + s_ident_run_idx + s_file_origin_idx + ".idXML";
      }
      else
      {
        out_fname = QFileInfo(rfi.out_basename.toQString()).completeBaseName().toStdString() + ".idXML";
      }

      String out = QDir::toNativeSeparators(output.append(QString("/")).append(out_fname.toQString())).toStdString();
      OPENMS_LOG_INFO << "Storing file: '" << out << "'." << std::endl;

      QDir dir(output_directory.toQString());
      FileHandler().storeIdentifications(out, rfc.prot_idents, rfc.pep_idents, {FileTypes::IDXML});
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPIDRipper tool;
  return tool.main(argc, argv);
}

/// @endcond
