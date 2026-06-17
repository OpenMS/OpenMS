// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Marc Sturm, Clemens Groepl, Lars Nilse, Chris Bielow $
// --------------------------------------------------------------------------

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileInfo.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FileInfo FileInfo
@brief Shows basic information about the data in an %OpenMS readable file.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; FileInfo &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format) </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none ; console or text file</td>
</tr>
</table>
</CENTER>

This tool can show basic information about the data in different file types, such as raw peak, featureXML and consensusXML files. It can
- show information about the data range of a file (m/z, RT, ion mobility, intensity)
- show a statistical summary for intensities, qualities, feature widths, precursor charges, activation methods
- show an overview of the metadata
- validate several XML formats against their XML schema
- check for corrupt data in a file (e.g., duplicate spectra)

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FileInfo.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FileInfo.html

In order to enrich the resulting data of your analysis pipeline or to quickly compare different outcomes of your pipeline you can invoke the aforementioned information of your input data and (intermediary) results.

The actual inspection logic lives in the reusable library class @ref OpenMS::FileInfo, which is also
available from pyOpenMS (it returns a structured result and can render the exact text/TSV shown here).
This tool is a thin command-line wrapper around that class.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFileInfo : public TOPPBase
{
public:
  TOPPFileInfo() : TOPPBase("FileInfo", "Shows basic information about the file, such as data ranges and file type.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    StringList in_types = { "mzData", "mzXML", "mzML", "sqMass", "dta", "dta2d", "mgf", "featureXML", "consensusXML", "idXML", "pepXML", "mzTab", "fid", "mzid", "trafoXML", "fasta", "pqp" };
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", in_types);
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content", false);
    setValidStrings_("in_type", in_types);
    registerOutputFile_("out", "<file>", "", "Optional output file. If left out, the output is written to the command line.", false);
    setValidFormats_("out", {"txt"});
    registerOutputFile_("out_tsv", "<file>", "", "Second optional output file. Tab separated flat text file.", false, true);
    setValidFormats_("out_tsv", {"tsv"});
    registerFlag_("m", "Show meta information about the whole experiment");
    registerFlag_("p", "Shows data processing information");
    registerFlag_("s", "Computes a five-number statistics of intensities, qualities, and widths");
    registerFlag_("d", "Show detailed listing of all spectra and chromatograms (peak files only)");
    registerFlag_("c", "Check for corrupt data in the file (peak files only)");
    registerFlag_("v", "Validate the file only (for mzML, mzData, mzXML, featureXML, idXML, consensusXML, pepXML)");
    registerFlag_("i", "Check whether a given mzML file contains valid indices (conforming to the indexedmzML standard)");
  }

  ExitCodes outputTo_(ostream& os, ostream& os_tsv)
  {
    const std::string in = getStringOption_("in");

    // Resolve the file type up-front so we can reproduce the historical CLI-specific
    // error handling (unknown type; '-i' restricted to mzML) before delegating to the library.
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));
    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = FileHandler::getType(in);
      writeDebug_(std::string("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }
    if (in_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }
    if (getFlag_("i") && in_type != FileTypes::MZML)
    {
      writeLogError_("Error: Can only validate indices for mzML files");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // Map the CLI flags onto the library options. All actual work (loading, computing
    // and rendering the human-readable / TSV report) lives in OpenMS::FileInfo now.
    FileInfo::Options opt;
    opt.forced_type   = in_type; // already resolved above
    opt.meta          = getFlag_("m") || !getStringOption_("out_tsv").empty();
    opt.processing    = getFlag_("p");
    opt.statistics    = getFlag_("s");
    opt.detailed      = getFlag_("d");
    opt.check_corrupt = getFlag_("c");
    opt.validate      = getFlag_("v");
    opt.check_index   = getFlag_("i");

    FileInfo fi;
    FileInfo::Result r = fi.run(in, opt);

    os << FileInfo::toText(r);
    os_tsv << FileInfo::toTSV(r);

    if (opt.check_index && r.validation.index_checked && !r.validation.index_valid)
    {
      return ILLEGAL_PARAMETERS;
    }
    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char **) override
  {
    std::string out = getStringOption_("out");
    std::string out_tsv = getStringOption_("out_tsv");

    ofstream os;
    ofstream os_tsv;
    boost::iostreams::filtering_ostream os_filt;
    boost::iostreams::filtering_ostream os_tsv_filt;


    if (out.empty())
    {
      os_filt.push(getGlobalLogInfo());
    }
    else
    {
      os.open(out);
      if (!os)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out);
      }
      os_filt.push(os);
    }

    if (out_tsv.empty())
    {
      os_tsv_filt.push(boost::iostreams::null_sink());
    }
    else
    {
      os_tsv.open(out_tsv);
      if (!os_tsv)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_tsv);
      }
      os_tsv_filt.push(os_tsv);
    }

    return outputTo_(os_filt, os_tsv_filt);
  }
};


int main(int argc, const char **argv)
{
  TOPPFileInfo tool;
  return tool.main(argc, argv);
}

/// @endcond
