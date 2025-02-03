// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <QFile>
#include <iomanip>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MzMLSplitter MzMLSplitter

@brief Splits an mzML file into multiple parts

This utility will split an input mzML file into @e N parts, with an approximately equal number of spectra and chromatograms in each part.
@e N is set by the parameter @p parts; optionally only spectra (parameter @p no_chrom) or only chromatograms (parameter @p no_spec) can be transferred to the output.

Alternatively to setting the number of parts directly, a target maximum file size for the parts can be specified (parameters @p size and @p unit).
The number of parts is then calculated by dividing the original file size by the target and rounding up.
Note that the resulting parts may actually be bigger than the target size (due to meta data that is included in every part) or
that more parts than necessary may be produced (if spectra or chromatograms are removed via @p no_spec/@p no_chrom).

This tool cannot be used as part of a TOPPAS workflow, because the number of output files is variable.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MzMLSplitter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MzMLSplitter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMzMLSplitter : public TOPPBase
{
public:
  TOPPMzMLSplitter() : TOPPBase("MzMLSplitter", "Splits an mzML file into multiple parts")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputPrefix_("out", "<prefix>", "", "Prefix for output files ('_part1of2.mzML' etc. will be appended; default: same as 'in' without the file extension)", false);
    registerIntOption_("parts", "<num>", 1, "Number of parts to split into (takes precedence over 'size' if set)", false);
    setMinInt_("parts", 1);
    registerIntOption_("size", "<num>", 0, "Approximate upper limit for resulting file sizes (in 'unit')", false);
    setMinInt_("size", 0);
    registerStringOption_("unit", "<choice>", "MB", "Unit for 'size' (base 1024)", false);
    setValidStrings_("unit", ListUtils::create<String>("KB,MB,GB"));
    // @TODO:
    // registerFlag_("precursor", "Make sure precursor spectra end up in the same part as their fragment spectra");
    registerFlag_("no_chrom", "Remove chromatograms, keep only spectra.");
    registerFlag_("no_spec", "Remove spectra, keep only chromatograms.");
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");

    if (out.empty())
    {
      out = FileHandler::stripExtension(in);
    }

    bool no_chrom = getFlag_("no_chrom"), no_spec = getFlag_("no_spec");
    if (no_chrom && no_spec)
    {
      writeLogError_("Error: 'no_chrom' and 'no_spec' cannot be used together");
      return ILLEGAL_PARAMETERS;
    }

    Size parts = getIntOption_("parts"), size = getIntOption_("size");
    if (parts == 1)
    {
      if (size == 0)
      {
        writeLogError_("Error: Higher value for parameter 'parts' or 'size' required");
        return ILLEGAL_PARAMETERS;
      }

      QFile mzml_file(in.toQString());
      // use float here to avoid too many decimals in output below:
      float total_size = mzml_file.size();
      String unit = getStringOption_("unit");
      if (unit == "KB")
        total_size /= 1024;
      else if (unit == "MB")
        total_size /= (1024 * 1024);
      else
        total_size /= (1024 * 1024 * 1024); // "GB"

      writeLogInfo_("File size: " + String(total_size) + " " + unit);
      parts = ceil(total_size / size);
    }
    writeLogInfo_("Splitting file into " + String(parts) + " parts...");

    PeakMap experiment;
    FileHandler().loadExperiment(in, experiment, {FileTypes::MZML});

    vector<MSSpectrum> spectra;
    vector<MSChromatogram> chromatograms;

    if (no_spec)
    {
      experiment.getSpectra().clear();
    }
    else
    {
      experiment.getSpectra().swap(spectra);
    }

    if (no_chrom)
    {
      experiment.getChromatograms().clear();
    }
    else
    {
      experiment.getChromatograms().swap(chromatograms);
    }

    writeLogInfo_("Total spectra: " + String(spectra.size()));
    writeLogInfo_("Total chromatograms: " + String(chromatograms.size()));

    Size spec_start = 0, chrom_start = 0;
    Size width = String(parts).size();
    for (Size counter = 1; counter <= parts; ++counter)
    {
      ostringstream out_name;
      out_name << out << "_part" << setw(width) << setfill('0') << counter << "of" << parts << ".mzML";
      PeakMap part = experiment;
      addDataProcessing_(part, getProcessingInfo_(DataProcessing::FILTERING));

      Size remaining = parts - counter + 1;
      Size n_spec = ceil((spectra.size() - spec_start) / double(remaining));
      if (n_spec > 0)
      {
        part.reserveSpaceSpectra(n_spec);
        for (Size i = spec_start; i < spec_start + n_spec; ++i)
        {
          part.addSpectrum(std::move(spectra[i]));
        }
      }
      spec_start += n_spec;

      Size n_chrom = ceil((chromatograms.size() - chrom_start) / double(remaining));
      if (n_chrom > 0)
      {
        part.reserveSpaceChromatograms(n_chrom);
        for (Size i = chrom_start; i < chrom_start + n_chrom; ++i)
        {
          part.addChromatogram(std::move(chromatograms[i]));
        }
      }
      chrom_start += n_chrom;

      writeLogInfo_("Part " + String(counter) + ": " + String(n_spec) + " spectra, " + String(n_chrom) + " chromatograms");
      FileHandler().storeExperiment(out_name.str(), part, {FileTypes::MZML});
    }

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPMzMLSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond
