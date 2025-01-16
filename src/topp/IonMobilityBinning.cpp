// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>

#include <iomanip>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IonMobilityBinning IonMobilityBinning

@brief Splits an mzML file with IonMobility frames into multiple mzML files by binning(merging) spectra by their IM values

Useful to convert IM data to a format that can be processed by tools that do not support IM data (e.g. FeatureFinderCentroided or SearchEngines).
The results of individual bins can be processed separately and then recombined afterwards.

To decide on the number of bins, try running @ref TOPP_FileInfo on the input file to get an idea of the range of IM values present.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IonMobilityBinning.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IonMobilityBinning.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIonMobilityBinning :
  public TOPPBase
{
public:

  TOPPIonMobilityBinning() :
    TOPPBase("IonMobilityBinning", "Splits an mzML file with IonMobility frames into multiple mzML files by binning(merging) spectra by their IM values")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (containing RT, IM, m/z, i.e. IM-frames).");
    setValidFormats_("in", {"mzML"});
    registerOutputPrefix_("out", "<directory>", "", "Path to the output directory to write the binned mzML files to.", true, false);
    registerIntOption_("bins", "<number>", 5, "Number of ion mobility bins to split the input file into", false, false);
    registerDoubleOption_("bin_extension_abs", "<number>", 0.0,
                          "Absolute extension of the bin in IM dimension (causes bins to overlap by 2x this value)", false, false);
    registerTOPPSubsection_("SpectraMerging", "Options for merging spectra within the same IM bin (from the same IM-frame)");
    registerDoubleOption_("SpectraMerging:mz_binning_width", "<number>", 0.01, "Width of the m/z bins", false, false);
    registerStringOption_("SpectraMerging:mz_binning_width_unit", "<unit>", "Da", "Unit of the m/z bin width", false, false);
    setValidStrings_("SpectraMerging:mz_binning_width_unit", {"Da", "ppm"});
    
  }

  ExitCodes main_(int, const char **) override
  {
    String input_file = getStringOption_("in");
    String out_prefix = getStringOption_("out");
    int bins = getIntOption_("bins");
    double bin_extension_abs = getDoubleOption_("bin_extension_abs");
    double mz_binning_width = getDoubleOption_("SpectraMerging:mz_binning_width");
    MZ_UNITS mz_binning_width_unit = getStringOption_("SpectraMerging:mz_binning_width_unit") == "Da" ? MZ_UNITS::DA : MZ_UNITS::PPM;

    PeakMap experiment;
    FileHandler().loadExperiment(input_file, experiment, {FileTypes::MZML});
 
    auto [mzML_bins, im_ranges] = IMDataConverter::splitExperimentByIonMobility(std::move(experiment), bins, bin_extension_abs, mz_binning_width, mz_binning_width_unit);
    
    const Size width = String(bins).size();
    for (Size counter = 0; counter < (Size)bins; ++counter)
    {
      ostringstream out_name;
      out_name << out_prefix << "_part" << setw(width) << setfill('0') << (1+counter) << "of" << bins << "_" << im_ranges[counter].getMin() << "-"
               << im_ranges[counter].getMax() << ".mzML ";
      addDataProcessing_(mzML_bins[counter], getProcessingInfo_(DataProcessing::ION_MOBILITY_BINNING));

      FileHandler().storeExperiment(out_name.str(), mzML_bins[counter], {FileTypes::MZML});
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIonMobilityBinning tool;
  return tool.main(argc, argv);
}

/// @endcond
