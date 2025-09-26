// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>

#include <iomanip>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IonMobilityBinning IonMobilityBinning

@brief Splits an mzML file with IonMobility frames into multiple mzML files by binning(merging) spectra by their IM values

This tool supports two modes:
- Regular ion mobility: Splits data into a user-defined number of bins
- FAIMS: Automatically splits data by the FAIMS compensation voltages (CVs) present in the file

For regular ion mobility data:
- Useful to convert IM data to a format that can be processed by tools that do not support IM data (e.g. FeatureFinderCentroided or SearchEngines)
- The results of individual bins can be processed separately and then recombined afterwards
- To decide on the number of bins, try running @ref TOPP_FileInfo on the input file to get an idea of the range of IM values present

For FAIMS data:
- Automatically detects FAIMS compensation voltages in the input file
- Creates one output file per unique CV value
- No binning parameters required as the splitting is based on the discrete CV values

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

  std::pair<std::vector<PeakMap>, Math::BinContainer> processFAIMSData_(PeakMap&& experiment)
  {
    set<double> faims_cvs = FAIMSHelper::getCompensationVoltages(experiment);
    auto mzML_bins = IMDataConverter::splitByFAIMSCV(std::move(experiment));
    Size n_bins = mzML_bins.size();
    
    Math::BinContainer im_ranges;
    for (Size i = 0; i < n_bins; ++i)
    {
      double faims_cv = *std::next(faims_cvs.begin(), i);
      im_ranges[i].setMax(faims_cv);
      im_ranges[i].setMin(faims_cv);
    }
    
    return {std::move(mzML_bins), std::move(im_ranges)};
  }

  void writeOutputFiles_(std::vector<PeakMap>& mzML_bins, 
    const Math::BinContainer& im_ranges,
    const String& out_prefix,
    Size n_bins)
  {
    const Size width = String(n_bins).size();
    for (Size i = 0; i < n_bins; ++i)
    {
      ostringstream out_name;
      out_name << out_prefix << "_part" 
      << setw(width) << setfill('0') << (1+i) 
      << "of" << n_bins << "_"
      << im_ranges[i].getMin() << "-"
      << im_ranges[i].getMax() << ".mzML";

      addDataProcessing_(mzML_bins[i], 
          getProcessingInfo_(DataProcessing::ION_MOBILITY_BINNING));
      FileHandler().storeExperiment(out_name.str(), mzML_bins[i], {FileTypes::MZML});
    }
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

    FileHandler().loadExperiment(input_file, experiment, {FileTypes::MZML}, log_type_);
 
    auto [mzML_bins, im_ranges] = IMDataConverter::splitExperimentByIonMobility(std::move(experiment), bins, bin_extension_abs, mz_binning_width, mz_binning_width_unit);
    
    const Size width = String(bins).size();
    for (Size counter = 0; counter < (Size)bins; ++counter)
    {
      ostringstream out_name;
      out_name << out_prefix << "_part" << setw(width) << setfill('0') << (1+counter) << "of" << bins << "_" << im_ranges[counter].getMin() << "-"
               << im_ranges[counter].getMax() << ".mzML ";
      addDataProcessing_(mzML_bins[counter], getProcessingInfo_(DataProcessing::ION_MOBILITY_BINNING));

    // Check if this is FAIMS data and get CVs
    auto cvs = FAIMSHelper::getCompensationVoltages(experiment);
    cvs.erase(IMTypes::DRIFTTIME_NOT_SET);  // remove invalid values

    Size n_bins{};
    std::vector<PeakMap> mzML_bins;
    Math::BinContainer im_ranges;
    if (!cvs.empty())
    {
      // TODO: do MS2 have FAIMS CVs annotated? if not how to keep them?
      // TOOD: currently this is fine as order of getCompensationVoltages is the same as splitByFAIMSCV. Safe?
      std::tie(mzML_bins, im_ranges) = processFAIMSData_(std::move(experiment));
      n_bins = mzML_bins.size();
    }
    else
    {
      // Regular IM data
      std::tie(mzML_bins, im_ranges) = IMDataConverter::splitExperimentByIonMobility(
          std::move(experiment),
          bins,
          bin_extension_abs,
          mz_binning_width,
          mz_binning_width_unit);
      n_bins = bins;
    }
    
    writeOutputFiles_(mzML_bins, im_ranges, out_prefix, n_bins);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIonMobilityBinning tool;
  return tool.main(argc, argv);
}

/// @endcond
