// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>
// TODO remove needed here for transform
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_NoiseFilterSGolay NoiseFilterSGolay

@brief  Executes a Savitzky Golay filter to reduce the noise in an MS experiment.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=4> &rarr; NoiseFilterSGolay &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes</td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_Resampler </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes</td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_BaselineFilter</td>
</tr>
</table>
</center>

The idea of the Savitzky Golay filter is to find filter-coefficients
that preserve higher moments, which means to approximate the underlying
function within the moving window by a polynomial of higher order
(typically quadratic or quartic) (see A. Savitzky and M. J. E. Golay,
''Smoothing and Differentiation of Data by Simplified Least Squares Procedures'').

@note The Savitzky Golay filter works only on uniform data (to generate equally spaced data use the @ref TOPP_Resampler tool).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_NoiseFilterSGolay.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_NoiseFilterSGolay.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPNoiseFilterSGolay :
  public TOPPBase
{
public:
  TOPPNoiseFilterSGolay() :
    TOPPBase("NoiseFilterSGolay", "Removes noise from profile spectra by using a Savitzky Golay filter (on uniform (equidistant) data).")
  {
  }

  /**
    @brief Helper class for the Low Memory Noise filtering
  */
  class NFSGolayMzMLConsumer :
    public MSDataWritingConsumer 
  {

  public:

    NFSGolayMzMLConsumer(const String& filename, const SavitzkyGolayFilter& sgf) :
      MSDataWritingConsumer(filename) 
    {
      sgf_ = sgf;
    }

    void processSpectrum_(MapType::SpectrumType& s) override
    {
      sgf_.filter(s);
    }

    void processChromatogram_(MapType::ChromatogramType& c) override 
    {
      sgf_.filter(c);
    }

  private:
    SavitzkyGolayFilter sgf_;
  };

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input raw data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output raw data file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerStringOption_("processOption", "<name>", "inmemory", "Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first", false, true);
    setValidStrings_("processOption", ListUtils::create<String>("inmemory,lowmemory"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return SavitzkyGolayFilter().getDefaults();
  }

  ExitCodes doLowMemAlgorithm(const SavitzkyGolayFilter& sgolay)
  {
    ///////////////////////////////////
    // Create the consumer object, add data processing
    ///////////////////////////////////
    NFSGolayMzMLConsumer sgolayConsumer(out, sgolay);
    sgolayConsumer.addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));

    ///////////////////////////////////
    // Create new MSDataReader and set our consumer
    ///////////////////////////////////
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.transform(in, &sgolayConsumer);

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    in = getStringOption_("in");
    out = getStringOption_("out");
    String process_option = getStringOption_("processOption");

    Param filter_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to filter", filter_param, 3);

    SavitzkyGolayFilter sgolay;
    sgolay.setLogType(log_type_);
    sgolay.setParameters(filter_param);

    if (process_option == "lowmemory")
    {
      return doLowMemAlgorithm(sgolay);
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    FileHandler mz_data_file;
    PeakMap exp;
    mz_data_file.loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    if (exp.empty() && exp.getChromatograms().empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    //check for peak type (profile data required)
    if (!exp.empty() && exp[0].getType(true) == SpectrumSettings::CENTROID)
    {
      writeLogWarn_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
    }

    //check if spectra are sorted
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (!exp[i].isSorted())
      {
        writeLogError_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //check if chromatograms are sorted
    for (Size i = 0; i < exp.getChromatograms().size(); ++i)
    {
      if (!exp.getChromatogram(i).isSorted())
      {
        writeLogError_("Error: Not all chromatograms are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    sgolay.filterExperiment(exp);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::SMOOTHING));

    mz_data_file.storeExperiment(out, exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

  String in;
  String out;

};


int main(int argc, const char ** argv)
{
  TOPPNoiseFilterSGolay tool;
  return tool.main(argc, argv);
}

/// @endcond
