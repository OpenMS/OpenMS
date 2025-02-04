// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
// TODO remove needed here for transform
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_NoiseFilterGaussian NoiseFilterGaussian

@brief  Executes a Gaussian filter to reduce the noise in an MS experiment.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=4> &rarr; NoiseFilterGaussian &rarr;</td>
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

The Gaussian filter is a peak area preserving low-pass filter and is characterized by narrow bandwidths,
sharp cutoffs, and low passband ripple.

@note The Gaussian filter works for uniform as well as for non-uniform data.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_NoiseFilterGaussian.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_NoiseFilterGaussian.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPNoiseFilterGaussian :
  public TOPPBase
{
public:
  TOPPNoiseFilterGaussian() :
    TOPPBase("NoiseFilterGaussian", "Removes noise from profile spectra by using Gaussian filter (on uniform as well as non-uniform data).")
  {
  }

  /**
    @brief Helper class for the Low Memory Noise filtering
  */
  class NFGaussMzMLConsumer :
    public MSDataWritingConsumer 
  {

  public:

    NFGaussMzMLConsumer(const String& filename, const GaussFilter& gf) :
      MSDataWritingConsumer(filename) 
    {
      gf_ = gf;
    }

    void processSpectrum_(MapType::SpectrumType& s) override
    {
      gf_.filter(s);
    }

    void processChromatogram_(MapType::ChromatogramType& c) override 
    {
      gf_.filter(c);
    }

  private:
    GaussFilter gf_;
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
    return GaussFilter().getDefaults();
  }

  ExitCodes doLowMemAlgorithm(const GaussFilter& gauss)
  {
    ///////////////////////////////////
    // Create the consumer object, add data processing
    ///////////////////////////////////
    NFGaussMzMLConsumer gaussConsumer(out, gauss);
    gaussConsumer.addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));

    ///////////////////////////////////
    // Create new MSDataReader and set our consumer
    ///////////////////////////////////
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.transform(in, &gaussConsumer);

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

    GaussFilter gauss;
    gauss.setLogType(log_type_);
    gauss.setParameters(filter_param);

    if (process_option == "lowmemory")
    {
      return doLowMemAlgorithm(gauss);
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

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
    try
    {
      gauss.filterExperiment(exp);
    }
    catch (Exception::IllegalArgument & e)
    {
      writeLogError_(String("Error: ") + e.what());
      return INCOMPATIBLE_INPUT_DATA;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::SMOOTHING));

    FileHandler().storeExperiment(out, exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

  String in;
  String out;
};


int main(int argc, const char ** argv)
{
  TOPPNoiseFilterGaussian tool;
  return tool.main(argc, argv);
}

/// @endcond
