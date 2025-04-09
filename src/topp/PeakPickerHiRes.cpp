// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
// TODO remove needed here for transform
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeakPickerHiRes PeakPickerHiRes

@brief A tool for peak detection in profile data. Executes the peak picking with @ref OpenMS::PeakPickerHiRes "high_res" algorithm.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=4> &rarr; PeakPickerHiRes &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_BaselineFilter </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=3> any tool operating on MS peak data @n (in mzML format)</td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterGaussian </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay </td>
</tr>
</table>
</center>

Reference:\n
Weisser <em>et al.</em>: <a href="https://doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

The conversion of the "raw" ion count data acquired
by the machine into peak lists for further processing
is usually called peak picking or centroiding. The choice of the algorithm
should mainly depend on the resolution of the data.
As the name implies, the @ref OpenMS::PeakPickerHiRes "high_res"
algorithm is fit for high resolution (Orbitrap or FTICR) data.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PeakPickerHiRes.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PeakPickerHiRes.html

For the parameters of the algorithm section see the algorithm documentation: @ref OpenMS::PeakPickerHiRes "PeakPickerHiRes"

Be aware that applying the algorithm to already picked data results in an error message and program exit or corrupted output data.
Advanced users may skip the check for already centroided data using the flag "-force" (useful e.g. if spectrum annotations in the data files are wrong).

In the following table you, can find example values of the most important algorithm parameters for
different instrument types. @n These parameters are not valid for all instruments of that type,
but can be used as a starting point for finding suitable parameters.
<table>
<tr BGCOLOR="#EBEBEB">
<td>&nbsp;</td>
<td><b>Q-TOF</b></td>
<td><b>LTQ Orbitrap</b></td>
</tr>
<tr>
<td BGCOLOR="#EBEBEB"><b>signal_to_noise</b></td>
<td>2</td>
<td>0</td>
</tr>
</table>
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPickerHiRes :
  public TOPPBase
{
public:
  TOPPPeakPickerHiRes() :
    TOPPBase("PeakPickerHiRes", "Finds mass spectrometric peaks in profile mass spectra.")
  {
  }

protected:

  /**
    @brief Helper class for the Low Memory peak-picking
  */
  class PPHiResMzMLConsumer :
    public MSDataWritingConsumer
  {

  public:

    PPHiResMzMLConsumer(String filename, const PeakPickerHiRes& pp) :
      MSDataWritingConsumer(std::move(filename)),
      ms_levels_(pp.getParameters().getValue("ms_levels").toIntVector())
    {
      pp_ = pp;
    }

    void processSpectrum_(MapType::SpectrumType& s) override
    {
      if (ms_levels_.empty()) //auto mode
      {
        if (s.getType() == SpectrumSettings::CENTROID)
        {
          return;
        }
      }
      else if (!ListUtils::contains(ms_levels_, s.getMSLevel()))
      {
        return;
      }

      MapType::SpectrumType sout;
      pp_.pick(s, sout);
      s = std::move(sout);
    }

    void processChromatogram_(MapType::ChromatogramType & c) override
    {
      MapType::ChromatogramType c_out;
      pp_.pick(c, c_out);
      c = std::move(c_out);
    }

  private:

    PeakPickerHiRes pp_;
    std::vector<Int> ms_levels_;
  };

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input profile data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output peak file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerStringOption_("processOption", "<name>", "inmemory", "Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first", false, true);
    setValidStrings_("processOption", ListUtils::create<String>("inmemory,lowmemory"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return PeakPickerHiRes().getDefaults();
  }

  ExitCodes doLowMemAlgorithm(const PeakPickerHiRes& pp)
  {
    ///////////////////////////////////
    // Create the consumer object, add data processing
    ///////////////////////////////////
    PPHiResMzMLConsumer pp_consumer(out, pp);
    pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));

    ///////////////////////////////////
    // Create new MSDataReader and set our consumer
    ///////////////////////////////////
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.transform(in, &pp_consumer);

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

    Param pepi_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes", pepi_param, 3);

    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pepi_param);

    if (process_option == "lowmemory")
    {
      return doLowMemAlgorithm(pp);
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    PeakMap ms_exp_raw;
    FileHandler().loadExperiment(in, ms_exp_raw, {FileTypes::MZML}, log_type_);

    if (ms_exp_raw.empty() && ms_exp_raw.getChromatograms().empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    //check if spectra are sorted
    for (Size i = 0; i < ms_exp_raw.size(); ++i)
    {
      if (!ms_exp_raw[i].isSorted())
      {
        writeLogError_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //check if chromatograms are sorted
    for (Size i = 0; i < ms_exp_raw.getChromatograms().size(); ++i)
    {
      if (!ms_exp_raw.getChromatogram(i).isSorted())
      {
        writeLogError_("Error: Not all chromatograms are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //-------------------------------------------------------------
    // pick
    //-------------------------------------------------------------
    PeakMap ms_exp_peaks;
    bool check_spectrum_type = !getFlag_("force");
    pp.pickExperiment(ms_exp_raw, ms_exp_peaks, check_spectrum_type);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    //annotate output with data processing info
    addDataProcessing_(ms_exp_peaks, getProcessingInfo_(DataProcessing::PEAK_PICKING));
    FileHandler().storeExperiment(out, ms_exp_peaks, {FileTypes::MZML});

    return EXECUTION_OK;
  }

  // parameters
  String in;
  String out;
};


int main(int argc, const char ** argv)
{
  TOPPPeakPickerHiRes tool;
  return tool.main(argc, argv);
}

/// @endcond
