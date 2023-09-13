// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeakPickerWavelet PeakPickerWavelet

@brief A tool for peak detection in profile data. Executes the peak picking with the algorithm described in described in Lange et al. (2006) Proc. PSB-06.
<CENTER>
<table>
    <tr>
        <th ALIGN = "center"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=4> &rarr; PeakPickerWavelet &rarr;</td>
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
</CENTER>

The conversion of the ''raw'' ion count data acquired
by the machine into peak lists for further processing
is usually called peak picking. The choice of the algorithm
should mainly depend on the resolution of the data.
As the name implies, the @ref OpenMS::PeakPickerHiRes "high_res"
algorithm is fit for high resolution data whereas in case
of low-resoluted data the @ref OpenMS::PeakPickerCWT "wavelet"
algorithm offers the ability to resolve highly convoluted
and asymmetric signals, separation of overlapping peaks
and nonlinear optimization.

@ref TOPP_example_signalprocessing_parameters is explained in the TOPP tutorial.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PeakPickerWavelet.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PeakPickerWavelet.html

For the parameters of the algorithm section see the algorithm documentation: @n
@ref OpenMS::PeakPickerCWT "PeakPickerCWT" @n

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
    <tr>
    <td BGCOLOR="#EBEBEB"><b>peak_width ("wavelet" only)</b></td>
        <td>0.1</td>
        <td>0.012</td>
    </tr>
</table>

In order to impove the results of the peak detection on low resolution data @ref TOPP_NoiseFilterSGolay or @ref TOPP_NoiseFilterGaussian and @ref TOPP_BaselineFilter can be applied.
For high resolution data this is not necessary.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPickerWavelet :
  public TOPPBase
{
public:
  TOPPPeakPickerWavelet() :
    TOPPBase("PeakPickerWavelet", "Finds mass spectrometric peaks in profile mass spectra.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input profile data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output peak file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    registerFlag_("write_peak_meta_data", "Write additional information about the picked peaks (maximal intensity, left and right area...) into the mzML-file. Attention: this can blow up files, since seven arrays are stored per spectrum!", true);

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return PeakPickerCWT().getDefaults();
  }

  ExitCodes main_(int, const char **) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    bool write_meta_data_arrays(getFlag_("write_peak_meta_data"));
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_exp_raw;
    mz_data_file.load(in, ms_exp_raw);

    if (ms_exp_raw.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    // check for peak type (profile data required)
    if (ms_exp_raw[0].getType(true) == SpectrumSettings::CENTROID)
    {
      writeLogWarn_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
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

    //-------------------------------------------------------------
    // pick
    //-------------------------------------------------------------
    PeakMap ms_exp_peaks;

    Param pepi_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to PeakPickerWavelet", pepi_param, 3);

    PeakPickerCWT pp;
    pp.setLogType(log_type_);
    pp.setParameters(pepi_param);
    try
    {
      pp.pickExperiment(ms_exp_raw, ms_exp_peaks);
    }
    catch (Exception::BaseException & e)
    {
      OPENMS_LOG_ERROR << "Exception caught: " << e.what() << "\n";
      return INTERNAL_ERROR;
    }
    if (!write_meta_data_arrays)
    {
      for (Size i = 0; i < ms_exp_peaks.size(); ++i)
      {
        ms_exp_peaks[i].getFloatDataArrays().clear();
      }
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(ms_exp_peaks, getProcessingInfo_(DataProcessing::PEAK_PICKING));

    mz_data_file.store(out, ms_exp_peaks);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPPeakPickerWavelet tool;
  return tool.main(argc, argv);
}

/// @endcond
