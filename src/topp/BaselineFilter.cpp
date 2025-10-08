// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_BaselineFilter BaselineFilter

@brief Executes the top-hat filter to remove the baseline of an MS experiment.

<CENTER>
<table>
    <tr>
        <th ALIGN = "center"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=2> &rarr; BaselineFilter &rarr;</td>
        <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay, @n @ref TOPP_NoiseFilterGaussian </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes @n (or ID engines on MS/MS data) </td>
    </tr>
</table>
</CENTER>

This nonlinear filter, known as the top-hat operator in morphological
mathematics (see Soille, ''Morphological Image Analysis''), is independent
of the underlying baseline shape.  It is able to detect an over brightness
even if the environment is not uniform.  The principle is based on the
subtraction of a signal from its opening (erosion followed by a dilation).
The size the structuring element (here a flat line) being conditioned by the
width of the lineament (in our case the maximum width of a mass
spectrometric peak) to be detected.

@note The top-hat filter works only on roughly uniform data!
      To generate equally-spaced data you can use the @ref TOPP_Resampler.

@note The length (given in Thomson) of the structuring element should be wider than the
maximum peak width in the raw data.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_BaselineFilter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_BaselineFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBaselineFilter :
  public TOPPBase
{
public:
  TOPPBaselineFilter() :
    TOPPBase("BaselineFilter", "Removes the baseline from profile spectra using a top-hat filter.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input raw data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output raw data file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    registerDoubleOption_("struc_elem_length", "<size>", 3, "Length of the structuring element (should be wider than maximal peak width - see documentation).", false);
    registerStringOption_("struc_elem_unit", "<unit>", "Thomson", "Unit of 'struc_elem_length' parameter.", false);
    setValidStrings_("struc_elem_unit", ListUtils::create<String>("Thomson,DataPoints"));
    registerStringOption_("method", "<string>", "tophat", "The name of the morphological filter to be applied. If you are unsure, use the default.", false);
    setValidStrings_("method", ListUtils::create<String>("identity,erosion,dilation,opening,closing,gradient,tophat,bothat,erosion_simple,dilation_simple"));
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap ms_exp;
    FileHandler().loadExperiment(in, ms_exp, {FileTypes::MZML}, log_type_);

    if (ms_exp.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    // check for peak type (raw data required)
    if (ms_exp[0].getType(true) == SpectrumSettings::CENTROID)
    {
      writeLogWarn_("Warning: OpenMS peak type estimation indicates that this is not raw data!");
    }

    //check if spectra are sorted
    for (Size i = 0; i < ms_exp.size(); ++i)
    {
      if (!ms_exp[i].isSorted())
      {
        writeLogError_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    MorphologicalFilter morph_filter;
    morph_filter.setLogType(log_type_);

    Param parameters;
    parameters.setValue("struc_elem_length", getDoubleOption_("struc_elem_length"));
    parameters.setValue("struc_elem_unit", getStringOption_("struc_elem_unit"));
    parameters.setValue("method", getStringOption_("method"));

    morph_filter.setParameters(parameters);
    morph_filter.filterExperiment(ms_exp);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(ms_exp, getProcessingInfo_(DataProcessing::BASELINE_REDUCTION));

    FileHandler().storeExperiment(out, ms_exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

};




int main(int argc, const char ** argv)
{
  TOPPBaselineFilter tool;
  return tool.main(argc, argv);
}

/// @endcond
