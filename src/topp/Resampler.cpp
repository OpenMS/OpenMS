// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>

#include <QtGui/QImage>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_Resampler Resampler

@brief Resampler can be used to transform an LC/MS map into a resampled map.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; Resampler &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay  </td>
</tr>
</table>
</CENTER>

When writing an peak file, all spectra are resampled with a new sampling
rate. The number of spectra does not change.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_Resampler.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_Resampler.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPResampler :
  public TOPPBase
{
public:
  TOPPResampler() :
    TOPPBase("Resampler",
             "Transforms an LC/MS map into a resampled map or a PNG image.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", {"mzML"});
    
    registerOutputFile_("out", "<file>", "", "output file in mzML format");
    setValidFormats_("out", {"mzML"});

    registerDoubleOption_("sampling_rate", "<rate>", 0.1,
                          "New sampling rate in m/z dimension (in Th unless ppm flag is set)", false);
    setMinFloat_("sampling_rate", 0.0);

    registerFlag_("ppm", "sampling_rate is given in ppm");
    registerFlag_("align_sampling", "Ensures that sampling is performed equally across the map (same raster is used for all spectra)");

    registerDoubleOption_("min_int_cutoff", "<min intensity>", -1.0,
                          "Intensity cutoff for peaks to be stored in output spectrum (only peaks above this cutoff will be stored, -1 means store all data)", false);

  }

  ExitCodes main_(int, const char **) override
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    double sampling_rate = getDoubleOption_("sampling_rate");
    double min_int_cutoff = getDoubleOption_("min_int_cutoff");
    bool align_sampling = getFlag_("align_sampling");
    bool ppm = getFlag_("ppm");
    PeakMap exp;
    exp.updateRanges();

    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    Param resampler_param;
    resampler_param.setValue("spacing", sampling_rate);
    resampler_param.setValue("ppm", ppm ? "true" : "false");

    LinearResamplerAlign lin_resampler; // LinearResampler does not know about ppm!
    lin_resampler.setParameters(resampler_param);
    if (!align_sampling)
    {
      // resample every scan
      for (Size i = 0; i < exp.size(); ++i)
      {
        lin_resampler.raster(exp[i]);
      }
    }
    else if(!exp.RangeRT::isEmpty())
    {
      // start with even position
      auto start_pos = floor(exp.getMinRT());

      // resample every scan
      for (Size i = 0; i < exp.size(); ++i)
      {
        lin_resampler.raster_align(exp[i], start_pos, exp.getMaxRT());
      }
    }

    if (min_int_cutoff >= 0.0)
    {
      ThresholdMower mow;
      Param p;
      p.setValue("threshold", min_int_cutoff);
      mow.setParameters(p);
      mow.filterPeakMap(exp);
    }

    //annotate output with data processing info
    addDataProcessing_(exp,
                       getProcessingInfo_(DataProcessing::DATA_PROCESSING));

    //store output
    FileHandler().storeExperiment(out, exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPResampler tool;
  return tool.main(argc, argv);
}

/// @endcond
