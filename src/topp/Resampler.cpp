// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>

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
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "",
                        "output file in mzML format");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

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
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, exp);

    Param resampler_param;
    resampler_param.setValue("spacing", sampling_rate);
    if (ppm) resampler_param.setValue("ppm", "true");
    else resampler_param.setValue("ppm", "false");

    if (!align_sampling)
    {
      LinearResampler lin_resampler;
      lin_resampler.setParameters(resampler_param);

      // resample every scan
      for (Size i = 0; i < exp.size(); ++i)
      {
        lin_resampler.raster(exp[i]);
      }
    }
    else
    {
      LinearResamplerAlign lin_resampler;
      lin_resampler.setParameters(resampler_param);

      bool start_pos_set = false;
      bool end_pos_set = false;
      double start_pos = 0.0;
      double end_pos = 0.0;
      // get max / min positions across whole map
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (!exp[i].empty() && (!start_pos_set || exp[i][0].getMZ() < start_pos) )
        {
          start_pos = exp[i][0].getMZ();
          start_pos_set = true;
        }
        if (!exp[i].empty() && (!end_pos_set || exp[i].back().getMZ() > end_pos) )
        {
          end_pos = exp[i].back().getMZ();
          end_pos_set = true;
        }
      }

      if (start_pos_set)
      {
        // start with even position
        start_pos = std::floor(start_pos);

        // resample every scan
        for (Size i = 0; i < exp.size(); ++i)
        {
          lin_resampler.raster_align(exp[i], start_pos, end_pos);
        }
      }
    }

    if (min_int_cutoff >= 0.0)
    {
      for (Size i = 0; i < exp.size(); ++i)
      {
        MSSpectrum tmp = exp[i];
        tmp.clear(false);
        for (Size j = 0; j < exp[i].size(); j++)
        {
          if (exp[i][j].getIntensity() > min_int_cutoff)
          {
            tmp.push_back(exp[i][j]);
          }
        }
        // swap
        exp[i] = tmp;
      }
    }

    //clear meta data because they are no longer meaningful
    exp.clearMetaDataArrays();

    //annotate output with data processing info
    addDataProcessing_(exp,
                       getProcessingInfo_(DataProcessing::DATA_PROCESSING));

    //store output
    f.store(out, exp);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPResampler tool;
  return tool.main(argc, argv);
}

/// @endcond
