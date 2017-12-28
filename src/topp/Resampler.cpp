// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ Resampler \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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
