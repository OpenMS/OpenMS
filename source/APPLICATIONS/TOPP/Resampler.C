// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/MISC/BilinearInterpolation.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <QtGui/QImage>

using namespace OpenMS;
using namespace OpenMS::Math;
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

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "",
                        "output file in mzML format");
    setValidFormats_("out", StringList::create("mzML"));

    registerDoubleOption_("sampling_rate", "<rate>", 0.1,
                          "New sampling rate in m/z dimension", false);
    setMinFloat_("sampling_rate", 0.0);

  }

  ExitCodes main_(int, const char **)
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    MSExperiment<> exp;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, exp);

    DoubleReal sampling_rate = getDoubleOption_("sampling_rate");

    LinearResampler lin_resampler;
    Param resampler_param;
    resampler_param.setValue("spacing", sampling_rate);
    lin_resampler.setParameters(resampler_param);

    // resample every scan
    for (Size i = 0; i < exp.size(); ++i)
    {
      lin_resampler.raster(exp[i]);
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
