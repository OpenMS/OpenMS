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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_MapAlignerSpectrum MapAlignerSpectrum

        @brief Corrects retention time distortions between maps by aligning spectra.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapAlignerSpectrum \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature finding algorithm) </td>
        </tr>
    </table>
</CENTER>

    This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

        All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

    The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types of data they can be applied to. Here, an experimental algorithm based on spectrum alignment is implemented. It is only applicable to peak maps (mzML format).
    For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmSpectrumAlignment "algorithm documentation".

    @see @ref TOPP_MapAlignerIdentification @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapRTTransformer

    Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested mostly with the "interpolated" model. The different available models are:
    - @ref OpenMS::TransformationModelLinear "linear": Linear model.
    - @ref OpenMS::TransformationModelInterpolated "interpolated": Smoothing spline (non-linear).

    The following parameters control the modeling of RT transformations (they can be set in the "model" section of the INI file):
    @htmlinclude OpenMS_MapAlignerSpectrumModel.parameters @n

    <B>The command line parameters of this tool are:</B> @n
    @verbinclude TOPP_MapAlignerSpectrum.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MapAlignerSpectrum.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerSpectrum :
  public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerSpectrum() :
    TOPPMapAlignerBase("MapAlignerSpectrum", "Corrects retention time distortions between maps by spectrum alignment.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    String formats = "mzML";
    // no support for a reference file yet:
    TOPPMapAlignerBase::registerOptionsAndFlags_(formats, REF_NONE);
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmSpectrumAlignment algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return getModelDefaults("interpolated");
    }
    return Param(); // shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes ret = checkParameters_();
    if (ret != EXECUTION_OK) return ret;

    MapAlignmentAlgorithmSpectrumAlignment algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    StringList ins = getStringList_("in");
    StringList outs = getStringList_("out");
    StringList trafos = getStringList_("trafo_out");
    Param model_params = getParam_().copy("model:", true);
    String model_type = model_params.getValue("type");
    model_params = model_params.copy(model_type + ":", true);
    std::vector<TransformationDescription> transformations;

    //-------------------------------------------------------------
    // perform peak alignment
    //-------------------------------------------------------------
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    // load input
    std::vector<PeakMap > peak_maps(ins.size());
    MzMLFile f;
    f.setLogType(log_type_);
    progresslogger.startProgress(0, ins.size(), "loading input files");
    for (Size i = 0; i < ins.size(); ++i)
    {
      progresslogger.setProgress(i);
      f.load(ins[i], peak_maps[i]);
    }
    progresslogger.endProgress();

    // try to align
    algorithm.align(peak_maps, transformations);
    if (model_type != "none")
    {
      for (vector<TransformationDescription>::iterator it = 
             transformations.begin(); it != transformations.end(); ++it)
      {
        it->fitModel(model_type, model_params);
      }
    }

    // write output
    progresslogger.startProgress(0, outs.size(), "applying RT transformations and writing output files");
    for (Size i = 0; i < outs.size(); ++i)
    {
      progresslogger.setProgress(i);

      MapAlignmentTransformer::transformRetentionTimes(peak_maps[i], 
                                                       transformations[i]);
      // annotate output with data processing info
      addDataProcessing_(peak_maps[i], 
                         getProcessingInfo_(DataProcessing::ALIGNMENT));

      f.store(outs[i], peak_maps[i]);
    }
    progresslogger.endProgress();

    if (!trafos.empty())
    {
      TransformationXMLFile trafo_file;
      for (Size i = 0; i < transformations.size(); ++i)
      {
        trafo_file.store(trafos[i], transformations[i]);
      }
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMapAlignerSpectrum tool;
  return tool.main(argc, argv);
}

/// @endcond
