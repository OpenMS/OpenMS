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
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <iostream> // for "cout"
#include <OpenMS/APPLICATIONS/TOPPBase.h>


#include <OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h> 

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_OpenSwathConfidenceScoring OpenSwathConfidenceScoring

    @brief Computes confidence scores for OpenSwath results.

    <CENTER>
        <table>
            <tr>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
                <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenSwathConfidenceScoring \f$ \longrightarrow \f$</td>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathFeatureXMLToTSV </td>
            </tr>
        </table>
    </CENTER>

    This is an implementation of the SRM scoring algorithm described in:

    Malmstroem, L.; Malmstroem, J.; Selevsek, N.; Rosenberger, G. & Aebersold, R.:\n
    <a href="http://dx.doi.org/10.1021/pr200844d">Automated workflow for large-scale selected reaction monitoring experiments.</a>\n
    J. Proteome Res., 2012, 11, 1644-1653

    It has been adapted for the scoring of OpenSwath results.

    The algorithm compares SRM/MRM features (peak groups) to assays and computes scores for the agreements. Every feature is compared not only to the "true" assay that was used to acquire the corresponding ion chromatograms, but also to a number (parameter @p decoys) of unrelated - but real - assays selected at random from the assay library (parameter @p lib). This serves to establish a background distribution of scores, against which the significance of the "true" score can be evaluated. The final confidence value of a feature is the local false discovery rate (FDR), calculated as the fraction of decoy assays that score higher than the "true" assay against the feature. In the output feature map, every feature is annotated with its local FDR in the meta value "local_FDR" (a "userParam" element in the featureXML), and its overall quality is set to "1 - local_FDR".

    The agreement of a feature and an assay is assessed based on the difference in retention time (RT) and on the deviation of relative transition intensities. The score @e S is computed using a binomial generalized linear model (GLM) of the form:

    @f[
    S = \frac{1}{1 + \exp(-(a + b \cdot \Delta_{RT}^2 + c \cdot d_{int}))}
    @f]

    The meanings of the model terms are as follows:

    @f$ \Delta_{RT} @f$: Observed retention times are first mapped to the scale of the assays (parameter @p trafo), then all RTs are scaled to the range 0 to 100 (based on the lowest/highest RT in the assay library). @f$ \Delta_{RT} @f$ is the absolute difference of the scaled RTs; note that this is squared in the scoring model.

    @f$ d_{int} @f$: To compute the intensity distance, the @e n (advanced
    parameter @p transitions) most intensive transitions of the feature are
    selected. For comparing against the "true" assay, the same transitions are
    considered; otherwise, the same number of most intensive transitions from
    the decoy assay. Transition intensities are scaled to a total of 1 per
    feature/assay and are ordered by the product (Q3) m/z value. Then the
    Manhattan distance of the intensity vectors is calculated (Malmstroem et
    al. used the RMSD instead, which has been replaced here to be independent
    of the number of transitions).

    @f$ a, b, c @f$: Model coefficients, stored in the advanced parameters @p GLM:intercept, @p GLM:delta_rt, and @p GLM:dist_int. The default values were estimated based on the training dataset used in the Malmstroem et al. study, reprocessed with the OpenSwath pipeline.

    In addition to the local FDRs, the scores of features against their "true" assays are recorded in the output - in the meta value "GLM_score" of the respective feature.


    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_OpenSwathConfidenceScoring.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_OpenSwathConfidenceScoring.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPOpenSwathConfidenceScoring :
  public TOPPBase
{
public:

  /// Constructor
  TOPPOpenSwathConfidenceScoring() :
    TOPPBase("OpenSwathConfidenceScoring", "Compute confidence scores for OpenSwath results")
  {
  }

  /// Docu in base class
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (OpenSwath results)");
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerInputFile_("lib", "<file>", "", "Assay library");
    setValidFormats_("lib", ListUtils::create<String>("traML"));
    registerOutputFile_("out", "<file>", "", 
                        "Output file (results with confidence scores)");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerInputFile_("trafo", "<file>", "", "Retention time transformation",
                       false);
    setValidFormats_("trafo", ListUtils::create<String>("trafoXML"));
    registerIntOption_("decoys", "<number>", 1000, "Number of decoy assays to "
                       "select from the library for every true assay (0 for "
                       "\"all\")", false);
    setMinInt_("decoys", 0);
    registerIntOption_("transitions", "<number>", 6, "Number of transitions "
                       "per feature to consider (highest intensities first; "
                       "0 for \"all\")", false);
    setMinInt_("transitions", 0);

    registerTOPPSubsection_("GLM",
                            "Parameters of the binomial GLM");
    registerDoubleOption_("GLM:intercept", "<value>", 3.87333466, 
                          "Intercept term", false, true);
    registerDoubleOption_("GLM:delta_rt", "<value>", -0.02898629, "Coefficient "
                          "of retention time difference", false, true);
    registerDoubleOption_("GLM:dist_int", "<value>", -7.75880768,
                          "Coefficient of intensity distance", false, true);
  }

  /// Docu in base class
  ExitCodes main_(int, const char**) override
  {
    TargetedExperiment library_; // assay library
    Size n_decoys_; // number of decoys to use (per feature/true assay)
    Size n_transitions_; // number of transitions to consider
    TransformationDescription rt_trafo_; /// RT transformation to map measured RTs to assay RTs

    LOG_DEBUG << "Reading parameters..." << endl;
    String in = getStringOption_("in");
    String lib = getStringOption_("lib");
    String out = getStringOption_("out");
    String trafo = getStringOption_("trafo");
    n_decoys_ = getIntOption_("decoys");
    n_transitions_ = getIntOption_("transitions");

    LOG_DEBUG << "Loading input files..." << endl;
    FeatureMap features;
    FeatureXMLFile().load(in, features);
    TraMLFile().load(lib, library_);

    if (trafo.empty())
    {
      LOG_WARN << "Warning: You have not supplied an RT transformation file "
               << "(parameter 'trafo'). You should be sure that the retention "
               << "times of your features ('in') and library ('lib') are on "
               << "the same scale." << endl;
    }
    else
    {
      TransformationXMLFile().load(trafo, rt_trafo_);
      if (rt_trafo_.getModelType() == "none") // fit a linear model now
      {
        rt_trafo_.fitModel("linear");
      }
    }

    ConfidenceScoring scoring(test_mode_);
    scoring.setLogType(log_type_);
    scoring.initialize(library_, n_decoys_, n_transitions_, rt_trafo_);
    scoring.initializeGlm(getDoubleOption_("GLM:intercept"), getDoubleOption_("GLM:delta_rt"), getDoubleOption_("GLM:dist_int"));
    scoring.scoreMap(features);

    LOG_DEBUG << "Storing results..." << endl;
    addDataProcessing_(features, 
                       getProcessingInfo_(DataProcessing::DATA_PROCESSING));
    FeatureXMLFile().store(out, features);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPOpenSwathConfidenceScoring t;
  return t.main(argc, argv);
}

/// @endcond
