// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h> 

#include <iostream> // for "cout"
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>

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
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; OpenSwathConfidenceScoring &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathFeatureXMLToTSV </td>
        </tr>
    </table>
</CENTER>

This is an implementation of the SRM scoring algorithm described in:

Malmstroem, L.; Malmstroem, J.; Selevsek, N.; Rosenberger, G. & Aebersold, R.:\n
<a href="https://doi.org/10.1021/pr200844d">Automated workflow for large-scale selected reaction monitoring experiments.</a>\n
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

    OPENMS_LOG_DEBUG << "Reading parameters..." << endl;
    String in = getStringOption_("in");
    String lib = getStringOption_("lib");
    String out = getStringOption_("out");
    String trafo = getStringOption_("trafo");
    n_decoys_ = getIntOption_("decoys");
    n_transitions_ = getIntOption_("transitions");

    OPENMS_LOG_DEBUG << "Loading input files..." << endl;
    FeatureMap features;
    FileHandler().loadFeatures(in, features, {FileTypes::FEATUREXML});
    FileHandler().loadTransitions(lib, library_, {FileTypes::TRAML});

    if (trafo.empty())
    {
      OPENMS_LOG_WARN << "Warning: You have not supplied an RT transformation file "
               << "(parameter 'trafo'). You should be sure that the retention "
               << "times of your features ('in') and library ('lib') are on "
               << "the same scale." << endl;
    }
    else
    {
      FileHandler().loadTransformations(trafo, rt_trafo_, true, {FileTypes::TRANSFORMATIONXML});
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

    OPENMS_LOG_DEBUG << "Storing results..." << endl;
    addDataProcessing_(features, 
                       getProcessingInfo_(DataProcessing::DATA_PROCESSING));
    FileHandler().storeFeatures(out, features, {FileTypes::FEATUREXML});

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPOpenSwathConfidenceScoring t;
  return t.main(argc, argv);
}

/// @endcond
