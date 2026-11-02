// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace Math; //PosteriorErrorProbabilityModel
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDPosteriorErrorProbability IDPosteriorErrorProbability

@brief  Tool to estimate the probability of peptide hits to be incorrectly assigned.

<CENTER>
<table>
<tr>
    <th ALIGN = "center"> potential predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=2> &rarr; IDPosteriorErrorProbability &rarr;</td>
    <th ALIGN = "center"> potential successor tools </td>
</tr>
<tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_CometAdapter (or other ID engines), @ref TOPP_PeptideIndexer </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate, @ref TOPP_Epifany </td>
</tr>
</table>
</CENTER>

@experimental This tool has not been tested thoroughly and might behave not as expected!

Posterior error probabilities (PEPs) are estimated separately for each search engine (and charge state, if @p split_charge is set).
The search engine scores are transformed such that higher scores are better (e.g. -log10 of e-values).
X! Tandem, Mascot, MyriMatch, OMSSA, MS-GF+, Comet, MSFragger, SpectraST, Sage, SimpleSearchEngine and ProSE have a dedicated score selection.
For other search engines, the main score of the identifications is used, unless it is already a probability (e.g. a q-value or PEP).

<B>Target-decoy based estimation (recommended)</B>

If the peptide hits are annotated as target or decoy (e.g. by @ref TOPP_PeptideIndexer after a search against a concatenated
target+decoy database), the decoys are used (parameter @p fit_algorithm:decoy_handling, default 'auto'):
- By default, the PEP is estimated non-parametrically from the local ratio of decoy and target top hits
  (target-decoy competition: each decoy top hit stands in for one incorrect target top hit). The decoy probability as a function of the score is
  estimated by isotonic regression, so the PEP decreases with the score. Only the ordering of the scores matters, so this works for any score type.
  This also works with very few decoys (or targets): even a single decoy carries information about the number of incorrect targets.
- With at least @p fit_algorithm:min_decoys but fewer than @p fit_algorithm:min_decoys_isotonic decoys, a semi-supervised mixture model
  is fitted instead, in which the decoys belong to the incorrect component.

Only top hits (one per spectrum) are used for fitting. Lower-ranked hits and decoys receive the PEP of a target hit with the same score.
Decoys therefore need to be kept in the data, i.e. run this tool before @ref TOPP_FalseDiscoveryRate (or use its option to keep decoys).
If a spectrum has several identifications (e.g. results of separate target and decoy searches of the same run), only its best top hit is used for fitting.
With @p split_charge, charge states with fewer than @p fit_algorithm:min_decoys decoys (or without targets, or whose model cannot be fitted) use the target-decoy data of all charge states.

<B>Mixture models</B>

Without target/decoy information (or with @p fit_algorithm:decoy_handling 'ignore'), a mixture model is fitted to all scores (or only the top hits, if @p top_hits_only is set):
By default an estimation is performed using the (inverse) Gumbel distribution for incorrectly assigned sequences
and a Gaussian distribution for correctly assigned sequences. The probabilities are calculated by using Bayes' law, similar to PeptideProphet.
Alternatively, a second Gaussian distribution can be used for incorrectly assigned sequences.
Without decoys, the mixture model assumes that the data contain correct assignments: for data with few or no correct assignments,
it might assign low PEPs to a few PSMs in the extreme upper tail of the scores. Use target/decoy information whenever possible.

<B>Plots</B>

In order to validate the computed probabilities an optional plot output can be generated.
There are two parameters for the plot:
The scores are plotted in the form of bins. Each bin represents a set of scores in a range of '(highest_score - smallest_score) / number_of_bins' (if all scores have positive values).
The midpoint of the bin is the mean of the scores it represents.
The parameter 'out_plot' should be used to give the plot a unique name. Two files are created. One with the binned scores and one with all steps of the estimation.
If target/decoy information is available, an additional plot will be generated with the score distributions of the target and decoy top hits and the estimated PEPs ('out_plot' must not be empty).
The plots are saved as a Gnuplot file. An attempt is made to call Gnuplot, which will create a PDF file containing all steps of the estimation. If this fails, the user has to run Gnuplot manually - or adjust the PATH environment such that Gnuplot can be found and retry.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDPosteriorErrorProbability.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDPosteriorErrorProbability.html

For the parameters of the algorithm section see the algorithms documentation: @n
@ref OpenMS::Math::PosteriorErrorProbabilityModel "fit_algorithm" @n

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDPosteriorErrorProbability :
  public TOPPBase
{
public:
  TOPPIDPosteriorErrorProbability() :
    TOPPBase("IDPosteriorErrorProbability", "Estimates probabilities for incorrectly assigned peptide sequences and a set of search engine scores using target-decoy information or a mixture model.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<std::string>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<std::string>("idXML"));
    registerOutputFile_("out_plot", "<file>", "", "txt file (if gnuplot is available, a corresponding PDF will be created as well.)", false);
    setValidFormats_("out_plot", ListUtils::create<std::string>("txt"));

    registerFlag_("split_charge", "The search engine scores are split by charge if this flag is set. Thus, for each charge state a new model will be computed.");
    registerFlag_("top_hits_only", "If set only the top hits of every PeptideIdentification will be used for fitting the mixture model without target/decoy information. (Target-decoy based estimates always use the top hits only.)");
    registerDoubleOption_("fdr_for_targets_smaller", "<value>", 0.05, "Deprecated and ignored. Target and decoy hits are determined by their target/decoy annotation.", false, true);
    registerFlag_("ignore_bad_data", "If set errors will be written but ignored. Useful for pipelines with many datasets where only a few are bad, but the pipeline should run through.");
    registerFlag_("prob_correct", "If set scores will be calculated as '1 - ErrorProbabilities' and can be interpreted as probabilities for correct identifications.");
    registerSubsection_("fit_algorithm", "Algorithm parameter subsection");
    addEmptyLine_();
  }

  //there is only one parameter at the moment
  Param getSubsectionDefaults_(const std::string& /*section*/) const override
  {
    Param p = PosteriorErrorProbabilityModel().getParameters();
    if (p.exists("out_plot"))
    { // hide from user -- we have a top-level param for that
      p.remove("out_plot");
    }
    else
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "INTERNAL ERROR: Param 'out_plot' was removed from fit-algorithm. Please update param handling internally!");
    }
    return p;
  }


  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    std::string inputfile_name = getStringOption_("in");
    std::string outputfile_name = getStringOption_("out");
    Param fit_algorithm = getParam_().copy("fit_algorithm:", true);
    fit_algorithm.setValue("out_plot", getStringOption_("out_plot")); // re-assemble full param (was moved to top-level)
    bool split_charge = getFlag_("split_charge");
    bool top_hits_only = getFlag_("top_hits_only");
    bool ignore_bad_data = getFlag_("ignore_bad_data");
    bool prob_correct = getFlag_("prob_correct");
    std::string outlier_handling = fit_algorithm.getValue("outlier_handling").toString();
    const std::string decoy_handling = fit_algorithm.getValue("decoy_handling").toString();
    const Size min_decoys = static_cast<Size>(int(fit_algorithm.getValue("min_decoys")));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    FileHandler file;
    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;
    file.loadIdentifications(inputfile_name, protein_ids, peptide_ids, {FileTypes::IDXML});
    PosteriorErrorProbabilityModel PEP_model;
    PEP_model.setParameters(fit_algorithm);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // map identifier "engine,charge" (if split_charge==true) or "engine"
    // to three extracted score vectors: [0] scores of all considered hits (for the unsupervised model),
    // [1]/[2] scores of target/decoy top hits (empty if there are no decoys).
    map<std::string, vector<vector<double> > > all_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(
     protein_ids,
     peptide_ids,
     split_charge,
     top_hits_only);

    // target/decoy top hits of all charge states (fallback for charge states with too few decoys)
    map<std::string, vector<vector<double> > > pooled_scores;
    if (split_charge && decoy_handling != "ignore")
    {
      pooled_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(protein_ids, peptide_ids, false, true);
    }

    // something to fit: enough scores for a mixture model, or target and decoy top hits (even a single one of each is informative)
    auto has_data = [](const std::pair<const std::string, vector<vector<double>>>& entry)
    {
      return entry.second[0].size() > 2 || (!entry.second[1].empty() && !entry.second[2].empty());
    };
    const bool any_data = std::any_of(all_scores.begin(), all_scores.end(), has_data) || std::any_of(pooled_scores.begin(), pooled_scores.end(), has_data);
    if (!any_data)
    {
      writeLogWarn_("No data collected. Check whether search engine is supported.");
      if (!ignore_bad_data)
      {
        return INPUT_FILE_EMPTY;
      }
    }

    std::string out_plot = fit_algorithm.getValue("out_plot").toString();
    StringUtils::trim(out_plot);

    std::set<std::pair<Size, Size>> updated_hits; // (index of identification, index of hit)
    for (auto & score : all_scores)
    {
      vector<std::string> engine_info;
      StringUtils::split(score.first, ',', engine_info);
      std::string engine = engine_info[0];
      Int charge = (engine_info.size() == 2) ? StringUtils::toInt32(engine_info[1]) : -1;
      const std::string model_name = engine + (split_charge ? " (charge " + StringUtils::toStr(charge) + ")" : std::string());

      if (split_charge)
      {
        // only adapt plot output if plot is requested (this badly violates the output rules and needs to change!)
        // one way to fix this: plot charges into a single file (no renaming of output file needed) - but this requires major code restructuring
        if (!out_plot.empty()) fit_algorithm.setValue("out_plot", out_plot + "_charge_" + StringUtils::toStr(charge));
        PEP_model.setParameters(fit_algorithm);
      }

      vector<double>* targets = &score.second[1];
      vector<double>* decoys = &score.second[2];

      // With split_charge, charge states with fewer than min_decoys decoys (or without targets) use the target-decoy data of all
      // charge states, if other charge states contribute top hits (targets or decoys). (Which estimator is used is decided by
      // the model, see decoy_handling.)
      auto pooled_it = pooled_scores.find(engine);
      const bool pooled_available = split_charge && pooled_it != pooled_scores.end() && !pooled_it->second[1].empty() && !pooled_it->second[2].empty();
      const bool pooled_has_more = pooled_available && pooled_it->second[1].size() + pooled_it->second[2].size() > targets->size() + decoys->size();
      // decoys present: target-decoy based estimate (without targets, nothing can be estimated; never model decoys without supervision)
      const bool use_decoys = (decoy_handling != "ignore") && (!decoys->empty() || pooled_available);
      bool use_pooled = false;
      auto switch_to_pooled = [&](const std::string& reason)
      {
        writeLogInfo_(reason + " Using the target-decoy model of all charge states for " + model_name + ".");
        targets = &pooled_it->second[1];
        decoys = &pooled_it->second[2];
        use_pooled = true;
      };
      if (use_decoys && pooled_has_more && (decoys->size() < min_decoys || targets->empty()))
      {
        switch_to_pooled("Only " + StringUtils::toStr(targets->size()) + " target and " + StringUtils::toStr(decoys->size()) + " decoy top hits for " + model_name + ".");
      }
      else if (!use_decoys && score.second[0].size() <= 2)
      {
        continue; // too few scores to fit a model (reported below)
      }

      // fit to score vector
      //TODO choose outlier handling based on search engine? If not set by user?
      //XTandem is prone to accumulation at min values/censoring
      //OMSSA is prone to outliers
      bool return_value = use_decoys ? PEP_model.fitTargetDecoy(*targets, *decoys, outlier_handling)
                                     : PEP_model.fit(score.second[0], outlier_handling);
      if (!return_value && use_decoys && pooled_available && !use_pooled)
      {
        switch_to_pooled("Unable to fit the data of " + model_name + ".");
        return_value = PEP_model.fitTargetDecoy(*targets, *decoys, outlier_handling);
      }

      if (!return_value)
      {
        writeLogWarn_("Unable to fit data. Algorithm did not run through for the following search engine: " + model_name);
        if (!ignore_bad_data)
        {
          return UNEXPECTED_RESULT;
        }
        continue;
      }

      std::string summary = "PEP model for " + model_name + ": " + PosteriorErrorProbabilityModel::modelTypeToString(PEP_model.getModelType());
      if (use_decoys)
      {
        summary += " (" + StringUtils::toStr(targets->size()) + " target and " + StringUtils::toStr(decoys->size()) + " decoy top hits)";
      }
      writeLogInfo_(summary);

      // plot target and decoy score distributions together with the estimated model
      if (!out_plot.empty() && use_decoys)
      {
        vector<double> target_copy(*targets), decoy_copy(*decoys);
        PEP_model.plotTargetDecoyEstimation(target_copy, decoy_copy);
      }

      bool unable_to_fit_data(true), data_might_not_be_well_fit(true);
      PosteriorErrorProbabilityModel::updateScores(
       PEP_model,
       engine,
       charge,
       prob_correct,
       split_charge,
       protein_ids,
       peptide_ids,
       unable_to_fit_data,
       data_might_not_be_well_fit,
       updated_hits);

      if (unable_to_fit_data)
      {
        writeLogWarn_(std::string("Unable to fit data for search engine: ") + model_name);
        if (!ignore_bad_data)
        {
          return UNEXPECTED_RESULT;
        }
      }
      else if (data_might_not_be_well_fit || (use_decoys && PEP_model.hasCalibrationWarning()))
      {
        writeLogWarn_(std::string("Data might not be well fitted for search engine: ") + model_name);
      }
    }
    // Unfortunately this cannot go into the algorithm since
    // you would overwrite some score types before they are extracted when you
    // do split_charge
    // Score types are kept consistent within each run: all identifications of runs with at least one updated hit (and of all runs
    // with -ignore_bad_data, so that pipelines can continue) are relabeled, and their hits without a PEP get a PEP of 1.
    std::set<std::string> updated_runs;
    for (const auto& updated : updated_hits)
    {
      updated_runs.insert(peptide_ids[updated.first].getIdentifier());
    }
    Size n_set_to_one = 0, n_kept = 0;
    for (Size index = 0; index < peptide_ids.size(); ++index)
    {
      PeptideIdentification& pep = peptide_ids[index];
      if (!ignore_bad_data && updated_runs.count(pep.getIdentifier()) == 0)
      {
        n_kept += pep.getHits().size();
        continue;
      }
      const std::string original_score_type = pep.getScoreType() + "_score";
      std::vector<PeptideHit>& hits = pep.getHits();
      for (Size hit_index = 0; hit_index < hits.size(); ++hit_index)
      {
        if (updated_hits.count({index, hit_index}) == 0)
        {
          hits[hit_index].setMetaValue(original_score_type, hits[hit_index].getScore());
          hits[hit_index].setScore(prob_correct ? 0.0 : 1.0);
          ++n_set_to_one;
        }
      }
      if (prob_correct)
      {
        pep.setScoreType("Posterior Probability");
        pep.setHigherScoreBetter(true);
      }
      else
      {
        pep.setScoreType("Posterior Error Probability");
        pep.setHigherScoreBetter(false);
      }
    }
    if (n_set_to_one > 0)
    {
      writeLogWarn_("Warning: " + StringUtils::toStr(n_set_to_one) + " hit(s) could not be assigned a PEP by a model (e.g. unusable score, or too few PSMs)"
                    " and were assigned a PEP of 1 (the original score is kept as meta value).");
    }
    if (n_kept > 0)
    {
      writeLogWarn_("Warning: " + StringUtils::toStr(n_kept) + " hit(s) of search runs without a PEP model keep their original scores and score types"
                    " (use -ignore_bad_data to assign them a PEP of 1).");
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    file.storeIdentifications(outputfile_name, protein_ids, peptide_ids, {FileTypes::IDXML});
    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPIDPosteriorErrorProbability tool;
  return tool.main(argc, argv);
}

/// @endcond
