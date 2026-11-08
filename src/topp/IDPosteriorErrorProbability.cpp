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
#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

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

Posterior error probabilities (PEPs) are estimated separately for each search engine (and charge state, if @p split_charge is set,
and MS run, if @p split_run is set).
The search engine scores are transformed such that higher scores are better (e.g. -log10 of e-values).
X! Tandem, Mascot, MyriMatch, OMSSA, MS-GF+, Comet, MSFragger, SpectraST, Sage, SimpleSearchEngine and ProSE have a dedicated score selection.
For other search engines, the main score of the identifications is used, unless it is already a probability (e.g. a q-value or PEP).

<B>Target-decoy based estimation (recommended)</B>

If the peptide hits are annotated as target or decoy (e.g. by @ref TOPP_PeptideIndexer after a search against a concatenated
target+decoy database), the decoys are used (parameter @p fit_algorithm:decoy_handling, default 'auto'):
the PEP is estimated non-parametrically from the local ratio of decoy and target top hits
(target-decoy competition: each decoy top hit stands in for one incorrect target top hit). The decoy probability as a function of the score is
estimated by isotonic regression, so the PEP decreases with the score. The estimate does not assume any shape of the score distributions
and only depends on the ordering of the scores (up to the interpolation within and between blocks), so it works for any score type.
This also works with very few decoys (or targets): even a single decoy carries information about the number of incorrect targets.

Only top hits (one per spectrum) are used for fitting. Lower-ranked target hits receive the PEP of a target top hit with the same score.
Decoys receive the PEP that they would get as targets (label swap), which makes decoys and incorrect targets exchangeable, as target-decoy
FDR estimates on the PEPs (e.g. @ref TOPP_FalseDiscoveryRate, or the protein FDR after protein inference) require.
Decoys therefore need to be kept in the data, i.e. run this tool before @ref TOPP_FalseDiscoveryRate (or use its option to keep decoys).
If a spectrum has several identifications (e.g. results of separate target and decoy searches of the same run), only its best top hit is used for fitting.
With @p split_charge, charge states with fewer than @p fit_algorithm:min_decoys decoys (or without targets, or whose model cannot be fitted,
or whose targets do not score better than their decoys) use the target-decoy data of all charge states.

<B>Several runs</B>

By default, all search runs of the input file that stem from the same search engine are modeled together (one PEP curve per search engine),
with target-decoy competition per spectrum. Score types must be comparable across the runs (same search engine and settings;
e-values are preferable to raw scores). One model for several runs is more accurate than a model per run if the runs are small
(e.g. fractions, or runs with few identifications): with few decoys, the PEPs of a run depend strongly on the chance positions of its decoys,
and the few incorrect targets that score above the best decoy of their run get low PEPs.
With @p split_run, the PEPs of each MS run are estimated separately, e.g. if the score distributions of the runs differ strongly.
An MS run is identified by its spectra file, so separate target and decoy searches of the same spectra file still compete per spectrum
(identifications whose spectra file is unknown are grouped by their identification run). With @p out_plot, the plot files of the MS runs
are numbered in the order of the runs in the input ('_run_1', '_run_2', ...).

<B>PSM or peptide level</B>

By default (@p fit_level 'psm'), the model is fitted to the PSMs. With @p fit_level 'peptide', it is fitted to the best PSM of each peptide
(modified sequence; with @p split_charge: of each peptide and charge state): for the target-decoy estimate, the best top hit of each peptide
after the competition per spectrum, as Percolator and mokapot do at the peptide level. Every PSM gets the PEP of this model at its score,
so the best PSM of each peptide carries the PEP of the peptide.
Use the peptide level for protein inference that takes the best PSM of each peptide (e.g. @ref TOPP_Epifany): the best of several PSMs is
selected by its score, which PSM-level PEPs do not account for, while the decoys go through the same selection at the peptide level.
For PSM-level uses of the PEPs (e.g. filtering PSMs by PEP, or the PEPs of all PSMs of a peptide), use the PSM level.
The parameter @p fit_algorithm:min_decoys counts decoy peptides at the peptide level, and @p fit_algorithm:min_mixture_scores the peptides
(their best scores) for the mixture model.

<B>Mixture models</B>

Without target/decoy information (or with @p fit_algorithm:decoy_handling 'ignore'), a mixture model is fitted to all scores (or only the top hits, if @p top_hits_only is set):
By default an estimation is performed using the (inverse) Gumbel distribution for incorrectly assigned sequences
and a Gaussian distribution for correctly assigned sequences. The probabilities are calculated by using Bayes' law, similar to PeptideProphet.
Alternatively, a second Gaussian distribution can be used for incorrectly assigned sequences.
Without decoys, the mixture model assumes that the data contain correct assignments: for data with few or no correct assignments,
it might assign low PEPs to a few PSMs in the extreme upper tail of the scores. On real data (Comet e-values of a yeast/UPS1 dataset,
validated with an entrapment database), the PEPs of the mixture model were about 1.5 times too high in total, and much too high for
confident PSMs, while the target-decoy estimate was within a few percent. Use target/decoy information whenever possible.
The mixture model needs at least @p fit_algorithm:min_mixture_scores scores (default 35, as PeptideProphet), also after outlier handling.
With @p split_charge, a charge state with fewer scores, or whose mixture model cannot be fitted, uses the mixture model of all charge states.
If there are still too few scores, the tool fails (with @p ignore_bad_data: the hits get a PEP of 1, the original score is kept as meta value).

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
    registerFlag_("split_run", "If set, the PEPs of each MS run (spectra file) are estimated separately. By default, one model is fitted to all runs of the input file (per search engine).");
    registerStringOption_("fit_level", "<level>", "psm", "Fit the model to the PSMs ('psm'), or to the best PSM of each peptide ('peptide'; for protein inference that uses the best PSM per peptide, e.g. Epifany). "
                          "Every PSM gets the PEP of the model at its score.", false);
    setValidStrings_("fit_level", ListUtils::create<std::string>("psm,peptide"));
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


  /// number of finite scores (the mixture model ignores non-finite ones)
  static Size finiteCount(const vector<double>& scores)
  {
    return static_cast<Size>(std::count_if(scores.begin(), scores.end(), [](double score) { return std::isfinite(score); }));
  }

  /**
    @brief groups the identifications by MS run

    An MS run is identified by its spectra file, so that identification runs of the same spectra file (e.g. separate target and decoy
    searches) belong to the same MS run. Identifications whose spectra file is unknown are grouped by their identification run.

    @return (name of the MS run, indices into @p peptide_ids) in the order of the first identification of each MS run
  */
  static std::vector<std::pair<std::string, std::vector<Size>>> groupByMSRun(const vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids)
  {
    IdentifierMSRunMapper ms_runs;
    try
    {
      ms_runs.create(protein_ids);
    }
    catch (Exception::InvalidValue&)
    {
      // several identification runs of the same spectra files: expected here; the mapping from runs to files is complete nevertheless
    }
    std::vector<std::pair<std::string, std::vector<Size>>> runs;
    std::map<std::string, Size> run_index; // name -> index into runs
    for (Size index = 0; index < peptide_ids.size(); ++index)
    {
      const PeptideIdentification& pep = peptide_ids[index];
      std::string name;
      try
      {
        ms_runs.validateMergeIndex(pep, index); // throws if the run has several spectra files, but the one of this identification is unknown
        const std::string ms_file = ms_runs.getPrimaryMSRunPath(pep);
        if (!ms_file.empty())
        {
          name = "spectra file '" + ms_file + "'";
        }
      }
      catch (Exception::BaseException&)
      {
        // spectra file unknown
      }
      if (name.empty())
      {
        name = "identification run '" + pep.getIdentifier() + "'";
      }
      auto [it, inserted] = run_index.emplace(name, runs.size());
      if (inserted)
      {
        runs.emplace_back(name, std::vector<Size>());
      }
      runs[it->second].second.push_back(index);
    }
    return runs;
  }

  /**
    @brief estimates the PEPs of the identifications in @p peptide_ids: one model per search engine (and charge state, with split_charge)

    @param[in,out] protein_ids the identification runs
    @param[in,out] peptide_ids the identifications; the hits that get a PEP are updated
    @param[in] run_name name of the MS run in log messages (empty: all runs of the input)
    @param[in] plot_suffix appended to the name of the plot files
    @param[in,out] any_data set to true if there are scores to fit
    @param[in,out] updated_hits the hits that got a PEP (index into @p peptide_ids, index of the hit) are added
    @return EXECUTION_OK, or the exit code of a failure (only without ignore_bad_data)
  */
  ExitCodes estimatePEPs_(vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids, const std::string& run_name,
                          const std::string& plot_suffix, bool& any_data, std::set<std::pair<Size, Size>>& updated_hits)
  {
    Param fit_algorithm = fit_algorithm_;
    const std::string out_plot = out_plot_.empty() ? std::string() : out_plot_ + plot_suffix;
    fit_algorithm.setValue("out_plot", out_plot);
    PosteriorErrorProbabilityModel PEP_model;
    PEP_model.setParameters(fit_algorithm);
    const std::string unit = peptide_level_ ? "peptides" : "top hits"; // what the target-decoy estimate counts
    const std::string run_suffix = run_name.empty() ? std::string() : " of " + run_name;

    // map identifier "engine,charge" (if split_charge==true) or "engine"
    // to three extracted score vectors: [0] scores of all considered hits (for the unsupervised model),
    // [1]/[2] scores of target/decoy top hits (empty if there are no decoys); with fit_level 'peptide', one score per peptide.
    map<std::string, vector<vector<double> > > all_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(
     protein_ids,
     peptide_ids,
     split_charge_,
     top_hits_only_,
     peptide_level_);

    // scores of all charge states (fallback for charge states with too few decoys, or too few scores for the mixture model)
    map<std::string, vector<vector<double> > > pooled_scores;
    if (split_charge_)
    {
      pooled_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(protein_ids, peptide_ids, false, top_hits_only_, peptide_level_);
    }

    // something to fit: scores for a mixture model, or target and decoy top hits (even a single one of each is informative)
    auto has_data = [](const std::pair<const std::string, vector<vector<double>>>& entry)
    {
      return !entry.second[0].empty() || (!entry.second[1].empty() && !entry.second[2].empty());
    };
    any_data = any_data || std::any_of(all_scores.begin(), all_scores.end(), has_data) || std::any_of(pooled_scores.begin(), pooled_scores.end(), has_data);

    for (auto & score : all_scores)
    {
      vector<std::string> engine_info;
      StringUtils::split(score.first, ',', engine_info);
      std::string engine = engine_info[0];
      Int charge = (engine_info.size() == 2) ? StringUtils::toInt32(engine_info[1]) : -1;
      const std::string model_name = engine + (split_charge_ ? " (charge " + StringUtils::toStr(charge) + ")" : std::string()) + run_suffix;

      if (split_charge_)
      {
        // only adapt plot output if plot is requested (this badly violates the output rules and needs to change!)
        // one way to fix this: plot charges into a single file (no renaming of output file needed) - but this requires major code restructuring
        if (!out_plot.empty()) fit_algorithm.setValue("out_plot", out_plot + "_charge_" + StringUtils::toStr(charge));
        PEP_model.setParameters(fit_algorithm);
      }

      vector<double>* targets = &score.second[1];
      vector<double>* decoys = &score.second[2];
      vector<double>* mixture_scores = &score.second[0];

      // With split_charge, charge states with fewer than min_decoys decoys (or without targets) use the target-decoy data of all
      // charge states, if other charge states contribute top hits (targets or decoys).
      auto pooled_it = pooled_scores.find(engine);
      const bool pooled_available = split_charge_ && pooled_it != pooled_scores.end() && !pooled_it->second[1].empty() && !pooled_it->second[2].empty();
      const bool pooled_has_more = pooled_available && pooled_it->second[1].size() + pooled_it->second[2].size() > targets->size() + decoys->size();
      // decoys present: target-decoy based estimate (without targets, nothing can be estimated; never model decoys without supervision)
      const bool use_decoys = (decoy_handling_ != "ignore") && (!decoys->empty() || pooled_available);
      bool use_pooled = false;
      auto switch_to_pooled = [&](const std::string& reason)
      {
        writeLogInfo_(reason + " Using the target-decoy model of all charge states for " + model_name + ".");
        targets = &pooled_it->second[1];
        decoys = &pooled_it->second[2];
        use_pooled = true;
      };
      // Without target-decoy data, charge states with too few scores for their own mixture model (or whose mixture model cannot be
      // fitted) use the mixture model of all charge states, if other charge states contribute scores.
      const bool pooled_mixture_available = split_charge_ && pooled_it != pooled_scores.end() && finiteCount(pooled_it->second[0]) > finiteCount(score.second[0]);
      auto switch_to_pooled_mixture = [&](const std::string& reason)
      {
        writeLogInfo_(reason + " Using the mixture model of all charge states (" + StringUtils::toStr(finiteCount(pooled_it->second[0])) + " scores) for " + model_name + ".");
        mixture_scores = &pooled_it->second[0];
      };
      if (use_decoys && pooled_has_more && (decoys->size() < min_decoys_ || targets->empty()))
      {
        switch_to_pooled("Only " + StringUtils::toStr(targets->size()) + " target and " + StringUtils::toStr(decoys->size()) + " decoy " + unit + " for " + model_name + ".");
      }
      else if (!use_decoys && finiteCount(*mixture_scores) < min_mixture_scores_)
      {
        if (pooled_mixture_available && finiteCount(pooled_it->second[0]) >= min_mixture_scores_)
        {
          switch_to_pooled_mixture("Only " + StringUtils::toStr(finiteCount(*mixture_scores)) + " scores for " + model_name + ".");
        }
        else
        {
          const bool with_pooled = split_charge_ && pooled_it != pooled_scores.end();
          const bool has_targets = !score.second[1].empty() || (with_pooled && !pooled_it->second[1].empty());
          const bool has_decoys = !score.second[2].empty() || (with_pooled && !pooled_it->second[2].empty());
          std::string advice;
          if (!has_decoys)
          {
            advice = "Search against a target+decoy database (DecoyDatabase, PeptideIndexer) to use the target-decoy estimate, which also works with few PSMs.";
          }
          else if (!has_targets)
          {
            advice = "The data contain decoy but no target top hits: annotate the target PSMs (meta value 'target_decoy', e.g. with PeptideIndexer)"
                     + std::string(decoy_handling_ == "ignore" ? " and use fit_algorithm:decoy_handling 'auto'" : "")
                     + " to use the target-decoy estimate, which also works with few PSMs.";
          }
          else // only reachable with decoy_handling 'ignore' (otherwise the target-decoy estimate is used)
          {
            advice = "The data contain decoys: with fit_algorithm:decoy_handling 'auto', the target-decoy estimate is used, which also works with few PSMs.";
          }
          writeLogWarn_("Only " + StringUtils::toStr(finiteCount(*mixture_scores)) + " scores for " + model_name
                        + (with_pooled ? " (all charge states: " + StringUtils::toStr(finiteCount(pooled_it->second[0])) + ")" : std::string())
                        + ", too few to fit the mixture model (fit_algorithm:min_mixture_scores: " + StringUtils::toStr(min_mixture_scores_) + "). " + advice);
          if (!ignore_bad_data_)
          {
            return UNEXPECTED_RESULT;
          }
          continue;
        }
      }

      // fit to score vector
      //TODO choose outlier handling based on search engine? If not set by user?
      //XTandem is prone to accumulation at min values/censoring
      //OMSSA is prone to outliers
      bool return_value = use_decoys ? PEP_model.fitTargetDecoy(*targets, *decoys, outlier_handling_)
                                     : PEP_model.fit(*mixture_scores, outlier_handling_);
      if (!return_value && use_decoys && pooled_available && !use_pooled)
      {
        switch_to_pooled("Unable to fit the data of " + model_name + ".");
        return_value = PEP_model.fitTargetDecoy(*targets, *decoys, outlier_handling_);
      }
      else if (return_value && use_decoys && pooled_has_more && !use_pooled && PEP_model.getTargetDecoyAUC() >= 0.0 && PEP_model.getTargetDecoyAUC() <= 0.5)
      {
        // a charge state whose targets do not score better than its decoys would get (almost) constant PEPs, for decoys as well
        switch_to_pooled("Target PSMs of " + model_name + " do not score better than decoy PSMs (AUC = " + StringUtils::toStr(PEP_model.getTargetDecoyAUC()) + ").");
        return_value = PEP_model.fitTargetDecoy(*targets, *decoys, outlier_handling_);
      }
      else if (!return_value && !use_decoys && pooled_mixture_available && mixture_scores == &score.second[0])
      {
        switch_to_pooled_mixture("Unable to fit the mixture model of " + model_name + ".");
        return_value = PEP_model.fit(*mixture_scores, outlier_handling_);
      }

      if (!return_value)
      {
        writeLogWarn_("Unable to fit data. Algorithm did not run through for the following search engine: " + model_name);
        if (!ignore_bad_data_)
        {
          return UNEXPECTED_RESULT;
        }
        continue;
      }

      std::string summary = "PEP model for " + model_name + ": " + PosteriorErrorProbabilityModel::modelTypeToString(PEP_model.getModelType());
      if (use_decoys)
      {
        summary += " (" + StringUtils::toStr(targets->size()) + " target and " + StringUtils::toStr(decoys->size()) + " decoy " + unit + ")";
      }
      else
      {
        summary += " (" + StringUtils::toStr(finiteCount(*mixture_scores)) + (peptide_level_ ? " peptide" : "") + " scores" + (mixture_scores != &score.second[0] ? " of all charge states" : "") + ")";
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
       prob_correct_,
       split_charge_,
       protein_ids,
       peptide_ids,
       unable_to_fit_data,
       data_might_not_be_well_fit,
       updated_hits);

      if (unable_to_fit_data)
      {
        writeLogWarn_(std::string("Unable to fit data for search engine: ") + model_name);
        if (!ignore_bad_data_)
        {
          return UNEXPECTED_RESULT;
        }
      }
      else if (data_might_not_be_well_fit || (use_decoys && PEP_model.hasCalibrationWarning()))
      {
        writeLogWarn_(std::string("Data might not be well fitted for search engine: ") + model_name);
      }
    }
    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    std::string inputfile_name = getStringOption_("in");
    std::string outputfile_name = getStringOption_("out");
    fit_algorithm_ = getParam_().copy("fit_algorithm:", true);
    out_plot_ = getStringOption_("out_plot");
    StringUtils::trim(out_plot_);
    split_charge_ = getFlag_("split_charge");
    const bool split_run = getFlag_("split_run");
    peptide_level_ = (getStringOption_("fit_level") == "peptide");
    top_hits_only_ = getFlag_("top_hits_only");
    ignore_bad_data_ = getFlag_("ignore_bad_data");
    prob_correct_ = getFlag_("prob_correct");
    outlier_handling_ = fit_algorithm_.getValue("outlier_handling").toString();
    decoy_handling_ = fit_algorithm_.getValue("decoy_handling").toString();
    min_decoys_ = static_cast<Size>(int(fit_algorithm_.getValue("min_decoys")));
    min_mixture_scores_ = static_cast<Size>(int(fit_algorithm_.getValue("min_mixture_scores")));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    FileHandler file;
    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;
    file.loadIdentifications(inputfile_name, protein_ids, peptide_ids, {FileTypes::IDXML});

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    bool any_data = false;
    std::set<std::pair<Size, Size>> updated_hits; // (index of identification, index of hit)
    if (split_run)
    {
      const std::vector<std::pair<std::string, std::vector<Size>>> ms_runs = groupByMSRun(protein_ids, peptide_ids);
      writeLogInfo_("Estimating the PEPs of " + StringUtils::toStr(ms_runs.size()) + " MS run(s) separately.");
      for (Size run = 0; run < ms_runs.size(); ++run)
      {
        const std::string& run_name = ms_runs[run].first;
        const std::vector<Size>& members = ms_runs[run].second;
        PeptideIdentificationList run_ids;
        run_ids.reserve(members.size());
        for (Size index : members)
        {
          run_ids.push_back(peptide_ids[index]);
        }
        // plots: one set of files per MS run, numbered in the order of the MS runs in the input
        const std::string plot_suffix = ms_runs.size() > 1 ? "_run_" + StringUtils::toStr(run + 1) : std::string();
        if (!out_plot_.empty() && !plot_suffix.empty())
        {
          writeLogInfo_("Plots of " + run_name + ": " + out_plot_ + plot_suffix);
        }
        std::set<std::pair<Size, Size>> run_updated_hits;
        const ExitCodes result = estimatePEPs_(protein_ids, run_ids, run_name, plot_suffix, any_data, run_updated_hits);
        if (result != EXECUTION_OK)
        {
          return result;
        }
        for (Size k = 0; k < members.size(); ++k)
        {
          peptide_ids[members[k]] = std::move(run_ids[k]);
        }
        for (const auto& updated : run_updated_hits)
        {
          updated_hits.emplace(members[updated.first], updated.second);
        }
      }
    }
    else
    {
      const ExitCodes result = estimatePEPs_(protein_ids, peptide_ids, std::string(), std::string(), any_data, updated_hits);
      if (result != EXECUTION_OK)
      {
        return result;
      }
    }
    if (!any_data)
    {
      writeLogWarn_("No data collected. Check whether search engine is supported.");
      if (!ignore_bad_data_)
      {
        return INPUT_FILE_EMPTY;
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
      if (!ignore_bad_data_ && updated_runs.count(pep.getIdentifier()) == 0)
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
          hits[hit_index].setScore(prob_correct_ ? 0.0 : 1.0);
          ++n_set_to_one;
        }
      }
      if (prob_correct_)
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
    // The search engine stays the original engine (the PEP model relies on it); the derived score records its own
    // producer in the runs it was assigned in, so runs from different engines share one definition.
    const std::string score_type = prob_correct_ ? "Posterior Probability" : "Posterior Error Probability";
    const std::string version = getProcessingInfo_(std::set<DataProcessing::ProcessingAction>{}).getSoftware().getVersion();
    for (auto& run : protein_ids)
    {
      if (ignore_bad_data_ || updated_runs.count(run.getIdentifier()) > 0) run.setScoreSoftware(score_type, toolName_(), version);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    file.storeIdentifications(outputfile_name, protein_ids, peptide_ids, {FileTypes::IDXML});
    return EXECUTION_OK;
  }

  // options (set in main_)
  Param fit_algorithm_; ///< parameters of the PEP model
  std::string out_plot_; ///< base name of the plot files (empty: no plots)
  bool split_charge_ = false;
  bool peptide_level_ = false;
  bool top_hits_only_ = false;
  bool ignore_bad_data_ = false;
  bool prob_correct_ = false;
  std::string outlier_handling_;
  std::string decoy_handling_;
  Size min_decoys_ = 0;
  Size min_mixture_scores_ = 0;
};


int main(int argc, const char** argv)
{
  TOPPIDPosteriorErrorProbability tool;
  return tool.main(argc, argv);
}

/// @endcond
