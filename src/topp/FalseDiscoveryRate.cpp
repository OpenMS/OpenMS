// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/ClassFDR.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_FalseDiscoveryRate FalseDiscoveryRate

@brief Tool to estimate the false discovery rate on peptide and protein level
<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; FalseDiscoveryRate &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_CometAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDFilter </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer </td>
        </tr>
    </table>
</CENTER>

This TOPP tool calculates the false discovery rate (FDR) for results of target-decoy searches. The FDR calculation can be performed for proteins and/or for peptides (more exactly, peptide spectrum matches).

The false discovery rate is defined as the number of false discoveries (decoy hits) divided by the number of false and correct discoveries (both target and decoy hits) with a score better than a given threshold.

@ref TOPP_PeptideIndexer must be applied to the search results (idXML file) to index the data and to annotate peptide and protein hits with their target/decoy status.

@note When no decoy hits were found you will get a warning like this:<br>
"FalseDiscoveryRate: #decoy sequences is zero! Setting all target sequences to q-value/FDR 0!"<br>
This should be a serious concern, since it indicates a possible problem with the target/decoy annotation step (@ref TOPP_PeptideIndexer), e.g. due to a misconfigured database.

@note FalseDiscoveryRate only annotates peptides and proteins with their FDR. By setting FDR:PSM or FDR:protein the maximum q-value (e.g., 0.05 corresponds to an FDR of 5%) can be controlled on the PSM and protein level.
Alternatively, FDR filtering can be performed in the @ref TOPP_IDFilter tool by setting score:pep and score:prot to the maximum q-value. After potential filtering, associations are
automatically updated and unreferenced proteins/peptides removed based on the advanced cleanup parameters.

<B>FDR within classes of PSMs</B>

With @p -PSM_class, PSMs are split into classes (by their modifications, or by the value of a peptide hit meta value) and
each target PSM also gets a q-value within its class (meta value "class_q-value"). Its PSM q-value becomes the larger of
the two, so @p FDR:PSM holds within every class, not only overall. Without this, a small class (e.g. PSMs with a rare
modification) can pass at a much higher error rate than the overall FDR, because its decoys are diluted among the
targets of the large class. See ClassFDR for the estimator. Only available for PSM-level FDR (@p -peptide false).

Since the class q-value can only raise a PSM's q-value, the class FDR only removes PSMs; it never adds any. It is
conservative in small classes and in classes whose targets are mostly false: there the null-fraction estimate is clipped
at 1, and PSMs of a class that is already calibrated are lost too (in a semi-specific search, about 40% of the
semi-specific PSMs at 1%). It cannot correct a large class whose own error rate is too high; that is up to the overall
FDR estimate.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FalseDiscoveryRate.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FalseDiscoveryRate.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFalseDiscoveryRate :
  public TOPPBase
{
public:
  TOPPFalseDiscoveryRate() :
    TOPPBase("FalseDiscoveryRate", "Estimates the false discovery rate on peptide and protein level using decoy searches.")
  {
  }

protected:

  Param getSubsectionDefaults_(const std::string& /*section*/) const override
  {
    return FalseDiscoveryRate().getDefaults();
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Identifications from searching a target-decoy database.");
    setValidFormats_("in", ListUtils::create<std::string>("idXML"));
    registerOutputFile_("out", "<file>", "", "Identifications with annotated FDR");
    setValidFormats_("out", ListUtils::create<std::string>("idXML"));
    registerStringOption_("PSM", "<FDR level>", "true", "Perform FDR calculation on PSM level", false);
    setValidStrings_("PSM", ListUtils::create<std::string>("true,false"));
    registerStringOption_("peptide", "<FDR level>", "false", "Perform FDR calculation on peptide level and annotates it as meta value\n(Note: if set, also calculates FDR/q-value on PSM level.)", false);
    setValidStrings_("peptide", ListUtils::create<std::string>("true,false"));
    registerStringOption_("PSM_peptide_base_score", "<score name or type>", "", "Set if you want to choose a different score than the last calculated main score for PSM or peptide level.", false);
    registerStringOption_("PSM_peptide_base_score_orientation", "<higher/lower>", "", "In case the score orientation cannot be inferred.", false, true);
    setValidStrings_("PSM_peptide_base_score_orientation", ListUtils::create<std::string>("higher_better, lower_better"));
    registerStringOption_("PSM_class", "<class>", "", "Also control the PSM FDR within classes of PSMs: 'modification' (by the variable modifications of a PSM), or the name of a peptide hit meta value whose values are the classes. A PSM's q-value becomes the larger of its overall and its class q-value, so this only removes PSMs: conservative in small or mostly false classes, and it cannot correct a large class that is itself miscalibrated. Empty: off.", false);
    registerStringOption_("protein", "<FDR level>", "true", "Perform FDR calculation on protein level", false);
    setValidStrings_("protein", ListUtils::create<std::string>("true,false"));
    registerStringOption_("proteingroup", "<FDR level>", "false", "Perform FDR calculation on (indist.) protein group level, too. Currently, this will enable protein FDR automatically (since internals need to be in-sync) but will affect the level at which it filters (if enabled).", false);
    setValidStrings_("proteingroup", ListUtils::create<std::string>("true,false"));

    registerStringOption_("protein_score", "<type>", "", "The protein score used to calculate the protein FDR. If empty, the main score is used.", false, true);
    auto ids = IDScoreSwitcherAlgorithm();
    setValidStrings_("protein_score", ids.getScoreNames()); // lists all scores (including PSM only scores)

    registerStringOption_("protein_base_score", "<score name or type>", "", "Set if you want to choose a different score than the last calculated main score for protein (group) level.", false);
    registerStringOption_("protein_base_score_orientation", "<higher/lower>", "", "Set if you want to choose a different score than the last calculated main score for protein (group) level.", false, true);
    setValidStrings_("protein_base_score_orientation", ListUtils::create<std::string>("higher_better, lower_better"));

    registerTOPPSubsection_("FDR", "FDR control");
    registerDoubleOption_("FDR:PSM", "<fraction>", 1, "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, disabled for 1)", false);
    setMinFloat_("FDR:PSM", 0);
    setMaxFloat_("FDR:PSM", 1);

    registerDoubleOption_("FDR:protein", "<fraction>", 1, "Filter proteins based on q-value (e.g., 0.05 = 5% FDR, disabled for 1)", false);
    setMinFloat_("FDR:protein", 0);
    setMaxFloat_("FDR:protein", 1);

    registerTOPPSubsection_("FDR:cleanup", "Cleanup references after FDR control");
    registerStringOption_("FDR:cleanup:remove_proteins_without_psms","<choice>", "true",
        "Remove proteins without PSMs (due to being decoy or below PSM FDR threshold).", false, true);
    setValidStrings_("FDR:cleanup:remove_proteins_without_psms", {"true","false"});
    registerStringOption_("FDR:cleanup:remove_psms_without_proteins","<choice>", "true",
        "Remove PSMs without proteins (due to being decoy or below protein FDR threshold).", false, true);
    setValidStrings_("FDR:cleanup:remove_psms_without_proteins", {"true","false"});
    registerStringOption_("FDR:cleanup:remove_spectra_without_psms","<choice>", "true",
        "Remove spectra without PSMs (due to being decoy or below protein FDR threshold)."
        " Caution: if remove_psms_without_proteins is false, protein level filtering does not propagate.", false, true);
    setValidStrings_("FDR:cleanup:remove_spectra_without_psms", {"true","false"});

    registerSubsection_("algorithm", "Parameter section for the FDR calculation algorithm");
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    Param alg_param = getParam_().copy("algorithm:", true);
    FalseDiscoveryRate fdr;

    fdr.setParameters(alg_param);
    writeDebug_("Parameters passed to FalseDiscoveryRate", alg_param, 3);

    // input/output files
    std::string in = getStringOption_("in");
    std::string out = getStringOption_("out");
    const double protein_fdr = getDoubleOption_("FDR:protein");
    const double psm_fdr = getDoubleOption_("FDR:PSM");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeptideIdentificationList pep_ids;
    vector<ProteinIdentification> prot_ids;

    FileHandler().loadIdentifications(in, prot_ids, pep_ids, {FileTypes::IDXML});

    Size n_prot_ids = prot_ids.size();
    Size n_prot_hits = IDFilter::countHits(prot_ids);
    Size n_pep_ids = pep_ids.size();
    Size n_pep_hits = IDFilter::countHits(pep_ids);

    bool filter_applied(false);

    try
    {
      bool groups = getStringOption_("proteingroup") == "true";
      if (getStringOption_("protein") == "true" || groups)
      {
        std::string protein_score = getStringOption_("protein_score");
        if (!protein_score.empty())
        {
          try 
          {
            IDScoreSwitcherAlgorithm::ScoreType score_type = IDScoreSwitcherAlgorithm::toScoreTypeEnum(protein_score);
            IDScoreSwitcherAlgorithm switcher;
            Size c = 0;
            switcher.switchToGeneralScoreType(prot_ids, score_type, c);
          }
          catch (Exception::MissingInformation& e)
          {
            IDScoreSwitcherAlgorithm switcher;
            auto params = switcher.getParameters();
            params.setValue("new_score", protein_score);
            params.setValue("new_score_orientation", getStringOption_("protein_base_score_orientation"));
            params.setValue("proteins", "true");
            switcher.setParameters(params);
            Size c = 0;
            for (auto& run : prot_ids)
            {
              switcher.switchScores(run,c);
            }
          }
        }

        for (auto& run : prot_ids)
        {
          if (!run.hasInferenceData() && !getFlag_("force"))
          {
            throw OpenMS::Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "It seems like protein inference was not yet performed."
              " Calculating Protein FDR is probably not meaningful. To override,"
              " use the force flag.");
          }
          else
          {
            fdr.applyBasic(run, groups);
            if (protein_fdr < 1)
            {
              if (groups)
              {
                OPENMS_LOG_INFO << "FDR control: Filtering protein groups..." << endl;
                IDFilter::filterGroupsByScore(run.getIndistinguishableProteins(), protein_fdr, run.isHigherScoreBetter());
                filter_applied = true;
              }
              else
              {
                OPENMS_LOG_INFO << "FDR control: Filtering proteins..." << endl;
                IDFilter::filterHitsByScore(run, protein_fdr);
                filter_applied = true;
              }
            }
          }
        }
      }

      bool peptide_level_fdr = getStringOption_("peptide") == "true";
      bool psm_level_fdr = getStringOption_("PSM") == "true";

      const std::string psm_class = getStringOption_("PSM_class");
      if (!psm_class.empty() && peptide_level_fdr)
      {
        OPENMS_LOG_ERROR << "Error: 'PSM_class' is only supported for PSM-level FDR ('peptide' false)." << endl;
        return ILLEGAL_PARAMETERS;
      }
      bool use_class_fdr = !psm_class.empty() && psm_level_fdr;
      if (use_class_fdr)
      {
        std::function<std::string(const PeptideHit&)> class_of;
        if (psm_class == "modification")
        {
          std::vector<std::string> fixed;
          for (const auto& run : prot_ids)
          {
            const auto& mods = run.getSearchParameters().fixed_modifications;
            fixed.insert(fixed.end(), mods.begin(), mods.end());
          }
          class_of = [fixed](const PeptideHit& hit) { return ClassFDR::classFromModifications(hit, fixed); };
        }
        else
        {
          class_of = [psm_class](const PeptideHit& hit) { return ClassFDR::classFromMetaValue(hit, psm_class); };
        }
        const auto classes = ClassFDR::annotateClassQValues(pep_ids, class_of, alg_param.getValue("use_all_hits").toBool());
        OPENMS_LOG_INFO << "PSM classes (" << psm_class << "):" << endl;
        for (const auto& c : classes)
        {
          OPENMS_LOG_INFO << "  " << c.name << ": " << c.targets << " targets, " << c.decoys << " decoys, "
                          << (c.method == ClassFDR::Method::POOLED_NULL ? "pooled-decoy q-values, pi0 " : "class target-decoy q-values");
          if (c.method == ClassFDR::Method::POOLED_NULL) OPENMS_LOG_INFO << c.pi0;
          OPENMS_LOG_INFO << endl;
        }
        if (classes.size() < 2)
        {
          OPENMS_LOG_INFO << "Only one PSM class: class FDR not applied." << endl;
          use_class_fdr = false;
        }
      }

      if (psm_level_fdr || peptide_level_fdr)
      {
        fdr.apply(pep_ids, peptide_level_fdr);
        if (use_class_fdr)
        {
          if (!pep_ids.empty() && pep_ids.front().getScoreType() != "q-value")
          {
            OPENMS_LOG_WARN << "Class FDR needs q-values ('no_qvalues' false): not applied." << endl;
          }
          else
          {
            const Size raised = ClassFDR::raiseToClassQValues(pep_ids);
            for (auto& id : pep_ids) id.sort();
            OPENMS_LOG_INFO << "Class FDR raised the q-value of " << raised << " PSMs." << endl;
          }
        }
        // TODO If no decoys are removed in the param settings, we shouldn't need cleanups
        //  but then all tests need to be changed since cleanup sorts.
        //if (alg_param.getValue("add_decoy_peptides").toBool())
        //{
        //  filter_applied = true;
        //}
        filter_applied = true;
        
        if (psm_fdr < 1)
        {
          filter_applied = true;
          OPENMS_LOG_INFO << "FDR control: Filtering PSMs..." << endl;
          IDFilter::filterHitsByScore(pep_ids, psm_fdr);
        }
      }
    }
    catch (Exception::MissingInformation& e)
    {
      OPENMS_LOG_FATAL_ERROR << "FalseDiscoveryRate failed due to missing information:\n"
      << e.what();
      return INCOMPATIBLE_INPUT_DATA;
    }

    if (filter_applied)
    {
      //remove_proteins_without_psms
      if (getStringOption_("FDR:cleanup:remove_proteins_without_psms") == "true")
      {
        IDFilter::removeUnreferencedProteins(prot_ids, pep_ids);
      }
      //remove_psms_without_proteins
      IDFilter::removeDanglingProteinReferences(pep_ids,
                                        prot_ids,
                                        getStringOption_("FDR:cleanup:remove_psms_without_proteins") == "true");
      //remove_spectra_without_psms
      if (getStringOption_("FDR:cleanup:remove_spectra_without_psms") == "true")
      {
        IDFilter::removeEmptyIdentifications(pep_ids);
      }

      // we want to keep "empty" protein ID runs because they contain search meta data

      // update protein groupings if necessary:
      for (auto& prot : prot_ids)
      {
        bool valid = IDFilter::updateProteinGroups(prot.getProteinGroups(),
                                                  prot.getHits());
        if (!valid)
        {
          OPENMS_LOG_WARN << "Warning: While updating protein groups, some prot_ids were removed from groups that are still present. "
                  << "The new grouping (especially the group probabilities) may not be completely valid any more." 
                  << endl;
        }

        valid = IDFilter::updateProteinGroups(
          prot.getIndistinguishableProteins(), prot.getHits());

        if (!valid)
        {
          OPENMS_LOG_WARN << "Warning: While updating indistinguishable prot_ids, some prot_ids were removed from groups that are still present. "
                  << "The new grouping (especially the group probabilities) may not be completely valid any more." 
                  << endl;
        }
      }
    }

    // some stats
    OPENMS_LOG_INFO << "Before filtering:\n"
             << n_prot_ids << " protein identification(s) with "
             << n_prot_hits << " protein hit(s),\n"
             << n_pep_ids << " peptide identification(s) with "
             << n_pep_hits << " pep_ids hit(s).\n"
             << "After filtering:\n"
             << prot_ids.size() << " protein identification(s) with "
             << IDFilter::countHits(prot_ids) << " protein hit(s),\n"
             << pep_ids.size() << " peptide identification(s) with "
             << IDFilter::countHits(pep_ids) << " pep_ids hit(s)." << endl;

    OPENMS_LOG_INFO << "Writing filtered output..." << endl;
    FileHandler().storeIdentifications(out, prot_ids, pep_ids, {FileTypes::IDXML});
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFalseDiscoveryRate tool;
  return tool.main(argc, argv);
}

/// @endcond
