// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>

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
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
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

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return FalseDiscoveryRate().getDefaults();
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Identifications from searching a target-decoy database.");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Identifications with annotated FDR");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerStringOption_("PSM", "<FDR level>", "true", "Perform FDR calculation on PSM level", false);
    setValidStrings_("PSM", ListUtils::create<String>("true,false"));
    registerStringOption_("peptide", "<FDR level>", "false", "Perform FDR calculation on peptide level and annotates it as meta value\n(Note: if set, also calculates FDR/q-value on PSM level.)", false);
    setValidStrings_("peptide", ListUtils::create<String>("true,false"));
    registerStringOption_("PSM_peptide_base_score", "<score name or type>", "", "Set if you want to choose a different score than the last calculated main score for PSM or peptide level.", false);
    registerStringOption_("PSM_peptide_base_score_orientation", "<higher/lower>", "", "In case the score orientation cannot be inferred.", false, true);
    setValidStrings_("PSM_peptide_base_score_orientation", ListUtils::create<String>("higher_better, lower_better"));
    registerStringOption_("protein", "<FDR level>", "true", "Perform FDR calculation on protein level", false);
    setValidStrings_("protein", ListUtils::create<String>("true,false"));
    registerStringOption_("proteingroup", "<FDR level>", "false", "Perform FDR calculation on (indist.) protein group level, too. Currently, this will enable protein FDR automatically (since internals need to be in-sync) but will affect the level at which it filters (if enabled).", false);
    setValidStrings_("proteingroup", ListUtils::create<String>("true,false"));

    registerStringOption_("protein_score", "<type>", "", "The protein score used to calculate the protein FDR. If empty, the main score is used.", false, true);
    auto ids = IDScoreSwitcherAlgorithm();
    setValidStrings_("protein_score", ids.getScoreTypeNames()); // lists all scores (including PSM only scores)

    registerStringOption_("protein_base_score", "<score name or type>", "", "Set if you want to choose a different score than the last calculated main score for protein (group) level.", false);
    registerStringOption_("protein_base_score_orientation", "<higher/lower>", "", "Set if you want to choose a different score than the last calculated main score for protein (group) level.", false, true);
    setValidStrings_("protein_base_score_orientation", ListUtils::create<String>("higher_better, lower_better"));

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
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const double protein_fdr = getDoubleOption_("FDR:protein");
    const double psm_fdr = getDoubleOption_("FDR:PSM");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    vector<PeptideIdentification> pep_ids;
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
        String protein_score = getStringOption_("protein_score");
        if (!protein_score.empty())
        {
          try 
          {
            IDScoreSwitcherAlgorithm::ScoreType score_type = IDScoreSwitcherAlgorithm::getScoreType(protein_score);
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

      if (psm_level_fdr || peptide_level_fdr)
      {
        fdr.apply(pep_ids, peptide_level_fdr);
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
      IDFilter::updateProteinReferences(pep_ids,
                                        prot_ids,
                                        getStringOption_("FDR:cleanup:remove_psms_without_proteins") == "true");
      //remove_spectra_without_psms
      if (getStringOption_("FDR:cleanup:remove_spectra_without_psms") == "true")
      {
        IDFilter::removeEmptyIdentifications(pep_ids);
      }

      IDFilter::updateHitRanks(prot_ids);
      IDFilter::updateHitRanks(pep_ids);

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
