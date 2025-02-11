// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <unordered_set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ConsensusID ConsensusID

@brief Computes a consensus from results of multiple peptide identification engines.

<CENTER>
<table>
    <tr>
        <th ALIGN = "center"> potential predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=4> &rarr; ConsensusID &rarr;</td>
        <th ALIGN = "center"> potential successor tools </td>
    </tr>
    <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_PeptideIndexer </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDFilter </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDMapper </td>
    </tr>
</table>
</CENTER>

<B>Reference:</B>

Nahnsen <em>et al.</em>: <a href="https://doi.org/10.1021/pr2002879">Probabilistic consensus scoring improves tandem mass spectrometry peptide identification</a> (J. Proteome Res., 2011, PMID: 21644507).

<B>Algorithms:</B>

ConsensusID offers several algorithms that can aggregate results from multiple peptide identification engines ("search engines") into consensus identifications - typically one per MS2 spectrum. This works especially well for search engines that provide more than one peptide hit per spectrum, i.e. that report not just the best hit, but also a list of runner-up candidates with corresponding scores.

The available algorithms are (see also @ref OpenMS::ConsensusIDAlgorithm and its subclasses):
@li @p PEPMatrix: Scoring based on posterior error probabilities (PEPs) and peptide sequence similarities. This algorithm uses a substitution matrix to score the similarity of sequences not listed by all search engines. It requires PEPs as the scores for all peptide hits.
@li @p PEPIons: Scoring based on posterior error probabilities (PEPs) and fragment ion similarities ("shared peak count"). This algorithm, too, requires PEPs as scores.
@li @p best: For each peptide ID, this uses the best score of any search engine as the consensus score. All peptide IDs must have the same score type.
@li @p worst: For each peptide ID, this uses the worst score of any search engine as the consensus score. All peptide IDs must have the same score type.
@li @p average: For each peptide ID, this uses the average score of all search engines as the consensus score. Again, all peptide IDs must have the same score type.
@li @p ranks: Calculates a consensus score based on the ranks of peptide IDs in the results of different search engines. The final score is in the range (0, 1], with 1 being the best score. The input peptide IDs do not need to have the same score type.

PEPs for search results can be calculated using the @ref TOPP_IDPosteriorErrorProbability tool, which supports a variety of search engines.

@note Important: All protein-level identification results will be lost by applying ConsensusID. (It is unclear how potentially conflicting protein-level results from different search engines should be combined.) If necessary, run the @ref TOPP_PeptideIndexer tool to add protein references for peptides again.

@note Peptides with different post-translational modifications (PTMs), or with different site localizations of the same PTMs, are treated as different peptides by all algorithms. However, a qualification applies for the @p PEPMatrix algorithm: The similarity scoring method used there can only take unmodified peptide sequences into account, so PTMs are ignored during that step. However, the PTMs are not removed from the peptides, and there will be separate results for differently-modified peptides.

<B>File types:</B>

Different input files types are supported:
@li idXML: A file containing multiple identification runs, typically from different search engines. Use @ref TOPP_IDMerger to merge individual idXML files from different search runs into one. During the ConsensusID analysis, the identification results will be grouped according to their originating MS2 spectra, based on retention time and precursor m/z information (see parameters @p rt_delta and @p mz_delta). One consensus identification will be generated for each group. With the per_spectrum flag you can also input multiple idXML files. A consensus will be made per combination of originating mzml file and spectrum_ref.
@li featureXML or consensusXML: Given (consensus) features annotated with peptide identifications from multiple search runs, one consensus identification is created for every annotated feature. Peptide identifications not assigned to features are not considered and will be removed. See @ref TOPP_IDMapper for the task of mapping peptide identifications to feature maps or consensus maps.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>Filtering:</B>

Generally, search results can be filtered according to various criteria using @ref TOPP_IDFilter before (or after) applying this tool. ConsensusID itself offers only a limited number of filtering options that are especially useful in its context (see the @p filter parameter section):
@li @p considered_hits: Limits the number of alternative peptide hits considered per spectrum/feature for each identification run. This helps to reduce runtime, especially for the @p PEPMatrix and @p PEPIons algorithms, which involve costly "all vs. all" comparisons of peptide hits.
@li @p min_support: This allows filtering of peptide hits based on agreement between search engines. Every peptide sequence in the analysis has been identified by at least one search run. This parameter defines which fraction (between 0 and 1) of the remaining search runs must "support" a peptide identification that should be kept. The meaning of "support" differs slightly between algorithms: For @p best, @p worst, @p average and @p rank, each search run supports peptides that it has also identified among its top @p considered_hits candidates. So @p min_support simply gives the fraction of additional search engines that must have identified a peptide. (For example, if there are three search runs, and only peptides identified by at least two of them should be kept, set @p min_support to 0.5.) For the similarity-based algorithms @p PEPMatrix and @p PEPIons, the "support" for a peptide is the average similarity of the most-similar peptide from each (other) search run. (In the context of the JPR publication, this is the average of the similarity scores used in the consensus score calculation for a peptide.)
@li @p count_empty: Typically not all search engines will provide results for all searched MS2 spectra. This parameter determines whether search runs that provided no results should be counted in the "support" calculation; by default, they are ignored.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ConsensusID.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ConsensusID.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPConsensusID :
  public TOPPBase
{
public:
  TOPPConsensusID() :
    TOPPBase("ConsensusID", "Computes a consensus of peptide identifications of several identification engines.")
  {
  }

protected:

  String algorithm_; // algorithm for consensus calculation (input parameter)
  bool keep_old_scores_;

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file(s)>", {}, "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML,featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML,featureXML,consensusXML"));

    addEmptyLine_();
    registerDoubleOption_("rt_delta", "<value>", 0.1, "[idXML input only] Maximum allowed retention time deviation between identifications belonging to the same spectrum.", false);
    setMinFloat_("rt_delta", 0.0);
    registerDoubleOption_("mz_delta", "<value>", 0.1, "[idXML input only] Maximum allowed precursor m/z deviation between identifications belonging to the same spectrum.", false);
    setMinFloat_("mz_delta", 0.0);

    registerFlag_("per_spectrum", "(only idXML) if set, mapping will be done based on exact matching of originating mzml file and spectrum_ref");

    // General algorithm parameters are defined in the abstract base class
    // "ConsensusIDAlgorithm", but we can't get them from there because we can't
    // instantiate the class. So we get those parameters from a subclass that
    // doesn't add any other parameters:
    registerTOPPSubsection_("filter", "Options for filtering peptide hits");
    registerFullParam_(ConsensusIDAlgorithmBest().getDefaults());

    registerStringOption_("algorithm", "<choice>", "PEPMatrix",
                          "Algorithm used for consensus scoring.\n"
                          "* PEPMatrix: Scoring based on posterior error probabilities (PEPs) and peptide sequence similarities (scored by a substitution matrix). Requires PEPs as scores.\n"
                          "* PEPIons: Scoring based on posterior error probabilities (PEPs) and fragment ion similarities ('shared peak count'). Requires PEPs as scores.\n"
                          "* best: For each peptide ID, use the best score of any search engine as the consensus score. Requires the same score type in all ID runs.\n"
                          "* worst: For each peptide ID, use the worst score of any search engine as the consensus score. Requires the same score type in all ID runs.\n"
                          "* average:  For each peptide ID, use the average score of all search engines as the consensus. Requires the same score type in all ID runs.\n"
                          "* ranks: Calculates a consensus score based on the ranks of peptide IDs in the results of different search engines. The final score is in the range (0, 1], with 1 being the best score. No requirements about score types.", false);
    setValidStrings_("algorithm", ListUtils::create<String>("PEPMatrix,PEPIons,best,worst,average,ranks"));

    // subsections appear in alphabetical (?) order, independent of the order
    // in which they were registered:
    registerSubsection_("PEPIons", "PEPIons algorithm parameters");
    registerSubsection_("PEPMatrix", "PEPMatrix algorithm parameters");

  }


  Param getSubsectionDefaults_(const String& section) const override
  {
    Param algo_params;
    if (section == "PEPMatrix")
    {
      algo_params = ConsensusIDAlgorithmPEPMatrix().getDefaults();
    }
    else // section == "PEPIons"
    {
      algo_params = ConsensusIDAlgorithmPEPIons().getDefaults();
    }
    // remove parameters defined in the base class (to avoid duplicates):
    algo_params.remove("filter:");
    return algo_params;
  }


  void setProteinIdentifications_(vector<ProteinIdentification>& prot_ids)
  {
    // modification params are necessary for further analysis tools (e.g. LuciPHOr2)
    set<String> fixed_mods_set;
    set<String> var_mods_set;
    StringList merged_spectra_data;
    String engine = prot_ids[0].getSearchEngine();
    String version = prot_ids[0].getSearchEngineVersion();

    // merge proteins
    unordered_set<String> seen_proteins;
    vector<ProteinHit> merged_proteins;
    for (Size i = 0; i < prot_ids.size(); ++i)
    {
      for (auto& hit : prot_ids[i].getHits())
      {
        const auto& iter_inserted = seen_proteins.emplace(hit.getAccession());
        if (iter_inserted.second)
        {
          merged_proteins.push_back(std::move(hit));
        }
      }
    }

    for (ProteinIdentification& prot : prot_ids)
    {
      ProteinIdentification::SearchParameters search_params(prot.getSearchParameters());
      std::copy(search_params.fixed_modifications.begin(), search_params.fixed_modifications.end(), std::inserter(fixed_mods_set, fixed_mods_set.end()));
      std::copy(search_params.variable_modifications.begin(), search_params.variable_modifications.end(), std::inserter(var_mods_set, var_mods_set.end()));
      StringList spectra_data;
      prot.getPrimaryMSRunPath(spectra_data);
      std::copy(spectra_data.begin(), spectra_data.end(), std::inserter(merged_spectra_data, merged_spectra_data.end()));
    }
    ProteinIdentification::SearchParameters search_params;
    std::vector<String> fixed_mods(fixed_mods_set.begin(), fixed_mods_set.end());
    std::vector<String> var_mods(var_mods_set.begin(), var_mods_set.end());
    search_params.fixed_modifications    = fixed_mods;
    search_params.variable_modifications = var_mods;


    prot_ids.clear();
    prot_ids.resize(1);
    prot_ids[0].getHits() = merged_proteins;
    prot_ids[0].setDateTime(DateTime::now());
    prot_ids[0].setSearchEngine("OpenMS/ConsensusID_" + algorithm_);
    prot_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    prot_ids[0].setSearchParameters(search_params);

    //TODO for completeness we could in the other algorithms, collect all search engines and put them here
    // or maybe put it in a DataProcessing step
    //TODO actually this only makes sense if there was only one search engine. (see the alternative
    // setProteinIdentificationSettings_)
    // best, worst, average can also be used on PEP scores for different search engines. IDPEP does not
    // overwrite the search engine (in contrast to PercolatorAdapter)
    if (algorithm_ == "best" || algorithm_ == "worst" || algorithm_ == "average")
    {
      prot_ids[0].setMetaValue("ConsensusIDBaseSearch", engine + String(":") + version);
    }

    // make file name entries unique
    std::sort(merged_spectra_data.begin(), merged_spectra_data.end());
    StringList::iterator last = std::unique(merged_spectra_data.begin(), merged_spectra_data.end());
    merged_spectra_data.erase(last, merged_spectra_data.end());
    prot_ids[0].setPrimaryMSRunPath(merged_spectra_data);
  }

  tuple<String, String, ProteinIdentification::SearchParameters> getOriginalSearchEngineSettings_(const ProteinIdentification& prot)
  {
    String engine = prot.getSearchEngine();
    const ProteinIdentification::SearchParameters& old_sp = prot.getSearchParameters();
    if (engine != "Percolator")
    {
      return std::tie(engine, prot.getSearchEngineVersion(), old_sp);
    }
    else
    {
      String original_SE = "Unknown";
      String original_SE_ver = "0.0";
      vector<String> mvkeys;

      old_sp.getKeys(mvkeys);
      for (const String & mvkey : mvkeys)
      {
        if (mvkey.hasPrefix("SE:"))
        {
          original_SE = mvkey.substr(3);
          original_SE_ver = old_sp.getMetaValue(mvkey);
          break; // multiSE percolator before consensusID not allowed; we take first only
        }
      }

      ProteinIdentification::SearchParameters sp{};
      for (const String & mvkey : mvkeys)
      {
        if (mvkey.hasPrefix(original_SE))
        {
          if (mvkey.hasSuffix("db"))
          {
            sp.db = old_sp.getMetaValue(mvkey);
          }
          else if (mvkey.hasSuffix("db_version"))
          {
            sp.db_version = old_sp.getMetaValue(mvkey);
          }
          else if (mvkey.hasSuffix("taxonomy"))
          {
            sp.taxonomy = old_sp.getMetaValue(mvkey);
          }
          else if (mvkey.hasSuffix("charges"))
          {
            sp.charges = old_sp.getMetaValue(mvkey);;
          }
          else if (mvkey.hasSuffix("fixed_modifications"))
          {
            const String& s = old_sp.getMetaValue(mvkey);
            s.split(',', sp.fixed_modifications);
          }
          else if (mvkey.hasSuffix("variable_modifications"))
          {
            const String& s = old_sp.getMetaValue(mvkey);
            s.split(',', sp.variable_modifications);
          }
          else if (mvkey.hasSuffix("missed_cleavages"))
          {
            sp.missed_cleavages = (UInt) old_sp.getMetaValue(mvkey);
          }
          else if (mvkey.hasSuffix("fragment_mass_tolerance"))
          {
            sp.fragment_mass_tolerance = (double) old_sp.getMetaValue(mvkey);
          }
          else if (mvkey.hasSuffix("fragment_mass_tolerance_ppm"))
          {
            sp.fragment_mass_tolerance_ppm = old_sp.getMetaValue(mvkey).toBool();
          }
          else if (mvkey.hasSuffix("precursor_mass_tolerance"))
          {
            sp.precursor_mass_tolerance = (double) old_sp.getMetaValue(mvkey);
          }
          else if (mvkey.hasSuffix("precursor_mass_tolerance_ppm"))
          {
            sp.precursor_mass_tolerance_ppm = old_sp.getMetaValue(mvkey).toBool();
          }
          else if (mvkey.hasSuffix("digestion_enzyme"))
          {
            Protease p = *(ProteaseDB::getInstance()->getEnzyme(old_sp.getMetaValue(mvkey)));
            sp.digestion_enzyme = p;
          }
          else if (mvkey.hasSuffix("enzyme_term_specificity"))
          {
            sp.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>((int) old_sp.getMetaValue(mvkey));
          }
        }
      }
      return std::tie(original_SE, original_SE_ver, sp);
    }
  }

  void setProteinIdentificationSettings_(ProteinIdentification& prot_id,
      vector<tuple<String, String, ProteinIdentification::SearchParameters>>& se_ver_settings,
      vector<tuple<String, String, vector<pair<String, String>>>>& rescore_ver_settings)
  {
    // modification params are necessary for further analysis tools (e.g. LuciPHOr2)
    set<String> fixed_mods_set;
    set<String> var_mods_set;
    set<EnzymaticDigestion::Specificity> specs;
    double prec_tol_ppm = 0.;
    double prec_tol_da = 0.;
    double frag_tol_ppm = 0.;
    double frag_tol_da = 0.;
    int min_chg = 10000;
    int max_chg = -10000;
    Size mc = 0;
    // we sort them to pick the same entries, no matter the order of the inputs
    set<String, std::greater<String>> enzymes;
    set<String, std::greater<String>> dbs;

    // use the first settings as basis (i.e. copy over db and enzyme and tolerance)
    // we assume that they are the same or similar
    ProteinIdentification::SearchParameters new_sp = get<2>(se_ver_settings[0]);

    // first check the rescoring procedure. Should at least be the same tool.
    // "" = IDPosteriorProbability. If parts were not rescored at all, they wont have a PEP annotated,
    // and the tool will fail in the next step (beginning of algorithm)
    // TODO maybe also consolidate/merge those settings. But they are currently only used for reporting.
    const auto& final_rescore_ver_setting = rescore_ver_settings[0];
    const String& final_rescore_algo = get<0>(final_rescore_ver_setting);
    const String& final_rescore_algo_version = get<1>(final_rescore_ver_setting);

    for (const auto& rescore_ver_setting : rescore_ver_settings)
    {
      if (get<0>(rescore_ver_setting) != final_rescore_algo
          || get<1>(rescore_ver_setting) != final_rescore_algo_version)
      {
        OPENMS_LOG_WARN << "Warning: Trying to use ConsensusID on searches with different rescoring algorithms. " +
                           get<0>(rescore_ver_setting) + " vs " + final_rescore_algo;
      }
    }
    if (!final_rescore_algo.empty()) new_sp.setMetaValue(final_rescore_algo, final_rescore_algo_version);
    for (const auto& s : get<2>(final_rescore_ver_setting))
    {
      // the metavalue names in s.first already contain the algorithm name. No need to prepend
      new_sp.setMetaValue(s.first, s.second);
    }

    bool allsamese = true;
    for (const auto& se_ver_setting : se_ver_settings)
    {
      allsamese = allsamese &&
          (get<0>(se_ver_setting) == get<0>(se_ver_settings[0]) &&
           get<1>(se_ver_setting) == get<1>(se_ver_settings[0]));

      const ProteinIdentification::SearchParameters& sp = get<2>(se_ver_setting);
      const String& SE = get<0>(se_ver_setting);
      new_sp.setMetaValue("SE:" + SE, get<1>(se_ver_setting));
      new_sp.setMetaValue(SE+":db",sp.db);
      new_sp.setMetaValue(SE+":db_version",sp.db_version);
      new_sp.setMetaValue(SE+":taxonomy",sp.taxonomy);
      new_sp.setMetaValue(SE+":charges",sp.charges);
      new_sp.setMetaValue(SE+":fixed_modifications",ListUtils::concatenate(sp.fixed_modifications, ","));
      new_sp.setMetaValue(SE+":variable_modifications",ListUtils::concatenate(sp.variable_modifications, ","));
      new_sp.setMetaValue(SE+":missed_cleavages",sp.missed_cleavages);
      new_sp.setMetaValue(SE+":fragment_mass_tolerance",sp.fragment_mass_tolerance);
      new_sp.setMetaValue(SE+":fragment_mass_tolerance_unit",sp.fragment_mass_tolerance_ppm ? "ppm" : "Da");
      new_sp.setMetaValue(SE+":precursor_mass_tolerance",sp.precursor_mass_tolerance);
      new_sp.setMetaValue(SE+":precursor_mass_tolerance_unit",sp.precursor_mass_tolerance_ppm  ? "ppm" : "Da");
      new_sp.setMetaValue(SE+":digestion_enzyme",sp.digestion_enzyme.getName());
      new_sp.setMetaValue(SE+":enzyme_term_specificity",EnzymaticDigestion::NamesOfSpecificity[sp.enzyme_term_specificity]);

      const auto& chg_pair = sp.getChargeRange();
      if (chg_pair.first != 0 && chg_pair.first < min_chg)
      {
        min_chg = chg_pair.first;
      }
      if (chg_pair.second != 0 && chg_pair.second > max_chg)
      {
        max_chg = chg_pair.second;
      }
      if (sp.missed_cleavages > mc )
      {
        mc = sp.missed_cleavages;
      }
      if (sp.fragment_mass_tolerance_ppm)
      {
        if (sp.fragment_mass_tolerance > frag_tol_ppm) frag_tol_ppm = sp.fragment_mass_tolerance;
      }
      else
      {
        if (sp.fragment_mass_tolerance > frag_tol_da) frag_tol_da = sp.fragment_mass_tolerance;
      }
      if (sp.precursor_mass_tolerance_ppm)
      {
        if (sp.precursor_mass_tolerance > prec_tol_ppm) prec_tol_ppm = sp.precursor_mass_tolerance;
      }
      else
      {
        if (sp.precursor_mass_tolerance > prec_tol_da) prec_tol_da = sp.precursor_mass_tolerance;
      }

      enzymes.insert(sp.digestion_enzyme.getName());
      dbs.insert(sp.db);
      specs.insert(sp.enzyme_term_specificity);

      std::copy(sp.fixed_modifications.begin(), sp.fixed_modifications.end(), std::inserter(fixed_mods_set, fixed_mods_set.end()));
      std::copy(sp.variable_modifications.begin(), sp.variable_modifications.end(), std::inserter(var_mods_set, var_mods_set.end()));
    }

    if (specs.find(EnzymaticDigestion::SPEC_NONE) != specs.end())
    {
      new_sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_NONE;
    }
    else if (specs.find(EnzymaticDigestion::SPEC_SEMI) != specs.end())
    {
      new_sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_SEMI;
    }
    else if (specs.find(EnzymaticDigestion::SPEC_NONTERM) != specs.end())
    {
      new_sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_NONTERM;
    }
    else if (specs.find(EnzymaticDigestion::SPEC_NOCTERM) != specs.end())
    {
      new_sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_NOCTERM;
    }
    else if (specs.find(EnzymaticDigestion::SPEC_FULL) != specs.end())
    {
      new_sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
    }

    std::vector<String> fixed_mods(fixed_mods_set.begin(), fixed_mods_set.end());
    std::vector<String> var_mods(var_mods_set.begin(), var_mods_set.end());
    new_sp.fixed_modifications    = fixed_mods;
    new_sp.variable_modifications = var_mods;

    String final_enz;
    for (const auto& enz : enzymes)
    {
      if (enz != "unknown_enzyme")
      {
        // Although the set should be sorted to start with the longest
        // versions, this extends "" to Trypsin and e.g. Trypsin to Trypsin/P
        if (enz.hasSubstring(final_enz))
        {
          final_enz = enz;
        }
        else if (!final_enz.hasSubstring(enz))
        {
          OPENMS_LOG_WARN << "Warning: Trying to use ConsensusID on searches with incompatible enzymes."
          " OpenMS officially supports only one enzyme per search. Using " + final_enz + " to (incompletely)"
          " represent the combined run. This might or might not lead to inconsistencies downstream.";
        }
      }
    }
    new_sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(final_enz);

    String final_db = *dbs.begin();
    String final_db_bn = final_db;
    final_db_bn.substitute("\\","/");
    final_db_bn = File::basename(final_db_bn);
    // we need to copy to substitute anyway
    for (auto db : dbs) // OMS_CODING_TEST_EXCLUDE
    {
      db.substitute("\\","/");
      if (File::basename(db) != final_db_bn)
      {
        OPENMS_LOG_WARN << "Warning: Trying to use ConsensusID on searches with different databases."
        " OpenMS officially supports only one database per search. Using " + final_db + " to (incompletely)"
        " represent the combined run. This might or might not lead to inconsistencies downstream.";
      }
    }

    new_sp.charges = String(min_chg) + "-" + String(max_chg);
    if (prec_tol_da > 0 && prec_tol_ppm > 0)
    {
      OPENMS_LOG_WARN << "Warning: Trying to use ConsensusID on searches with incompatible "
      "precursor tolerance units. Using Da for the combined run.";
    }
    if (prec_tol_da > 0)
    {
      new_sp.precursor_mass_tolerance = prec_tol_da;
      new_sp.precursor_mass_tolerance_ppm = false;
    }
    else
    {
      new_sp.precursor_mass_tolerance = prec_tol_ppm;
      new_sp.precursor_mass_tolerance_ppm = true;
    }
    if (frag_tol_da > 0 && frag_tol_ppm > 0)
    {
      OPENMS_LOG_WARN << "Warning: Trying to use ConsensusID on searches with incompatible "
      "fragment tolerance units. Using Da for the combined run.";
    }
    if (frag_tol_da > 0)
    {
      new_sp.fragment_mass_tolerance = frag_tol_da;
      new_sp.fragment_mass_tolerance_ppm = false;
    }
    else
    {
      new_sp.fragment_mass_tolerance = frag_tol_ppm;
      new_sp.fragment_mass_tolerance_ppm = true;
    }

    new_sp.missed_cleavages = mc;

    prot_id.setDateTime(DateTime::now());
    prot_id.setSearchEngine("OpenMS/ConsensusID_" + algorithm_);
    prot_id.setSearchEngineVersion(VersionInfo::getVersion());
    prot_id.setSearchParameters(new_sp);

    //TODO for completeness we could in the other algorithms, collect all search engines and put them here
    // or maybe put it in a DataProcessing step
    if (allsamese)
    {
      prot_id.setMetaValue("ConsensusIDBaseSearch", get<0>(se_ver_settings[0]) + String(":") + get<1>(se_ver_settings[0]));
    }
  }


  template <typename MapType>
  void processFeatureOrConsensusMap_(MapType& input_map,
                                     ConsensusIDAlgorithm* consensus)
  {
    // Problem with feature data: IDs from multiple spectra may be attached to
    // a (consensus) feature, so we may have multiple IDs from the same search
    // engine. This means that we can't just use the number of search runs as
    // our "baseline" for the number of identifications (parameter
    // "number_of_runs")! To work around this, we multiply the number of
    // different ID runs with the max. number of times we see the same ID run
    // in the annotations of a feature.

    map<String, String> runid_to_se;
    map<String, Size> id_mapping; // mapping: run ID -> index
    Size number_of_runs = input_map.getProteinIdentifications().size();
    for (Size i = 0; i < number_of_runs; ++i)
    {
      const auto& prot = input_map.getProteinIdentifications()[i];
      id_mapping[prot.getIdentifier()] = i;
      if (keep_old_scores_)
      {
        runid_to_se[prot.getIdentifier()] = prot.getOriginalSearchEngineName();
      }
    }

    // compute consensus:
    for (typename MapType::Iterator map_it = input_map.begin();
         map_it != input_map.end(); ++map_it)
    {
      vector<PeptideIdentification>& ids = map_it->getPeptideIdentifications();
      vector<Size> times_seen(number_of_runs);
      for (PeptideIdentification& pep : ids)
      {
        ++times_seen[id_mapping[pep.getIdentifier()]];
      }
      Size n_repeats = *max_element(times_seen.begin(), times_seen.end());

      consensus->apply(ids, runid_to_se, number_of_runs * n_repeats);
    }

    // create new identification run:
    setProteinIdentifications_(input_map.getProteinIdentifications());
    // remove outdated information (protein references will be broken):
    input_map.getUnassignedPeptideIdentifications().clear();
  }


  ExitCodes main_(int, const char**) override
  {
    StringList in = getStringList_("in");
    FileTypes::Type in_type = FileHandler::getType(in[0]);
    String out = getStringOption_("out");
    double rt_delta = getDoubleOption_("rt_delta");
    double mz_delta = getDoubleOption_("mz_delta");
    keep_old_scores_ = getFlag_("filter:keep_old_scores");

    //----------------------------------------------------------------
    // set up ConsensusID
    //----------------------------------------------------------------
    ConsensusIDAlgorithm* consensus;
    // general algorithm parameters:
    Param algo_params = ConsensusIDAlgorithmBest().getDefaults();
    algorithm_ = getStringOption_("algorithm");
    if (algorithm_ == "PEPMatrix")
    {
      consensus = new ConsensusIDAlgorithmPEPMatrix();
      // add algorithm-specific parameters:
      algo_params.merge(getParam_().copy("PEPMatrix:", true));
    }
    else if (algorithm_ == "PEPIons")
    {
      consensus = new ConsensusIDAlgorithmPEPIons();
      // add algorithm-specific parameters:
      algo_params.merge(getParam_().copy("PEPIons:", true));
    }
    else if (algorithm_ == "best")
    {
      consensus = new ConsensusIDAlgorithmBest();
    }
    else if (algorithm_ == "worst")
    {
      consensus = new ConsensusIDAlgorithmWorst();
    }
    else if (algorithm_ == "average")
    {
      consensus = new ConsensusIDAlgorithmAverage();
    }
    else // algorithm_ == "ranks"
    {
      consensus = new ConsensusIDAlgorithmRanks();
    }
    algo_params.update(getParam_(), false, OpenMS_Log_debug); // update general params.
    consensus->setParameters(algo_params);

    //----------------------------------------------------------------
    // idXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::IDXML)
    {
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      if (getFlag_("per_spectrum"))
      {
        map<String, unordered_map<String, vector<PeptideIdentification>>> grouping_per_file;
        map<String, unordered_set<String>> seen_proteins_per_file;
        map<String, Size> runid_to_old_run_idx;
        map<String, String> runid_to_old_se;
        // the values (new_run_idx) in mzml_to_new_run_idx correspond to the indices in mzml_to_sesettings
        map<String, Size> mzml_to_new_run_idx;
        vector<vector<tuple<String, String, ProteinIdentification::SearchParameters>>> mzml_to_sesettings;
        vector<vector<tuple<String, String, vector<pair<String,String>>>>> mzml_to_rescoresettings;

        for (const auto& infile : in)
        {
          vector<ProteinIdentification> tmp_prot_ids;
          vector<PeptideIdentification> tmp_pep_ids;
          FileHandler().loadIdentifications(infile, tmp_prot_ids, tmp_pep_ids, {FileTypes::IDXML});
          Size idx(0);
          for (const auto& prot : tmp_prot_ids)
          {
            runid_to_old_run_idx[prot.getIdentifier()] = idx++;
            if (keep_old_scores_)
            {
              runid_to_old_se[prot.getIdentifier()] = prot.getOriginalSearchEngineName();
            }
            StringList original_files;
            prot.getPrimaryMSRunPath(original_files);
            for (auto& f : original_files)
            {
              std::replace( f.begin(), f.end(), '\\', '/');
              f = FileHandler::stripExtension(File::basename(f)); // some SE adapters write full paths, some may use raw
            }
            if (original_files.size() != 1)
            {
              //TODO in theory you could also compare the whole StringList (if you want to consensusID
              // a whole "merge" of multiple ID files (e.g. fractions)
              throw Exception::InvalidValue(
                  __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  "Currently only ID runs on exactly one mzML file are supported. "
                  "Run " + prot.getIdentifier() + " contains no 'file_origin' UserParam (report issue) or has multiple entries in it (avoid merging).", String(original_files.size()));
            }
            String original_file = original_files[0];
            auto iter_inserted = seen_proteins_per_file.emplace(original_file, unordered_set<String>{});
            const auto se_ver_settings = getOriginalSearchEngineSettings_(prot);
            tuple<String, String, vector<pair<String,String>>> rescore_ver_settings{"","",vector<pair<String,String>>()};
            //TODO find a way to get/check IDPEP params.
            if (prot.getSearchEngine() == "Percolator")
            {
              get<0>(rescore_ver_settings) = prot.getSearchEngine();
              get<1>(rescore_ver_settings) = prot.getSearchEngineVersion();
              const auto& sp = prot.getSearchParameters();
              vector<String> mvkeys;
              sp.getKeys(mvkeys);
              for (const String & mvkey : mvkeys)
              {
                if (mvkey.hasPrefix("Percolator:"))
                {
                  // we do not cut the tool (here Percolator) prefix since we will use it as is
                  // in the new params
                  get<2>(rescore_ver_settings).emplace_back(mvkey, sp.getMetaValue(mvkey));
                }
              }
            }

            if (iter_inserted.second)
            {
              mzml_to_new_run_idx[original_file] = prot_ids.size();
              mzml_to_sesettings.emplace_back(vector<tuple<String, String, ProteinIdentification::SearchParameters>>{});
              mzml_to_sesettings.back().emplace_back(se_ver_settings);
              mzml_to_rescoresettings.emplace_back(vector<tuple<String, String, vector<pair<String,String>>>>{});
              mzml_to_rescoresettings.back().emplace_back(rescore_ver_settings);
              prot_ids.emplace_back(ProteinIdentification());
              prot_ids.back().setIdentifier("ConsensusID for " + original_file);
            }
            else
            {
              mzml_to_sesettings[mzml_to_new_run_idx[original_file]].emplace_back(se_ver_settings);
              mzml_to_rescoresettings[mzml_to_new_run_idx[original_file]].emplace_back(rescore_ver_settings);
            }
            for (auto& hit : prot.getHits())
            {
              auto acciter_inserted = iter_inserted.first->second.emplace(hit.getAccession());
              if (acciter_inserted.second)
              {
                prot_ids[mzml_to_new_run_idx[original_file]].getHits().emplace_back(std::move(hit));
              }
            }
          }

          for (auto& pep_id : tmp_pep_ids)
          {
            StringList original_files;
            const ProteinIdentification& old = tmp_prot_ids[runid_to_old_run_idx[pep_id.getIdentifier()]];
            old.getPrimaryMSRunPath(original_files); // the size should have been checked during the loop over proteins
            for (auto& f : original_files)
            {
              std::replace( f.begin(), f.end(), '\\', '/');
              f = FileHandler::stripExtension(File::basename(f)); // some SE adapters write full paths, some may use raw
            }
            String original_file = original_files[0];
            auto iter_inserted = grouping_per_file.emplace(original_file, unordered_map<String,vector<PeptideIdentification>>{});
            if (pep_id.metaValueExists("spectrum_reference"))
            {
              String nativeID = pep_id.getSpectrumReference();
              auto nativeid_iter_inserted = iter_inserted.first->second.emplace(nativeID, vector<PeptideIdentification>{});
              nativeid_iter_inserted.first->second.emplace_back(std::move(pep_id));
            }
          }
        }
        for (auto& file_ref_peps : grouping_per_file)
        {
          Size new_run_id = mzml_to_new_run_idx[file_ref_peps.first];
          ProteinIdentification& to_put = prot_ids[new_run_id];
          // Note: we assume that at least one of the inputs had mzML as an extension
          // we could keep track of it but IMHO we should not allow raw there at all (just complicates things)
          to_put.setPrimaryMSRunPath({file_ref_peps.first + ".mzML"});
          setProteinIdentificationSettings_(to_put, mzml_to_sesettings[new_run_id], mzml_to_rescoresettings[new_run_id]);
          for (const auto& ref_peps : file_ref_peps.second)
          {
            vector<PeptideIdentification> peps = ref_peps.second;
            if (peps.empty())
            {
              continue; //sth went wrong. skip
            }
            double mz = peps[0].getMZ();
            double rt = peps[0].getRT();
            // has to have a ref, save it, since apply might modify everything
            String ref = peps[0].getSpectrumReference();
            consensus->apply(peps, runid_to_old_se, mzml_to_sesettings[new_run_id].size());
            for (auto& p : peps)
            {
              p.setIdentifier(to_put.getIdentifier());
              p.setMZ(mz);
              p.setRT(rt);
              p.setSpectrumReference( ref);
              //TODO copy other meta values from the originals? They need to be collected
              // in the algorithm subclasses though first
              pep_ids.emplace_back(std::move(p));
            }
          }
        }
      }
      else // link spectra by RT and mz proximity and do ConsensusID on their PSMs
      {
        if (in.size() != 1)
        {
          OPENMS_LOG_FATAL_ERROR << "ConsensusID on idXML without the --per_spectrum flag, expects a single idXML file."
          "Please merge the files with IDMerger using its default settings." << std::endl;
        }
        // note: this requires a single merged idXML file.
        FileHandler().loadIdentifications(in[0], prot_ids, pep_ids, {FileTypes::IDXML});

        if (prot_ids.size() == 1)
        {
          OPENMS_LOG_FATAL_ERROR << "ConsensusID on idXML without the --per_spectrum flag expects a merged idXML file"
          "with multiple runs. Only one run found in the first file." << std::endl;
        }

        // merge peptide IDs by precursor position - this is equivalent to a
        // feature linking problem (peptide IDs from different ID runs <->
        // features from different maps), so we bring the data into a format
        // suitable for a feature grouping algorithm:
        vector<FeatureMap> maps(prot_ids.size());
        map<String, String> runid_to_se;
        map<String, Size> id_mapping; // mapping: run ID -> index (of feature map)
        for (Size i = 0; i < prot_ids.size(); ++i)
        {
          id_mapping[prot_ids[i].getIdentifier()] = i;
          if (keep_old_scores_)
          {
            runid_to_se[prot_ids[i].getIdentifier()] = prot_ids[i].getOriginalSearchEngineName();
          }
        }

        for (PeptideIdentification& pep : pep_ids)
        {
          String run_id = pep.getIdentifier();
          if (!pep.hasRT() || !pep.hasMZ())
          {
            OPENMS_LOG_FATAL_ERROR << "Peptide ID without RT and/or m/z information found in identification run '" + run_id + "'.\nMake sure that this information is included for all IDs when generating/converting search results. Aborting!" << endl;
            return INCOMPATIBLE_INPUT_DATA;
          }

          Feature feature;
          feature.setRT(pep.getRT());
          feature.setMZ(pep.getMZ());
          feature.getPeptideIdentifications().push_back(pep);
          maps[id_mapping[run_id]].push_back(feature);
        }
        // precondition for "FeatureGroupingAlgorithmQT::group":
        for (FeatureMap& map : maps)
        {
          map.updateRanges();
        }

        FeatureGroupingAlgorithmQT linker;
        Param linker_params = linker.getDefaults();
        linker_params.setValue("use_identifications", "false");
        linker_params.setValue("ignore_charge", "true");
        linker_params.setValue("distance_RT:max_difference", rt_delta);
        linker_params.setValue("distance_MZ:max_difference", mz_delta);
        linker_params.setValue("distance_MZ:unit", "Da");
        linker.setParameters(linker_params);

        ConsensusMap grouping;
        linker.group(maps, grouping);

        Size old_size = prot_ids.size();
        // create new identification run
        setProteinIdentifications_(prot_ids);

        // compute consensus
        pep_ids.clear();
        for (auto& cfeature : grouping)
        {
          auto& ids = cfeature.getPeptideIdentifications();
          consensus->apply(ids, runid_to_se, old_size);

          if (!ids.empty())
          {
            PeptideIdentification& pep_id = ids[0];
            // hits may be empty due to filtering (parameter "min_support");
            // in that case skip to avoid a warning from "IDXMLFile::store":
            if (!pep_id.getHits().empty())
            {
              pep_id.setIdentifier(prot_ids[0].getIdentifier());
              pep_id.setRT(cfeature.getRT());
              pep_id.setMZ(cfeature.getMZ());
              pep_ids.push_back(pep_id);
            }
          }
        }
      }
      // store consensus
      FileHandler().storeIdentifications(out, prot_ids, pep_ids, {FileTypes::IDXML});
    }

    //----------------------------------------------------------------
    // featureXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap map;
      FileHandler().loadFeatures(in[0], map, {FileTypes::FEATUREXML});

      processFeatureOrConsensusMap_(map, consensus);

      FileHandler().storeFeatures(out, map, {FileTypes::FEATUREXML});
    }

    //----------------------------------------------------------------
    // consensusXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::CONSENSUSXML)
    {
      ConsensusMap map;
      FileHandler().loadConsensusFeatures(in[0], map, {FileTypes::CONSENSUSXML});

      processFeatureOrConsensusMap_(map, consensus);

      FileHandler().storeConsensusFeatures(out, map, {FileTypes::CONSENSUSXML});
    }

    delete consensus;

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPConsensusID tool;
  return tool.main(argc, argv);
}

/// @endcond
