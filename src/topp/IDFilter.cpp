// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <limits>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDFilter IDFilter

@brief Filters peptide/protein identification results by different criteria.
<CENTER>
 <table>
  <tr>
   <th ALIGN = "center"> potential predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=5> &rarr; IDFilter &rarr;</td>
   <th ALIGN = "center"> potential successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapterOnline (or other ID engines) </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinInference </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier (for spectral counting) </td>
  </tr>
 </table>
</CENTER>

This tool is used to filter the identifications found by a peptide/protein identification engine like Mascot.
Different filters can be applied.
To enable any of them, just change their default value.
All active filters are applied in order.

Most filtering options should be straight-forward - see the documentation of the different parameters.
For some filters that warrent further discussion, see below.

<b>Score filters</b> (@p score:peptide, @p score:protein):

Peptide or protein hits with scores at least as good as the given cut-off are retained by the filter; hits with worse scores are removed.
Whether scores should be higher or lower than the cut-off depends on the type/orientation of the score.

The score that was most recently set by a processing step is considered for filtering.
For example, it could be a Mascot score (if MascotAdapterOnline was applied) or an FDR (if FalseDiscoveryRate was applied), etc.
@ref TOPP_IDScoreSwitcher is useful to switch to a particular score before filtering.

<b>Protein accession filters</b> (@p whitelist:proteins, @p whitelist:protein_accessions, @p blacklist:proteins, @p blacklist:protein_accessions):

These filters retain only peptide and protein hits that @e do (whitelist) or <em>do not</em> (blacklist) match any of the proteins from a given set.
This set of proteins can be given through a FASTA file (<tt>...:proteins</tt>) or as a list of accessions (<tt>...:protein_accessions</tt>).

Note that even in the case of a FASTA file, matching is only done by protein accession, not by sequence.
If necessary, use @ref TOPP_PeptideIndexer to generate protein references for peptide hits via sequence look-up.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDFilter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFilter :
  public TOPPBase
{
public:
  TOPPIDFilter() :
    TOPPBase("IDFilter", "Filters results from protein or peptide identification engines based on different criteria.")
  {

  }

protected:

  void registerOptionsAndFlags_() override
  {
    vector<String> all_mods;
    StringList all_enzymes;
    StringList specificity;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    specificity.assign(EnzymaticDigestion::NamesOfSpecificity, EnzymaticDigestion::NamesOfSpecificity + 3); //only allow none,semi,full for now

    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", {"idXML","consensusXML"});
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", {"idXML","consensusXML"});

    registerTOPPSubsection_("precursor", "Filtering by precursor attributes (RT, m/z, charge, length)");
    registerStringOption_("precursor:rt", "[min]:[max]", ":", "Retention time range to extract.", false);
    registerStringOption_("precursor:mz", "[min]:[max]", ":", "Mass-to-charge range to extract.", false);
    registerStringOption_("precursor:length", "[min]:[max]", ":", "Keep only peptide hits with a sequence length in this range.", false);
    registerStringOption_("precursor:charge", "[min]:[max]", ":", "Keep only peptide hits with charge states in this range.", false);

    registerTOPPSubsection_("score", "Filtering by peptide/protein score.");
    auto ids = IDScoreSwitcherAlgorithm();
    registerDoubleOption_("score:psm", "<score>", NAN, "The score which should be reached by a peptide hit to be kept. (use 'NAN' to disable this filter)", false);
    registerDoubleOption_("score:peptide", "<score>", NAN, "The score which should be reached by a peptide hit to be kept.  (use 'NAN' to disable this filter)", false);
    registerStringOption_("score:type_peptide", "<type>", "", "Score used for filtering. If empty, the main score is used.", false, true);
    setValidStrings_("score:type_peptide", ids.getScoreTypeNames());

    registerDoubleOption_("score:protein", "<score>", NAN, "The score which should be reached by a protein hit to be kept. All proteins are filtered based on their singleton scores irrespective of grouping. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides. (use 'NAN' to disable this filter)", false);
    registerStringOption_("score:type_protein", "<type>", "", "The type of the score which should be reached by a protein hit to be kept. If empty, the most recently set score is used.", false, true);
    setValidStrings_("score:type_protein", ids.getScoreTypeNames());

    registerDoubleOption_("score:proteingroup", "<score>", NAN, "The score which should be reached by a protein group to be kept. Performs group level score filtering (including groups of single proteins). Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides. (use 'NAN' to disable this filter)", false);

    registerTOPPSubsection_("whitelist", "Filtering by whitelisting (only peptides/proteins from a given set can pass)");
    registerInputFile_("whitelist:proteins", "<file>", "", "Filename of a FASTA file containing protein sequences.\n"
                                                           "All peptides that are not referencing a protein in this file are removed.\n"
                                                           "All proteins whose accessions are not present in this file are removed.", false);
    setValidFormats_("whitelist:proteins", {"fasta"});
    registerStringList_("whitelist:protein_accessions", "<accessions>", vector<String>(), "All peptides that do not reference at least one of the provided protein accession are removed.\nOnly proteins of the provided list are retained.", false);
    registerInputFile_("whitelist:peptides", "<file>", "", "Only peptides with the same sequence and modification assignment as any peptide in this file are kept. Use with 'whitelist:ignore_modifications' to only compare by sequence.\n", false);
    setValidFormats_("whitelist:peptides", {"idXML"});
    registerFlag_("whitelist:ignore_modifications", "Compare whitelisted peptides by sequence only.", true);
    registerStringList_("whitelist:modifications", "<selection>", vector<String>(), "Keep only peptides with sequences that contain (any of) the selected modification(s)", false, true);
    setValidStrings_("whitelist:modifications", all_mods);

    registerTOPPSubsection_("blacklist", "Filtering by blacklisting (only peptides/proteins NOT present in a given set can pass)");
    registerInputFile_("blacklist:proteins", "<file>", "", "Filename of a FASTA file containing protein sequences.\n"
                                                           "All peptides that are referencing a protein in this file are removed.\n"
                                                           "All proteins whose accessions are present in this file are removed.", false);
    setValidFormats_("blacklist:proteins", {"fasta"});
    registerStringList_("blacklist:protein_accessions", "<accessions>", vector<String>(), "All peptides that reference at least one of the provided protein accession are removed.\nOnly proteins not in the provided list are retained.", false);
    registerInputFile_("blacklist:peptides", "<file>", "", "Peptides with the same sequence and modification assignment as any peptide in this file are filtered out. Use with 'blacklist:ignore_modifications' to only compare by sequence.\n", false);
    setValidFormats_("blacklist:peptides", {"idXML"});
    registerFlag_("blacklist:ignore_modifications", "Compare blacklisted peptides by sequence only.", true);
    registerStringList_("blacklist:modifications", "<selection>", vector<String>(), "Remove all peptides with sequences that contain (any of) the selected modification(s)", false, true);
    setValidStrings_("blacklist:modifications", all_mods);
    registerStringOption_("blacklist:RegEx",  "<selection>", "", "Remove all peptides with (unmodified) sequences matched by the RegEx e.g. [BJXZ] removes ambiguous peptides.", false, true);
    registerTOPPSubsection_("in_silico_digestion", "This filter option removes peptide hits which are not in the list of in silico peptides generated by the rules specified below");
    registerInputFile_("in_silico_digestion:fasta", "<file>", "", "fasta protein sequence database.", false);
    setValidFormats_("in_silico_digestion:fasta", {"fasta"});
    registerStringOption_("in_silico_digestion:enzyme", "<enzyme>", "Trypsin", "enzyme used for the digestion of the sample",false, true);
    setValidStrings_("in_silico_digestion:enzyme", all_enzymes);
    registerStringOption_("in_silico_digestion:specificity", "<specificity>", specificity[EnzymaticDigestion::SPEC_FULL], "Specificity of the filter", false, true);
    setValidStrings_("in_silico_digestion:specificity", specificity);
    registerIntOption_("in_silico_digestion:missed_cleavages", "<integer>", -1,
                       "range of allowed missed cleavages in the peptide sequences\n"
                       "By default missed cleavages are ignored", false, true);
    setMinInt_("in_silico_digestion:missed_cleavages", -1);
    registerFlag_("in_silico_digestion:methionine_cleavage", "Allow methionine cleavage at the N-terminus of the protein.", true);

    registerTOPPSubsection_("missed_cleavages", "This filter option removes peptide hits which do not confirm with the allowed missed cleavages specified below.");
    registerStringOption_("missed_cleavages:number_of_missed_cleavages", "[min]:[max]", ":",
                          "range of allowed missed cleavages in the peptide sequences.\n"
                          "For example: 0:1 -> peptides with two or more missed cleavages will be removed,\n"
                          "0:0 -> peptides with any missed cleavages will be removed", false);
    registerStringOption_("missed_cleavages:enzyme", "<enzyme>", "Trypsin", "enzyme used for the digestion of the sample", false, true);
    setValidStrings_("missed_cleavages:enzyme", all_enzymes);

    registerTOPPSubsection_("rt", "Filtering by RT predicted by 'RTPredict'");
    registerDoubleOption_("rt:p_value", "<float>", 0.0, "Retention time filtering by the p-value predicted by RTPredict.", false, true);
    registerDoubleOption_("rt:p_value_1st_dim", "<float>", 0.0, "Retention time filtering by the p-value predicted by RTPredict for first dimension.", false, true);
    setMinFloat_("rt:p_value", 0);
    setMaxFloat_("rt:p_value", 1);
    setMinFloat_("rt:p_value_1st_dim", 0);
    setMaxFloat_("rt:p_value_1st_dim", 1);

    registerTOPPSubsection_("mz", "Filtering by mass error");
    registerDoubleOption_("mz:error", "<float>", -1, "Filtering by deviation to theoretical mass (disabled for negative values).", false, true);
    registerStringOption_("mz:unit", "<String>", "ppm", "Absolute or relative error.", false, true);
    setValidStrings_("mz:unit", ListUtils::create<String>("Da,ppm"));

    registerTOPPSubsection_("best", "Filtering best hits per spectrum (for peptides) or from proteins");
    registerIntOption_("best:n_spectra", "<integer>", 0, "Keep only the 'n' best spectra (i.e., PeptideIdentifications) (for n > 0). A spectrum is considered better if it has a higher scoring peptide hit than the other spectrum.", false);
    setMinInt_("best:n_spectra", 0);
    registerIntOption_("best:n_peptide_hits", "<integer>", 0, "Keep only the 'n' highest scoring peptide hits per spectrum (for n > 0).", false);
    setMinInt_("best:n_peptide_hits", 0);
    registerStringOption_("best:spectrum_per_peptide", "<String>", "false", "Keep one spectrum per peptide. Value determines if same sequence but different charges or modifications are treated as separate peptides or the same peptide. (default: false = filter disabled).", false);
    setValidStrings_("best:spectrum_per_peptide", {"false", "sequence", "sequence+charge", "sequence+modification", "sequence+charge+modification"});    
    registerIntOption_("best:n_protein_hits", "<integer>", 0, "Keep only the 'n' highest scoring protein hits (for n > 0).", false);
    setMinInt_("best:n_protein_hits", 0);
    registerFlag_("best:strict", "Keep only the highest scoring peptide hit.\n"
                                 "Similar to n_peptide_hits=1, but if there are ties between two or more highest scoring hits, none are kept.");
    registerStringOption_("best:n_to_m_peptide_hits", "[min]:[max]", ":", "Peptide hit rank range to extracts", false, true);


    registerFlag_("var_mods", "Keep only peptide hits with variable modifications (as defined in the 'SearchParameters' section of the input file).", false);

    registerFlag_("remove_duplicate_psm", "Removes duplicated PSMs per spectrum and retains the one with the higher score.", true);
    registerFlag_("remove_shared_peptides", "Only peptides matching exactly one protein are kept. Remember that isoforms count as different proteins!");
    registerFlag_("keep_unreferenced_protein_hits", "Proteins not referenced by a peptide are retained in the IDs.");
    registerFlag_("remove_decoys", "Remove proteins according to the information in the user parameters. Usually used in combination with 'delete_unreferenced_peptide_hits'.");
    registerFlag_("delete_unreferenced_peptide_hits", "Peptides not referenced by any protein are deleted in the IDs. Usually used in combination with 'score:protein' or 'thresh:prot'.");

    registerStringList_("remove_peptide_hits_by_metavalue", "<name> 'lt|eq|gt|ne' <value>", StringList(), "Expects a 3-tuple (=3 entries in the list), i.e. <name> 'lt|eq|gt|ne' <value>; the first is the name of meta value, followed by the comparison operator (equal, less, greater, not equal) and the value to compare to. All comparisons are done after converting the given value to the corresponding data value type of the meta value (for lists, this simply compares length, not content!)!", false, true);
  }


  ExitCodes main_(int, const char**) override
  {
    const String tmp_feature_id_metaval_ = "tmp_feature_id";
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    //only used for cxml
    ConsensusMap cmap;
    unordered_map<UInt64, ConsensusFeature*> id_to_featureref;

    const auto& infiletype = FileHandler::getType(inputfile_name);
    if (infiletype == FileTypes::IDXML)
    {
      FileHandler().loadIdentifications(inputfile_name, proteins, peptides, {FileTypes::IDXML});
    }
    else if (infiletype == FileTypes::CONSENSUSXML)
    {
      FileHandler().loadConsensusFeatures(inputfile_name, cmap, {FileTypes::CONSENSUSXML});
      for (auto& f : cmap)
      {
        UInt64 id = f.getUniqueId();
        id_to_featureref[id] = &f;
        for (auto& p : f.getPeptideIdentifications())
        {
          p.setMetaValue(tmp_feature_id_metaval_, String(id));
          peptides.push_back(std::move(p));
          //if ((UInt64)peptides.back().getMetaValue(tmp_feature_id_metaval_) != id) std::cout << "WHAT THE FUCK" << std::endl;
        }
        f.getPeptideIdentifications().clear();
      }
      auto& unassigned = cmap.getUnassignedPeptideIdentifications();
      peptides.reserve(peptides.size() + unassigned.size());
      std::move(std::begin(unassigned),std::end(unassigned),
          std::back_inserter(peptides));
      unassigned.clear();

      std::swap(proteins, cmap.getProteinIdentifications());
    }


    Size n_prot_ids = proteins.size();
    Size n_prot_hits = IDFilter::countHits(proteins);
    Size n_pep_ids = peptides.size();
    Size n_pep_hits = IDFilter::countHits(peptides);

    // handle remove_meta
    StringList meta_info = getStringList_("remove_peptide_hits_by_metavalue");
    bool remove_meta_enabled = (!meta_info.empty());
    if (remove_meta_enabled && meta_info.size() != 3)
    {
      writeLogError_("Param 'remove_peptide_hits_by_metavalue' has invalid number of arguments. Expected 3, got " + String(meta_info.size()) + ". Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    if (remove_meta_enabled && !(meta_info[1] == "lt" || meta_info[1] == "eq" || meta_info[1] == "gt" || meta_info[1] == "ne"))
    {
      writeLogError_("Param 'remove_peptide_hits_by_metavalue' has invalid second argument. Expected one of 'lt', 'eq', 'gt' or 'ne'. Got '" + meta_info[1] + "'. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // Filtering peptide identification according to set criteria

    double rt_high = numeric_limits<double>::infinity(), rt_low = -rt_high;
    if (parseRange_(getStringOption_("precursor:rt"), rt_low, rt_high))
    {
      OPENMS_LOG_INFO << "Filtering peptide IDs by precursor RT..." << endl;
      IDFilter::filterPeptidesByRT(peptides, rt_low, rt_high);
    }

    double mz_high = numeric_limits<double>::infinity(), mz_low = -mz_high;
    if (parseRange_(getStringOption_("precursor:mz"), mz_low, mz_high))
    {
      OPENMS_LOG_INFO << "Filtering peptide IDs by precursor m/z...";
      IDFilter::filterPeptidesByMZ(peptides, mz_low, mz_high);
    }


    // Filtering peptide hits according to set criteria

    if (getFlag_("remove_duplicate_psm"))
    {
      OPENMS_LOG_INFO << "Removing duplicated psms..." << endl;
      IDFilter::removeDuplicatePeptideHits(peptides);
    }

    if (getFlag_("remove_shared_peptides"))
    {
      OPENMS_LOG_INFO << "Filtering peptides by unique match to a protein..." << endl;
      IDFilter::keepUniquePeptidesPerProtein(peptides);
    }

    double pred_rt_pv = getDoubleOption_("rt:p_value");
    if (pred_rt_pv > 0)
    {
      OPENMS_LOG_INFO << "Filtering by RT prediction p-value..." << endl;
      IDFilter::filterPeptidesByRTPredictPValue(
        peptides, "predicted_RT_p_value", pred_rt_pv);
    }

    double pred_rt_pv_1d = getDoubleOption_("rt:p_value_1st_dim");
    if (pred_rt_pv_1d > 0)
    {
      OPENMS_LOG_INFO << "Filtering by RT prediction p-value (first dim.)..." << endl;
      IDFilter::filterPeptidesByRTPredictPValue(
        peptides, "predicted_RT_p_value_first_dim", pred_rt_pv_1d);
    }

    String whitelist_fasta = getStringOption_("whitelist:proteins").trim();
    if (!whitelist_fasta.empty())
    {
      OPENMS_LOG_INFO << "Filtering by protein whitelisting (FASTA input)..." << endl;
      // load protein accessions from FASTA file:
      vector<FASTAFile::FASTAEntry> fasta;
      FASTAFile().load(whitelist_fasta, fasta);
      set<String> accessions;
      for (vector<FASTAFile::FASTAEntry>::iterator it = fasta.begin();
           it != fasta.end(); ++it)
      {
        accessions.insert(it->identifier);
      }
      IDFilter::keepHitsMatchingProteins(peptides, accessions);
      IDFilter::keepHitsMatchingProteins(proteins, accessions);
    }

    vector<String> whitelist_accessions =
      getStringList_("whitelist:protein_accessions");
    if (!whitelist_accessions.empty())
    {
      OPENMS_LOG_INFO << "Filtering by protein whitelisting (accessions input)..."
               << endl;
      set<String> accessions(whitelist_accessions.begin(),
                             whitelist_accessions.end());
      IDFilter::keepHitsMatchingProteins(peptides, accessions);
      IDFilter::keepHitsMatchingProteins(proteins, accessions);
    }

    String whitelist_peptides = getStringOption_("whitelist:peptides").trim();
    if (!whitelist_peptides.empty())
    {
      OPENMS_LOG_INFO << "Filtering by inclusion peptide whitelisting..." << endl;
      vector<PeptideIdentification> inclusion_peptides;
      vector<ProteinIdentification> inclusion_proteins; // ignored
      FileHandler().loadIdentifications(whitelist_peptides, inclusion_proteins,
                       inclusion_peptides, {FileTypes::IDXML});
      bool ignore_mods = getFlag_("whitelist:ignore_modifications");
      IDFilter::keepPeptidesWithMatchingSequences(peptides, inclusion_peptides,
                                                  ignore_mods);
    }

    vector<String> whitelist_mods = getStringList_("whitelist:modifications");
    if (!whitelist_mods.empty())
    {
      OPENMS_LOG_INFO << "Filtering peptide IDs by modification whitelisting..."
               << endl;
      set<String> good_mods(whitelist_mods.begin(), whitelist_mods.end());
      IDFilter::keepPeptidesWithMatchingModifications(peptides, good_mods);
    }

    String blacklist_fasta = getStringOption_("blacklist:proteins").trim();
    if (!blacklist_fasta.empty())
    {
      OPENMS_LOG_INFO << "Filtering by protein blacklisting (FASTA input)..." << endl;
      // load protein accessions from FASTA file:
      vector<FASTAFile::FASTAEntry> fasta;
      FASTAFile().load(blacklist_fasta, fasta);
      set<String> accessions;
      for (FASTAFile::FASTAEntry& ft : fasta)
      {
        accessions.insert(ft.identifier);
      }
      IDFilter::removeHitsMatchingProteins(peptides, accessions);
      IDFilter::removeHitsMatchingProteins(proteins, accessions);
    }

    vector<String> blacklist_accessions =
      getStringList_("blacklist:protein_accessions");
    if (!blacklist_accessions.empty())
    {
      OPENMS_LOG_INFO << "Filtering by protein blacklisting (accessions input)..."
               << endl;
      set<String> accessions(blacklist_accessions.begin(),
                             blacklist_accessions.end());
      IDFilter::removeHitsMatchingProteins(peptides, accessions);
      IDFilter::removeHitsMatchingProteins(proteins, accessions);
    }

    String blacklist_peptides = getStringOption_("blacklist:peptides").trim();
    if (!blacklist_peptides.empty())
    {
      OPENMS_LOG_INFO << "Filtering by exclusion peptide blacklisting..." << endl;
      vector<PeptideIdentification> exclusion_peptides;
      vector<ProteinIdentification> exclusion_proteins; // ignored
      FileHandler().loadIdentifications(blacklist_peptides, exclusion_proteins,
                       exclusion_peptides, {FileTypes::IDXML});
      bool ignore_mods = getFlag_("blacklist:ignore_modifications");
      IDFilter::removePeptidesWithMatchingSequences(
        peptides, exclusion_peptides, ignore_mods);
    }

    vector<String> blacklist_mods = getStringList_("blacklist:modifications");
    if (!blacklist_mods.empty())
    {
      OPENMS_LOG_INFO << "Filtering peptide IDs by modification blacklisting..."
               << endl;
      set<String> bad_mods(blacklist_mods.begin(), blacklist_mods.end());
      IDFilter::removePeptidesWithMatchingModifications(peptides, bad_mods);
    }

    String blacklist_regex = getStringOption_("blacklist:RegEx");
    if (!blacklist_regex.empty())
    {
      IDFilter::removePeptidesWithMatchingRegEx(peptides, blacklist_regex);
    }

    if (getFlag_("best:strict"))
    {
      OPENMS_LOG_INFO << "Filtering by best peptide hits..." << endl;
      IDFilter::keepBestPeptideHits(peptides, true);
    }


    Int min_length = 0, max_length = 0;
    if (parseRange_(getStringOption_("precursor:length"), min_length, max_length))
    {
      OPENMS_LOG_INFO << "Filtering by peptide length..." << endl;
      if ((min_length < 0) || (max_length < 0))
      {
        OPENMS_LOG_ERROR << "Fatal error: negative values are not allowed for parameter 'precursor:length'" << endl;
        return ILLEGAL_PARAMETERS;
      }
      IDFilter::filterPeptidesByLength(peptides, Size(min_length),
                                       Size(max_length));
    }

    // Filter by digestion enzyme product

    String protein_fasta = getStringOption_("in_silico_digestion:fasta").trim();
    if (!protein_fasta.empty())
    {
      OPENMS_LOG_INFO << "Filtering peptides by digested protein (FASTA input)..." << endl;
      // load protein accessions from FASTA file:
      vector<FASTAFile::FASTAEntry> fasta;
      FASTAFile().load(protein_fasta, fasta);

      // Configure Enzymatic digestion
      ProteaseDigestion digestion;
      String enzyme = getStringOption_("in_silico_digestion:enzyme").trim();
      if (!enzyme.empty())
      {
        digestion.setEnzyme(enzyme);
      }

      String specificity = getStringOption_("in_silico_digestion:specificity").trim();
      if (!specificity.empty())
      {
        digestion.setSpecificity(digestion.getSpecificityByName(specificity));
      }

      Int missed_cleavages = getIntOption_("in_silico_digestion:missed_cleavages");
      bool ignore_missed_cleavages = true;
      if (missed_cleavages > -1)
      {
        ignore_missed_cleavages = false;
        if (digestion.getSpecificity() == EnzymaticDigestion::SPEC_FULL)
        {
          OPENMS_LOG_WARN << "Specificity not full, missed_cleavages option is redundant" << endl;
        }
        digestion.setMissedCleavages(missed_cleavages);
      }

      bool methionine_cleavage = false;
      if (getFlag_("in_silico_digestion:methionine_cleavage"))
      {
        methionine_cleavage = true;
      }

      // Build the digest filter function
      IDFilter::DigestionFilter filter(fasta,
                                       digestion,
                                       ignore_missed_cleavages,
                                       methionine_cleavage);
      // Filter peptides
      filter.filterPeptideEvidences(peptides);
    }

    // Filter peptide hits by missing cleavages

    Int min_cleavages, max_cleavages;
    min_cleavages = max_cleavages = IDFilter::PeptideDigestionFilter::disabledValue();

    if (parseRange_(getStringOption_("missed_cleavages:number_of_missed_cleavages"), min_cleavages, max_cleavages))
    {
      // Configure Enzymatic digestion
      ProteaseDigestion digestion;
      String enzyme = getStringOption_("missed_cleavages:enzyme");
      if (!enzyme.empty())
      {
        digestion.setEnzyme(enzyme);
      }

      OPENMS_LOG_INFO << "Filtering peptide hits by their missed cleavages count with enzyme " << digestion.getEnzymeName() << "..." << endl;

      // Build the digest filter function
      IDFilter::PeptideDigestionFilter filter(digestion, min_cleavages, max_cleavages);

      // Filter peptide hits
      for (auto& peptide : peptides)
      {
        filter.filterPeptideSequences(peptide.getHits());
      }
    }

    if (getFlag_("var_mods"))
    {
      OPENMS_LOG_INFO << "Filtering for variable modifications..." << endl;
      // gather possible variable modifications from search parameters:
      set<String> var_mods;
      for (ProteinIdentification& prot : proteins)
      {
        const ProteinIdentification::SearchParameters& params =
          prot.getSearchParameters();
        for (vector<String>::const_iterator mod_it =
               params.variable_modifications.begin(); mod_it !=
               params.variable_modifications.end(); ++mod_it)
        {
          var_mods.insert(*mod_it);
        }
      }
      IDFilter::keepPeptidesWithMatchingModifications(peptides, var_mods);
    }

    double psm_score = getDoubleOption_("score:psm");
    if (!std::isnan(psm_score))
    {
      OPENMS_LOG_INFO << "Filtering by PSM score (better than " << psm_score << ")..." << endl;
      IDFilter::filterHitsByScore(peptides, psm_score);
    }
    else
    {
      OPENMS_LOG_INFO << "No 'score:psm' threshold set. Not filtering by peptide score." << endl;
    }

    double pep_score = getDoubleOption_("score:peptide");
    String score_type = getStringOption_("score:type_peptide");

    if (!std::isnan(pep_score))
    {
      OPENMS_LOG_INFO << "Filtering by peptide score (better than " << pep_score << ")..." << endl;

      if (!score_type.empty())
      {
        IDScoreSwitcherAlgorithm::ScoreType score_type_enum = IDScoreSwitcherAlgorithm::getScoreType(score_type);
        IDFilter::filterHitsByScore(proteins, pep_score, score_type_enum);
      }
      else
      {
        IDFilter::filterHitsByScore(peptides, pep_score);
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No 'score:peptide' threshold set. Not filtering by peptide score." << endl;
    }


    Int min_charge = numeric_limits<Int>::min(), max_charge =
      numeric_limits<Int>::max();
    if (parseRange_(getStringOption_("precursor:charge"), min_charge, max_charge))
    {
      OPENMS_LOG_INFO << "Filtering by peptide charge..." << endl;
      IDFilter::filterPeptidesByCharge(peptides, min_charge, max_charge);
    }

    const Size best_n_spectra = getIntOption_("best:n_spectra");
    if (best_n_spectra > 0)
    {
      OPENMS_LOG_INFO << "Filtering by best n spectra..." << endl;
      IDFilter::keepNBestSpectra(peptides, best_n_spectra);
    }

    Size best_n_pep = getIntOption_("best:n_peptide_hits");
    if (best_n_pep > 0)
    {
      OPENMS_LOG_INFO << "Filtering by best n peptide hits..." << endl;
      IDFilter::keepNBestHits(peptides, best_n_pep);
    }

    String spectrum_per_peptide = getStringOption_("best:spectrum_per_peptide");
    if (spectrum_per_peptide != "false")
    {
      OPENMS_LOG_INFO << "Keeping best spectrum per " << spectrum_per_peptide << endl;
      if (spectrum_per_peptide == "sequence") // group by sequence and return best spectrum (->smallest number of spectra)
      {
        IDFilter::keepBestPerPeptide(peptides, true, true, 1);
      } else if (spectrum_per_peptide == "sequence+modification")
      {
        IDFilter::keepBestPerPeptide(peptides, false, true, 1);
      } else if (spectrum_per_peptide == "sequence+charge")
      {
        IDFilter::keepBestPerPeptide(peptides, true, false, 1);
      } else if (spectrum_per_peptide == "sequence+charge+modification") // group by sequence, modificationm, charge combination and return best spectrum (->largest number of spectra)
      {
        IDFilter::keepBestPerPeptide(peptides, false, false, 1);
      }
    }

    Int min_rank = 0, max_rank = 0;
    if (parseRange_(getStringOption_("best:n_to_m_peptide_hits"), min_rank,
                    max_rank))
    {
      OPENMS_LOG_INFO << "Filtering by peptide hit ranks..." << endl;
      if ((min_rank < 0) || (max_rank < 0))
      {
        OPENMS_LOG_ERROR << "Fatal error: negative values are not allowed for parameter 'best:n_to_m_peptide_hits'" << endl;
        return ILLEGAL_PARAMETERS;
      }
      IDFilter::filterHitsByRank(peptides, Size(min_rank), Size(max_rank));
    }

    double mz_error = getDoubleOption_("mz:error");
    if (mz_error > 0)
    {
      OPENMS_LOG_INFO << "Filtering by mass error..." << endl;
      bool unit_ppm = (getStringOption_("mz:unit") == "ppm");
      IDFilter::filterPeptidesByMZError(peptides, mz_error, unit_ppm);
    }


    // Filtering protein identifications according to set criteria
    double prot_score = getDoubleOption_("score:protein");
    String score_type_prot = getStringOption_("score:type_protein");

    if (!std::isnan(prot_score))
    {
      OPENMS_LOG_INFO << "Filtering by protein score  (better than " << prot_score << ") ..." << endl;
      if (!score_type_prot.empty())
      {
        IDScoreSwitcherAlgorithm::ScoreType score_type_prot_enum = IDScoreSwitcherAlgorithm::getScoreType(score_type_prot);
        IDFilter::filterHitsByScore(proteins, prot_score, score_type_prot_enum);
      }
      else
      {
        IDFilter::filterHitsByScore(proteins, prot_score);
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No 'score:protein' threshold set. Not filtering by protein score." << endl;
    }


    Size best_n_prot = getIntOption_("best:n_protein_hits");
    if (best_n_prot > 0)
    {
      OPENMS_LOG_INFO << "Filtering by best n protein hits (" << best_n_prot << ") ... " << endl;
      IDFilter::keepNBestHits(proteins, best_n_prot);
    }

    if (getFlag_("remove_decoys"))
    {
      OPENMS_LOG_INFO << "Removing decoy hits..." << endl;
      IDFilter::removeDecoyHits(peptides);
      IDFilter::removeDecoyHits(proteins);
    }

    // Filtering protein identifications according to set criteria
    double prot_grp_score = getDoubleOption_("score:proteingroup");
    if (!std::isnan(prot_grp_score))
    {
      for (auto& proteinid : proteins)
      {
        OPENMS_LOG_INFO << "Filtering by protein group score..." << endl;
        IDFilter::filterGroupsByScore(proteinid.getIndistinguishableProteins(), prot_grp_score, proteinid.isHigherScoreBetter());
        IDFilter::filterGroupsByScore(proteinid.getProteinGroups(), prot_grp_score, proteinid.isHigherScoreBetter());
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No 'score:proteingroup' threshold set. Not filtering by protein group score." << endl;
    }


    // remove peptide hits with meta values:
    if (remove_meta_enabled)
    {
      auto checkMVs = [this, &meta_info](PeptideHit& ph)->bool
      {
        if (!ph.metaValueExists(meta_info[0])) return true; // not having the meta value means passing the test
        DataValue v_data = ph.getMetaValue(meta_info[0]);
        DataValue v_user;
        switch (v_data.valueType())
        {
          case DataValue::STRING_VALUE : v_user = String(meta_info[2]); break;
          case DataValue::INT_VALUE : v_user = String(meta_info[2]).toInt(); break;
          case DataValue::DOUBLE_VALUE : v_user = String(meta_info[2]).toDouble(); break;
          case DataValue::STRING_LIST : v_user = (StringList)ListUtils::create<String>(meta_info[2]); break;
          case DataValue::INT_LIST : v_user = ListUtils::create<Int>(meta_info[2]); break;
          case DataValue::DOUBLE_LIST : v_user = ListUtils::create<double>(meta_info[2]); break;
          case DataValue::EMPTY_VALUE : v_user = DataValue::EMPTY; break;
          default: throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unkown!"); break;
        }

        if (meta_info[1] == "lt")
        {
          return !(v_data < v_user);
        }
        else if (meta_info[1] == "eq")
        {
          return !(v_data == v_user);
        }
        else if (meta_info[1] == "gt")
        {
          return !(v_data > v_user);
        }
        else if (meta_info[1] == "ne")
        {
          return (v_data == v_user);
        }
        else
        {
          writeLogError_("Internal Error. Meta value filtering got invalid comparison operator ('" + meta_info[1] + "'), which should have been caught before! Aborting!");
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Illegal meta value filtering operator!");
        }
      }; // of lambda

      for (auto & pid : peptides)
      {
        vector<PeptideHit>& phs = pid.getHits();
        phs.erase(remove_if(phs.begin(), phs.end(), checkMVs), phs.end());
      }
    }

    // Clean-up:

    // propagate filter from PSM level to protein level
    if (!getFlag_("keep_unreferenced_protein_hits"))
    {
      OPENMS_LOG_INFO << "Removing unreferenced protein hits..." << endl;
      IDFilter::removeUnreferencedProteins(proteins, peptides);
    }

    // propagate filter from protein level to protein group level
    for (ProteinIdentification& prot : proteins)
    {
      bool valid = IDFilter::updateProteinGroups(prot.getProteinGroups(),
                                                 prot.getHits());
      if (!valid)
      {
        OPENMS_LOG_WARN << "Warning: While updating protein groups, some proteins were removed from groups that are still present. The new grouping (especially the group probabilities) may not be completely valid any more." << endl;
      }

      valid = IDFilter::updateProteinGroups(
          prot.getIndistinguishableProteins(), prot.getHits());
      if (!valid)
      {
        OPENMS_LOG_WARN << "Warning: While updating indistinguishable proteins, some proteins were removed from groups that are still present. The new grouping (especially the group probabilities) may not be completely valid any more." << endl;
      }

      if (!std::isnan(prot_grp_score))
      {
        // Pass potential filtering on group level down to proteins
        IDFilter::removeUngroupedProteins(prot.getIndistinguishableProteins(), prot.getHits());
      }
    }

    // remove non-existant protein references from peptides (and optionally:
    // remove peptides with no proteins):
    bool rm_pep = getFlag_("delete_unreferenced_peptide_hits");
    if (rm_pep)
    {
      OPENMS_LOG_INFO << "Removing peptide hits without protein references..." << endl;
    }
    IDFilter::updateProteinReferences(peptides, proteins, rm_pep);

    IDFilter::removeEmptyIdentifications(peptides);
    // we want to keep "empty" protein IDs because they contain search meta data

    IDFilter::updateHitRanks(proteins);
    IDFilter::updateHitRanks(peptides);

    // some stats
    OPENMS_LOG_INFO << "Before filtering:\n"
             << n_prot_ids << " identification runs with "
             << n_prot_hits << " proteins,\n"
             << n_pep_ids << " spectra identified with "
             << n_pep_hits << " spectrum matches.\n"
             << "After filtering:\n"
             << proteins.size() << " identification runs with "
             << IDFilter::countHits(proteins) << " proteins,\n"
             << peptides.size() << " spectra identified with "
             << IDFilter::countHits(peptides) << " spectrum matches." << endl;

    if (infiletype == FileTypes::IDXML)
    {
      FileHandler().storeIdentifications(outputfile_name, proteins, peptides, {FileTypes::IDXML});
    }
    else if (infiletype == FileTypes::CONSENSUSXML)
    {
      for (auto& p : peptides)
      {
        if (p.metaValueExists(tmp_feature_id_metaval_))
        {
          UInt64 id(0);
          std::stringstream(p.getMetaValue(tmp_feature_id_metaval_).toString()) >> id;
          p.removeMetaValue(tmp_feature_id_metaval_);
          auto& f = id_to_featureref[id];
          f->getPeptideIdentifications().push_back(std::move(p));
        }
        else
        {
          cmap.getUnassignedPeptideIdentifications().push_back(std::move(p));
        }
      }
      peptides.clear();
      std::swap(proteins, cmap.getProteinIdentifications());
      proteins.clear();
      FileHandler().storeConsensusFeatures(outputfile_name, cmap, {FileTypes::CONSENSUSXML});
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDFilter tool;

  return tool.main(argc, argv);
}

/// @endcond
