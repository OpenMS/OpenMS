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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
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
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=5> \f$ \longrightarrow \f$ IDFilter \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
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

 <b>Score filters</b> (@p score:pep, @p score:prot, @p thresh:pep, @p thresh:prot):

 Peptide or protein hits with scores at least as good as the given cut-off are retained by the filter; hits with worse scores are removed.
 Whether scores should be higher or lower than the cut-off depends on the type/orientation of the score.

 The score that was most recently set by a processing step is considered for filtering.
 For example, it could be a Mascot score (if MascotAdapterOnline was applied) or an FDR (if FalseDiscoveryRate was applied), etc.
 @ref UTILS_IDScoreSwitcher is useful to switch to a particular score before filtering.

 An example to illustrate the significance threshold filters (@p thresh:pep, @p thresh:prot):
 Assume a peptide hit has a score of 30, the significance threshold is 40, and higher scores are better. Then the hit will be kept if the cut-off value is set to 0.75 or lower, and removed for higher cut-offs.

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
    specificity.assign(EnzymaticDigestion::NamesOfSpecificity, EnzymaticDigestion::NamesOfSpecificity + EnzymaticDigestion::SIZE_OF_SPECIFICITY);

    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerTOPPSubsection_("precursor", "Filtering by precursor RT or m/z");
    registerStringOption_("precursor:rt", "[min]:[max]", ":", "Retention time range to extract.", false);
    registerStringOption_("precursor:mz", "[min]:[max]", ":", "Mass-to-charge range to extract.", false);

    registerTOPPSubsection_("score", "Filtering by peptide/protein score.");
    registerDoubleOption_("score:pep", "<score>", 0, "The score which should be reached by a peptide hit to be kept.", false);
    registerDoubleOption_("score:prot", "<score>", 0, "The score which should be reached by a protein hit to be kept. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides.", false);
    registerTOPPSubsection_("thresh", "Filtering by significance threshold");
    registerDoubleOption_("thresh:pep", "<fraction>", 0.0, "Keep a peptide hit only if its score is above this fraction of the peptide significance threshold.", false, true);
    registerDoubleOption_("thresh:prot", "<fraction>", 0.0, "Keep a protein hit only if its score is above this fraction of the protein significance threshold. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides.", false, true);

    registerTOPPSubsection_("whitelist", "Filtering by whitelisting (only peptides/proteins from a given set can pass)");
    registerInputFile_("whitelist:proteins", "<file>", "", "Filename of a FASTA file containing protein sequences.\n"
                                                           "All peptides that are not referencing a protein in this file are removed.\n"
                                                           "All proteins whose accessions are not present in this file are removed.", false);
    setValidFormats_("whitelist:proteins", ListUtils::create<String>("fasta"));
    registerStringList_("whitelist:protein_accessions", "<accessions>", vector<String>(), "All peptides that do not reference at least one of the provided protein accession are removed.\nOnly proteins of the provided list are retained.", false);
    registerInputFile_("whitelist:peptides", "<file>", "", "Only peptides with the same sequence and modification assignment as any peptide in this file are kept. Use with 'whitelist:ignore_modifications' to only compare by sequence.\n", false);
    setValidFormats_("whitelist:peptides", ListUtils::create<String>("idXML"));
    registerFlag_("whitelist:ignore_modifications", "Compare whitelisted peptides by sequence only.", false);
    registerStringList_("whitelist:modifications", "<selection>", vector<String>(), "Keep only peptides with sequences that contain (any of) the selected modification(s)", false);
    setValidStrings_("whitelist:modifications", all_mods);

    registerTOPPSubsection_("blacklist", "Filtering by blacklisting (only peptides/proteins NOT present in a given set can pass)");
    registerInputFile_("blacklist:proteins", "<file>", "", "Filename of a FASTA file containing protein sequences.\n"
                                                           "All peptides that are referencing a protein in this file are removed.\n"
                                                           "All proteins whose accessions are present in this file are removed.", false);
    setValidFormats_("blacklist:proteins", ListUtils::create<String>("fasta"));
    registerStringList_("blacklist:protein_accessions", "<accessions>", vector<String>(), "All peptides that reference at least one of the provided protein accession are removed.\nOnly proteins not in the provided list are retained.", false);
    registerInputFile_("blacklist:peptides", "<file>", "", "Peptides with the same sequence and modification assignment as any peptide in this file are filtered out. Use with 'blacklist:ignore_modifications' to only compare by sequence.\n", false);
    setValidFormats_("blacklist:peptides", ListUtils::create<String>("idXML"));
    registerFlag_("blacklist:ignore_modifications", "Compare blacklisted peptides by sequence only.", false);
    registerStringList_("blacklist:modifications", "<selection>", vector<String>(), "Remove all peptides with sequences that contain (any of) the selected modification(s)", false);
    setValidStrings_("blacklist:modifications", all_mods);

    registerTOPPSubsection_("in_silico_digestion", "This filter option removes peptide hits which are not in the list of in silico peptides generated by the rules specified below");
    registerInputFile_("in_silico_digestion:fasta", "<file>", "", "fasta protein sequence database.", false);
    setValidFormats_("in_silico_digestion:fasta", ListUtils::create<String>("fasta"));
    registerStringOption_("in_silico_digestion:enzyme", "<enzyme>", "Trypsin", "enzyme used for the digestion of the sample",false);
    setValidStrings_("in_silico_digestion:enzyme", all_enzymes);
    registerStringOption_("in_silico_digestion:specificity", "<specificity>", specificity[EnzymaticDigestion::SPEC_FULL], "Specificity of the filter", false);
    setValidStrings_("in_silico_digestion:specificity", specificity);
    registerIntOption_("in_silico_digestion:missed_cleavages", "<integer>", -1, 
                       "range of allowed missed cleavages in the peptide sequences\n"
                       "By default missed cleavages are ignored", false);
    setMinInt_("in_silico_digestion:missed_cleavages", -1);
    registerFlag_("in_silico_digestion:methionine_cleavage", "Allow methionine cleavage at the N-terminus of the protein.", false);

    registerTOPPSubsection_("missed_cleavages", "This filter option removes peptide hits which do not confirm with the allowed missed cleavages specified below.");
    registerStringOption_("missed_cleavages:number_of_missed_cleavages", "[min]:[max]", ":",
                          "range of allowed missed cleavages in the peptide sequences.\n"
                          "For example: 0:1 -> peptides with two or more missed cleavages will be removed,\n"
                          "0:0 -> peptides with any missed cleavages will be removed", false);
    registerStringOption_("missed_cleavages:enzyme", "<enzyme>", "Trypsin", "enzyme used for the digestion of the sample", false);
    setValidStrings_("missed_cleavages:enzyme", all_enzymes);

    registerTOPPSubsection_("rt", "Filtering by RT predicted by 'RTPredict'");
    registerDoubleOption_("rt:p_value", "<float>", 0.0, "Retention time filtering by the p-value predicted by RTPredict.", false, true);
    registerDoubleOption_("rt:p_value_1st_dim", "<float>", 0.0, "Retention time filtering by the p-value predicted by RTPredict for first dimension.", false, true);
    setMinFloat_("rt:p_value", 0);
    setMaxFloat_("rt:p_value", 1);
    setMinFloat_("rt:p_value_1st_dim", 0);
    setMaxFloat_("rt:p_value_1st_dim", 1);

    registerTOPPSubsection_("mz", "Filtering by mass error");
    registerDoubleOption_("mz:error", "<float>", -1, "Filtering by deviation to theoretical mass (disabled for negative values).", false);
    registerStringOption_("mz:unit", "<String>", "ppm", "Absolute or relative error.", false);
    setValidStrings_("mz:unit", ListUtils::create<String>("Da,ppm"));

    registerTOPPSubsection_("best", "Filtering best hits per spectrum (for peptides) or from proteins");
    registerIntOption_("best:n_peptide_hits", "<integer>", 0, "Keep only the 'n' highest scoring peptide hits per spectrum (for n > 0).", false);
    setMinInt_("best:n_peptide_hits", 0);
    registerIntOption_("best:n_protein_hits", "<integer>", 0, "Keep only the 'n' highest scoring protein hits (for n > 0).", false);
    setMinInt_("best:n_protein_hits", 0);
    registerFlag_("best:strict", "Keep only the highest scoring peptide hit.\n"
                                 "Similar to n_peptide_hits=1, but if there are ties between two or more highest scoring hits, none are kept.");
    registerStringOption_("best:n_to_m_peptide_hits", "[min]:[max]", ":", "Peptide hit rank range to extracts", false, true);

    registerStringOption_("length", "[min]:[max]", ":", "Keep only peptide hits with a sequence length in this range.", false);

    registerStringOption_("charge", "[min]:[max]", ":", "Keep only peptide hits with charge states in this range.", false);

    registerFlag_("var_mods", "Keep only peptide hits with variable modifications (as defined in the 'SearchParameters' section of the input file).", false);

    registerFlag_("unique", "If a peptide hit occurs more than once per peptide ID, only one instance is kept.");
    registerFlag_("unique_per_protein", "Only peptides matching exactly one protein are kept. Remember that isoforms count as different proteins!");
    registerFlag_("keep_unreferenced_protein_hits", "Proteins not referenced by a peptide are retained in the IDs.");
    registerFlag_("remove_decoys", "Remove proteins according to the information in the user parameters. Usually used in combination with 'delete_unreferenced_peptide_hits'.");
    registerFlag_("delete_unreferenced_peptide_hits", "Peptides not referenced by any protein are deleted in the IDs. Usually used in combination with 'score:prot' or 'thresh:prot'.");

  }


  ExitCodes main_(int, const char**) override
  {
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    IdXMLFile().load(inputfile_name, proteins, peptides);

    Size n_prot_ids = proteins.size();
    Size n_prot_hits = IDFilter::countHits(proteins);
    Size n_pep_ids = peptides.size();
    Size n_pep_hits = IDFilter::countHits(peptides);


    // Filtering peptide identification according to set criteria

    double rt_high = numeric_limits<double>::infinity(), rt_low = -rt_high;
    if (parseRange_(getStringOption_("precursor:rt"), rt_low, rt_high))
    {
      LOG_INFO << "Filtering peptide IDs by precursor RT..." << endl;
      IDFilter::filterPeptidesByRT(peptides, rt_low, rt_high);
    }

    double mz_high = numeric_limits<double>::infinity(), mz_low = -mz_high;
    if (parseRange_(getStringOption_("precursor:mz"), mz_low, mz_high))
    {
      LOG_INFO << "Filtering peptide IDs by precursor m/z...";
      IDFilter::filterPeptidesByMZ(peptides, mz_low, mz_high);
    }


    // Filtering peptide hits according to set criteria

    if (getFlag_("unique"))
    {
      LOG_INFO << "Removing duplicate peptide hits..." << endl;
      IDFilter::removeDuplicatePeptideHits(peptides);
    }

    if (getFlag_("unique_per_protein"))
    {
      LOG_INFO << "Filtering peptides by unique match to a protein..." << endl;
      IDFilter::keepUniquePeptidesPerProtein(peptides);
    }

    double peptide_significance = getDoubleOption_("thresh:pep");
    if (peptide_significance > 0)
    {
      LOG_INFO << "Filtering by peptide significance threshold..." << endl;
      IDFilter::filterHitsBySignificance(peptides, peptide_significance);
    }

    double pred_rt_pv = getDoubleOption_("rt:p_value");
    if (pred_rt_pv > 0)
    {
      LOG_INFO << "Filtering by RT prediction p-value..." << endl;
      IDFilter::filterPeptidesByRTPredictPValue(
        peptides, "predicted_RT_p_value", pred_rt_pv);
    }

    double pred_rt_pv_1d = getDoubleOption_("rt:p_value_1st_dim");
    if (pred_rt_pv_1d > 0)
    {
      LOG_INFO << "Filtering by RT prediction p-value (first dim.)..." << endl;
      IDFilter::filterPeptidesByRTPredictPValue(
        peptides, "predicted_RT_p_value_first_dim", pred_rt_pv_1d);
    }

    String whitelist_fasta = getStringOption_("whitelist:proteins").trim();
    if (!whitelist_fasta.empty())
    {
      LOG_INFO << "Filtering by protein whitelisting (FASTA input)..." << endl;
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
      LOG_INFO << "Filtering by protein whitelisting (accessions input)..."
               << endl;
      set<String> accessions(whitelist_accessions.begin(),
                             whitelist_accessions.end());
      IDFilter::keepHitsMatchingProteins(peptides, accessions);
      IDFilter::keepHitsMatchingProteins(proteins, accessions);
    }

    String whitelist_peptides = getStringOption_("whitelist:peptides").trim();
    if (!whitelist_peptides.empty())
    {
      LOG_INFO << "Filtering by inclusion peptide whitelisting..." << endl;
      vector<PeptideIdentification> inclusion_peptides;
      vector<ProteinIdentification> inclusion_proteins; // ignored
      IdXMLFile().load(whitelist_peptides, inclusion_proteins,
                       inclusion_peptides);
      bool ignore_mods = getFlag_("whitelist:ignore_modifications");
      IDFilter::keepPeptidesWithMatchingSequences(peptides, inclusion_peptides,
                                                  ignore_mods);
    }

    vector<String> whitelist_mods = getStringList_("whitelist:modifications");
    if (!whitelist_mods.empty())
    {
      LOG_INFO << "Filtering peptide IDs by modification whitelisting..."
               << endl;
      set<String> good_mods(whitelist_mods.begin(), whitelist_mods.end());
      IDFilter::keepPeptidesWithMatchingModifications(peptides, good_mods);
    }

    String blacklist_fasta = getStringOption_("blacklist:proteins").trim();
    if (!blacklist_fasta.empty())
    {
      LOG_INFO << "Filtering by protein blacklisting (FASTA input)..." << endl;
      // load protein accessions from FASTA file:
      vector<FASTAFile::FASTAEntry> fasta;
      FASTAFile().load(blacklist_fasta, fasta);
      set<String> accessions;
      for (vector<FASTAFile::FASTAEntry>::iterator it = fasta.begin();
           it != fasta.end(); ++it)
      {
        accessions.insert(it->identifier);
      }
      IDFilter::removeHitsMatchingProteins(peptides, accessions);
      IDFilter::removeHitsMatchingProteins(proteins, accessions);
    }

    vector<String> blacklist_accessions =
      getStringList_("blacklist:protein_accessions");
    if (!blacklist_accessions.empty())
    {
      LOG_INFO << "Filtering by protein blacklisting (accessions input)..."
               << endl;
      set<String> accessions(blacklist_accessions.begin(),
                             blacklist_accessions.end());
      IDFilter::removeHitsMatchingProteins(peptides, accessions);
      IDFilter::removeHitsMatchingProteins(proteins, accessions);
    }

    String blacklist_peptides = getStringOption_("blacklist:peptides").trim();
    if (!blacklist_peptides.empty())
    {
      LOG_INFO << "Filtering by exclusion peptide blacklisting..." << endl;
      vector<PeptideIdentification> exclusion_peptides;
      vector<ProteinIdentification> exclusion_proteins; // ignored
      IdXMLFile().load(blacklist_peptides, exclusion_proteins,
                       exclusion_peptides);
      bool ignore_mods = getFlag_("blacklist:ignore_modifications");
      IDFilter::removePeptidesWithMatchingSequences(
        peptides, exclusion_peptides, ignore_mods);
    }

    vector<String> blacklist_mods = getStringList_("blacklist:modifications");
    if (!blacklist_mods.empty())
    {
      LOG_INFO << "Filtering peptide IDs by modification blacklisting..."
               << endl;
      set<String> bad_mods(blacklist_mods.begin(), blacklist_mods.end());
      IDFilter::removePeptidesWithMatchingModifications(peptides, bad_mods);
    }


    if (getFlag_("best:strict"))
    {
      LOG_INFO << "Filtering by best peptide hits..." << endl;
      IDFilter::keepBestPeptideHits(peptides, true);
    }


    Int min_length = 0, max_length = 0;
    if (parseRange_(getStringOption_("length"), min_length, max_length))
    {
      LOG_INFO << "Filtering by peptide length..." << endl;
      if ((min_length < 0) || (max_length < 0))
      {
        LOG_ERROR << "Fatal error: negative values are not allowed for parameter 'length'" << endl;
        return ILLEGAL_PARAMETERS;
      }
      IDFilter::filterPeptidesByLength(peptides, Size(min_length),
                                       Size(max_length));
    }

    // Filter by digestion enzyme product

    String protein_fasta = getStringOption_("in_silico_digestion:fasta").trim();
    if (!protein_fasta.empty())
    {
      LOG_INFO << "Filtering peptides by digested protein (FASTA input)..." << endl;
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
          LOG_WARN << "Specificity not full, missed_cleavages option is redundant" << endl;
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

      LOG_INFO << "Filtering peptide hits by their missed cleavages count with enzyme " << digestion.getEnzymeName() << "..." << endl;

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
      LOG_INFO << "Filtering for variable modifications..." << endl;
      // gather possible variable modifications from search parameters:
      set<String> var_mods;
      for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
           prot_it != proteins.end(); ++prot_it)
      {
        const ProteinIdentification::SearchParameters& params =
          prot_it->getSearchParameters();
        for (vector<String>::const_iterator mod_it =
               params.variable_modifications.begin(); mod_it !=
               params.variable_modifications.end(); ++mod_it)
        {
          var_mods.insert(*mod_it);
        }
      }
      IDFilter::keepPeptidesWithMatchingModifications(peptides, var_mods);
    }

    double pep_score = getDoubleOption_("score:pep");
    // @TODO: what if 0 is a reasonable cut-off for some score?
    if (pep_score != 0)
    {
      LOG_INFO << "Filtering by peptide score..." << endl;
      IDFilter::filterHitsByScore(peptides, pep_score);
    }

    Int min_charge = numeric_limits<Int>::min(), max_charge =
      numeric_limits<Int>::max();
    if (parseRange_(getStringOption_("charge"), min_charge, max_charge))
    {
      LOG_INFO << "Filtering by peptide charge..." << endl;
      IDFilter::filterPeptidesByCharge(peptides, min_charge, max_charge);
    }

    Size best_n_pep = getIntOption_("best:n_peptide_hits");
    if (best_n_pep > 0)
    {
      LOG_INFO << "Filtering by best n peptide hits..." << endl;
      IDFilter::keepNBestHits(peptides, best_n_pep);
    }

    Int min_rank = 0, max_rank = 0;
    if (parseRange_(getStringOption_("best:n_to_m_peptide_hits"), min_rank,
                    max_rank))
    {
      LOG_INFO << "Filtering by peptide hit ranks..." << endl;
      if ((min_rank < 0) || (max_rank < 0))
      {
        LOG_ERROR << "Fatal error: negative values are not allowed for parameter 'best:n_to_m_peptide_hits'" << endl;
        return ILLEGAL_PARAMETERS;
      }
      IDFilter::filterHitsByRank(peptides, Size(min_rank), Size(max_rank));
    }

    double mz_error = getDoubleOption_("mz:error");
    if (mz_error > 0)
    {
      LOG_INFO << "Filtering by mass error..." << endl;
      bool unit_ppm = (getStringOption_("mz:unit") == "ppm");
      IDFilter::filterPeptidesByMZError(peptides, mz_error, unit_ppm);
    }


    // Filtering protein identifications according to set criteria

    double protein_significance = getDoubleOption_("thresh:prot");
    if (protein_significance > 0)
    {
      LOG_INFO << "Filtering by protein significance threshold..." << endl;
      IDFilter::filterHitsBySignificance(proteins, protein_significance);
    }

    double prot_score = getDoubleOption_("score:prot");
    // @TODO: what if 0 is a reasonable cut-off for some score?
    if (prot_score != 0)
    {
      LOG_INFO << "Filtering by protein score..." << endl;
      IDFilter::filterHitsByScore(proteins, prot_score);
    }

    Size best_n_prot = getIntOption_("best:n_protein_hits");
    if (best_n_prot > 0)
    {
      LOG_INFO << "Filtering by best n protein hits..." << endl;
      IDFilter::keepNBestHits(proteins, best_n_prot);
    }

    if (getFlag_("remove_decoys"))
    {
      LOG_INFO << "Removing decoy hits..." << endl;
      IDFilter::removeDecoyHits(peptides);
      IDFilter::removeDecoyHits(proteins);
    }


    // Clean-up:

    if (!getFlag_("keep_unreferenced_protein_hits"))
    {
      LOG_INFO << "Removing unreferenced protein hits..." << endl;
      IDFilter::removeUnreferencedProteins(proteins, peptides);
    }

    IDFilter::updateHitRanks(proteins);
    IDFilter::updateHitRanks(peptides);

    // remove non-existant protein references from peptides (and optionally:
    // remove peptides with no proteins):
    bool rm_pep = getFlag_("delete_unreferenced_peptide_hits");
    if (rm_pep) LOG_INFO << "Removing peptide hits without protein references..." << endl;
    IDFilter::updateProteinReferences(peptides, proteins, rm_pep);

    IDFilter::removeEmptyIdentifications(peptides);
    // we want to keep "empty" protein IDs because they contain search meta data

    // update protein groupings if necessary:
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      bool valid = IDFilter::updateProteinGroups(prot_it->getProteinGroups(),
                                                 prot_it->getHits());
      if (!valid)
      {
        LOG_WARN << "Warning: While updating protein groups, some proteins were removed from groups that are still present. The new grouping (especially the group probabilities) may not be completely valid any more." << endl;
      }

      valid = IDFilter::updateProteinGroups(
        prot_it->getIndistinguishableProteins(), prot_it->getHits());
      if (!valid)
      {
        LOG_WARN << "Warning: While updating indistinguishable proteins, some proteins were removed from groups that are still present. The new grouping (especially the group probabilities) may not be completely valid any more." << endl;
      }
    }

    // some stats
    LOG_INFO << "Before filtering:\n"
             << n_prot_ids << " protein identification(s) with "
             << n_prot_hits << " protein hit(s),\n"
             << n_pep_ids << " peptide identification(s) with "
             << n_pep_hits << " peptides hit(s).\n"
             << "After filtering:\n"
             << proteins.size() << " protein identification(s) with "
             << IDFilter::countHits(proteins) << " protein hit(s),\n"
             << peptides.size() << " peptide identification(s) with "
             << IDFilter::countHits(peptides) << " peptides hit(s)." << endl;

    IdXMLFile().store(outputfile_name, proteins, peptides);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDFilter tool;

  return tool.main(argc, argv);
}

/// @endcond
