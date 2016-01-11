// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <limits>
#include <cmath>
#include <set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_IDFilter IDFilter

 @brief Filters protein identification engine results by different criteria.
<CENTER>
 <table>
  <tr>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=5> \f$ \longrightarrow \f$ IDFilter \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinInference </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
  </tr>
 </table>
</CENTER>

 This tool is used to filter the identifications found by
 a peptide/protein identification tool like Mascot. Different filters can be applied:

 To enable any of the filters, just change their default value.
 All active filters will be applied in order.

 <ul>

  <li>
   <b>precursor:rt</b>:<br> Precursor RT range for the peptide identification to be kept.
  </li>
  <li>
   <b>precursor:mz</b>:<br> Precursor m/z range for the peptide identification to be kept.
  </li>
  <li>
   <b>score:pep</b>:<br> The score a peptide hit should have to be kept.
  </li>
  <li>
   <b>score:prot</b>:<br> The score a protein hit should have to be kept.
  </li>
  <li>
   <b>thresh:pep</b>:<br> The fraction of the significance threshold that should
   be reached by a peptide hit to be kept. If for example a peptide
   has score 30 and the significance threshold is 40, the
   peptide will only be kept by the filter if the significance
   threshold fraction is set to 0.75 or lower.
  </li>
  <li>
   <b>thresh:prot</b>:<br> This parameter
   behaves in the same way as the peptide significance threshold
   fraction parameter. The only difference is that it is used
   to filter protein hits.
  </li>
  <li>
   <b>whitelist:proteins</b>:<br> If you know which proteins
   are in the measured sample you can specify a FASTA file
   which contains the protein sequences of those proteins. All
   peptides which are not a substring of a protein contained
   in the sequences file will be filtered out. The filtering is based on the
      protein identifiers attached to the peptide hits. Protein Hits not matching
      any FASTA protein are also removed.<br>
      If you want filtering using the sequence alone, then use the flag @em WhiteList:by_seq_only.
  </li>
  <li>
   <b>blacklist:peptides</b>:<br> For this option you specify an idXML file.
   All peptides that are present in both files (in-file and exclusion peptides
   file) will be dropped. Protein Hits are not affected.
  </li>
  <li><b>rt</b>:<br> To filter identifications according to their
   predicted retention times you have to set 'rt:p_value' and/or 'rt:p_value_1st_dim' larger than 0, depending which RT
      dimension you want to filter.
   This filter can only be applied to idXML files produced by @ref TOPP_RTPredict.
  </li>
  <li>
   <b>best:n_peptide_hits</b>:<br> Only the best n peptide hits of a spectrum are kept. If two hits have the same score, their order is random.
  </li>
  <li>
   <b>best:n_protein_hits</b>:<br> Only the best n protein hits of a spectrum are kept. If two hits have the same score, their order is random.
  </li>
  <li>
   <b>best:strict</b>:<br> Only the best hit of a spectrum is kept.
   If there is more than one hit for a spectrum with the maximum score, then
   none of the hits will be kept. This is similar to n_peptide_hits=1, but if there are two or more highest scoring hits, none are kept.
  </li>
 </ul>

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

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerTOPPSubsection_("precursor", "Filtering by precursor RT or m/z");
    registerStringOption_("precursor:rt", "[min]:[max]", ":", "Retention time range to extract.", false);
    registerStringOption_("precursor:mz", "[min]:[max]", ":", "Mass-to-charge range to extract.", false);

    registerTOPPSubsection_("score", "Filtering by peptide/protein score. To enable any of the filters below, just change their default value. All active filters will be applied in order.");
    registerDoubleOption_("score:pep", "<score>", 0, "The score which should be reached by a peptide hit to be kept. The score is dependent on the most recent(!) preprocessing - it could be Mascot scores (if a MascotAdapter was applied before), or an FDR (if FalseDiscoveryRate was applied before), etc.", false);
    registerDoubleOption_("score:prot", "<score>", 0, "The score which should be reached by a protein hit to be kept. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides.", false);
    registerTOPPSubsection_("thresh", "Filtering by significance threshold");
    registerDoubleOption_("thresh:pep", "<fraction>", 0.0, "Keep a peptide hit only if its score is above this fraction of the peptide significance threshold.", false);
    registerDoubleOption_("thresh:prot", "<fraction>", 0.0, "Keep a protein hit only if its score is above this fraction of the protein significance threshold. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides.", false);

    registerTOPPSubsection_("whitelist", "Filtering by whitelisting (only instances also present in a whitelist file can pass)");
    registerInputFile_("whitelist:proteins", "<file>", "", "Filename of a FASTA file containing protein sequences.\n"
                                                           "All peptides that are not referencing a protein in the file are removed.\n"
                                                           "All proteins whose accession is not present in this file are removed.", false);
    setValidFormats_("whitelist:proteins", ListUtils::create<String>("fasta"));

    registerStringList_("whitelist:protein_accessions", "<accessions>", ListUtils::create<String>(""), "All peptides that are not referencing at least one of the provided protein accession are removed.\nOnly proteins of the provided list are retained.", false);

    registerTOPPSubsection_("blacklist", "Filtering by blacklisting (only instances not present in a blacklist file can pass)");
    registerInputFile_("blacklist:peptides", "<file>", "", "Peptides having the same sequence and modification assignment as any peptide in this file will be filtered out. Use with blacklist:ignore_modifications flag to only compare by sequence.\n", false);
    setValidFormats_("blacklist:peptides", ListUtils::create<String>("idXML"));
    registerFlag_("blacklist:ignore_modifications", "Compare blacklisted peptides by sequence only.\n", false);

    registerTOPPSubsection_("rt", "Filtering by RT predicted by 'RTPredict'");
    registerDoubleOption_("rt:p_value", "<float>", 0.0, "Retention time filtering by the p-value predicted by RTPredict.", false);
    registerDoubleOption_("rt:p_value_1st_dim", "<float>", 0.0, "Retention time filtering by the p-value predicted by RTPredict for first dimension.", false);
    setMinFloat_("rt:p_value", 0);
    setMaxFloat_("rt:p_value", 1);
    setMinFloat_("rt:p_value_1st_dim", 0);
    setMaxFloat_("rt:p_value_1st_dim", 1);

    registerTOPPSubsection_("mz", "Filtering by mass error");
    registerDoubleOption_("mz:error", "<float>", -1, "Filtering by deviation to theoretical mass (disabled for negative values).", false);
    registerStringOption_("mz:unit", "<String>", "ppm", "Absolute or relative error.", false);
    setValidStrings_("mz:unit", ListUtils::create<String>("Da,ppm"));

    registerTOPPSubsection_("best", "Filtering best hits per spectrum (for peptides) or from proteins");
    registerIntOption_("best:n_peptide_hits", "<integer>", 0, "Keep only the 'n' highest scoring peptide hits per spectrum (for n>0).", false);
    setMinInt_("best:n_peptide_hits", 0);
    registerIntOption_("best:n_protein_hits", "<integer>", 0, "Keep only the 'n' highest scoring protein hits (for n > 0).", false);
    setMinInt_("best:n_protein_hits", 0);
    registerFlag_("best:strict", "Keep only the highest scoring peptide hit.\n"
                                 "Similar to n_peptide_hits=1, but if there are two or more highest scoring hits, none are kept.");
    registerStringOption_("best:n_to_m_peptide_hits", "[min]:[max]", ":", "peptide hit rank range to extracts", false, true);

    registerStringOption_("length", "[min]:[max]", ":", "Keep only peptide hits with a sequence length in this range.", false);

    registerStringOption_("charge", "[min]:[max]", ":", "Keep only peptide hits with charge states in this range.", false);

    registerFlag_("var_mods", "Keep only peptide hits with variable modifications (fixed modifications from SearchParameters will be ignored).", false);

    registerFlag_("unique", "If a peptide hit occurs more than once per peptide ID, only one instance is kept.");
    registerFlag_("unique_per_protein", "Only peptides matching exactly one protein are kept. Remember that isoforms count as different proteins!");
    registerFlag_("keep_unreferenced_protein_hits", "Proteins not referenced by a peptide are retained in the ids.");
    registerFlag_("remove_decoys", "Remove proteins according to the information in the user parameters. Usually used in combination with 'delete_unreferenced_peptide_hits'.");
    registerFlag_("delete_unreferenced_peptide_hits", "Peptides not referenced by any protein are deleted in the IDs. Usually used in combination with 'score:prot' or 'thresh:prot'.");

  }


  ExitCodes main_(int, const char**)
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
      LOG_INFO << "Filtering by removing duplicate peptide hits..." << endl;
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

    String fasta_name = getStringOption_("whitelist:proteins").trim();
    if (!fasta_name.empty())
    {
      LOG_INFO << "Filtering by protein whitelisting (FASTA input)..." << endl;
      // load protein accessions from FASTA file:
      vector<FASTAFile::FASTAEntry> fasta;
      FASTAFile().load(fasta_name, fasta);
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

    String blacklist_name = getStringOption_("blacklist:peptides").trim();
    if (!blacklist_name.empty())
    {
      LOG_INFO << "Filtering by exclusion peptide blacklisting..." << endl;
      vector<PeptideIdentification> exclusion_peptides;
      vector<ProteinIdentification> exclusion_proteins; // ignored
      IdXMLFile().load(blacklist_name, exclusion_proteins, exclusion_peptides);
      bool ignore_mods = getFlag_("blacklist:ignore_modifications");
      IDFilter::removePeptidesWithMatchingSequences(
        peptides, exclusion_peptides, ignore_mods);
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

    // vector<String> mods = getStringList_("modifications");
    // if (!mods.empty())
    // {
    //   set<String> good_mods(mods.begin(), mods.end());
    //   IDFilter::keepPeptidesWithMatchingModifications(peptides, good_mods);
    // }

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
      LOG_INFO << "Filtering by best n peptide hits...\n" << endl;
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
      LOG_INFO << "Filtering by best n protein hits...\n" << endl;
      IDFilter::keepNBestHits(proteins, best_n_prot);
    }
        
    if (getFlag_("remove_decoys"))
    {
      IDFilter::removeDecoyHits(peptides);
      IDFilter::removeDecoyHits(proteins);
    }


    // Clean-up:

    if (!getFlag_("keep_unreferenced_protein_hits"))
    {
      LOG_INFO << "Filtering unreferenced protein hits..." << endl;
      IDFilter::removeUnreferencedProteins(proteins, peptides);
    }

    IDFilter::updateHitRanks(proteins);
    IDFilter::updateHitRanks(peptides);

    // remove non-existant protein references from peptides (and optionally:
    // remove peptides with no proteins):
    bool rm_pep = getFlag_("delete_unreferenced_peptide_hits");
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
