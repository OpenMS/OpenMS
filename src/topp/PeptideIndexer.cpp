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
// $Authors: Andreas Bertsch, Chris Bielow, Knut Reinert $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_PeptideIndexer PeptideIndexer

    @brief Refreshes the protein references for all peptide hits from an idXML file and adds target/decoy information.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ PeptideIndexer \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
        </tr>
    </table>
</CENTER>

  All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy", depending on whether the protein accession contains the decoy pattern (parameter @p decoy_string) as a suffix or prefix, respectively (see parameter @p prefix). For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins, only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool. (For FDR calculations, "target+decoy" peptide hits count as target hits.)

  @note Make sure that your protein names in the database contain a correctly formatted decoy string. This can be ensured by using @ref UTILS_DecoyDatabase.
        If the decoy identifier is not recognized successfully all proteins will be assumed to stem from the target-part of the query.<br>
        E.g., "sw|P33354_REV|YEHR_ECOLI Uncharacterized lipop..." is <b>invalid</b>, since the tool has no knowledge of how SwissProt entries are build up.
        A correct identifier could be "rev_sw|P33354|YEHR_ECOLI Uncharacterized li ..." or "sw|P33354|YEHR_ECOLI_rev Uncharacterized li", depending on whether you are
        using prefix annotation or not.<br>
        This tool will also report some helpful target/decoy statistics when it is done.

  PeptideIndexer supports relative database filenames, which (when not found in the current working directory) are looked up in the directories specified 
  by @p OpenMS.ini:id_db_dir (see @subpage TOPP_advanced).

  By default this tool will fail if an unmatched peptide occurs, i.e. if the database does not contain the corresponding protein.
  You can force it to return successfully in this case by using the flag @p allow_unmatched.

  Some search engines (such as Mascot) will replace ambiguous amino acids ('B', 'Z', and 'X') in the protein database with unambiguous amino acids in the reported peptides, e.g. exchange 'X' with 'H'.
  This will cause such peptides to not be found by exactly matching their sequences to the database. However, we can recover these cases by using tolerant search (done automatically).

  Two search modes are available:
    - exact: [default mode] Peptide sequences require exact match in protein database.
             If at least one protein hit is found, no tolerant search is used for this peptide.
             If no protein for this peptide can be found, tolerant matching is automatically used for this peptide.
    - tolerant:
             Allow ambiguous amino acids in protein sequence, e.g., 'M' in peptide will match 'X' in protein.
             This mode might yield more protein hits for some peptides (those that contain ambiguous amino acids).
             Tolerant search also allows for real sequence mismatches (see 'mismatches_max'), in case you want to find related proteins which
             might be the origin of a peptide if it had a SNP for example. Runtime increase is moderate when allowing a single
             mismatch, but rises drastically for two or more.

  The exact mode is much faster (about 10 times) and consumes less memory (about 2.5 times),
  but might fail to report a few protein hits with ambiguous amino acids for some peptides. Usually these proteins are putative, however.
  The exact mode also supports usage of multiple threads (@p threads option) to speed up computation even further, at the cost of some memory.
  If tolerant searching needs to be done for unassigned peptides,
  the latter will consume the major share of the runtime.
  Independent of whether exact or tolerant search is used, we require ambiguous amino acids in peptide sequences to match exactly in the protein DB (i.e. 'X' in a peptide only matches 'X' in the database).

  Leucine/Isoleucine:
  Further complications can arise due to the presence of the isobaric amino acids isoleucine ('I') and leucine ('L') in protein sequences. 
  Since the two have the exact same chemical composition and mass, they generally cannot be distinguished by mass spectrometry. 
  If a peptide containing 'I' was reported as a match for a spectrum, a peptide containing 'L' instead would be an equally good match (and vice versa). 
  To account for this inherent ambiguity, setting the flag @p IL_equivalent causes 'I' and 'L' to be considered as indistinguishable.@n
  For example, if the sequence "PEPTIDE" (matching "Protein1") was identified as a search hit, 
  but the database additionally contained "PEPTLDE" (matching "Protein2"), running PeptideIndexer with the @p IL_equivalent option would 
  report both "Protein1" and "Protein2" as accessions for "PEPTIDE".
  (This is independent of the error-tolerant search controlled by @p full_tolerant_search and @p aaa_max.)

  Enzyme specificity:
  Once a peptide sequence is found in a protein sequence, this does <b>not</b> imply that the hit is valid! This is where enzyme specificity comes into play.
  By default, we demand that the peptide is fully tryptic (i.e. the enzyme parameter is set to "trypsin" and specificity is "full").
  So unless the peptide coincides with C- and/or N-terminus of the protein, the peptide's cleavage pattern should fulfill the trypsin cleavage rule [KR][^P].
  We make one exception for peptides starting at the second amino acid of a protein if the first amino acid of that protein is methionine (M), 
  which is usually cleaved off in vivo. For example, the two peptides AAAR and MAAAR would both match a protein starting with MAAAR.
  You can relax the requirements further by choosing <tt>semi-tryptic</tt> (only one of two "internal" termini must match requirements) 
  or <tt>none</tt> (essentially allowing all hits, no matter their context).

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PeptideIndexer.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_PeptideIndexer.html
*/

struct PeptideProteinMatchInformation
{
  /// index of the protein the peptide is contained in
  OpenMS::Size protein_index;

  /// the amino acid after the peptide in the protein
  char AABefore;

  /// the amino acid before the peptide in the protein
  char AAAfter;

  /// the position of the peptide in the protein
  OpenMS::Int position;

  bool operator<(const PeptideProteinMatchInformation& other) const
  {
    if (protein_index != other.protein_index)
    {
      return protein_index < other.protein_index;
    }
    else if (position != other.position)
    {
      return position < other.position;
    }
    else if (AABefore != other.AABefore)
    {
      return AABefore < other.AABefore;
    }
    else if (AAAfter != other.AAAfter)
    {
      return AAAfter < other.AAAfter;
    }
    return false;
  }

  bool operator==(const PeptideProteinMatchInformation& other) const
  {
    return protein_index == other.protein_index &&
           position == other.position &&
           AABefore == other.AABefore &&
           AAAfter == other.AAAfter;
  }

};

namespace seqan
{

  struct FoundProteinFunctor
  {
public:
    typedef OpenMS::Map<OpenMS::Size, std::set<PeptideProteinMatchInformation> > MapType;

    /// peptide index --> protein indices
    MapType pep_to_prot;

    /// number of accepted hits (passing addHit() constraints)
    OpenMS::Size filter_passed;

    /// number of rejected hits (not passing addHit())
    OpenMS::Size filter_rejected;

private:
    EnzymaticDigestion enzyme_;

public:
    explicit FoundProteinFunctor(const EnzymaticDigestion& enzyme) :
      pep_to_prot(), filter_passed(0), filter_rejected(0), enzyme_(enzyme)
    {
    }

    template <typename TIter1, typename TIter2>
    void operator()(const TIter1& iter_pep, const TIter2& iter_prot)
    {
      // the peptide sequence (will not change)
      const OpenMS::String tmp_pep(begin(representative(iter_pep)),
                                   end(representative(iter_pep)));

      // remember mapping of proteins to peptides and vice versa
      const OpenMS::Size count_occ = countOccurrences(iter_pep);
      for (OpenMS::Size i_pep = 0; i_pep < count_occ; ++i_pep)
      {
        const OpenMS::Size idx_pep = getOccurrences(iter_pep)[i_pep].i1;
        const OpenMS::Size count_occ_prot = countOccurrences(iter_prot);
        for (OpenMS::Size i_prot = 0; i_prot < count_occ_prot; ++i_prot)
        {
          const seqan::Pair<int> prot_occ = getOccurrences(iter_prot)[i_prot];
          // the protein sequence (will change for every Occurrence -- hitting
          // multiple proteins)
          const OpenMS::String tmp_prot(
            begin(indexText(container(iter_prot))[getSeqNo(prot_occ)]),
            end(indexText(container(iter_prot))[getSeqNo(prot_occ)]));
          // check if hit is valid and add (if valid)
          addHit(idx_pep, prot_occ.i1, tmp_pep, tmp_prot,
                 getSeqOffset(prot_occ));
        }
      }
    }

    void addHit(OpenMS::Size idx_pep, OpenMS::Size idx_prot,
                const OpenMS::String& seq_pep, const OpenMS::String& protein,
                OpenMS::Size position)
    {
      if (enzyme_.isValidProduct(AASequence::fromString(protein), position,
                                 seq_pep.length()))
      {
        PeptideProteinMatchInformation match;
        match.protein_index = idx_prot;
        match.position = position;
        match.AABefore = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : protein[position - 1];
        match.AAAfter = (position + seq_pep.length() >= protein.size()) ? PeptideEvidence::C_TERMINAL_AA : protein[position + seq_pep.length()];
        pep_to_prot[idx_pep].insert(match);
        ++filter_passed;
      }
      else
      {
        // LOG_WARN << "Peptide " << seq_pep << " is not a valid hit to protein "
        //          << protein << " at position " << position << std::endl;
        ++filter_rejected;
      }
    }

    bool operator==(const FoundProteinFunctor& rhs) const
    {
      if (pep_to_prot.size() != rhs.pep_to_prot.size())
      {
        LOG_ERROR << "Size " << pep_to_prot.size() << " "
                  << rhs.pep_to_prot.size() << std::endl;
        return false;
      }

      MapType::const_iterator it1 = pep_to_prot.begin();
      MapType::const_iterator it2 = rhs.pep_to_prot.begin();
      while (it1 != pep_to_prot.end())
      {
        if (it1->first != it2->first)
        {
          LOG_ERROR << "Index of " << it1->first << " " << it2->first
                    << std::endl;
          return false;
        }
        if (it1->second.size() != it2->second.size())
        {
          LOG_ERROR << "Size of " << it1->first << " " << it1->second.size()
                    << "--" << it2->second.size() << std::endl;
          return false;
        }
        if (!equal(it1->second.begin(), it1->second.end(), it2->second.begin()))
        {
          LOG_ERROR << "not equal set for " << it1->first << std::endl;
          return false;
        }
        ++it1;
        ++it2;
      }
      return true;
    }

    bool operator!=(const FoundProteinFunctor& rhs) const
    {
      return !(*this == rhs);
    }

  };


  // saving some memory for the SA
  template <>
  struct SAValue<Index<StringSet<Peptide>, IndexWotd<> > >
  {
    typedef Pair<unsigned> Type;
  };

  template <typename T = void>
  struct EquivalenceClassAA_
  {
    static unsigned const VALUE[24];
  };
  template <typename T>
  unsigned const EquivalenceClassAA_<T>::VALUE[24] =
  {
    1, // 0 Ala Alanine (A)
    2, // 1 Arg Arginine (R)
    4, // 2 Asn Asparagine (N)
    8, // 3 Asp Aspartic Acid (D)
    16, // 4 Cys Cystine (C)
    32, // 5 Gln Glutamine (Q)
    64, // 6 Glu Glutamic Acid (E)
    128, // 7 Gly Glycine (G)
    256, // 8 His Histidine (H)
    512, // 9 Ile Isoleucine (I)
    1024, // 10 Leu Leucine (L)
    2048, // 11 Lys Lysine (K)
    4096, // 12 Met Methionine (M)
    8192, // 13 Phe Phenylalanine (F)
    16384, // 14 Pro Proline (P)
    32768, // 15 Ser Serine (S)
    65536, // 16 Thr Threonine (T)
    131072, // 17 Trp Tryptophan (W)
    262144, // 18 Tyr Tyrosine (Y)
    524288, // 19 Val Valine (V)
    // ambiguous AA's
    4+8, //  Aspartic Acid (D), Asparagine(N) == (B)
    32+64, // Glutamic Acid(E), Glutamine(Q) == (Z)
     static_cast<unsigned>(-1), // 22 Unknown (matches ALL)
    static_cast<unsigned>(-1), // 23 Terminator (dummy)
  };


  template <bool enumerateA, bool enumerateB, typename TOnFoundFunctor,
            typename TTreeIteratorA, typename TIterPosA,
            typename TTreeIteratorB, typename TIterPosB, typename TErrors>
  inline void _approximateAminoAcidTreeSearch(TOnFoundFunctor& onFoundFunctor,
                                              TTreeIteratorA iterA,
                                              TIterPosA iterPosA,
                                              TTreeIteratorB iterB_,
                                              TIterPosB iterPosB,
                                              TErrors errorsLeft, // usually 0 for our case
                                              TErrors classErrorsLeft) // ambiguous AA's allowed
  {
    if (enumerateA && !goDown(iterA)) return;

    if (enumerateB && !goDown(iterB_)) return;

    do
    {
      TTreeIteratorB iterB = iterB_;
      do
      {
        TErrors e = errorsLeft;
        TErrors ec = classErrorsLeft;
        TIterPosA ipA = iterPosA;
        TIterPosB ipB = iterPosB;

        while (true)
        {
          if (ipA == repLength(iterA))
          {
            if (isLeaf(iterA))
            {
              onFoundFunctor(iterA, iterB);
            }
            else if (ipB == repLength(iterB) && !isLeaf(iterB))
            {
              _approximateAminoAcidTreeSearch<true, true>(
                onFoundFunctor, iterA, ipA, iterB, ipB, e, ec);
            }
            else
            {
              _approximateAminoAcidTreeSearch<true, false>(
                onFoundFunctor, iterA, ipA, iterB, ipB, e, ec);
            }
            break;
          }
          else
          {
            if (ipB == repLength(iterB))
            {
              if (!isLeaf(iterB))
              {
                _approximateAminoAcidTreeSearch<false, true>(
                  onFoundFunctor, iterA, ipA, iterB, ipB, e, ec);
              }
              break;
            }
          }

          if (_charComparator(representative(iterA)[ipA],
                              representative(iterB)[ipB],
                              EquivalenceClassAA_<char>::VALUE))
          {
            // matched (including character classes) - look at ambiguous AA in
            // PROTEIN tree (peptide tree is not considered!)
            const char x_prot = representative(iterB)[ipB];
            if ((x_prot == 'X') || (x_prot == 'B') || (x_prot == 'Z'))
            {
              if (ec == 0) break;
              --ec; // decrease class error tokens
            }

            // dealing with 'X' in peptide sequence: only match exactly 'X' in
            // proteinDB, not just any representative of 'X' (this is how X!
            // Tandem would report results)
            const char x_pep = representative(iterA)[ipA];
            if ((x_pep == 'X') || (x_pep == 'B') || (x_pep == 'Z'))
            {
              if (x_pep != x_prot) break;
            }
          }
          else
          { // real mismatches (e.g., when checking for SNP's)
            if (e == 0) break;
            --e;
          }

          ++ipA;
          ++ipB;
        }
      }
      while (enumerateB && goRight(iterB));
    }
    while (enumerateA && goRight(iterA));
  }

  template <typename TEquivalenceTable>
  inline bool _charComparator(AminoAcid charA, AminoAcid charB,
                              TEquivalenceTable equivalence)
  {
    const unsigned a_index = ordValue(charA);
    const unsigned b_index = ordValue(charB);
    return (equivalence[a_index] & equivalence[b_index]) != 0;
  }

}


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeptideIndexer :
  public TOPPBase
{
public:
  TOPPPeptideIndexer() :
    TOPPBase("PeptideIndexer",
             "Refreshes the protein references for all peptide hits.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input idXML file containing the identifications.");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerInputFile_("fasta", "<file>", "", "Input sequence database in FASTA format. Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "Output idXML file.");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerStringOption_("decoy_string", "<string>", "_rev", "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.", false);
    registerFlag_("prefix", "If set, protein accessions in the database contain 'decoy_string' as prefix.");
    registerStringOption_("missing_decoy_action", "<action>", "error", "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message)", false);
    setValidStrings_("missing_decoy_action", ListUtils::create<String>("error,warn"));

    registerTOPPSubsection_("enzyme", "The enzyme determines valid cleavage sites; cleavage specificity determines to what extent validity is enforced.");

    registerStringOption_("enzyme:name", "", EnzymaticDigestion::NamesOfEnzymes[0], "Enzyme which determines valid cleavage sites - e.g. trypsin cleaves after lysine (K) or arginine (R), but not before proline (P).", false);
    StringList enzymes;
    enzymes.assign(EnzymaticDigestion::NamesOfEnzymes, EnzymaticDigestion::NamesOfEnzymes + EnzymaticDigestion::SIZE_OF_ENZYMES);
    setValidStrings_("enzyme:name", enzymes);

    registerStringOption_("enzyme:specificity", "", EnzymaticDigestion::NamesOfSpecificity[0], "Specificity of the enzyme."
                          "\n  '" + EnzymaticDigestion::NamesOfSpecificity[0] + "': both internal cleavage sites must match."
                          "\n  '" + EnzymaticDigestion::NamesOfSpecificity[1] + "': one of two internal cleavage sites must match."
                          "\n  '" + EnzymaticDigestion::NamesOfSpecificity[2] + "': allow all peptide hits no matter their context. Therefore, the enzyme chosen does not play a role here", false);
    StringList spec;
    spec.assign(EnzymaticDigestion::NamesOfSpecificity, EnzymaticDigestion::NamesOfSpecificity + EnzymaticDigestion::SIZE_OF_SPECIFICITY);
    setValidStrings_("enzyme:specificity", spec);

    registerFlag_("write_protein_sequence", "If set, the protein sequences are stored as well.");
    registerFlag_("write_protein_description", "If set, the protein description is stored as well.");
    registerFlag_("keep_unreferenced_proteins", "If set, protein hits which are not referenced by any peptide are kept.");
    registerFlag_("allow_unmatched", "If set, unmatched peptide sequences are allowed. By default (i.e. if this flag is not set) the program terminates with an error on unmatched peptides.");
    registerFlag_("full_tolerant_search", "If set, all peptide sequences are matched using tolerant search. Thus potentially more proteins (containing ambiguous amino acids) are associated. This is much slower!");
    registerIntOption_("aaa_max", "<number>", 4, "[tolerant search only] Maximal number of ambiguous amino acids (AAA) allowed when matching to a protein database with AAA's. AAA's are 'B', 'Z' and 'X'", false);
    setMinInt_("aaa_max", 0);
    registerIntOption_("mismatches_max", "<number>", 0, "[tolerant search only] Maximal number of real mismatches (will be used after checking for ambiguous AA's (see 'aaa_max' option). In general this param should only be changed if you want to look for other potential origins of a peptide which might have unknown SNPs or alike.", false);
    setMinInt_("mismatches_max", 0);
    registerFlag_("IL_equivalent", "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable)");
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    bool write_protein_sequence = getFlag_("write_protein_sequence");
    bool write_protein_description = getFlag_("write_protein_description");
    bool keep_unreferenced_proteins = getFlag_("keep_unreferenced_proteins");
    bool allow_unmatched = getFlag_("allow_unmatched");
    bool il_equivalent = getFlag_("IL_equivalent");

    String decoy_string = getStringOption_("decoy_string");
    bool prefix = getFlag_("prefix");

    String db_name = getStringOption_("fasta");
    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }

    EnzymaticDigestion enzyme;
    enzyme.setEnzyme(enzyme.getEnzymeByName(getStringOption_("enzyme:name")));
    enzyme.setSpecificity(enzyme.getSpecificityByName(getStringOption_("enzyme:specificity")));


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    // we stream the Fasta file
    vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(db_name, proteins);

    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;

    IdXMLFile().load(in, prot_ids, pep_ids);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    if (proteins.size() == 0) // we do not allow an empty database
    {
      LOG_ERROR << "Error: An empty FASTA file was provided. Mapping makes no sense. Aborting..." << std::endl;
      return INPUT_FILE_EMPTY;
    }

    if (pep_ids.size() == 0) // Aho-Corasick requires non-empty input
    {
      LOG_WARN << "Warning: An empty idXML file was provided. Output will be empty as well." << std::endl;
      if (!getFlag_("keep_unreferenced_proteins"))
      {
        prot_ids.clear();
      }
      IdXMLFile().store(out, prot_ids, pep_ids);
      return EXECUTION_OK;
    }

    writeDebug_("Collecting peptides...", 1);

    seqan::FoundProteinFunctor func(enzyme); // stores the matches (need to survive local scope which follows)
    Map<String, Size> acc_to_prot; // build map: accessions to FASTA protein index

    { // new scope - forget data after search

      /**
       BUILD Protein DB
      */
      seqan::StringSet<seqan::Peptide> prot_DB;

      bool has_DB_duplicates(false);

      for (Size i = 0; i != proteins.size(); ++i)
      {
        String seq = proteins[i].sequence.remove('*');
        if (il_equivalent)
        { // convert  L to I; warning: do not use 'J', since Seqan does not know about it and will convert 'J' to 'X'
          seq.substitute('L', 'I');
        }

        
        String acc = proteins[i].identifier;
        // check for duplicate proteins
        if (acc_to_prot.has(acc))
        {
          LOG_WARN << "PeptideIndexer: Warning, protein identifiers should be unique to a database. Identifier '" << acc << "' found multiple times.\n";
          has_DB_duplicates = true;
          // check if sequence is identical
          const seqan::Peptide& tmp_prot = prot_DB[acc_to_prot[acc]];
          if (String(begin(tmp_prot), end(tmp_prot)) != seq)
          {
            LOG_ERROR << "PeptideIndexer: protein identifier '" << acc << "' found multiple times with different sequences" << (il_equivalent ? " (I/L substituted)" : "") 
                      << ":\n" << tmp_prot << "\nvs.\n" << seq << "\n! Please fix the database and run PeptideIndexer again!" << std::endl;
            return INPUT_FILE_CORRUPT;
          }
          // remove duplicate sequence from 'proteins', since 'prot_DB' and 'proteins' need to correspond 1:1 (later indexing depends on it)
          // The other option would be to allow two identical entries, but later on, only the last one will be reported (making the first protein an orphan; implementation details below)
          // Thus, the only safe option is to remove the duplicate from 'proteins' and not to add it to 'prot_DB'
          proteins.erase(proteins.begin()+i);
          // try this index again in the next loop (--i is save since this condition is met only when i>0)
          --i;
        } 
        else
        {
          // extend protein DB
          seqan::appendValue(prot_DB, seq.c_str());
          acc_to_prot[acc] = i;
        }
        
      }
      // make sure the warnings above are printed to screen
      if (has_DB_duplicates) LOG_WARN << std::endl;

      /**
        BUILD Peptide DB
      */
      seqan::StringSet<seqan::Peptide> pep_DB;
      for (vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
      {
        //String run_id = it1->getIdentifier();
        vector<PeptideHit> hits = it1->getHits();
        for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
        {
          String seq = it2->getSequence().toUnmodifiedString().remove('*');
          if (il_equivalent)
          { // convert  L to I; warning: do not use 'J', since Seqan does not know about it and will convert 'J' to 'X'
            seq.substitute('L', 'I');
          }
          appendValue(pep_DB, seq.c_str());
        }
      }

      writeLog_(String("Mapping ") + length(pep_DB) + " peptides to " + length(prot_DB) + " proteins.");

      /** first, try Aho Corasick (fast) -- using exact matching only */
      bool SA_only = getFlag_("full_tolerant_search");
      UInt max_mismatches = getIntOption_("mismatches_max");
      
      if (!SA_only && (max_mismatches > 0))
      { // this combination is not allowed, and we want the user to make a conscious decision about it
        LOG_ERROR << "Error: Exact matching in combination with #mismatches > 0 is not allowed.\n"
                  << "       Either use full tolerant search ('full_tolerant_search') or set 'mismatches_max' back to '0'."
                  << "       Aborting..." << std::endl;
        return ILLEGAL_PARAMETERS;
      };

      if (!SA_only)
      {
        StopWatch sw;
        sw.start();
        SignedSize protDB_length = (SignedSize) length(prot_DB);
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          seqan::Pattern<seqan::StringSet<seqan::Peptide>, seqan::AhoCorasick> pattern(pep_DB);
          seqan::FoundProteinFunctor func_threads(enzyme);
          writeDebug_("Finding peptide/protein matches ...", 1);

#pragma omp for
          for (SignedSize i = 0; i < protDB_length; ++i)
          {
            seqan::Finder<seqan::Peptide> finder(prot_DB[i]);
            while (find(finder, pattern))
            {
              //seqan::appendValue(pat_hits, seqan::Pair<Size, Size>(position(pattern), position(finder)));

              //func_threads.pep_to_prot[position(pattern)].insert(i);
              // String(seqan::String<char, seqan::CStyle>(prot_DB[i])), position(finder))
              // target.assign(begin(source, Standard()), end(source, Standard()));
              const seqan::Peptide& tmp_pep = pep_DB[position(pattern)];
              const seqan::Peptide& tmp_prot = prot_DB[i];

              func_threads.addHit(position(pattern), i, String(begin(tmp_pep), end(tmp_pep)), String(begin(tmp_prot), end(tmp_prot)), position(finder));
            }
          }

          // join results again
#ifdef _OPENMP
#pragma omp critical(PeptideIndexer_joinAC)
#endif
          {
            func.filter_passed += func_threads.filter_passed;
            func.filter_rejected += func_threads.filter_rejected;
            for (seqan::FoundProteinFunctor::MapType::const_iterator it = func_threads.pep_to_prot.begin(); it != func_threads.pep_to_prot.end(); ++it)
            {
              func.pep_to_prot[it->first].insert(func_threads.pep_to_prot[it->first].begin(), func_threads.pep_to_prot[it->first].end());
            }

          }
        } // end parallel

        sw.stop();

        writeLog_(String("Aho-Corasick done. Found ") + func.filter_passed + " hits in " + func.pep_to_prot.size() + " of " + length(pep_DB) + " peptides (time: " + sw.getClockTime() + " s (wall), " + sw.getCPUTime() + " s (CPU)).");
      }

      /// check if every peptide was found:
      if (func.pep_to_prot.size() != length(pep_DB))
      {
        // search using SA, which supports mismatches (introduced by resolving ambiguous AA's by e.g. Mascot) -- expensive!
        writeLog_(String("Using suffix array to find ambiguous matches..."));

        // search peptides which remained unidentified during Aho-Corasick (might be all if 'full_tolerant_search' is enabled)
        seqan::StringSet<seqan::Peptide> pep_DB_SA;
        Map<Size, Size> missed_pep;
        for (Size p = 0; p < length(pep_DB); ++p)
        {
          if (!func.pep_to_prot.has(p))
          {
            missed_pep[length(pep_DB_SA)] = p;
            appendValue(pep_DB_SA, pep_DB[p]);
          }
        }

        writeLog_(String("... for ") + length(pep_DB_SA) + " peptide(s).");

        seqan::FoundProteinFunctor func_SA(enzyme);

        typedef seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::IndexWotd<> > TIndex;
        TIndex prot_Index(prot_DB);
        TIndex pep_Index(pep_DB_SA);

        // use only full peptides in Suffix Array
        const Size length_SA = length(pep_DB_SA);
        resize(indexSA(pep_Index), length_SA);
        for (Size i = 0; i < length_SA; ++i)
        {
          indexSA(pep_Index)[i].i1 = (unsigned)i;
          indexSA(pep_Index)[i].i2 = 0;
        }

        typedef seqan::Iterator<TIndex, seqan::TopDown<seqan::PreorderEmptyEdges> >::Type TTreeIter;

        //seqan::open(indexSA(prot_Index), "c:\\tmp\\prot_Index.sa");
        //seqan::open(indexDir(prot_Index), "c:\\tmp\\prot_Index.dir");

        TTreeIter prot_Iter(prot_Index);
        TTreeIter pep_Iter(pep_Index);

        UInt max_aaa = getIntOption_("aaa_max");
        seqan::_approximateAminoAcidTreeSearch<true, true>(func_SA, pep_Iter, 0u, prot_Iter, 0u, max_mismatches, max_aaa);

        // augment results with SA hits
        func.filter_passed += func_SA.filter_passed;
        func.filter_rejected += func_SA.filter_rejected;
        for (seqan::FoundProteinFunctor::MapType::const_iterator it = func_SA.pep_to_prot.begin(); it != func_SA.pep_to_prot.end(); ++it)
        {
          func.pep_to_prot[missed_pep[it->first]] = it->second;
        }

      }

    } // end local scope

    // write some stats
    LOG_INFO << "Peptide hits passing enzyme filter: " << func.filter_passed << "\n"
             << "     ... rejected by enzyme filter: " << func.filter_rejected << std::endl;

    /* do mapping */
    writeDebug_("Reindexing peptide/protein matches...", 1);

    /// index existing proteins
    Map<String, Size> runid_to_runidx; // identifier to index
    for (Size run_idx = 0; run_idx < prot_ids.size(); ++run_idx)
    {
      runid_to_runidx[prot_ids[run_idx].getIdentifier()] = run_idx;
    }

    /// store target/decoy status of proteins
    Map<String, bool> protein_is_decoy; // accession -> is decoy?

    /// for peptides --> proteins
    Size stats_matched_unique(0);
    Size stats_matched_multi(0);
    Size stats_unmatched(0);
    Size stats_count_m_t(0);
    Size stats_count_m_d(0);
    Size stats_count_m_td(0);
    Map<Size, set<Size> > runidx_to_protidx; // in which protID do appear which proteins (according to mapped peptides)

    Size pep_idx(0);
    for (vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
    {
      // which ProteinIdentification does the peptide belong to?
      Size run_idx = runid_to_runidx[it1->getIdentifier()];

      vector<PeptideHit> hits = it1->getHits();

      for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
      {
        // clear protein accessions
        it2->setPeptideEvidences(vector<PeptideEvidence>());

        // add new protein references
        for (set<PeptideProteinMatchInformation>::const_iterator it_i = func.pep_to_prot[pep_idx].begin();
             it_i != func.pep_to_prot[pep_idx].end();
             ++it_i)
        {
          const String& accession = proteins[it_i->protein_index].identifier;
          PeptideEvidence pe;
          pe.setProteinAccession(accession);
          pe.setStart(it_i->position);
          pe.setEnd(it_i->position + it2->getSequence().size() - 1);
          pe.setAABefore(it_i->AABefore);
          pe.setAAAfter(it_i->AAAfter);
          it2->addPeptideEvidence(pe);

          runidx_to_protidx[run_idx].insert(it_i->protein_index); // fill protein hits

          if (!protein_is_decoy.has(accession))
          {
            protein_is_decoy[accession] = (prefix && accession.hasPrefix(decoy_string)) || (!prefix && accession.hasSuffix(decoy_string));
          }
        }

        ///
        /// is this a decoy hit?
        ///
        bool matches_target(false);
        bool matches_decoy(false);

        set<String> protein_accessions = it2->extractProteinAccessions();
        for (set<String>::const_iterator it = protein_accessions.begin(); it != protein_accessions.end(); ++it)
        {
          if (protein_is_decoy[*it])
          {
            matches_decoy = true;
          }
          else
          {
            matches_target = true;
          }
          // this is rare in practice, so the test may not really save time:
          // if (matches_decoy && matches_target)
          // {
          //   break; // no need to check remaining accessions
          // }
        }
        String target_decoy;
        if (matches_decoy && matches_target)
        {
          target_decoy = "target+decoy";
          ++stats_count_m_td;
        }
        else if (matches_target)
        {
          target_decoy = "target";
          ++stats_count_m_t;
        }
        else if (matches_decoy)
        {
          target_decoy = "decoy";
          ++stats_count_m_d;
        }
        it2->setMetaValue("target_decoy", target_decoy);

        if (protein_accessions.size() == 1)
        {
          it2->setMetaValue("protein_references", "unique");
          ++stats_matched_unique;
        }
        else if (protein_accessions.size() > 1)
        {
          it2->setMetaValue("protein_references", "non-unique");
          ++stats_matched_multi;
        }
        else
        {
          it2->setMetaValue("protein_references", "unmatched");
          ++stats_unmatched;
          if (stats_unmatched < 5) LOG_INFO << "Unmatched peptide: " << it2->getSequence() << "\n";
          else if (stats_unmatched == 5) LOG_INFO << "Unmatched peptide: ...\n";
        }

        ++pep_idx; // next hit
      }
      it1->setHits(hits);
    }

    LOG_INFO << "Statistics of peptides (target/decoy):\n";
    LOG_INFO << "  match to target DB only: " << stats_count_m_t << "\n";
    LOG_INFO << "  match to decoy DB only : " << stats_count_m_d << "\n";
    LOG_INFO << "  match to both          : " << stats_count_m_td << "\n";

    LOG_INFO << "Statistics of peptides (mapping to proteins):\n";
    LOG_INFO << "  no match (to 0 protein)         : " << stats_unmatched << "\n";
    LOG_INFO << "  unique match (to 1 protein)     : " << stats_matched_unique << "\n";
    LOG_INFO << "  non-unique match (to >1 protein): " << stats_matched_multi << std::endl;


    /// exit if no peptides were matched to decoy
    if ((stats_count_m_d + stats_count_m_td) == 0)
    {
      String msg("No peptides were matched to the decoy portion of the database! Did you provide the correct concatenated database? Are your 'decoy_string' (=" + getStringOption_("decoy_string") + ") and 'prefix' (=" + String(getFlag_("prefix")) + ") settings correct?");
      if (getStringOption_("missing_decoy_action") == "error")
      {
        LOG_ERROR << "Error: " << msg << "\nSet 'missing_decoy_action' to 'warn' if you are sure this is ok!\nAborting ..." << std::endl;
        return UNEXPECTED_RESULT;
      }
      else
      {
        LOG_WARN << "Warn: " << msg << "\nSet 'missing_decoy_action' to 'error' if you want to elevate this to an error!" << std::endl;
      }
    }

    /// for proteins --> peptides

    Int stats_new_proteins(0);
    Int stats_orphaned_proteins(0);

    // all peptides contain the correct protein hit references, now update the protein hits
    for (Size run_idx = 0; run_idx < prot_ids.size(); ++run_idx)
    {
      set<Size> masterset = runidx_to_protidx[run_idx]; // all found protein matches

      vector<ProteinHit> new_protein_hits;
      // go through existing hits and update (do not create from anew, as there might be other information (score, rank, etc.) which
      // we want to preserve
      for (vector<ProteinHit>::iterator p_hit = prot_ids[run_idx].getHits().begin(); p_hit != prot_ids[run_idx].getHits().end(); ++p_hit)
      {
        const String& acc = p_hit->getAccession();
        if (acc_to_prot.has(acc) // accession needs to exist in new FASTA file
            && masterset.find(acc_to_prot[acc]) != masterset.end())
        { // this accession was there already
          String seq;
          if (write_protein_sequence) 
          {
            seq = proteins[acc_to_prot[acc]].sequence;
          }
          p_hit->setSequence(seq);
          
          if (write_protein_description)
          {
            const String& description = proteins[acc_to_prot[acc]].description;
            //std::cout << "Description = " << description << "\n";
            p_hit->setDescription(description);
          }
          
          new_protein_hits.push_back(*p_hit);
          masterset.erase(acc_to_prot[acc]); // remove from master (at the end only new proteins remain)
        }
        else // old hit is orphaned
        {
          ++stats_orphaned_proteins;
          if (keep_unreferenced_proteins) new_protein_hits.push_back(*p_hit);
        }
      }

      // add remaining new hits
      for (set<Size>::const_iterator it = masterset.begin();
           it != masterset.end();
           ++it)
      {
        ProteinHit hit;
        hit.setAccession(proteins[*it].identifier);
        if (write_protein_sequence)
        {
          hit.setSequence(proteins[*it].sequence);
        }
        
        if (write_protein_description)
        {
          //std::cout << "Description = " << proteins[*it].description << "\n";
          hit.setDescription(proteins[*it].description);
        }
        
        new_protein_hits.push_back(hit);
        ++stats_new_proteins;
      }

      prot_ids[run_idx].setHits(new_protein_hits);
    }

    // annotate target/decoy status of proteins:
    for (vector<ProteinIdentification>::iterator id_it = prot_ids.begin(); id_it != prot_ids.end(); ++id_it)
    {
      for (vector<ProteinHit>::iterator hit_it = id_it->getHits().begin(); hit_it != id_it->getHits().end(); ++hit_it)
      {
        hit_it->setMetaValue("target_decoy", (protein_is_decoy[hit_it->getAccession()] ? "decoy" : "target"));
      }
    }

    LOG_INFO << "Statistics of proteins:\n";
    LOG_INFO << "  new proteins: " << stats_new_proteins << "\n";
    LOG_INFO << "  orphaned proteins: " << stats_orphaned_proteins << (keep_unreferenced_proteins ? " (all kept)" : " (all removed)") << "\n";

    writeDebug_("Reindexing finished!", 1);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    IdXMLFile().store(out, prot_ids, pep_ids);

    if ((!allow_unmatched) && (stats_unmatched > 0))
    {
      LOG_WARN << "PeptideIndexer found unmatched peptides, which could not be associated to a protein.\n"
               << "Potential solutions:\n"
               << "   - check your FASTA database for completeness\n"
               << "   - set 'enzyme:specificity' to match the identification parameters of the search engine\n"
               << "   - some engines (e.g. X! Tandem) employ loose cutting rules generating non-tryptic peptides;\n"
               << "     if you trust them, disable enzyme specificity\n"
               << "   - increase 'aaa_max' to allow more ambiguous amino acids\n"
               << "   - as a last resort: use the 'allow_unmatched' option to accept unmatched peptides\n"
               << "     (note that unmatched peptides cannot be used for FDR calculation or quantification)\n";

      LOG_WARN << "Result files were written, but PeptideIndexer will exit with an error code." << std::endl;
      return UNEXPECTED_RESULT;
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPPeptideIndexer tool;
  return tool.main(argc, argv);
}

/// @endcond
