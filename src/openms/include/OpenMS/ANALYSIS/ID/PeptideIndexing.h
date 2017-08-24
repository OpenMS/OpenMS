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
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_PEPTIDEINDEXING_H
#define OPENMS_ANALYSIS_ID_PEPTIDEINDEXING_H


#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <fstream>

namespace OpenMS
{

/**
  @brief Refreshes the protein references for all peptide hits in a vector of PeptideIdentifications and adds target/decoy information.

  All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy", 
  depending on whether the protein accession contains the decoy pattern (parameter @p decoy_string) as a suffix or prefix, respectively (see parameter @p prefix).
  For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins,
  only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool.
  (For FDR calculations, "target+decoy" peptide hits count as target hits.)

  @note Make sure that your protein names in the database contain a correctly formatted decoy string. This can be ensured by using @ref UTILS_DecoyDatabase.
        If the decoy identifier is not recognized successfully all proteins will be assumed to stem from the target-part of the query.<br>
        E.g., "sw|P33354_DECOY|YEHR_ECOLI Uncharacterized lipop..." is <b>invalid</b>, since the tool has no knowledge of how SwissProt entries are build up.
        A correct identifier could be "DECOY_sw|P33354|YEHR_ECOLI Uncharacterized li ..." or "sw|P33354|YEHR_ECOLI_DECOY Uncharacterized li", depending on whether you are
        using prefix or suffix annotation.<br>
        This tool will also report some helpful target/decoy statistics when it is done.

  By default this tool will fail if an unmatched peptide occurs, i.e. if the database does not contain the corresponding protein.
  You can force it to return successfully in this case by using the flag @p allow_unmatched.

  Search engines (such as Mascot) will replace ambiguous amino acids ('B', 'J', 'Z' and 'X') in the protein database with unambiguous amino acids in the reported peptides, e.g. exchange 'X' with 'H'.
  This will cause such peptides to not be found by exactly matching their sequences to the protein database.
  However, we can recover these cases by using tolerant search for ambiguous amino acids in the protein sequence. This is done by default with up to four amino acids
  per peptide hit. If you only want exact matches, set @p aaa_max to zero (but expect that unmatched peptides might occur)!

  Leucine/Isoleucine:
  Further complications can arise due to the presence of the isobaric amino acids isoleucine ('I') and leucine ('L') in protein sequences.
  Since the two have the exact same chemical composition and mass, they generally cannot be distinguished by mass spectrometry.
  If a peptide containing 'I' was reported as a match for a spectrum, a peptide containing 'L' instead would be an equally good match (and vice versa).
  To account for this inherent ambiguity, setting the flag @p IL_equivalent causes 'I' and 'L' to be considered as indistinguishable.@n
  For example, if the sequence "PEPTIDE" (matching "Protein1") was identified as a search hit,
  but the database additionally contained "PEPTLDE" (matching "Protein2"), running PeptideIndexer with the @p IL_equivalent option would
  report both "Protein1" and "Protein2" as accessions for "PEPTIDE".
  (This is independent of ambiguous matching via @p aaa_max.)
  Additionally, setting this flag will convert all 'J's in any protein sequence to 'I'. This way, no tolerant search is required for 'J' (but is still possible for all
  the other ambiguous amino acids).
  If @p write_protein_sequences is requested and @p IL_equivalent is set as well, both the I/L-version and unmodified protein sequences need to be stored internally.
  This requires some extra memory, roughly equivalent to the size of the FASTA database file itself.

  Enzyme specificity:
  Once a peptide sequence is found in a protein sequence, this does <b>not</b> imply that the hit is valid! This is where enzyme specificity comes into play.
  By default, we demand that the peptide is fully tryptic (i.e. the enzyme parameter is set to "trypsin" and specificity is "full").
  So unless the peptide coincides with C- and/or N-terminus of the protein, the peptide's cleavage pattern should fulfill the trypsin cleavage rule [KR][^P].
  
  We make two exceptions to the specificity constraints:
  1) for peptides starting at the second or third position of a protein are still considered N-terminally specific,
  since the residues can be cleaved off in vivo; X!Tandem reports these peptides. For example, the two peptides ABAR and LABAR would both match a protein starting with MLABAR.
  2) adventitious cleavage at Asp|Pro (Aspartate/D | Proline/P) is allowed for all enzymes (as supported by X!Tandem), i.e. counts as a proper cleavage site (see http://www.thegpm.org/tandem/release.html).
  
  You can relax the requirements further by choosing <tt>semi-tryptic</tt> (only one of two "internal" termini must match requirements)
  or <tt>none</tt> (essentially allowing all hits, no matter their context). These settings should not be used (due to high risk of reporting false positives),
  unless the search engine was instructed to search peptides in the same way.
  
  Threading:
  This tool support multiple threads (@p threads option) to speed up computation, at the cost of little extra memory.
*/

 class OPENMS_DLLAPI PeptideIndexing :
    public DefaultParamHandler, public ProgressLogger
  {
public:

    /// Exit codes
    enum ExitCodes
    {
      EXECUTION_OK,
      DATABASE_EMPTY,
      PEPTIDE_IDS_EMPTY,
      DATABASE_CONTAINS_MULTIPLES,
      ILLEGAL_PARAMETERS,
      UNEXPECTED_RESULT
    };

    /// Default constructor
    PeptideIndexing();

    /// Default destructor
    virtual ~PeptideIndexing();

    /// main method of PeptideIndexing
    ExitCodes run(std::vector<FASTAFile::FASTAEntry>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

protected:
    struct PeptideProteinMatchInformation
    {
      /// index of the protein the peptide is contained in
      OpenMS::Size protein_index;

      /// the position of the peptide in the protein
      OpenMS::Int position;

      /// the amino acid after the peptide in the protein
      char AABefore;

      /// the amino acid before the peptide in the protein
      char AAAfter;

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
    struct FoundProteinFunctor
    {
    public:
      typedef std::map<OpenMS::Size, std::set<PeptideProteinMatchInformation> > MapType;

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

      void merge(FoundProteinFunctor& other)
      {
        if (pep_to_prot.empty())
        { // first merge is easy
          pep_to_prot.swap(other.pep_to_prot);
        }
        else
        {
          for (FoundProteinFunctor::MapType::const_iterator it = other.pep_to_prot.begin(); it != other.pep_to_prot.end(); ++it)
          { // augment set
            this->pep_to_prot[it->first].insert(other.pep_to_prot[it->first].begin(), other.pep_to_prot[it->first].end());
          }
        }
        // cheap members
        this->filter_passed += other.filter_passed;
        this->filter_rejected += other.filter_rejected;
      }

      void addHit(const OpenMS::Size idx_pep,
        const OpenMS::Size idx_prot,
        const OpenMS::Size len_pep,
        const OpenMS::String& seq_prot,
        OpenMS::Int position)
      {
        if (enzyme_.isValidProduct(seq_prot, position, len_pep, true))
        {
          PeptideProteinMatchInformation match;
          match.protein_index = idx_prot;
          match.position = position;
          match.AABefore = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : seq_prot[position - 1];
          match.AAAfter = (position + len_pep >= seq_prot.size()) ? PeptideEvidence::C_TERMINAL_AA : seq_prot[position + len_pep];
          pep_to_prot[idx_pep].insert(match);
          ++filter_passed;
          DEBUG_ONLY std::cerr << "Hit: " << len_pep << " (peplen) with hit to protein " << seq_prot << " at position " << position << std::endl;
        }
        else
        {
          //std::cerr << "REJECTED Peptide " << seq_pep << " with hit to protein "
          //  << seq_prot << " at position " << position << std::endl;
          ++filter_rejected;
        }
      }

    };

    inline void addHits_(AhoCorasickAmbiguous& fuzzyAC, const AhoCorasickAmbiguous::FuzzyACPattern& pattern, const AhoCorasickAmbiguous::PeptideDB& pep_DB, const String& prot, const String& full_prot, SignedSize i, Int offset, FoundProteinFunctor& func_threads) const
    {
      fuzzyAC.setProtein(prot);
      while (fuzzyAC.findNext(pattern))
      {
        const seqan::Peptide& tmp_pep = pep_DB[fuzzyAC.getHitDBIndex()];
        func_threads.addHit(fuzzyAC.getHitDBIndex(), i, length(tmp_pep), full_prot, fuzzyAC.getHitProteinPosition() + offset);
      }

    }

    virtual void updateMembers_();

    String decoy_string_;
    bool prefix_;
    String missing_decoy_action_;
    String enzyme_name_;
    String enzyme_specificity_;

    bool write_protein_sequence_;
    bool write_protein_description_;
    bool keep_unreferenced_proteins_;
    bool allow_unmatched_;
    bool IL_equivalent_;

    Int aaa_max_;

  };
}

#endif // OPENMS_ANALYSIS_ID_PEPTIDEINDEXING_H
