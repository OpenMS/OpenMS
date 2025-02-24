// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtilsSimple.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

namespace OpenMS
{

/**
  @brief Refreshes the protein references for all peptide hits in a vector of PeptideIdentifications and adds target/decoy information.

  All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy", 
  depending on whether the protein accession contains the decoy pattern (parameter @p decoy_string) as a suffix or prefix, respectively (see parameter @p prefix).
  For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins,
  only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool.
  (For FDR calculations, "target+decoy" peptide hits count as target hits.)

  @note Make sure that your protein names in the database contain a correctly formatted decoy string. This can be ensured by using @ref TOPP_DecoyDatabase.
        If the decoy identifier is not recognized successfully, all proteins will be assumed to stem from the target-part of the query.<br>
        E.g., "sw|P33354_DECOY|YEHR_ECOLI Uncharacterized lipop..." is <b>invalid</b>, since the tool has no knowledge of how SwissProt entries are build up.
        A correct identifier could be "DECOY_sw|P33354|YEHR_ECOLI Uncharacterized li ..." or "sw|P33354|YEHR_ECOLI_DECOY Uncharacterized li", depending on whether you are
        using prefix or suffix annotation.<br>
  
  Some helpful target/decoy statistics will be reported when done.

  By default this tool will fail if an unmatched peptide occurs, i.e. if the database does not contain the corresponding protein.
  You can force it to return successfully in this case by setting @p '-unmatched_action' to accept or even remove those hits.

  Search engines (such as X!Tandem) will replace ambiguous amino acids ('B', 'J', 'Z' and 'X') in the protein database with unambiguous amino 
  acids in the reported peptides, e.g. exchange 'X' with 'H'.
  This will cause such peptides to not be found by exactly matching their sequences to the protein database.
  However, we can recover these cases by using tolerant search for ambiguous amino acids in the protein sequence. This is done by default with
  up to three amino acids per peptide hit. If you only want exact matches, set @p aaa_max to zero (but expect that unmatched peptides might occur)!

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
  By default, the enzyme and the specificity used during search is derived from metadata in the idXML files ('auto' setting).
  
  We make two exceptions to any specificity constraints:
  1) for peptides starting at the second or third position of a protein are still considered N-terminally specific,
  since the residues can be cleaved off in vivo; X!Tandem reports these peptides. For example, the two peptides ABAR and LABAR would both match a protein starting with MLABAR.
  2) adventitious cleavage at Asp|Pro (Aspartate/D | Proline/P) is allowed for all enzymes (as supported by X!Tandem), i.e. counts as a proper cleavage site (see http://www.thegpm.org/tandem/release.html).
  
  You can relax the requirements further by choosing <tt>semi-tryptic</tt> (only one of two "internal" termini must match requirements)
  or <tt>none</tt> (essentially allowing all hits, no matter their context). These settings should not be used (due to high risk of reporting false positives),
  unless the search engine was instructed to search peptides in the same way (but then the default 'auto' setting will do the correct thing).

  X!Tandem treats any occurrence of 'X' as stop codon (and thus as cleavage site). The resulting peptide will be non- or semi-tryptic.
  Those hits will not be matched and need to be removed using @p '-unmatched_action' (do not use termini specificity to cheat around it! It adds more false hits!).
  
  The FASTA file should not contain duplicate protein accessions (since accessions are not validated) if a correct unique-matching annotation is important (target/decoy annotation is still correct).

  Threading:
  This tool support multiple threads (@p threads option) to speed up computation, at the cost of little extra memory.

*/

 class OPENMS_DLLAPI PeptideIndexing :
    public DefaultParamHandler, public ProgressLogger
  {
public:
    /// name of enzyme/specificity which signals that the enzyme/specificity should be taken from meta information
    static char const* const AUTO_MODE; /* = 'auto' */ 

    /// Exit codes
    enum ExitCodes
    {
      EXECUTION_OK,
      DATABASE_EMPTY,
      PEPTIDE_IDS_EMPTY,
      ILLEGAL_PARAMETERS,
      UNEXPECTED_RESULT
    };

    /// Action to take when peptide hits could not be matched
    enum class Unmatched
    {
      IS_ERROR,  ///< throws an error (and returns no results)
      WARN,      ///< skips annotation with target/decoy but returns with 'success'
      REMOVE,    ///< removes unmatched hits entirely and returns with 'success'
      SIZE_OF_UNMATCHED
    };
    static const std::array<std::string, (Size)Unmatched::SIZE_OF_UNMATCHED> names_of_unmatched;
    
    enum class MissingDecoy
    {
      IS_ERROR,
      WARN,
      SILENT,
      SIZE_OF_MISSING_DECOY
    };
    static const std::array<std::string, (Size)MissingDecoy::SIZE_OF_MISSING_DECOY> names_of_missing_decoy;

    /// Default constructor
    PeptideIndexing();

    /// Default destructor
    ~PeptideIndexing() override;

    /// forward for old interface and pyOpenMS; use other run() methods for more control
    ExitCodes run(std::vector<FASTAFile::FASTAEntry>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    /**
    @brief Re-index peptide identifications honoring enzyme cutting rules, ambiguous amino acids and target/decoy hits.
    
    Template parameter 'T' can be either TFI_File or TFI_Vector. If the data is already available, use TFI_Vector and pass the vector.
    If the data is still in a FASTA file and its not needed afterwards for additional processing, use TFI_File and pass the filename.

    PeptideIndexer refreshes target/decoy information and mapping of peptides to proteins.
    The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool. (For FDR calculations, "target+decoy" peptide hits count as target hits.)

    PeptideIndexer allows for ambiguous amino acids (B|J|Z|X) in the protein database, but not in the peptide sequences. 
    For the latter only I/L can be treated as equivalent (see 'IL_equivalent' flag), but 'J' is not allowed.
  
    Enzyme cutting rules and partial specificity can be specified.

    Resulting protein hits appear in the order of the FASTA file, except for orphaned proteins, which will appear first with an empty target_decoy metavalue.
    Duplicate protein accessions & sequences will not raise a warning, but create multiple hits (PeptideIndexer scans over the FASTA file once for efficiency
    reasons, and thus might not see all accessions & sequences at once).
      
    All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". 
    For proteins the possible values are "target" and "decoy", depending on whether the protein accession contains the decoy pattern (parameter @p decoy_string) 
    as a suffix or prefix, respectively (see parameter @p prefix). 
  
    Peptide hits are annotated with metavalue 'protein_references', and if matched to at least one protein also with metavalue 'target_decoy'.
    The possible values for 'target_decoy' are "target", "decoy" and "target+decoy", 
    depending on whether the peptide sequence is found only in target proteins, only in decoy proteins, or in both. The metavalue is not present, if the peptide is unmatched.
  
    Runtime: PeptideIndexer is usually very fast (loading and storing the data takes the most time) and search speed can be further improved (linearly), but using more threads. 
    Avoid allowing too many (>=4) ambiguous amino acids if your database contains long stretches of 'X' (exponential search space).

    @param proteins A list of proteins -- either read piecewise from a FASTA file or as existing vector of FASTAEntries.
    @param prot_ids Resulting protein identifications associated to pep_ids (will be re-written completely)
    @param pep_ids Peptide identifications which should be search within @p proteins and then linked to @p prot_ids
    @return Exit status codes.

    */
    ExitCodes run(FASTAContainer<TFI_File>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    /// Same as run() with TFI_File, but for proteins which are already in memory
    ExitCodes run(FASTAContainer<TFI_Vector>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    /// Which string is used to determine if a protein is a decoy or not
    const String& getDecoyString() const;

    /// Is the decoy string position a prefix or suffix?
    bool isPrefix() const;

 protected:
    void updateMembers_() override;

    template<typename T> ExitCodes run_(FASTAContainer<T>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    String decoy_string_{};
    bool prefix_{ false };
    MissingDecoy missing_decoy_action_ = MissingDecoy::IS_ERROR;
    String enzyme_name_{};
    String enzyme_specificity_{};

    bool write_protein_sequence_{ false };
    bool write_protein_description_{ false };
    bool keep_unreferenced_proteins_{ false };
    Unmatched unmatched_action_ = Unmatched::IS_ERROR;
    bool IL_equivalent_{ false };
    bool allow_nterm_protein_cleavage_{ true };

    Int aaa_max_{0};
    Int mm_max_{0};
 };
}

