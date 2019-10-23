// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once


#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <atomic>
#include <algorithm>
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
  
  Some helpful target/decoy statistics will be reported when done.

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
  
  The FASTA file should not contain duplicate protein accessions (since accessions are not validated) if a correct unique-matching annotation is important (target/decoy annotation is still correct).

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
      ILLEGAL_PARAMETERS,
      UNEXPECTED_RESULT
    };

    /// Default constructor
    PeptideIndexing();

    /// Default destructor
    ~PeptideIndexing() override;


     /// forward for old interface and pyOpenMS; use run<T>() for more control
    inline ExitCodes run(std::vector<FASTAFile::FASTAEntry>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
    {
      FASTAContainer<TFI_Vector> protein_container(proteins);
      return run<TFI_Vector>(protein_container, prot_ids, pep_ids);
    }

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
    template<typename T>
    ExitCodes run(FASTAContainer<T>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
    {
      // no decoy string provided? try to deduce from data
      if (decoy_string_.empty())
      {
        auto r = DecoyHelper::findDecoyString(proteins);
        proteins.reset();
        if (!r.success)
        {
          r.is_prefix = true;
          r.name = "DECOY_";
          OPENMS_LOG_WARN << "Unable to determine decoy string automatically (not enough decoys were detected)! Using default " << (r.is_prefix ? "prefix" : "suffix") << " decoy string '" << r.name << "'\n"
                          << "If you think that this is incorrect, please provide a decoy_string and its position manually!" << std::endl;
        }
        prefix_ = r.is_prefix;
        decoy_string_ = r.name;
        // decoy string and position was extracted successfully
        OPENMS_LOG_INFO << "Using " << (prefix_ ? "prefix" : "suffix") << " decoy string '" << decoy_string_ << "'" << std::endl;
      }

      //---------------------------------------------------------------
      // parsing parameters, correcting xtandem and MSGFPlus parameters
      //---------------------------------------------------------------
      ProteaseDigestion enzyme;
      enzyme.setEnzyme(enzyme_name_);
      enzyme.setSpecificity(enzyme.getSpecificityByName(enzyme_specificity_));

      bool xtandem_fix_parameters = true, msgfplus_fix_parameters = true;

      // specificity is none or semi? don't automate xtandem 
      if (enzyme.getSpecificity() == EnzymaticDigestion::SPEC_SEMI ||
          enzyme.getSpecificity() == EnzymaticDigestion::SPEC_NONE) 
      {
        xtandem_fix_parameters = false;
      }

      // determine if search engine is solely xtandem or MSGFPlus
      for (const auto& prot_id : prot_ids)
      {
        String search_engine = prot_id.getSearchEngine();
        StringUtils::toUpper(search_engine);
        if (search_engine != "XTANDEM") { xtandem_fix_parameters = false; }
        if (!(search_engine == "MSGFPLUS" || search_engine == "MS-GF+")) { msgfplus_fix_parameters = false; }
      }

      // solely MSGFPlus -> Trypsin/P as enzyme
      if (msgfplus_fix_parameters && enzyme.getEnzymeName() == "Trypsin")
      {
        OPENMS_LOG_WARN << "MSGFPlus detected but enzyme cutting rules were set to Trypsin. Correcting to Trypsin/P to copy with special cutting rule in MSGFPlus." << std::endl;
        enzyme.setEnzyme("Trypsin/P");
      }

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
      // cache the first proteins
      const size_t PROTEIN_CACHE_SIZE = 4e5; // 400k should be enough for most DB's and is not too hard on memory either (~200 MB FASTA)

      this->startProgress(0, 1, "Load first DB chunk");
      proteins.cacheChunk(PROTEIN_CACHE_SIZE);
      this->endProgress();

      if (proteins.empty()) // we do not allow an empty database
      {
        OPENMS_LOG_ERROR << "Error: An empty database was provided. Mapping makes no sense. Aborting..." << std::endl;
        return DATABASE_EMPTY;
      }

      if (pep_ids.empty()) // Aho-Corasick requires non-empty input; but we allow this case, since the TOPP tool should not crash when encountering a bad raw file (with no PSMs)
      {
        OPENMS_LOG_WARN << "Warning: An empty set of peptide identifications was provided. Output will be empty as well." << std::endl;
        if (!keep_unreferenced_proteins_)
        {
          // delete only protein hits, not whole ID runs incl. meta data:
          for (std::vector<ProteinIdentification>::iterator it = prot_ids.begin();
            it != prot_ids.end(); ++it)
          {
            it->getHits().clear();
          }
        }
        return PEPTIDE_IDS_EMPTY;
      }

      FoundProteinFunctor func(enzyme, xtandem_fix_parameters); // store the matches
      Map<String, Size> acc_to_prot; // map: accessions --> FASTA protein index
      std::vector<bool> protein_is_decoy; // protein index -> is decoy?
      std::vector<std::string> protein_accessions; // protein index -> accession

      bool invalid_protein_sequence = false; // check for proteins with modifications, i.e. '[' or '(', and throw an exception

      { // new scope - forget data after search
      
        /*
        BUILD Peptide DB
        */
        bool has_illegal_AAs(false);
        AhoCorasickAmbiguous::PeptideDB pep_DB;
        for (std::vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
        {
          //String run_id = it1->getIdentifier();
          const std::vector<PeptideHit>& hits = it1->getHits();
          for (std::vector<PeptideHit>::const_iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
          {
            //
            // Warning:
            // do not skip over peptides here, since the results are iterated in the same way
            //
            String seq = it2->getSequence().toUnmodifiedString().remove('*'); // make a copy, i.e. do NOT change the peptide sequence!
            if (seqan::isAmbiguous(seqan::AAString(seq.c_str())))
            { // do not quit here, to show the user all sequences .. only quit after loop
              OPENMS_LOG_ERROR << "Peptide sequence '" << it2->getSequence() << "' contains one or more ambiguous amino acids (B|J|Z|X).\n";
              has_illegal_AAs = true;
            }
            if (IL_equivalent_) // convert L to I;
            {
              seq.substitute('L', 'I');
            }
            appendValue(pep_DB, seq.c_str());
          }
        }
        if (has_illegal_AAs)
        {
          OPENMS_LOG_ERROR << "One or more peptides contained illegal amino acids. This is not allowed!"
                    << "\nPlease either remove the peptide or replace it with one of the unambiguous ones (while allowing for ambiguous AA's to match the protein)." << std::endl;;
        }

        OPENMS_LOG_INFO << "Mapping " << length(pep_DB) << " peptides to " << (proteins.size() == PROTEIN_CACHE_SIZE ? "? (unknown number of)" : String(proteins.size()))  << " proteins." << std::endl;

        if (length(pep_DB) == 0)
        { // Aho-Corasick will crash if given empty needles as input
          OPENMS_LOG_WARN << "Warning: Peptide identifications have no hits inside! Output will be empty as well." << std::endl;
          return PEPTIDE_IDS_EMPTY;
        }

        /*
           Aho Corasick (fast)
        */
        OPENMS_LOG_INFO << "Searching with up to " << aaa_max_ << " ambiguous amino acid(s) and " << mm_max_ << " mismatch(es)!" << std::endl;
        SysInfo::MemUsage mu;
        OPENMS_LOG_INFO << "Building trie ...";
        StopWatch s;
        s.start();
        AhoCorasickAmbiguous::FuzzyACPattern pattern;
        AhoCorasickAmbiguous::initPattern(pep_DB, aaa_max_, mm_max_, pattern);
        s.stop();
        OPENMS_LOG_INFO << " done (" << int(s.getClockTime()) << "s)" << std::endl;
        s.reset();

        uint16_t count_j_proteins(0);
        bool has_active_data = true; // becomes false if end of FASTA file is reached
        const std::string jumpX(aaa_max_ + mm_max_ + 1, 'X'); // jump over stretches of 'X' which cost a lot of time; +1 because  AXXA is a valid hit for aaa_max == 2 (cannot split it)
        // use very large target value for progress if DB size is unknown (did not fit into first chunk)
        this->startProgress(0, proteins.size() == PROTEIN_CACHE_SIZE ? std::numeric_limits<SignedSize>::max() : proteins.size(), "Aho-Corasick");
        std::atomic<int> progress_prots(0);
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          FoundProteinFunctor func_threads(enzyme, xtandem_fix_parameters);
          Map<String, Size> acc_to_prot_thread; // map: accessions --> FASTA protein index
          AhoCorasickAmbiguous fuzzyAC;
          String prot;

          while (true) 
          {
            #pragma omp barrier // all threads need to be here, since we are about to swap protein data
            #pragma omp single
            {
              DEBUG_ONLY std::cerr << " activating cache ...\n";
              has_active_data = proteins.activateCache(); // swap in last cache
              protein_accessions.resize(proteins.getChunkOffset() + proteins.chunkSize());
            } // implicit barrier here
            
            if (!has_active_data) break; // leave while-loop
            SignedSize prot_count = (SignedSize)proteins.chunkSize();

            #pragma omp master
            {
              DEBUG_ONLY std::cerr << "Filling Protein Cache ...";
              proteins.cacheChunk(PROTEIN_CACHE_SIZE);
              protein_is_decoy.resize(proteins.getChunkOffset() + prot_count);
              for (SignedSize i = 0; i < prot_count; ++i)
              { // do this in master only, to avoid false sharing
                const String& seq = proteins.chunkAt(i).identifier;
                protein_is_decoy[i + proteins.getChunkOffset()] = (prefix_ ? seq.hasPrefix(decoy_string_) : seq.hasSuffix(decoy_string_));
              }
              DEBUG_ONLY std::cerr << " done" << std::endl;
            }
            DEBUG_ONLY std::cerr << " starting for loop \n";
            // search all peptides in each protein
            #pragma omp for schedule(dynamic, 100) nowait
            for (SignedSize i = 0; i < prot_count; ++i)
            {
              ++progress_prots; // atomic
              if (omp_get_thread_num() == 0)
              {
                this->setProgress(progress_prots);
              }

              prot = proteins.chunkAt(i).sequence;
              prot.remove('*');

              // check for invalid sequences with modifications
              if (prot.has('[') || prot.has('('))
              { 
                 invalid_protein_sequence = true; // not omp-critical because its write-only
                 // we cannot throw an exception here, since we'd need to catch it within the parallel region
              }
              
              // convert  L/J to I; also replace 'J' in proteins
              if (IL_equivalent_)
              {
                prot.substitute('L', 'I');
                prot.substitute('J', 'I');
              }
              else
              { // warn if 'J' is found (it eats into aaa_max)
                if (prot.has('J'))
                {
                 #pragma omp atomic
                 ++count_j_proteins;
                }
              }

              Size prot_idx = i + proteins.getChunkOffset();
              
              // test if protein was a hit
              Size hits_total = func_threads.filter_passed + func_threads.filter_rejected;

              // check if there are stretches of 'X'
              if (prot.has('X'))
              {
                // create chunks of the protein (splitting it at stretches of 'X..X') and feed them to AC one by one
                size_t offset = -1, start = 0;
                while ((offset = prot.find(jumpX, offset + 1)) != std::string::npos)
                {
                  //std::cout << "found X..X at " << offset << " in protein " << proteins[i].identifier << "\n";
                  addHits_(fuzzyAC, pattern, pep_DB, prot.substr(start, offset + jumpX.size() - start), prot, prot_idx, (int)start, func_threads);
                  // skip ahead while we encounter more X...
                  while (offset + jumpX.size() < prot.size() && prot[offset + jumpX.size()] == 'X') ++offset;
                  start = offset;
                  //std::cout << "  new start: " << start << "\n";
                }
                // last chunk
                if (start < prot.size())
                {
                  addHits_(fuzzyAC, pattern, pep_DB, prot.substr(start), prot, prot_idx, (int)start, func_threads);
                }
              }
              else
              {
                addHits_(fuzzyAC, pattern, pep_DB, prot, prot, prot_idx, 0, func_threads);
              }
              // was protein found?
              if (hits_total < func_threads.filter_passed + func_threads.filter_rejected)
              {
                protein_accessions[prot_idx] = proteins.chunkAt(i).identifier;
                acc_to_prot_thread[protein_accessions[prot_idx]] = prot_idx;
              }
            } // end parallel FOR

            // join results again
            DEBUG_ONLY std::cerr << " critical now \n";
            #ifdef _OPENMP
            #pragma omp critical(PeptideIndexer_joinAC)
            #endif
            {
              s.start();
              // hits
              func.merge(func_threads);
              // accession -> index
              acc_to_prot.insert(acc_to_prot_thread.begin(), acc_to_prot_thread.end());
              acc_to_prot_thread.clear();
              s.stop();
            } // OMP end critical
          } // end readChunk
        } // OMP end parallel
        this->endProgress();
        std::cout << "Merge took: " << s.toString() << "\n";
        mu.after();
        std::cout << mu.delta("Aho-Corasick") << "\n\n";

        OPENMS_LOG_INFO << "\nAho-Corasick done:\n  found " << func.filter_passed << " hits for " << func.pep_to_prot.size() << " of " << length(pep_DB) << " peptides.\n";

        // write some stats
        OPENMS_LOG_INFO << "Peptide hits passing enzyme filter: " << func.filter_passed << "\n"
                 << "     ... rejected by enzyme filter: " << func.filter_rejected << std::endl;

        if (count_j_proteins)
        {
          OPENMS_LOG_WARN << "PeptideIndexer found " << count_j_proteins << " protein sequences in your database containing the amino acid 'J'."
            << "To match 'J' in a protein, an ambiguous amino acid placeholder for I/L will be used.\n"
            << "This costs runtime and eats into the 'aaa_max' limit, leaving less opportunity for B/Z/X matches.\n"
            << "If you want 'J' to be treated as unambiguous, enable '-IL_equivalent'!" << std::endl;
        }

      } // end local scope

      //
      //   do mapping 
      //
      // index existing proteins
      Map<String, Size> runid_to_runidx; // identifier to index
      for (Size run_idx = 0; run_idx < prot_ids.size(); ++run_idx)
      {
        runid_to_runidx[prot_ids[run_idx].getIdentifier()] = run_idx;
      }
      
      // for peptides --> proteins
      Size stats_matched_unique(0);
      Size stats_matched_multi(0);
      Size stats_unmatched(0);    // no match to DB
      Size stats_count_m_t(0);    // match to Target DB
      Size stats_count_m_d(0);    // match to Decoy DB
      Size stats_count_m_td(0);   // match to T+D DB

      Map<Size, std::set<Size> > runidx_to_protidx; // in which protID do appear which proteins (according to mapped peptides)

      Size pep_idx(0);
      for (std::vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
      {
        // which ProteinIdentification does the peptide belong to?
        Size run_idx = runid_to_runidx[it1->getIdentifier()];

        std::vector<PeptideHit>& hits = it1->getHits();

        for (std::vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
        {
          // clear protein accessions
          it2->setPeptideEvidences(std::vector<PeptideEvidence>());
          
          //
          // is this a decoy hit?
          //
          bool matches_target(false);
          bool matches_decoy(false);

          std::set<Size> prot_indices; /// protein hits of this peptide
          // add new protein references
          for (std::set<PeptideProteinMatchInformation>::const_iterator it_i = func.pep_to_prot[pep_idx].begin();
            it_i != func.pep_to_prot[pep_idx].end(); ++it_i)
          {
            prot_indices.insert(it_i->protein_index);
            const String& accession = protein_accessions[it_i->protein_index];
            PeptideEvidence pe(accession, it_i->position, it_i->position + (int)it2->getSequence().size() - 1, it_i->AABefore, it_i->AAAfter);
            it2->addPeptideEvidence(pe);

            runidx_to_protidx[run_idx].insert(it_i->protein_index); // fill protein hits

            if (protein_is_decoy[it_i->protein_index])
            {
              matches_decoy = true;
            }
            else
            {
              matches_target = true;
            }
          }

          if (matches_decoy && matches_target)
          {
            it2->setMetaValue("target_decoy", "target+decoy");
            ++stats_count_m_td;
          }
          else if (matches_target)
          {
            it2->setMetaValue("target_decoy", "target");
            ++stats_count_m_t;
          }
          else if (matches_decoy)
          {
            it2->setMetaValue("target_decoy", "decoy");
            ++stats_count_m_d;
          } // else: could match to no protein (i.e. both are false)
          //else ... // not required (handled below; see stats_unmatched);

          if (prot_indices.size() == 1)
          {
            it2->setMetaValue("protein_references", "unique");
            ++stats_matched_unique;
          }
          else if (prot_indices.size() > 1)
          {
            it2->setMetaValue("protein_references", "non-unique");
            ++stats_matched_multi;
          }
          else
          {
            it2->setMetaValue("protein_references", "unmatched");
            ++stats_unmatched;
            if (stats_unmatched < 15) OPENMS_LOG_INFO << "Unmatched peptide: " << it2->getSequence() << "\n";
            else if (stats_unmatched == 15) OPENMS_LOG_INFO << "Unmatched peptide: ...\n";
          }

          ++pep_idx; // next hit
        }

      }

      Size total_peptides = stats_count_m_t + stats_count_m_d + stats_count_m_td + stats_unmatched;
      OPENMS_LOG_INFO << "-----------------------------------\n";
      OPENMS_LOG_INFO << "Peptide statistics\n";
      OPENMS_LOG_INFO << "\n";
      OPENMS_LOG_INFO << "  unmatched                : " << stats_unmatched << " (" << stats_unmatched * 100 / total_peptides << " %)\n";
      OPENMS_LOG_INFO << "  target/decoy:\n";
      OPENMS_LOG_INFO << "    match to target DB only: " << stats_count_m_t << " (" << stats_count_m_t * 100 / total_peptides << " %)\n";
      OPENMS_LOG_INFO << "    match to decoy DB only : " << stats_count_m_d << " (" << stats_count_m_d * 100 / total_peptides << " %)\n";
      OPENMS_LOG_INFO << "    match to both          : " << stats_count_m_td << " (" << stats_count_m_td * 100 / total_peptides << " %)\n";
      OPENMS_LOG_INFO << "\n";
      OPENMS_LOG_INFO << "  mapping to proteins:\n";
      OPENMS_LOG_INFO << "    no match (to 0 protein)         : " << stats_unmatched << "\n";
      OPENMS_LOG_INFO << "    unique match (to 1 protein)     : " << stats_matched_unique << "\n";
      OPENMS_LOG_INFO << "    non-unique match (to >1 protein): " << stats_matched_multi << std::endl;

      /// for proteins --> peptides
      Size stats_matched_proteins(0), stats_matched_new_proteins(0), stats_orphaned_proteins(0), stats_proteins_target(0), stats_proteins_decoy(0);

      // all peptides contain the correct protein hit references, now update the protein hits
      for (Size run_idx = 0; run_idx < prot_ids.size(); ++run_idx)
      {
        std::set<Size> masterset = runidx_to_protidx[run_idx]; // all protein matches from above

        std::vector<ProteinHit>& phits = prot_ids[run_idx].getHits();
        {
          // go through existing protein hits and count orphaned proteins (with no peptide hits)
          std::vector<ProteinHit> orphaned_hits;
          for (std::vector<ProteinHit>::iterator p_hit = phits.begin(); p_hit != phits.end(); ++p_hit)
          {
            const String& acc = p_hit->getAccession();
            if (!acc_to_prot.has(acc)) // acc_to_prot only contains found proteins from current run
            { // old hit is orphaned
              ++stats_orphaned_proteins;
              if (keep_unreferenced_proteins_)
              {
                p_hit->setMetaValue("target_decoy", "");
                orphaned_hits.push_back(*p_hit);
              }
            }
          }
          // only keep orphaned hits (if any)
          phits = orphaned_hits;
        }

        // add new protein hits
        FASTAFile::FASTAEntry fe;
        phits.reserve(phits.size() + masterset.size());
        for (std::set<Size>::const_iterator it = masterset.begin(); it != masterset.end(); ++it)
        {
          ProteinHit hit;
          hit.setAccession(protein_accessions[*it]);
          
          if (write_protein_sequence_ || write_protein_description_)
          {
            proteins.readAt(fe, *it);
            if (write_protein_sequence_)
            {
              hit.setSequence(fe.sequence);
            } // no else, since sequence is empty by default
            if (write_protein_description_)
            {
              hit.setDescription(fe.description);
            } // no else, since description is empty by default
          }
          if (protein_is_decoy[*it])
          {
            hit.setMetaValue("target_decoy", "decoy");
            ++stats_proteins_decoy;
          }
          else
          {
            hit.setMetaValue("target_decoy", "target");
            ++stats_proteins_target;
          }
          phits.push_back(hit);
          ++stats_matched_new_proteins;
        }
        stats_matched_proteins += phits.size();
      }


      OPENMS_LOG_INFO << "-----------------------------------\n";
      OPENMS_LOG_INFO << "Protein statistics\n";
      OPENMS_LOG_INFO << "\n";
      OPENMS_LOG_INFO << "  total proteins searched: " << proteins.size() << "\n";
      OPENMS_LOG_INFO << "  matched proteins       : " << stats_matched_proteins << " (" << stats_matched_new_proteins << " new)\n";
      if (stats_matched_proteins)
      { // prevent Division-by-0 Exception
        OPENMS_LOG_INFO << "  matched target proteins: " << stats_proteins_target << " (" << stats_proteins_target * 100 / stats_matched_proteins << " %)\n";
        OPENMS_LOG_INFO << "  matched decoy proteins : " << stats_proteins_decoy << " (" << stats_proteins_decoy * 100 / stats_matched_proteins << " %)\n";
      }
      OPENMS_LOG_INFO << "  orphaned proteins      : " << stats_orphaned_proteins << (keep_unreferenced_proteins_ ? " (all kept)" : " (all removed)\n");
      OPENMS_LOG_INFO << "-----------------------------------" << std::endl;


      /// exit if no peptides were matched to decoy
      bool has_error = false;

      if (invalid_protein_sequence)
      {
        OPENMS_LOG_ERROR << "Error: One or more protein sequences contained the characters '[' or '(', which are illegal in protein sequences."
                 << "\nPeptide hits might be masked by these characters (which usually indicate presence of modifications).\n";
        has_error = true;
      }

      if ((stats_count_m_d + stats_count_m_td) == 0)
      {
        String msg("No peptides were matched to the decoy portion of the database! Did you provide the correct concatenated database? Are your 'decoy_string' (=" + String(decoy_string_) + ") and 'decoy_string_position' (=" + String(param_.getValue("decoy_string_position")) + ") settings correct?");
        if (missing_decoy_action_ == "error")
        {
          OPENMS_LOG_ERROR << "Error: " << msg << "\nSet 'missing_decoy_action' to 'warn' if you are sure this is ok!\nAborting ..." << std::endl;
          has_error = true;
        }
        else if (missing_decoy_action_ == "warn")
        {
          OPENMS_LOG_WARN << "Warn: " << msg << "\nSet 'missing_decoy_action' to 'error' if you want to elevate this to an error!" << std::endl;
        }
        else // silent
        {
        }
      }

      if ((!allow_unmatched_) && (stats_unmatched > 0))
      {
        OPENMS_LOG_ERROR << "PeptideIndexer found unmatched peptides, which could not be associated to a protein.\n"
                  << "Potential solutions:\n"
                  << "   - check your FASTA database for completeness\n"
                  << "   - set 'enzyme:specificity' to match the identification parameters of the search engine\n"
                  << "   - some engines (e.g. X! Tandem) employ loose cutting rules generating non-tryptic peptides;\n"
                  << "     if you trust them, disable enzyme specificity\n"
                  << "   - increase 'aaa_max' to allow more ambiguous amino acids\n"
                  << "   - as a last resort: use the 'allow_unmatched' option to accept unmatched peptides\n"
                  << "     (note that unmatched peptides cannot be used for FDR calculation or quantification)\n";
        has_error = true;
      }

      if (has_error)
      {
        OPENMS_LOG_ERROR << "Result files will be written, but PeptideIndexer will exit with an error code." << std::endl;
        return UNEXPECTED_RESULT;
      }
      return EXECUTION_OK;
    }

     const String& getDecoyString() const;

     bool isPrefix() const;

 protected:

    struct PeptideProteinMatchInformation
    {
      OpenMS::Size protein_index; //< index of the protein the peptide is contained in
      OpenMS::Int position; //< the position of the peptide in the protein
      char AABefore; //< the amino acid after the peptide in the protein
      char AAAfter; //< the amino acid before the peptide in the protein

      const std::tuple<const Size&, const Int&, const char&, const char&> tie() const
      {
        return std::tie(protein_index, position, AABefore, AAAfter);
      }
      bool operator<(const PeptideProteinMatchInformation& other) const
      {
        return tie() < other.tie();
      }
      bool operator==(const PeptideProteinMatchInformation& other) const
      {
        return tie() == other.tie();
      }
    };

    struct FoundProteinFunctor
    {
    public:
      typedef std::map<OpenMS::Size, std::set<PeptideProteinMatchInformation> > MapType;
      MapType pep_to_prot; //< peptide index --> protein indices
      OpenMS::Size filter_passed; //< number of accepted hits (passing addHit() constraints)
      OpenMS::Size filter_rejected; //< number of rejected hits (not passing addHit())

    private:
      ProteaseDigestion enzyme_;
      bool xtandem_; //< are we checking xtandem cleavage rules?

    public:
      explicit FoundProteinFunctor(const ProteaseDigestion& enzyme, bool xtandem) :
        pep_to_prot(), filter_passed(0), filter_rejected(0), enzyme_(enzyme), xtandem_(xtandem)
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
          other.pep_to_prot.clear();
        }
        // cheap members
        this->filter_passed += other.filter_passed;
        other.filter_passed = 0;
        this->filter_rejected += other.filter_rejected;
        other.filter_rejected = 0;
      }

      void addHit(const OpenMS::Size idx_pep,
        const OpenMS::Size idx_prot,
        const OpenMS::Size len_pep,
        const OpenMS::String& seq_prot,
        OpenMS::Int position)
      {
        if (enzyme_.isValidProduct(seq_prot, position, len_pep, true, true, xtandem_))
        {
          PeptideProteinMatchInformation match;
          match.protein_index = idx_prot;
          match.position = position;
          match.AABefore = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : seq_prot[position - 1];
          match.AAAfter = (position + len_pep >= seq_prot.size()) ? PeptideEvidence::C_TERMINAL_AA : seq_prot[position + len_pep];
          pep_to_prot[idx_pep].insert(match);
          ++filter_passed;
        }
        else
        {
          //std::cerr << "REJECTED Peptide " << seq_pep << " with hit to protein "
          //  << seq_prot << " at position " << position << std::endl;
          ++filter_rejected;
        }
      }

    };

    inline void addHits_(AhoCorasickAmbiguous& fuzzyAC, const AhoCorasickAmbiguous::FuzzyACPattern& pattern, const AhoCorasickAmbiguous::PeptideDB& pep_DB, const String& prot, const String& full_prot, SignedSize idx_prot, Int offset, FoundProteinFunctor& func_threads) const
    {
      fuzzyAC.setProtein(prot);
      while (fuzzyAC.findNext(pattern))
      {
        const seqan::Peptide& tmp_pep = pep_DB[fuzzyAC.getHitDBIndex()];
        func_threads.addHit(fuzzyAC.getHitDBIndex(), idx_prot, length(tmp_pep), full_prot, fuzzyAC.getHitProteinPosition() + offset);
      }
    }

    void updateMembers_() override;

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
    Int mm_max_;
 };
}

