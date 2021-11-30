// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>

using namespace OpenMS;
using namespace std;


  char const* const PeptideIndexing::AUTO_MODE = "auto";
  const std::array<std::string, (Size)PeptideIndexing::Unmatched::SIZE_OF_UNMATCHED> PeptideIndexing::names_of_unmatched = { "error", "warn", "remove" };
  const std::array<std::string, (Size)PeptideIndexing::MissingDecoy::SIZE_OF_MISSING_DECOY> PeptideIndexing::names_of_missing_decoy = { "error" , "warn" , "silent" };


  // internal datastructure to store match information (not exported)
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

  // internal functor (not exported)
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
      //TODO we could read and double-check missed cleavages as well
      if (enzyme_.isValidProduct(seq_prot, position, len_pep, true, true, xtandem_))
      {
        PeptideProteinMatchInformation match
        {
          idx_prot,
          position,
          (position == 0) ? PeptideEvidence::N_TERMINAL_AA : seq_prot[position - 1],
          (position + len_pep >= seq_prot.size()) ?
                          PeptideEvidence::C_TERMINAL_AA :
                          seq_prot[position + len_pep]
        };
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


  // free function (not exproted) used to add hits
  void addHits_(AhoCorasickAmbiguous& fuzzyAC, const AhoCorasickAmbiguous::FuzzyACPattern& pattern, const AhoCorasickAmbiguous::PeptideDB& pep_DB, const String& prot, const String& full_prot, SignedSize idx_prot, Int offset, FoundProteinFunctor& func_threads)
  {
    fuzzyAC.setProtein(prot);
    while (fuzzyAC.findNext(pattern))
    {
      const seqan::Peptide& tmp_pep = pep_DB[fuzzyAC.getHitDBIndex()];
      func_threads.addHit(fuzzyAC.getHitDBIndex(), idx_prot, length(tmp_pep), full_prot, fuzzyAC.getHitProteinPosition() + offset);
    }
  }

  PeptideIndexing::PeptideIndexing()
    : DefaultParamHandler("PeptideIndexing")
  {

    defaults_.setValue("decoy_string", "", "String that was appended (or prefixed - see 'decoy_string_position' flag below) to the accessions in the protein database to indicate decoy proteins. If empty (default), it's determined automatically (checking for common terms, both as prefix and suffix).");

    defaults_.setValue("decoy_string_position", "prefix", "Is the 'decoy_string' prepended (prefix) or appended (suffix) to the protein accession? (ignored if decoy_string is empty)");
    defaults_.setValidStrings("decoy_string_position", { "prefix", "suffix" });

    defaults_.setValue("missing_decoy_action", names_of_missing_decoy[(Size)MissingDecoy::IS_ERROR], "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message), 'silent' (no action is taken, not even a warning)");
    defaults_.setValidStrings("missing_decoy_action", std::vector<std::string>(names_of_missing_decoy.begin(), names_of_missing_decoy.end()));

    defaults_.setValue("enzyme:name", AUTO_MODE, "Enzyme which determines valid cleavage sites - e.g. trypsin cleaves after lysine (K) or arginine (R), but not before proline (P). Default: deduce from input");

    StringList enzymes{};
    ProteaseDB::getInstance()->getAllNames(enzymes);
    enzymes.emplace(enzymes.begin(), AUTO_MODE); // make it the first item

    defaults_.setValidStrings("enzyme:name", ListUtils::create<std::string>(enzymes));

    defaults_.setValue("enzyme:specificity", AUTO_MODE, "Specificity of the enzyme. Default: deduce from input."
      "\n  '" + EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL] + "': both internal cleavage sites must match."
      "\n  '" + EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_SEMI] + "': one of two internal cleavage sites must match."
      "\n  '" + EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_NONE] + "': allow all peptide hits no matter their context (enzyme is irrelevant).");

    defaults_.setValidStrings("enzyme:specificity", {AUTO_MODE,
                                                     EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL],
                                                     EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_SEMI],
                                                     EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_NONE]});

    defaults_.setValue("write_protein_sequence", "false", "If set, the protein sequences are stored as well.");
    defaults_.setValidStrings("write_protein_sequence", { "true", "false" });

    defaults_.setValue("write_protein_description", "false", "If set, the protein description is stored as well.");
    defaults_.setValidStrings("write_protein_description", { "true", "false" });

    defaults_.setValue("keep_unreferenced_proteins", "false", "If set, protein hits which are not referenced by any peptide are kept.");
    defaults_.setValidStrings("keep_unreferenced_proteins", { "true", "false" });

    defaults_.setValue("unmatched_action", names_of_unmatched[(Size)Unmatched::IS_ERROR], "If peptide sequences cannot be matched to any protein: 1) raise an error; 2) warn (unmatched PepHits will miss target/decoy annotation with downstream problems); 3) remove the hit.");
    defaults_.setValidStrings("unmatched_action", std::vector<std::string>(names_of_unmatched.begin(), names_of_unmatched.end()));

    defaults_.setValue("aaa_max", 3, "Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are B, J, Z and X!");
    defaults_.setMinInt("aaa_max", 0);
    defaults_.setMaxInt("aaa_max", 10);
    
    defaults_.setValue("mismatches_max", 0, "Maximal number of mismatched (mm) amino acids allowed when matching to a protein database."
                                            " The required runtime is exponential in the number of mm's; apply with care."
                                            " MM's are allowed in addition to AAA's.");
    defaults_.setMinInt("mismatches_max", 0);
    defaults_.setMaxInt("mismatches_max", 10);

    defaults_.setValue("IL_equivalent", "false", "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable). Also occurrences of 'J' will be treated as 'I' thus avoiding ambiguous matching.");
    defaults_.setValidStrings("IL_equivalent", { "true", "false" });

    defaultsToParam_();
  }

  PeptideIndexing::~PeptideIndexing() = default;


  void PeptideIndexing::updateMembers_()
  {
    decoy_string_ = param_.getValue("decoy_string").toString();
    prefix_ = (param_.getValue("decoy_string_position") == "prefix" ? true : false);
    missing_decoy_action_ = (MissingDecoy)Helpers::indexOf(names_of_missing_decoy, param_.getValue("missing_decoy_action"));
    enzyme_name_ = param_.getValue("enzyme:name").toString();
    enzyme_specificity_ = param_.getValue("enzyme:specificity").toString();

    write_protein_sequence_ = param_.getValue("write_protein_sequence").toBool();
    write_protein_description_ = param_.getValue("write_protein_description").toBool();
    keep_unreferenced_proteins_ = param_.getValue("keep_unreferenced_proteins").toBool();
    unmatched_action_ = (Unmatched)Helpers::indexOf(names_of_unmatched, param_.getValue("unmatched_action"));
    IL_equivalent_ = param_.getValue("IL_equivalent").toBool();
    aaa_max_ = static_cast<Int>(param_.getValue("aaa_max"));
    mm_max_ = static_cast<Int>(param_.getValue("mismatches_max"));
  }

const String &PeptideIndexing::getDecoyString() const
{
  return decoy_string_;
}

bool PeptideIndexing::isPrefix() const
{
  return prefix_;
}

template<typename T>
PeptideIndexing::ExitCodes PeptideIndexing::run(FASTAContainer<T>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
{
  if ((enzyme_name_ == "Chymotrypsin" || enzyme_name_ == "Chymotrypsin/P" || enzyme_name_ == "TrypChymo")
    && IL_equivalent_)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "The used enzyme " + enzyme_name_ + "differentiates between I and L, therefore the IL_equivalent option cannot be used.");
  }
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
  if (!enzyme_name_.empty() && (enzyme_name_.compare(AUTO_MODE) != 0))
  { // use param (not empty, not 'auto')
    enzyme.setEnzyme(enzyme_name_);
  }
  else if (!prot_ids.empty() && prot_ids[0].getSearchParameters().digestion_enzyme.getName() != "unknown_enzyme")
  { // take from meta (this assumes all runs used the same enzyme)
    OPENMS_LOG_INFO << "Info: using '" << prot_ids[0].getSearchParameters().digestion_enzyme.getName() << "' as enzyme (obtained from idXML) for digestion." << std::endl;
    enzyme.setEnzyme(&prot_ids[0].getSearchParameters().digestion_enzyme);
  }
  else
  { // fallback
    OPENMS_LOG_WARN << "Warning: Enzyme name neither given nor deduceable from input. Defaulting to Trypsin!" << std::endl;
    enzyme.setEnzyme("Trypsin");
  } 

  bool xtandem_fix_parameters = false;
  bool msgfplus_fix_parameters = false;

  // determine if at least one search engine was xtandem or MSGFPlus to enable special rules
  for (const auto& prot_id : prot_ids)
  {
    String search_engine = prot_id.getOriginalSearchEngineName();
    StringUtils::toUpper(search_engine);
    OPENMS_LOG_INFO << "Peptide identification engine: " << search_engine << std::endl;
    if (search_engine == "XTANDEM" || prot_id.getSearchParameters().metaValueExists("SE:XTandem")) { xtandem_fix_parameters = true; }
    if (search_engine == "MS-GF+" || search_engine == "MSGFPLUS" || prot_id.getSearchParameters().metaValueExists("SE:MS-GF+")) { msgfplus_fix_parameters = true; }
  }

  // including MSGFPlus -> Trypsin/P as enzyme
  if (msgfplus_fix_parameters && enzyme.getEnzymeName() == "Trypsin")
  {
    OPENMS_LOG_WARN << "MSGFPlus detected but enzyme cutting rules were set to Trypsin. Correcting to Trypsin/P to cope with special cutting rule in MSGFPlus." << std::endl;
    enzyme.setEnzyme("Trypsin/P");
  }

  OPENMS_LOG_INFO << "Enzyme: " << enzyme.getEnzymeName() << std::endl;

  if (!enzyme_specificity_.empty() && (enzyme_specificity_.compare(AUTO_MODE) != 0))
  { // use param (not empty and not 'auto')
    enzyme.setSpecificity(ProteaseDigestion::getSpecificityByName(enzyme_specificity_));
  }
  else if (!prot_ids.empty() && prot_ids[0].getSearchParameters().enzyme_term_specificity != ProteaseDigestion::SPEC_UNKNOWN)
  { // deduce from data ('auto')
    enzyme.setSpecificity(prot_ids[0].getSearchParameters().enzyme_term_specificity);
    OPENMS_LOG_INFO << "Info: using '" << EnzymaticDigestion::NamesOfSpecificity[prot_ids[0].getSearchParameters().enzyme_term_specificity] << "' as enzyme specificity (obtained from idXML) for digestion." << std::endl;
  }
  else
  { // fallback
    OPENMS_LOG_WARN << "Warning: Enzyme specificity neither given nor present in the input file. Defaulting to 'full'!" << std::endl;
    enzyme.setSpecificity(ProteaseDigestion::SPEC_FULL);
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

    for (std::vector<PeptideHit>::iterator it_hit = hits.begin(); it_hit != hits.end(); /* no increase here! we might need to skip it; see below */)
    {
      // clear protein accessions
      it_hit->setPeptideEvidences(std::vector<PeptideEvidence>());
      
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
        PeptideEvidence pe(accession, it_i->position, it_i->position + (int)it_hit->getSequence().size() - 1, it_i->AABefore, it_i->AAAfter);
        it_hit->addPeptideEvidence(pe);

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
      ++pep_idx; // next hit

      if (matches_decoy && matches_target)
      {
        it_hit->setMetaValue("target_decoy", "target+decoy");
        ++stats_count_m_td;
      }
      else if (matches_target)
      {
        it_hit->setMetaValue("target_decoy", "target");
        ++stats_count_m_t;
      }
      else if (matches_decoy)
      {
        it_hit->setMetaValue("target_decoy", "decoy");
        ++stats_count_m_d;
      } // else: could match to no protein (i.e. both are false)
      //else ... // not required (handled below; see stats_unmatched);

      if (prot_indices.size() == 1)
      {
        it_hit->setMetaValue("protein_references", "unique");
        ++stats_matched_unique;
      }
      else if (prot_indices.size() > 1)
      {
        it_hit->setMetaValue("protein_references", "non-unique");
        ++stats_matched_multi;
      }
      else
      {
        ++stats_unmatched;
        if (stats_unmatched < 15) OPENMS_LOG_INFO << "Unmatched peptide: " << it_hit->getSequence() << "\n";
        else if (stats_unmatched == 15) OPENMS_LOG_INFO << "Unmatched peptide: ...\n";
        if (unmatched_action_ == Unmatched::REMOVE)
        {
          it_hit = hits.erase(it_hit);
          continue; // already points to the next hit
        }
        else
        {
          it_hit->setMetaValue("protein_references", "unmatched");
        }
      }

      ++it_hit; // next hit
    } // all hits

  } // next PepID

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
    String msg("No peptides were matched to the decoy portion of the database! Did you provide the correct concatenated database? Are your 'decoy_string' (=" + decoy_string_ + ") and 'decoy_string_position' (=" + std::string(param_.getValue("decoy_string_position")) + ") settings correct?");
    if (missing_decoy_action_ == MissingDecoy::IS_ERROR)
    {
      OPENMS_LOG_ERROR << "Error: " << msg << "\nSet 'missing_decoy_action' to 'warn' if you are sure this is ok!\nAborting ..." << std::endl;
      has_error = true;
    }
    else if (missing_decoy_action_ == MissingDecoy::WARN)
    {
      OPENMS_LOG_WARN << "Warn: " << msg << "\nSet 'missing_decoy_action' to 'error' if you want to elevate this to an error!" << std::endl;
    }
    else // silent
    {
    }
  }

  if (stats_unmatched > 0)
  {
    OPENMS_LOG_ERROR << "PeptideIndexer found unmatched peptides, which could not be associated to a protein.\n";
    if (unmatched_action_ == Unmatched::IS_ERROR)
    {
      OPENMS_LOG_ERROR
        << "Potential solutions:\n"
        << "   - check your FASTA database is identical to the search DB (or use 'auto')\n"
        << "   - set 'enzyme:specificity' and 'enzyme:name' to 'auto' to match the parameters of the search engine\n"
        << "   - increase 'aaa_max' to allow more ambiguous amino acids\n"
        << "   - as a last resort: use the 'unmatched_action' option to accept or even remove unmatched peptides\n"
        << "     (note that unmatched peptides cannot be used for FDR calculation or quantification)\n";
      has_error = true;
    }
    else if (unmatched_action_ == Unmatched::WARN)
    {
      OPENMS_LOG_ERROR << "  Warning: " << stats_unmatched << " unmatched hits have been found, but were not removed!\n"
        << "These are not annotated with target/decoy information and might lead to issues with downstream tools (such as FDR).\n"
        << "Switch to '" << names_of_unmatched[(Size)Unmatched::REMOVE] << "' if you want to avoid these problems.\n";
    }
    else if (unmatched_action_ == Unmatched::REMOVE)
    {
      OPENMS_LOG_ERROR << "  Warning: " << stats_unmatched <<" unmatched hits have been removed!\n"
                        << "Make sure that these hits are actually a violation of the cutting rules by inspecting the database!\n";
      if (xtandem_fix_parameters) OPENMS_LOG_ERROR << "Since the results are from X!Tandem, this is probably ok (check anyways).\n";
    }
    else
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
  }
  

  if (has_error)
  {
    OPENMS_LOG_ERROR << "Result files will be written, but PeptideIndexer will exit with an error code." << std::endl;
    return UNEXPECTED_RESULT;
  }
  return EXECUTION_OK;
}

// specialize templates here so they get instantiated only once (speeds up compile time and reduces memory usage)
template OPENMS_DLLAPI PeptideIndexing::ExitCodes PeptideIndexing::run<class TFI_File>(FASTAContainer<TFI_File>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);
template OPENMS_DLLAPI PeptideIndexing::ExitCodes PeptideIndexing::run<class TFI_Vector>(FASTAContainer<TFI_Vector>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

PeptideIndexing::ExitCodes PeptideIndexing::run(std::vector<FASTAFile::FASTAEntry>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
{
  FASTAContainer<TFI_Vector> protein_container(proteins);
  return run<TFI_Vector>(protein_container, prot_ids, pep_ids);
}

/// @endcond
