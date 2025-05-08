// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <atomic>
#include <map>
#include <array>


#ifdef _OPENMP 
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;


  char const* const PeptideIndexing::AUTO_MODE = "auto";
  const std::array<std::string, (Size)PeptideIndexing::Unmatched::SIZE_OF_UNMATCHED> PeptideIndexing::names_of_unmatched = { "error", "warn", "remove" };
  const std::array<std::string, (Size)PeptideIndexing::MissingDecoy::SIZE_OF_MISSING_DECOY> PeptideIndexing::names_of_missing_decoy = { "error" , "warn" , "silent" };


  // internal data structure to store match information (not exported)
  struct PeptideProteinMatchInformation
  {
    Hit::T peptide_index; ///< index of the peptide
    Hit::T protein_index; ///< index of the protein the peptide is contained in
    Hit::T position;      ///< the position of the peptide in the protein
    char AABefore; //< the amino acid after the peptide in the protein
    char AAAfter; //< the amino acid before the peptide in the protein

    PeptideProteinMatchInformation(Hit::T pep_index, Hit::T prot_index, Hit::T pep_pos, char aa_before, char aa_after) :
        peptide_index(pep_index), protein_index(prot_index), position(pep_pos), AABefore(aa_before), AAAfter(aa_after) 
    {
      static_assert(sizeof(Hit::T) == 4); // make sure we do not waste huge amounts of memory; there are instances where we find 10^9 hits, which amounts to 14GB of data!
    }

    const std::tuple<const Hit::T&, const Hit::T&, const Hit::T&, const char&, const char&> tie() const
    { // sorting in the order peptide_index, then protein_index, then the rest is paramount for the code below to work!
      return std::tie(peptide_index, protein_index, position, AABefore, AAAfter);
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
    using MapType = std::vector< PeptideProteinMatchInformation >;
    MapType pep_to_prot; ///< peptide index --> protein indices as flat vector
    Size filter_passed{}; ///< number of accepted hits (passing addHit() constraints)
    Size filter_rejected{}; ///< number of rejected hits (not passing addHit())

    ProteaseDigestion enzyme;

  private:
    bool xtandem_; //< are we checking XTandem cleavage rules?

  public:
    explicit FoundProteinFunctor(const ProteaseDigestion& enzyme, bool xtandem) :
      enzyme(enzyme), xtandem_(xtandem)
    {
    }

    void merge(FoundProteinFunctor& other)
    {
      if (pep_to_prot.empty())
      { // first merge is cheap
        pep_to_prot = std::move(other.pep_to_prot);
      }
      else
      {
        pep_to_prot.insert(pep_to_prot.end(), other.pep_to_prot.cbegin(), other.pep_to_prot.cend());
        other.pep_to_prot.clear();
      }
      
      // cheap members
      this->filter_passed += other.filter_passed;
      other.filter_passed = 0;
      this->filter_rejected += other.filter_rejected;
      other.filter_rejected = 0;
    }

    bool validate(const String& seq_prot, const Int position, const Int len_pep, const bool allow_nterm_protein_cleavage)
    {
      const bool ignore_missed_cleavages = true;
      return enzyme.isValidProduct(seq_prot, position, len_pep, ignore_missed_cleavages, allow_nterm_protein_cleavage, xtandem_);
    }

    void addHit(
      const bool is_valid,
      const Hit::T idx_pep,
      const Hit::T idx_prot,
      const Hit::T len_pep,
      const String& seq_prot,
      const Hit::T position)
    {
      //TODO we could read and double-check missed cleavages as well
      if (is_valid)
      {
        // append a PeptideProteinMatchInformation
        pep_to_prot.emplace_back(
          idx_pep,
          idx_prot,
          position,
          (position == 0) ? PeptideEvidence::N_TERMINAL_AA : seq_prot[position - 1],
          (position + len_pep >= seq_prot.size()) ? PeptideEvidence::C_TERMINAL_AA : seq_prot[position + len_pep]
        );
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


  // free function (not exported) used to add hits
  void search(ACTrie& trie, ACTrieState& state, const String& prot, const String& full_prot, size_t prot_offset, Hit::T idx_prot, 
              FoundProteinFunctor& func_threads, const bool allow_nterm_protein_cleavage)
  {
    state.setQuery(prot);
    trie.getAllHits(state);
    
    if (state.hits.empty()) 
    {
      return; // avoid expensive tokenize()
    }

    // by design in AC, hits are ordered by occurrence in protein: duplicate peptides will be consecutive.
    Hit old{};
    bool last_valid = false;
    for (const auto& hit : state.hits)
    {
      Hit::T pos = Hit::T(hit.query_pos + prot_offset);
      if (!(hit.needle_length == old.needle_length && hit.query_pos == old.query_pos))
      { // new stretch in protein. Validate if cutting site
        last_valid = func_threads.validate(full_prot, pos, hit.needle_length, allow_nterm_protein_cleavage);
      }
      func_threads.addHit(last_valid, hit.needle_index, idx_prot, hit.needle_length, full_prot, pos);
      old = hit;
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

    defaults_.setValue("aaa_max", 3, "Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are 'B', 'J', 'Z' and 'X'.");
    defaults_.setMinInt("aaa_max", 0);
    defaults_.setMaxInt("aaa_max", 10);
    
    defaults_.setValue("mismatches_max", 0, "Maximal number of mismatched (mm) amino acids allowed when matching to a protein database."
                                            " The required runtime is exponential in the number of mm's; apply with care."
                                            " MM's are allowed in addition to AAA's.");
    defaults_.setMinInt("mismatches_max", 0);
    defaults_.setMaxInt("mismatches_max", 10);

    defaults_.setValue("IL_equivalent", "false", "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable). Also occurrences of 'J' will be treated as 'I' thus avoiding ambiguous matching.");
    defaults_.setValidStrings("IL_equivalent", { "true", "false" });

    defaults_.setValue("allow_nterm_protein_cleavage", "true", "Allow the protein N-terminus amino acid to clip.");
    defaults_.setValidStrings("allow_nterm_protein_cleavage", { "true", "false" });

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
    allow_nterm_protein_cleavage_ = param_.getValue("allow_nterm_protein_cleavage").toBool();
  }

PeptideIndexing::ExitCodes PeptideIndexing::run(std::vector<FASTAFile::FASTAEntry>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
{
  FASTAContainer<TFI_Vector> protein_container(proteins);
  return run_<TFI_Vector>(protein_container, prot_ids, pep_ids);
}

PeptideIndexing::ExitCodes PeptideIndexing::run(FASTAContainer<TFI_File>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
{
  return run_<TFI_File>(proteins, prot_ids, pep_ids);
}

PeptideIndexing::ExitCodes PeptideIndexing::run(FASTAContainer<TFI_Vector>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
{
  return run_<TFI_Vector>(proteins, prot_ids, pep_ids);
}

const String& PeptideIndexing::getDecoyString() const
{
  return decoy_string_;
}

bool PeptideIndexing::isPrefix() const
{
  return prefix_;
}

template<typename T>
PeptideIndexing::ExitCodes PeptideIndexing::run_(FASTAContainer<T>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
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
  // parsing parameters, correcting XTandem and MSGFPlus parameters
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
  { // fall-back
    OPENMS_LOG_WARN << "Warning: Enzyme name neither given nor deducible from input. Defaulting to Trypsin!" << std::endl;
    enzyme.setEnzyme("Trypsin");
  } 

  bool xtandem_fix_parameters = false;
  bool msgfplus_fix_parameters = false;

  // determine if at least one search engine was XTandem or MSGFPlus to enable special rules
  for (const auto& prot_id : prot_ids)
  {
    String search_engine = prot_id.getOriginalSearchEngineName();
    StringUtils::toUpper(search_engine);
    OPENMS_LOG_INFO << "Peptide identification engine: " << search_engine << std::endl;
    if (search_engine == "XTANDEM" || prot_id.getSearchParameters().metaValueExists("SE:XTandem")) { xtandem_fix_parameters = true; }
    if (search_engine == "MS-GF+" || search_engine == "MSGFPLUS" || prot_id.getSearchParameters().metaValueExists("SE:MS-GF+")) { msgfplus_fix_parameters = true; }
  }

  if (xtandem_fix_parameters)
  {
    OPENMS_LOG_WARN << "X!Tandem detected. Allowing random Asp/Pro cleavage." << std::endl;
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
  { // fall-back
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
  std::map<String, Size> acc_to_prot; // map: accessions --> FASTA protein index
  std::vector<bool> protein_is_decoy; // protein index -> is decoy?
  std::vector<std::string> protein_accessions; // protein index -> accession

  bool invalid_protein_sequence = false; // check for proteins with modifications, i.e. '[' or '(', and throw an exception

  { // new scope - forget data after search
    /*
        Aho Corasick (fast)
    */
    ACTrie ac_trie(aaa_max_, mm_max_);
    SysInfo::MemUsage mu;
    OPENMS_LOG_INFO << "Building trie ...";
    StopWatch s;
    s.start();
    bool peptide_has_X {false}; // if any peptide contains an 'X', we switch off protein-X splitting (see below)
    for (const auto& pep : pep_ids)
    {
      for (const auto& hit : pep.getHits())
      {
        //
        // Warning:
        // do not skip over peptides here, since the results are iterated in the same way
        //
        String seq = hit.getSequence().toUnmodifiedString().remove('*'); // make a copy, i.e. do NOT change the peptide sequence!
        if (IL_equivalent_)                                              // convert L to I;
        {
          seq.substitute('L', 'I');
        }
        peptide_has_X |= seq.has('X');
        ac_trie.addNeedle(seq);
      }
    }
    s.stop();
    OPENMS_LOG_INFO << " done (" << int(s.getClockTime()) << "s)" << std::endl;
    if (ac_trie.getNeedleCount() == 0)
    { // Aho-Corasick will crash if given empty needles as input
      OPENMS_LOG_WARN << "Warning: Peptide identifications have no hits inside! Output will be empty as well." << std::endl;
      return PEPTIDE_IDS_EMPTY;
    }
    s.start();
    OPENMS_LOG_INFO << "Compressing trie to BFS format ..." << std::endl;
    ac_trie.compressTrie();
    s.stop();
    OPENMS_LOG_INFO << " done (" << int(s.getClockTime()) << "s)" << std::endl;
    s.reset();
    OPENMS_LOG_INFO << "Mapping " << ac_trie.getNeedleCount() << " peptides to " << (proteins.size() == PROTEIN_CACHE_SIZE ? "? (unknown number of)" : String(proteins.size())) << " proteins."
                    << std::endl;

    OPENMS_LOG_INFO << "Searching with up to " << aaa_max_ << " ambiguous amino acid(s) and " << mm_max_ << " mismatch(es)!" << std::endl;

    uint16_t count_j_proteins(0);
    bool has_active_data = true; // becomes false if end of FASTA file is reached
    const std::string jumpX(aaa_max_ + mm_max_ + 1, 'X'); // jump over stretches of 'X' which cost a lot of time; +1 because AXXA is a valid hit for aaa_max == 2 (cannot split it)
    // use very large target value for progress if DB size is unknown (did not fit into first chunk)
    this->startProgress(0, proteins.size() == PROTEIN_CACHE_SIZE ? std::numeric_limits<SignedSize>::max() : proteins.size(), "Aho-Corasick");
    std::atomic<int> progress_prots(0);
    
    #pragma omp parallel
    {
      FoundProteinFunctor func_threads(enzyme, xtandem_fix_parameters);
      std::map<String, Size> acc_to_prot_thread; // map: accessions --> FASTA protein index
      ACTrieState ac_state;
      String prot;

      while (true) 
      {
        #pragma omp barrier // all threads need to be here, since we are about to swap protein data
        #pragma omp single
        {
          has_active_data = proteins.activateCache(); // swap in last cache
          protein_accessions.resize(proteins.getChunkOffset() + proteins.chunkSize());
        } // implicit barrier here
        
        if (!has_active_data) break; // leave while-loop
        SignedSize prot_count = (SignedSize)proteins.chunkSize();

        #pragma omp master
        {
          proteins.cacheChunk(PROTEIN_CACHE_SIZE);
          protein_is_decoy.resize(proteins.getChunkOffset() + prot_count);
          for (SignedSize i = 0; i < prot_count; ++i)
          { // do this in master only, to avoid false sharing
            const String& seq = proteins.chunkAt(i).identifier;
            protein_is_decoy[i + proteins.getChunkOffset()] = (prefix_ ? seq.hasPrefix(decoy_string_) : seq.hasSuffix(decoy_string_));
          }
        }

        // search all peptides in each protein
        #pragma omp for schedule(dynamic, 100) nowait
        for (SignedSize i = 0; i < prot_count; ++i)
        {
          ++progress_prots; // atomic
          #ifdef _OPENMP // without OMP, we always set progress
          if (omp_get_thread_num() == 0)
          #endif
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

          const Hit::T prot_idx = Hit::T(i + proteins.getChunkOffset());
          
          // grab #hits before searching protein; we know its a hit if this number changes
          const Size hits_total = func_threads.filter_passed + func_threads.filter_rejected;

          // check if there are stretches of 'X' in the protein, but not in the peptide
          if (!peptide_has_X && prot.has('X'))
          {
            // create chunks of the protein (splitting it at stretches of 'X..X') and feed them to AC one by one
            size_t offset = -1, start = 0;
            while ((offset = prot.find(jumpX, offset + 1)) != std::string::npos)
            {
              //std::cout << "found X..X at " << offset << " in protein " << proteins[i].identifier << "\n";
              search(ac_trie, ac_state, prot.substr(start, offset + jumpX.size() - start), prot, start, prot_idx, func_threads,
                     allow_nterm_protein_cleavage_);
              // skip ahead while we encounter more X...
              while (offset + jumpX.size() < prot.size() && prot[offset + jumpX.size()] == 'X') ++offset;
              start = offset;
              //std::cout << "  new start: " << start << "\n";
            }
            // last chunk
            if (start < prot.size())
            {
              search(ac_trie, ac_state, prot.substr(start), prot, start, prot_idx, func_threads, allow_nterm_protein_cleavage_);
            }
          }
          else // search the whole protein at once
          {
            search(ac_trie, ac_state, prot, prot, 0, prot_idx, func_threads, allow_nterm_protein_cleavage_);
          }
          // was protein found?
          if (hits_total < func_threads.filter_passed + func_threads.filter_rejected)
          {
            protein_accessions[prot_idx] = proteins.chunkAt(i).identifier;
            acc_to_prot_thread[protein_accessions[prot_idx]] = prot_idx;
          }
        } // end parallel FOR

        // join results
        #pragma omp critical(PeptideIndexer_joinAC)
        {
          s.start();
          // hits
          func.merge(func_threads);
          // sort hits by peptide index
          std::sort(func.pep_to_prot.begin(), func.pep_to_prot.end());
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
    
    {
      // count number of peptides found
      // the vector 'pep_to_prot' is sorted by peptide_index, and then by protein_index 
      size_t found_peptide_count{0};
      Hit::T last_peptide_idx = -1;
      for (const auto& hit : func.pep_to_prot)
      {
        if (hit.peptide_index != last_peptide_idx)
        {
          last_peptide_idx = hit.peptide_index;
          ++found_peptide_count;
        }
      }

      OPENMS_LOG_INFO << "\nAho-Corasick done:\n  found " << func.filter_passed << " hits for " << found_peptide_count << " of " << ac_trie.getNeedleCount() << " peptides.\n";
    }
    
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
  std::map<String, Size> runid_to_runidx; // identifier to index
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

  std::map<Size, std::set<Size> > runidx_to_protidx; // in which protID do appear which proteins (according to mapped peptides)

  Size pep_idx(0);
  Size func_hits_idx(0); ///< current position in func.pep_to_prot[] which has a stretch of matches for current pep_idx
  const Size func_hits_size = func.pep_to_prot.size(); 
  for (std::vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
  {
    // which ProteinIdentification does the peptide belong to?
    Size run_idx = runid_to_runidx[it1->getIdentifier()];

    std::vector<PeptideHit>& hits = it1->getHits();

    for (std::vector<PeptideHit>::iterator it_hit = hits.begin(); it_hit != hits.end(); /* no increase here! we might need to erase it; see below */)
    {
      // clear protein accessions
      it_hit->setPeptideEvidences(std::vector<PeptideEvidence>());
      
      //
      // is this a decoy hit?
      //
      bool matches_target(false);
      bool matches_decoy(false);

      size_t prot_count_of_current_pep {0}; ///< protein hits of this peptide
      Hit::T last_prot_index = -1;
      // add new protein references
      // the vector 'pep_to_prot' is sorted by peptide_index, and then by protein_index 
      while ((func_hits_idx < func_hits_size) && (func.pep_to_prot[func_hits_idx].peptide_index == pep_idx))
      {
        const auto& pe = func.pep_to_prot[func_hits_idx];

        if (last_prot_index != pe.protein_index)
        { // span of hits for a new protein starts;
          last_prot_index = pe.protein_index;
          ++prot_count_of_current_pep;
        }

        const String& accession = protein_accessions[pe.protein_index];
        it_hit->addPeptideEvidence(PeptideEvidence(accession, pe.position, pe.position + (int)it_hit->getSequence().size() - 1, pe.AABefore, pe.AAAfter));

        runidx_to_protidx[run_idx].insert(pe.protein_index); // fill protein hits

        if (protein_is_decoy[pe.protein_index])
        {
          matches_decoy = true;
        }
        else
        {
          matches_target = true;
        }
        ++func_hits_idx;
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

      if (prot_count_of_current_pep == 1)
      {
        it_hit->setMetaValue("protein_references", "unique");
        ++stats_matched_unique;
      }
      else if (prot_count_of_current_pep > 1)
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
    const auto& masterset = runidx_to_protidx[run_idx]; // all protein matches from above

    std::vector<ProteinHit>& phits = prot_ids[run_idx].getHits();
    {
      // go through existing protein hits and count orphaned proteins (with no peptide hits)
      std::vector<ProteinHit> orphaned_hits;
      for (std::vector<ProteinHit>::iterator p_hit = phits.begin(); p_hit != phits.end(); ++p_hit)
      {
        const String& acc = p_hit->getAccession();
        if (acc_to_prot.find(acc) == acc_to_prot.end()) // acc_to_prot only contains found proteins from current run
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

/// @endcond
