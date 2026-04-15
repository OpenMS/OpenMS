// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <bit>
#include <cmath>
#include <functional>
#include <mutex>
#include <boost/sort/sort.hpp>


using namespace std;



namespace OpenMS
{

  // Static member definitions
  std::array<double, 128> FragmentIndex::residue_mass_table_{};
  std::once_flag FragmentIndex::mass_table_once_flag_;
  FragmentIndex::IonOffsets FragmentIndex::ion_offsets_{};

  void FragmentIndex::initResidueMassTable_()
  {
    std::call_once(mass_table_once_flag_, []() {
      residue_mass_table_.fill(0.0);
      const ResidueDB* rdb = ResidueDB::getInstance();
      for (char c = 'A'; c <= 'Z'; ++c)
      {
        const Residue* r = rdb->getResidue(static_cast<unsigned char>(c));
        if (r != nullptr)
        {
          residue_mass_table_[static_cast<size_t>(c)] = r->getMonoWeight(Residue::Internal);
        }
      }

      // Precompute ion-type offsets
      ion_offsets_.b_offset = Residue::getInternalToBIon().getMonoWeight();
      ion_offsets_.y_offset = Residue::getInternalToYIon().getMonoWeight();
      ion_offsets_.a_offset = Residue::getInternalToAIon().getMonoWeight();
      ion_offsets_.c_offset = Residue::getInternalToCIon().getMonoWeight();
      ion_offsets_.x_offset = Residue::getInternalToXIon().getMonoWeight();
      ion_offsets_.z_offset = Residue::getInternalToZIon().getMonoWeight();
    });
  }

  void FragmentIndex::initModificationTables_()
  {
    if (mod_tables_initialized_) return;

    fixed_mod_deltas_.fill(0.0);
    fixed_mod_ptrs_.fill(nullptr);
    for (auto& v : variable_mod_table_) v.clear();
    variable_nterm_mods_.clear();
    variable_cterm_mods_.clear();
    fixed_nterm_delta_ = 0.0;
    fixed_cterm_delta_ = 0.0;
    fixed_nterm_mod_ptr_ = nullptr;
    fixed_cterm_mod_ptr_ = nullptr;

    // Build fixed modification lookup
    if (!modifications_fixed_.empty())
    {
      auto fixed_map = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
      for (const auto& [mod_ptr, residue_ptr] : fixed_map.val)
      {
        auto term_spec = mod_ptr->getTermSpecificity();
        double delta = mod_ptr->getDiffMonoMass();
        char origin = mod_ptr->getOrigin();

        if (origin == 'X' || origin == '.') // terminal-only, no specific AA
        {
          if (term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::PROTEIN_N_TERM)
          {
            fixed_nterm_delta_ = delta;
            fixed_nterm_mod_ptr_ = mod_ptr;
          }
          else if (term_spec == ResidueModification::C_TERM || term_spec == ResidueModification::PROTEIN_C_TERM)
          {
            fixed_cterm_delta_ = delta;
            fixed_cterm_mod_ptr_ = mod_ptr;
          }
        }
        else
        {
          // Residue-specific fixed mod (e.g., Carbamidomethyl on C)
          // For ANYWHERE: applies at all matching positions
          // For N_TERM/C_TERM: only at terminal positions (handled during enumeration)
          if (term_spec == ResidueModification::ANYWHERE)
          {
            fixed_mod_deltas_[static_cast<unsigned char>(origin)] = delta;
            fixed_mod_ptrs_[static_cast<unsigned char>(origin)] = mod_ptr;
          }
          else if (term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::PROTEIN_N_TERM)
          {
            fixed_nterm_delta_ = delta;
            fixed_nterm_mod_ptr_ = mod_ptr;
          }
          else if (term_spec == ResidueModification::C_TERM || term_spec == ResidueModification::PROTEIN_C_TERM)
          {
            fixed_cterm_delta_ = delta;
            fixed_cterm_mod_ptr_ = mod_ptr;
          }
        }
      }
    }

    // Build variable modification lookup
    if (!modifications_variable_.empty())
    {
      auto var_map = ModifiedPeptideGenerator::getModifications(modifications_variable_);
      for (const auto& [mod_ptr, residue_ptr] : var_map.val)
      {
        auto term_spec = mod_ptr->getTermSpecificity();
        double delta = mod_ptr->getDiffMonoMass();
        char origin = mod_ptr->getOrigin();
        VarModEntry entry{delta, mod_ptr, term_spec};

        if (origin == 'X' || origin == '.')
        {
          // Pure terminal mod (no specific AA)
          if (term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::PROTEIN_N_TERM)
          {
            variable_nterm_mods_.push_back(entry);
          }
          else if (term_spec == ResidueModification::C_TERM || term_spec == ResidueModification::PROTEIN_C_TERM)
          {
            variable_cterm_mods_.push_back(entry);
          }
        }
        else
        {
          // Residue-specific variable mod — add to the table for this AA
          variable_mod_table_[static_cast<unsigned char>(origin)].push_back(entry);
        }
      }
    }

    mod_tables_initialized_ = true;
  }

  size_t FragmentIndex::buildModSlots_(const char* sequence, size_t seq_len, ModSlot* out_slots,
                                       bool is_protein_nterm, bool is_protein_cterm) const
  {
    size_t n_slots = 0;

    // 1. Pure N-terminal variable mods (not residue-specific, origin='X')
    for (const auto& entry : variable_nterm_mods_)
    {
      // PROTEIN_N_TERM: only for peptides at protein start
      // N_TERM: for any peptide's N-terminus
      if (entry.term_spec == ResidueModification::PROTEIN_N_TERM && !is_protein_nterm) continue;
      if (n_slots >= MAX_MOD_SLOTS) break;
      out_slots[n_slots++] = {ModSlot::NTERM_SLOT, entry.delta_mass, entry.mod_ptr};
    }

    // 2. Per-residue variable mods, left-to-right
    for (size_t i = 0; i < seq_len; ++i)
    {
      unsigned char aa = static_cast<unsigned char>(sequence[i]);
      const auto& var_mods = variable_mod_table_[aa];
      for (const auto& entry : var_mods)
      {
        if (n_slots >= MAX_MOD_SLOTS) break;
        // ANYWHERE: any position
        // N_TERM: peptide N-term (position 0)
        // PROTEIN_N_TERM: only position 0 AND peptide is at protein start
        // C_TERM: peptide C-term (last position)
        // PROTEIN_C_TERM: only last position AND peptide is at protein end
        bool applies = false;
        if (entry.term_spec == ResidueModification::ANYWHERE)
        {
          applies = true;
        }
        else if (entry.term_spec == ResidueModification::N_TERM && i == 0)
        {
          applies = true;
        }
        else if (entry.term_spec == ResidueModification::PROTEIN_N_TERM && i == 0 && is_protein_nterm)
        {
          applies = true;
        }
        else if (entry.term_spec == ResidueModification::C_TERM && i == seq_len - 1)
        {
          applies = true;
        }
        else if (entry.term_spec == ResidueModification::PROTEIN_C_TERM && i == seq_len - 1 && is_protein_cterm)
        {
          applies = true;
        }
        if (applies)
        {
          out_slots[n_slots++] = {static_cast<uint16_t>(i), entry.delta_mass, entry.mod_ptr};
        }
      }
    }

    // 3. Pure C-terminal variable mods (not residue-specific, origin='X')
    for (const auto& entry : variable_cterm_mods_)
    {
      if (entry.term_spec == ResidueModification::PROTEIN_C_TERM && !is_protein_cterm) continue;
      if (n_slots >= MAX_MOD_SLOTS) break;
      out_slots[n_slots++] = {ModSlot::CTERM_SLOT, entry.delta_mass, entry.mod_ptr};
    }

    return n_slots;
  }

  void FragmentIndex::generateFragmentsLightweight_(
    std::vector<Fragment>& fragments,
    const char* sequence,
    size_t seq_len,
    UInt32 peptide_idx,
    double n_term_mod_mass,
    double c_term_mod_mass,
    const double* residue_mod_masses) const
  {
    const double proton = Constants::PROTON_MASS_U;
    const auto& table = residue_mass_table_;

    // Generate prefix ions (b, a, c) - left to right cumulative sum
    // Fragment charge is always 1 for the index (matching original TSG call)
    // Ion index is 1-based: i=0 produces ion 1 (b1/a1/c1), so skip when (i+1) <= min_ion_index_
    if (add_b_ions_ || add_a_ions_ || add_c_ions_)
    {
      {
        constexpr int z = 1;
        double base_mass = proton * z + n_term_mod_mass;
        double cumulative = base_mass;

        for (size_t i = 0; i + 1 < seq_len; ++i)
        {
          double res_mass = table[static_cast<unsigned char>(sequence[i])];
          if (residue_mod_masses) res_mass += residue_mod_masses[i];
          cumulative += res_mass;

          if (i + 1 <= min_ion_index_) continue; // skip ions below min_ion_index

          if (add_b_ions_)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.b_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_a_ions_)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.a_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_c_ions_)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.c_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
        }
      }
    }

    // Generate suffix ions (y, x, z) - right to left cumulative sum
    // Suffix ion index: first iteration produces y1, second y2, etc.
    if (add_y_ions_ || add_x_ions_ || add_z_ions_)
    {
      {
        constexpr int z = 1;
        double base_mass = proton * z + c_term_mod_mass;
        double cumulative = base_mass;
        size_t suffix_ion_num = 0;

        for (size_t j = seq_len; j > 1; --j)
        {
          double res_mass = table[static_cast<unsigned char>(sequence[j - 1])];
          if (residue_mod_masses) res_mass += residue_mod_masses[j - 1];
          cumulative += res_mass;
          ++suffix_ion_num;

          if (suffix_ion_num <= min_ion_index_) continue; // skip ions below min_ion_index

          if (add_y_ions_)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.y_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_x_ions_)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.x_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_z_ions_)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.z_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
        }
      }
    }
  }

#ifdef DEBUG_FRAGMENT_INDEX
    static void print_slice(const std::vector<FragmentIndex::Fragment>& slice, size_t low, size_t high)
    {
      cout << "Slice: ";
      for(size_t i = low; i <= high; i++)
      {
        cout << slice[i].fragment_mz_ << " ";
      }
      cout << endl;
    }


  void FragmentIndex::addSpecialPeptide( OpenMS::AASequence& peptide, Size source_idx)
  {
    float temp_mono = peptide.getMonoWeight();
    fi_peptides_.push_back({AASequence(std::move(peptide)), source_idx,temp_mono});
  }
#endif

  void FragmentIndex::clear()
  {
    fi_fragments_.clear();
    fi_peptides_.clear();
    bucket_min_mz_.clear();
    is_build_ = false;
    mod_tables_initialized_ = false;
  }

  AASequence FragmentIndex::reconstructModifiedSequence(
    const Peptide& peptide,
    const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const
  {
    const string& protein_seq = fasta_entries[peptide.protein_idx].sequence;
    AASequence seq = AASequence::fromString(protein_seq.substr(peptide.sequence_.first, peptide.sequence_.second));

    const bool has_modifications = !(modifications_fixed_.empty() && modifications_variable_.empty());
    if (!has_modifications) return seq;

    // Apply fixed modifications at each residue
    for (size_t i = 0; i < seq.size(); ++i)
    {
      unsigned char aa = static_cast<unsigned char>(protein_seq[peptide.sequence_.first + i]);
      if (fixed_mod_ptrs_[aa] != nullptr)
      {
        seq.setModification(i, fixed_mod_ptrs_[aa]);
      }
    }
    // Apply fixed terminal modifications
    if (fixed_nterm_mod_ptr_ != nullptr)
    {
      seq.setNTerminalModification(fixed_nterm_mod_ptr_);
    }
    if (fixed_cterm_mod_ptr_ != nullptr)
    {
      seq.setCTerminalModification(fixed_cterm_mod_ptr_);
    }

    // Apply variable modifications from bitmask
    if (peptide.mod_bitmask_ != 0)
    {
      const char* seq_ptr = protein_seq.c_str() + peptide.sequence_.first;
      size_t seq_len = peptide.sequence_.second;
      bool is_prot_nterm = (peptide.sequence_.first == 0);
      bool is_prot_cterm = (peptide.sequence_.first + seq_len == protein_seq.size());
      ModSlot slots[MAX_MOD_SLOTS];
      size_t n_slots = buildModSlots_(seq_ptr, seq_len, slots, is_prot_nterm, is_prot_cterm);

      for (size_t s = 0; s < n_slots; ++s)
      {
        if (!(peptide.mod_bitmask_ & (1u << s))) continue;

        if (slots[s].position == ModSlot::NTERM_SLOT)
        {
          seq.setNTerminalModification(slots[s].mod_ptr);
        }
        else if (slots[s].position == ModSlot::CTERM_SLOT)
        {
          seq.setCTerminalModification(slots[s].mod_ptr);
        }
        else
        {
          seq.setModification(slots[s].position, slots[s].mod_ptr);
        }
      }
    }

    return seq;
  }


  /// Compute precursor m/z at charge 1 (M+H)+ directly from amino acid chars.
  /// Formula: (sum_of_internal_masses + H2O + proton) / 1
  static float computePrecursorMzFromChars(const char* seq, size_t len, const std::array<double, 128>& table)
  {
    // M+H = sum(internal masses) + H2O + proton
    static const double water = Residue::getInternalToFull().getMonoWeight(); // 18.0105646834
    double mass = water + Constants::PROTON_MASS_U;
    for (size_t i = 0; i < len; ++i)
    {
      mass += table[static_cast<unsigned char>(seq[i])];
    }
    return static_cast<float>(mass);
  }

  void FragmentIndex::generatePeptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
      initResidueMassTable_();

      const bool has_modifications = !(modifications_fixed_.empty() && modifications_variable_.empty());
      const bool has_variable_mods = !modifications_variable_.empty();

      if (has_modifications)
      {
        initModificationTables_();
      }

      size_t skipped_peptides = 0;

      ProteaseDigestion digestor;
      digestor.setEnzyme(digestion_enzyme_);
      digestor.setMissedCleavages(missed_cleavages_);
      digestor.setSpecificity(enzyme_specificity_);

      OPENMS_LOG_INFO << "Generating peptides..." << std::endl;

      // Per-thread peptide vectors to avoid omp critical
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
#else
      const int num_threads = 1;
#endif
      vector<vector<Peptide>> thread_peptides(num_threads);
      const size_t est_per_thread = (fasta_entries.size() * 5) / num_threads + 1;
      for (int t = 0; t < num_threads; ++t)
        thread_peptides[t].reserve(est_per_thread);

      vector<pair<size_t, size_t>> digested_peptides;
      #pragma omp parallel for private(digested_peptides)
      for (SignedSize protein_idx = 0; protein_idx < (SignedSize)fasta_entries.size(); ++protein_idx)
      {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        digested_peptides.clear();
        const FASTAFile::FASTAEntry& protein = fasta_entries[protein_idx];
        digestor.digestUnmodified(StringView(protein.sequence), digested_peptides, peptide_min_length_, peptide_max_length_);

        for (const pair<size_t, size_t>& digested_peptide : digested_peptides)
        {
          // skip peptides containing unknown or ambiguous AA codes (X, B, Z)
          {
            const auto sub = protein.sequence.substr(digested_peptide.first, digested_peptide.second);
            if (sub.find_first_of("XBZ") != string::npos)
            {
              #pragma omp atomic
              skipped_peptides++;
              continue;
            }
          }

          const char* seq_ptr = protein.sequence.c_str() + digested_peptide.first;
          size_t seq_len = digested_peptide.second;

          // Compute base precursor mass from lookup table (includes fixed mod deltas)
          static const double water = Residue::getInternalToFull().getMonoWeight();
          double base_mass = water + Constants::PROTON_MASS_U + fixed_nterm_delta_ + fixed_cterm_delta_;
          for (size_t i = 0; i < seq_len; ++i)
          {
            base_mass += residue_mass_table_[static_cast<unsigned char>(seq_ptr[i])]
                       + fixed_mod_deltas_[static_cast<unsigned char>(seq_ptr[i])];
          }

          if (has_variable_mods)
          {
            // Bitmask-based variable modification enumeration
            bool is_prot_nterm = (digested_peptide.first == 0);
            bool is_prot_cterm = (digested_peptide.first + seq_len == protein.sequence.size());
            ModSlot slots[MAX_MOD_SLOTS];
            size_t n_slots = buildModSlots_(seq_ptr, seq_len, slots, is_prot_nterm, is_prot_cterm);

            if (n_slots == 0)
            {
              // No variable mod sites on this peptide — just the fixed-mod version
              float mz = static_cast<float>(base_mass);
              if (peptide_min_mass_ <= mz && mz <= peptide_max_mass_)
              {
                thread_peptides[tid].emplace_back(static_cast<UInt32>(protein_idx), uint32_t(0),
                  std::make_pair(static_cast<uint16_t>(digested_peptide.first),
                                 static_cast<uint16_t>(seq_len)), mz);
              }
            }
            else
            {
              // Pre-compute which slots share a position (conflict groups)
              // Build position-to-slot mapping for conflict detection
              // Two slots conflict if they map to the same residue position
              // (mutually exclusive: at most one variable mod per position)
              uint32_t conflict_mask[MAX_MOD_SLOTS] = {};
              for (size_t a = 0; a < n_slots; ++a)
              {
                for (size_t b = a + 1; b < n_slots; ++b)
                {
                  if (slots[a].position == slots[b].position)
                  {
                    conflict_mask[a] |= (1u << b);
                    conflict_mask[b] |= (1u << a);
                  }
                }
              }

              // Enumerate all valid bitmask subsets
              uint32_t max_bitmask = (1u << n_slots);
              for (uint32_t bitmask = 0; bitmask < max_bitmask; ++bitmask)
              {
                // Check max variable mods constraint
                unsigned int popcount = std::popcount(bitmask);
                if (popcount > max_variable_mods_per_peptide_) continue;

                // Check position conflicts: no two set bits can map to the same position
                bool conflict = false;
                for (size_t s = 0; s < n_slots && !conflict; ++s)
                {
                  if ((bitmask & (1u << s)) && (bitmask & conflict_mask[s] & ~(1u << s)))
                  {
                    conflict = true;
                  }
                }
                if (conflict) continue;

                // Compute variant precursor mass
                double variant_mass = base_mass;
                for (size_t s = 0; s < n_slots; ++s)
                {
                  if (bitmask & (1u << s))
                  {
                    variant_mass += slots[s].delta_mass;
                  }
                }

                float mz = static_cast<float>(variant_mass);
                if (peptide_min_mass_ <= mz && mz <= peptide_max_mass_)
                {
                  thread_peptides[tid].emplace_back(static_cast<UInt32>(protein_idx), bitmask,
                    std::make_pair(static_cast<uint16_t>(digested_peptide.first),
                                   static_cast<uint16_t>(seq_len)), mz);
                }
              }
            }
          }
          else if (has_modifications)
          {
            // Fixed mods only — no variable mods to enumerate
            float mz = static_cast<float>(base_mass);
            if (peptide_min_mass_ <= mz && mz <= peptide_max_mass_)
            {
              thread_peptides[tid].emplace_back(static_cast<UInt32>(protein_idx), uint32_t(0),
                std::make_pair(static_cast<uint16_t>(digested_peptide.first),
                               static_cast<uint16_t>(seq_len)), mz);
            }
          }
          else
          {
            // No modifications at all
            float unmodified_mz = static_cast<float>(base_mass);
            if (peptide_min_mass_ <= unmodified_mz && unmodified_mz <= peptide_max_mass_)
            {
              thread_peptides[tid].emplace_back(static_cast<UInt32>(protein_idx), uint32_t(0),
                std::make_pair(static_cast<uint16_t>(digested_peptide.first),
                               static_cast<uint16_t>(seq_len)), unmodified_mz);
            }
          }
        }
      }
      if (skipped_peptides > 0)
      {
        OPENMS_LOG_WARN << skipped_peptides << " peptides skipped due to unknown or ambiguous AA (X/B/Z)\n";
      }

      // Merge per-thread peptide vectors
      size_t total_peptides = 0;
      for (int t = 0; t < num_threads; ++t) total_peptides += thread_peptides[t].size();
      fi_peptides_.reserve(total_peptides);
      for (int t = 0; t < num_threads; ++t)
      {
        fi_peptides_.insert(fi_peptides_.end(), thread_peptides[t].begin(), thread_peptides[t].end());
        vector<Peptide>().swap(thread_peptides[t]);
      }

      OPENMS_LOG_INFO << "Sorting peptides..." << std::endl;
      sort(fi_peptides_.begin(), fi_peptides_.end(), [](const Peptide& a, const Peptide& b)
           {
        return std::tie(a.precursor_mz_, a.protein_idx) < std::tie(b.precursor_mz_, b.protein_idx);
           });
      OPENMS_LOG_INFO << "done." << std::endl;
  }

  void FragmentIndex::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
      /// generate all Peptides (also initializes residue mass table and mod tables)
      generatePeptides(fasta_entries);

      const bool has_modifications = !(modifications_fixed_.empty() && modifications_variable_.empty());

      OPENMS_LOG_INFO << "Generating fragments..." << std::endl;

      // Per-thread fragment vectors to avoid omp critical serialization
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
#else
      const int num_threads = 1;
#endif
      const size_t est_per_thread = (fi_peptides_.size() * 2 * peptide_min_length_) / num_threads + 1;
      vector<vector<Fragment>> thread_fragments(num_threads);
      for (int t = 0; t < num_threads; ++t)
      {
        thread_fragments[t].reserve(est_per_thread);
      }

      // Unified fragment generation path for all cases.
      // For modified peptides: reconstruct per-residue deltas from bitmask + mod tables.
      // No AASequence construction, no ModifiedPeptideGenerator.
      #pragma omp parallel for
      for (SignedSize peptide_idx = 0; peptide_idx < (SignedSize)fi_peptides_.size(); peptide_idx++)
      {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        const Peptide& pep = fi_peptides_[peptide_idx];
        const char* seq_ptr = fasta_entries[pep.protein_idx].sequence.c_str() + pep.sequence_.first;
        size_t seq_len = pep.sequence_.second;

        if (!has_modifications || pep.mod_bitmask_ == 0)
        {
          // No modifications or unmodified variant: just fixed mod deltas (if any)
          if (!has_modifications)
          {
            generateFragmentsLightweight_(
              thread_fragments[tid], seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
              0.0, 0.0, nullptr);
          }
          else
          {
            // Fixed mods only — build delta array from fixed_mod_deltas_
            vector<double> mod_masses(seq_len, 0.0);
            bool has_residue_mods = false;
            for (size_t i = 0; i < seq_len; ++i)
            {
              double delta = fixed_mod_deltas_[static_cast<unsigned char>(seq_ptr[i])];
              if (delta != 0.0)
              {
                mod_masses[i] = delta;
                has_residue_mods = true;
              }
            }
            generateFragmentsLightweight_(
              thread_fragments[tid], seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
              fixed_nterm_delta_, fixed_cterm_delta_,
              has_residue_mods ? mod_masses.data() : nullptr);
          }
        }
        else
        {
          // Variable modifications active: reconstruct delta array from bitmask
          vector<double> mod_masses(seq_len, 0.0);
          bool has_residue_mods = false;
          double n_term_mod = fixed_nterm_delta_;
          double c_term_mod = fixed_cterm_delta_;

          // Apply fixed mod deltas first
          for (size_t i = 0; i < seq_len; ++i)
          {
            double delta = fixed_mod_deltas_[static_cast<unsigned char>(seq_ptr[i])];
            if (delta != 0.0)
            {
              mod_masses[i] = delta;
              has_residue_mods = true;
            }
          }

          // Rebuild mod slots and apply variable mods from bitmask
          const string& prot_seq = fasta_entries[pep.protein_idx].sequence;
          bool is_prot_nterm = (pep.sequence_.first == 0);
          bool is_prot_cterm = (pep.sequence_.first + seq_len == prot_seq.size());
          ModSlot slots[MAX_MOD_SLOTS];
          size_t n_slots = buildModSlots_(seq_ptr, seq_len, slots, is_prot_nterm, is_prot_cterm);
          for (size_t s = 0; s < n_slots; ++s)
          {
            if (!(pep.mod_bitmask_ & (1u << s))) continue;

            if (slots[s].position == ModSlot::NTERM_SLOT)
            {
              n_term_mod += slots[s].delta_mass;
            }
            else if (slots[s].position == ModSlot::CTERM_SLOT)
            {
              c_term_mod += slots[s].delta_mass;
            }
            else
            {
              mod_masses[slots[s].position] += slots[s].delta_mass;
              has_residue_mods = true;
            }
          }

          generateFragmentsLightweight_(
            thread_fragments[tid], seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
            n_term_mod, c_term_mod,
            has_residue_mods ? mod_masses.data() : nullptr);
        }
      }

      // Merge per-thread vectors into global fragment array
      size_t total_fragments = 0;
      for (int t = 0; t < num_threads; ++t) total_fragments += thread_fragments[t].size();
      fi_fragments_.reserve(total_fragments);
      for (int t = 0; t < num_threads; ++t)
      {
        fi_fragments_.insert(fi_fragments_.end(), thread_fragments[t].begin(), thread_fragments[t].end());
        vector<Fragment>().swap(thread_fragments[t]);
      }

      OPENMS_LOG_INFO << "Sorting fragments..." << std::endl;

      /// 1.) First all Fragments are sorted by their own mass (parallel via Boost.Sort)
      boost::sort::block_indirect_sort(fi_fragments_.begin(), fi_fragments_.end(), [](const Fragment& a, const Fragment& b)
      {
        return std::tie(a.fragment_mz_, a.peptide_idx_) < std::tie(b.fragment_mz_, b.peptide_idx_);
      });

      // Empty database (no peptide passed length / mass / motif filters): nothing to bucket.
      // Mark as built and return — guards against the OMP loop below dividing by zero
      // when bucketsize_ becomes 0. This is a real risk for immunopeptidomics FASTAs that
      // contain entries shorter than peptide:min_size.
      if (fi_fragments_.empty())
      {
        bucketsize_ = 1; // keep non-zero to preserve bucket-walking loop invariants
        OPENMS_LOG_INFO << "[FragmentIndex] No fragments generated — index is empty." << std::endl;
        is_build_ = true;
        return;
      }

      /// Calculate the bucket size
      bucketsize_ = sqrt(fi_fragments_.size()); //Todo: MSFragger uses a different approach, which might be better
      OPENMS_LOG_INFO << "Creating DB with bucket_size " << bucketsize_ << endl;

      /// 2.) next sort after precursor mass and save the min_mz of each bucket
      #pragma omp parallel for
      for (SignedSize i = 0; i < (SignedSize)fi_fragments_.size(); i += bucketsize_)
      {

        #pragma omp critical
        bucket_min_mz_.emplace_back(fi_fragments_[i].fragment_mz_);

        auto bucket_start = fi_fragments_.begin() + i;
        auto bucket_end = (i + bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : bucket_start + bucketsize_;

//TODO: is this thread safe????
        sort(bucket_start, bucket_end, [](const Fragment& a, const Fragment& b) {
          return a.peptide_idx_ < b.peptide_idx_; // we don´t need a tie, because the idx are unique
        });
      }
      OPENMS_LOG_INFO << "Sorting by bucket min m/z:" << bucketsize_ << endl;
      //Resort in case the parallelization block above messed something up TODO: check if this can happen
      std::sort( bucket_min_mz_.begin(), bucket_min_mz_.end());
      is_build_ = true;
      OPENMS_LOG_INFO << "Fragment index built!" << endl;
  }

  std::pair<size_t, size_t> FragmentIndex::getPeptidesInMassWindow(float precursor_mass,
                                                                   const std::pair<float, float>& window) const
  {
    // Defensive: a reversed window (first > second) yields an empty half-open range.
    // Under normal computeMassWindow_() usage this never happens (lower is always <= 0,
    // upper is always >= 0), but if a caller builds a window by hand, we must not let
    // (second - first) underflow size_t downstream in searchDifferentPrecursorRanges.
    if (window.first > window.second)
    {
      return {0u, 0u};
    }

    auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                    precursor_mass + window.first,
                                    [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
    auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                     precursor_mass + window.second,
                                     [](float b, const Peptide& a) { return b < a.precursor_mz_; });
    return std::make_pair(std::distance(fi_peptides_.begin(), left_it),
                          std::distance(fi_peptides_.begin(), right_it));
  }

  std::pair<float, float> FragmentIndex::computeMassWindow_(float precursor_mass) const
  {
    if (precursor_mass_tolerance_unit_ppm_)
    {
      const float lo = -Math::ppmToMass<float>(static_cast<float>(precursor_mass_tolerance_lower_),
                                               precursor_mass);
      const float hi =  Math::ppmToMass<float>(static_cast<float>(precursor_mass_tolerance_upper_),
                                               precursor_mass);
      return {lo, hi};
    }
    return {-static_cast<float>(precursor_mass_tolerance_lower_),
             static_cast<float>(precursor_mass_tolerance_upper_)};
  }

  vector<FragmentIndex::Hit> FragmentIndex::query(const OpenMS::Peak1D& peak,
                                                  const pair<size_t, size_t>& peptide_idx_range,
                                                  uint16_t peak_charge)
  {
      float adjusted_mass = peak.getMZ() * (float)peak_charge -((peak_charge-1) * Constants::PROTON_MASS_U);

      float frag_tol = fragment_mz_tolerance_unit_ppm_ ? Math::ppmToMass(fragment_mz_tolerance_, adjusted_mass) : fragment_mz_tolerance_;

      auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), adjusted_mass - frag_tol);
      auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), adjusted_mass + frag_tol);

      if (left_it != bucket_min_mz_.begin()) --left_it;

      auto in_range_buckets = make_pair(std::distance(bucket_min_mz_.begin(), left_it), std::distance(bucket_min_mz_.begin(), right_it));

      vector<FragmentIndex::Hit> hits;
      hits.reserve(peptide_idx_range.second - peptide_idx_range.first);


      for (UInt32 j = in_range_buckets.first; j < in_range_buckets.second; j++)
      {
        auto slice_begin = fi_fragments_.begin() + (j*bucketsize_);
        auto slice_end = ((j+1) * bucketsize_) >= fi_fragments_.size() ? fi_fragments_.end() : (fi_fragments_.begin() + ((j+1) * bucketsize_)) ;

        auto left_iter = std::lower_bound(slice_begin, slice_end, peptide_idx_range.first, [](Fragment a, UInt32 b) { return a.peptide_idx_ < b;} );

        while (left_iter != slice_end) // sequential scan
        {
          // peptide_idx_range is half-open [first, second) — stop BEFORE index second.
          if (left_iter->peptide_idx_ >= peptide_idx_range.second) break;

          if ((adjusted_mass >= left_iter->fragment_mz_ - frag_tol ) && adjusted_mass <= (left_iter->fragment_mz_+ frag_tol))
          {

            hits.emplace_back(left_iter->peptide_idx_, left_iter->fragment_mz_);
            #ifdef DEBUG_FRAGMENT_INDEX
            if (left_iter->peptide_idx_ < peptide_idx_range.first || left_iter->peptide_idx_ >= peptide_idx_range.second)
              OPENMS_LOG_WARN << "idx out of range" << endl;
            #endif
          }
          ++left_iter;
        }
      }

      return hits;
  }

  void FragmentIndex::queryPeaks(SpectrumMatchesTopN& candidates, const MSSpectrum& spectrum,
                                const std::pair<size_t, size_t>& candidates_range,
                                const int16_t isotope_error,
                                const uint16_t precursor_charge)
  {


      for (const Peak1D& peak : spectrum)
      {
        vector<Hit> query_hits;
        uint16_t actual_max = std::min(precursor_charge, max_fragment_charge_);
        for (uint16_t fragment_charge = 1; fragment_charge <= actual_max; fragment_charge++)
        {
          query_hits = query(peak, candidates_range, fragment_charge);

          for (const auto& hit : query_hits)
          {
            {
              size_t idx = hit.peptide_idx - candidates_range.first;

              auto& source = candidates.hits_[idx];
              if (source.num_matched_ == 0)
              {
                source.precursor_charge_ = precursor_charge;
                source.peptide_idx_ = hit.peptide_idx;
                source.isotope_error_ = isotope_error;
              }
              ++source.num_matched_;
            }
          }
        }
      }
  }

  void FragmentIndex::trimHits(OpenMS::FragmentIndex::SpectrumMatchesTopN& init_hits) const
  {
      if (init_hits.hits_.size() > max_processed_hits_)
      {
        std::partial_sort(init_hits.hits_.begin(), init_hits.hits_.begin() + max_processed_hits_, init_hits.hits_.end(), [](const SpectrumMatch& a,const SpectrumMatch& b){
          if (a.num_matched_ != b.num_matched_)
          {
            return a.num_matched_ > b.num_matched_;
          }
          else
          {
            // Prefer isotope_error close to 0: abs(isotope_error), then isotope_error, then precursor_charge
            const auto abs_iso_a = a.isotope_error_ < 0 ? -a.isotope_error_ : a.isotope_error_;
            const auto abs_iso_b = b.isotope_error_ < 0 ? -b.isotope_error_ : b.isotope_error_;
            if (abs_iso_a != abs_iso_b) return abs_iso_a < abs_iso_b;
            if (a.isotope_error_ != b.isotope_error_) return a.isotope_error_ < b.isotope_error_;
            return a.precursor_charge_ < b.precursor_charge_;
          }
        });

        init_hits.hits_.resize(max_processed_hits_);
      }
      else
      {
        std::sort(init_hits.hits_.begin(), init_hits.hits_.end(), [](const SpectrumMatch& a, const SpectrumMatch& b) {
          if (a.num_matched_ != b.num_matched_)
          {
            return a.num_matched_ > b.num_matched_;
          }
          else
          {
            // Prefer isotope_error close to 0: abs(isotope_error), then isotope_error, then precursor_charge
            const auto abs_iso_a = a.isotope_error_ < 0 ? -a.isotope_error_ : a.isotope_error_;
            const auto abs_iso_b = b.isotope_error_ < 0 ? -b.isotope_error_ : b.isotope_error_;
            if (abs_iso_a != abs_iso_b) return abs_iso_a < abs_iso_b;
            if (a.isotope_error_ != b.isotope_error_) return a.isotope_error_ < b.isotope_error_;
            return a.precursor_charge_ < b.precursor_charge_;
          }
        });
      }
      if (init_hits.hits_.size() > 0  )
      {
        if (init_hits.hits_[0].num_matched_ < min_matched_peaks_)
          init_hits.hits_.resize(0);
      }


      for (auto hit_iter = init_hits.hits_.rbegin(); hit_iter != init_hits.hits_.rend(); ++hit_iter)
      {
        if (hit_iter->num_matched_ >= min_matched_peaks_)           // search for the first element that should be included
        {
          init_hits.hits_.resize(init_hits.hits_.size() - (distance(init_hits.hits_.rbegin(), hit_iter)));
          break;
        }
      }
      /* alternative code
       * auto it_zero = std::lower_bound(init_hits.hits_.begin(), init_hits.hits_.end(), min_matched_peaks_ , [](const SpectrumMatch& sm, uint32_t b){
return sm.num_matched_ > b;
});

if (it_zero != init_hits.hits_.end() && it_zero->num_matched_ == 0)
{
init_hits.hits_.erase(it_zero, init_hits.hits_.end());
}
       * */

  }

  void FragmentIndex::searchDifferentPrecursorRanges(const MSSpectrum& spectrum,
                                                     float precursor_mass,
                                                     SpectrumMatchesTopN& sms,
                                                     uint16_t charge)
  {
    // Open mode absorbs isotope shifts into the wide window — no per-isotope iteration.
    const bool open_mode = isOpenSearchMode_();
    const int16_t iso_lo = open_mode ? 0 : min_isotope_error_;
    const int16_t iso_hi = open_mode ? 0 : max_isotope_error_;

    for (int16_t isotope_error = iso_lo; isotope_error <= iso_hi; ++isotope_error)
    {
      const float shifted_mass = precursor_mass
        + static_cast<float>(isotope_error) * static_cast<float>(Constants::C13C12_MASSDIFF_U);

      const auto window = computeMassWindow_(shifted_mass);

      SpectrumMatchesTopN candidates_iso_error;
      auto candidates_range = getPeptidesInMassWindow(shifted_mass, window);
      // candidates_range is half-open [first, second) — size the hits vector exactly for
      // (second - first) entries. queryPeaks indexes via (peptide_idx - first), and the
      // loop in FragmentIndex::query() stops strictly before peptide_idx == second.
      candidates_iso_error.hits_.resize(candidates_range.second - candidates_range.first);

      queryPeaks(candidates_iso_error, spectrum, candidates_range, isotope_error, charge);

      sms += candidates_iso_error;
    }
  }

  void FragmentIndex::querySpectrum(const OpenMS::MSSpectrum& spectrum,
                                    OpenMS::FragmentIndex::SpectrumMatchesTopN& sms)
  {
      if (!isBuild())
      {
        OPENMS_LOG_WARN << "FragmentIndex not yet build \n";
        return;
      }

      if (spectrum.empty() || (spectrum.getMSLevel() != 2))
      {
        return;
      }

      const auto& precursor = spectrum.getPrecursors();
      if (precursor.size() != 1)
      {
        OPENMS_LOG_WARN << "Number of precursors is not equal 1 \n";
        return;
      }

      // two posible modes. Precursor has a charge or we test all possible charges
      vector<size_t> charges;
      //cout << "precursor charge = " << precursor[0].getCharge() << endl;
      if (precursor[0].getCharge())
      {
        //cout << "precursor charge found" << endl;
        charges.push_back(precursor[0].getCharge());
      }
      else
      {
        for (uint16_t i = min_precursor_charge_; i <= max_precursor_charge_; i++)
        {
          charges.push_back(i);
        }
      }
      // loop over all PRECURSOR-charges



      for (uint16_t charge : charges)
      {
        SpectrumMatchesTopN candidates_charge;
        float mz;
        mz = (float)precursor[0].getMZ() * charge - ((charge-1) * Constants::PROTON_MASS_U);
        searchDifferentPrecursorRanges(spectrum, mz, candidates_charge, charge);

        sms += candidates_charge;
      }
      trimHits(sms);
  }




  FragmentIndex::FragmentIndex() : DefaultParamHandler("FragmentIndex")
  {
    defaults_.setValue("ions:add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("ions:add_y_ions", {"true","false"});
    
    defaults_.setValue("ions:add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("ions:add_b_ions", {"true","false"});
    
    defaults_.setValue("ions:add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("ions:add_a_ions", {"true","false"});
    
    defaults_.setValue("ions:add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("ions:add_c_ions", {"true","false"});
    
    defaults_.setValue("ions:add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("ions:add_x_ions", {"true","false"});
    
    defaults_.setValue("ions:add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("ions:add_z_ions", {"true","false"});
    defaults_.setSectionDescription("ions", "Theoretical ion series toggles");


    defaults_.setValue("precursor:mass_tolerance_lower", 20.0,
                       "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                       "is [-lower, +upper] around the precursor). "
                       "When strongly asymmetric, also review precursor:isotope_error_min.");
    defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
    defaults_.setValue("precursor:mass_tolerance_upper", 20.0,
                       "Upper-side precursor-mass tolerance (positive magnitude).");
    defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", {"ppm", "Da"});

    defaults_.setValue("fragment:mass_tolerance", 10.0, "Fragment mass tolerance");
    std::vector<std::string> fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");
    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    defaults_.setValue("precursor:min_charge", 2, "min precursor charge");
    defaults_.setValue("precursor:max_charge", 5, "max precursor charge");

    defaults_.setValue("fragment:min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment:max_mz", 2000, "Maximal fragment mz for database");
    defaults_.setValue("fragment:min_ion_index", 2, "Ions with index less than or equal to this value are not added to the fragment index (use 0 to include all ions; 2 skips b1/b2/y1/y2). Low-index ions are often noisy and unreliable.");
    defaults_.setMinInt("fragment:min_ion_index", 0);

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    defaults_.setValue("enzyme", "Trypsin", "Enzyme for digestion");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));


    defaults_.setValue("peptide:missed_cleavages", 1, "Missed cleavages for digestion");
    defaults_.setValue("peptide:enzyme_specificity", "full",
      "Enzyme cleavage specificity required for both peptide termini.\n"
      "  'full' : both termini must be enzyme-specific (canonical, e.g. tryptic).\n"
      "  'semi' : only one terminus needs to be enzyme-specific (semi-tryptic).\n"
      "  'none' : no enzyme constraint at either terminus; every substring of length\n"
      "           [min_size, max_size] is enumerated. This is the canonical setting for\n"
      "           immunopeptidomics (e.g. HLA peptides 8..12mers). For very large search\n"
      "           spaces consider tightening 'peptide:min_size'/'peptide:max_size'.");
    defaults_.setValidStrings("peptide:enzyme_specificity", {"full", "semi", "none"});
    defaults_.setValue("peptide:min_size", 7, "Minimal peptide length for database");
    defaults_.setValue("peptide:max_size", 40, "Maximal peptide length for database");

    defaults_.setValue("peptide:min_mass", 100, "Minimal peptide mass for database");
    defaults_.setValue("peptide:max_mass", 9000, "Maximal peptide mass for database"); //Todo: set unlimited option


    is_build_ = false; // TODO: remove this and build on construction

    //Search-related params

    defaults_.setValue("fragment:min_matched_ions", 5, "Minimal number of matched ions to report a PSM");
    defaults_.setValue("precursor:isotope_error_min", -1, "Minimum allowed precursor isotope error");
    defaults_.setValue("precursor:isotope_error_max", 1, "Maximum allowed precursor isotope error");
    
    defaults_.setValue("fragment:max_charge", 2, "max fragment charge");
    defaults_.setValue("scoring:max_candidates_per_spectrum", 50, "The number of initial hits for which we calculate a score");
    defaults_.setSectionDescription("scoring", "Search/Scoring Limits");

    //defaults from the searchEngine that are not needed for this class, but otherwise we would generate a warning
    defaults_.setValue("decoys", "false", "Should decoys be generated?");
    defaults_.setValidStrings("decoys", {"true","false"} );
    defaults_.setValue("annotate:PSM",  std::vector<std::string>{"ALL"}, "Annotations added to each PSM.");
    defaults_.setValidStrings("annotate:PSM",
                              std::vector<std::string>{
                                "ALL",
                                Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM,
                                Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
                                Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
                                Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION}
    );
    defaults_.setValue("report:top_hits", 1, "Maximum number of top scoring hits per spectrum that are reported.");
    defaults_.setSectionDescription("report", "Reporting Options");
    defaults_.setValue("peptide:motif", "", "If set, only peptides that contain this motif (provided as RegEx) will be considered.");
    defaults_.setSectionDescription("peptide", "Peptide Options");

    IntList isotopes = {0, 1};
    defaults_.setValue("precursor:isotopes", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)");

    defaultsToParam_();
}

  void FragmentIndex::updateMembers_()
  {
    add_b_ions_ = param_.getValue("ions:add_b_ions").toBool();
    add_y_ions_ = param_.getValue("ions:add_y_ions").toBool();
    add_a_ions_ = param_.getValue("ions:add_a_ions").toBool();
    add_c_ions_ = param_.getValue("ions:add_c_ions").toBool();
    add_x_ions_ = param_.getValue("ions:add_x_ions").toBool();
    add_z_ions_ = param_.getValue("ions:add_z_ions").toBool();
    digestion_enzyme_ = param_.getValue("enzyme").toString();
    enzyme_specificity_ = EnzymaticDigestion::getSpecificityByName(
      param_.getValue("peptide:enzyme_specificity").toString());
    missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_min_mass_ = param_.getValue("peptide:min_mass");
    peptide_max_mass_ = param_.getValue("peptide:max_mass");
    peptide_min_length_ = param_.getValue("peptide:min_size");
    peptide_max_length_ = param_.getValue("peptide:max_size");
    fragment_min_mz_ = param_.getValue("fragment:min_mz");
    fragment_max_mz_ = param_.getValue("fragment:max_mz");
    min_ion_index_ = param_.getValue("fragment:min_ion_index");
    
    precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
    precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
    precursor_mass_tolerance_unit_ppm_ = param_.getValue("precursor:mass_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_ = param_.getValue("fragment:mass_tolerance");
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment:mass_tolerance_unit").toString() == "ppm";

    // Validation — setMinFloat(0.0) rejects negatives via checkDefaults_, but NaN/+inf slip past
    // (NaN < 0 is false). NaN would break lower_bound's strict-weak-ordering.
    if (!std::isfinite(precursor_mass_tolerance_lower_) ||
        !std::isfinite(precursor_mass_tolerance_upper_))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "precursor:mass_tolerance_lower and mass_tolerance_upper must be finite");
    }
    if (precursor_mass_tolerance_lower_ + precursor_mass_tolerance_upper_ <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "precursor window has zero width (lower + upper must be > 0)");
    }

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    max_variable_mods_per_peptide_ = param_.getValue("modifications:variable_max_per_peptide");

    min_matched_peaks_ = param_.getValue("fragment:min_matched_ions");
    min_isotope_error_ = param_.getValue("precursor:isotope_error_min");
    max_isotope_error_ = param_.getValue("precursor:isotope_error_max");
    min_precursor_charge_ = param_.getValue("precursor:min_charge");
    max_precursor_charge_ = param_.getValue("precursor:max_charge");
    max_fragment_charge_ = param_.getValue("fragment:max_charge");
    max_processed_hits_ = param_.getValue("scoring:max_candidates_per_spectrum");

    if (isOpenSearchMode_())
    {
      OPENMS_LOG_WARN << "[FragmentIndex] Open-search mode auto-triggered: window [-"
                      << precursor_mass_tolerance_lower_ << ", +"
                      << precursor_mass_tolerance_upper_ << "] "
                      << (precursor_mass_tolerance_unit_ppm_ ? "ppm" : "Da")
                      << " exceeds threshold. Isotope-error iteration collapses to [0, 0]."
                      << std::endl;
    }
  }
 
  bool FragmentIndex::isBuild() const
  {
    return is_build_;
  }

  const vector<FragmentIndex::Peptide>& FragmentIndex::getPeptides() const
  {
    return fi_peptides_;
  }

}
