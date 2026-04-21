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
#include <unordered_map>
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

  std::vector<double> FragmentIndex::computeSnesSigmaDeltaSet_(bool include_prot_nterm_mods,
                                                                bool include_prot_cterm_mods) const
  {
    // Precondition: initModificationTables_() has been called.
    // updateMembers_() guarantees this: it resets mod_tables_initialized_ and
    // calls initModificationTables_() at the end, so any setParameters() call
    // will have populated the tables before this helper is invoked.

    // Collect all per-mod deltas that should participate in the enumeration.
    // Respect term-specificity flags: PROTEIN_N_TERM and PROTEIN_C_TERM mods
    // are gated by the caller.
    std::vector<double> eligible_deltas;
    auto collect = [&](const std::vector<VarModEntry>& entries, bool /*residue_bound*/)
    {
      for (const auto& e : entries)
      {
        if (e.term_spec == ResidueModification::PROTEIN_N_TERM && !include_prot_nterm_mods) continue;
        if (e.term_spec == ResidueModification::PROTEIN_C_TERM && !include_prot_cterm_mods) continue;
        eligible_deltas.push_back(e.delta_mass);
      }
    };
    collect(variable_nterm_mods_, /*residue_bound=*/false);
    collect(variable_cterm_mods_, /*residue_bound=*/false);
    for (const auto& per_aa : variable_mod_table_)
    {
      collect(per_aa, /*residue_bound=*/true);
    }

    // Dedup eligible_deltas: multiple mods sharing the same mass shift (e.g.
    // Deamidated(N) and Deamidated(Q), both +0.984016 Da) produce identical
    // BFS paths. Collapsing them to a single representative halves BFS work.
    std::sort(eligible_deltas.begin(), eligible_deltas.end());
    eligible_deltas.erase(
        std::unique(eligible_deltas.begin(), eligible_deltas.end(),
                    [](double a, double b) { return std::abs(a - b) < 1e-6; }),
        eligible_deltas.end());

    // Enumerate multisets of size 0..max_per_peptide with replacement from
    // eligible_deltas. Store unique Σ values within a 1e-6 Da tolerance
    // (absorbs FP error across ~16 summed deltas in double precision).
    std::vector<double> result;
    result.push_back(0.0);

    if (eligible_deltas.empty() || max_variable_mods_per_peptide_ == 0)
    {
      return result;
    }

    // BFS: at level m, we have all Σ values reachable with exactly m mods.
    // We iterate m = 1..max_per_peptide, extending each level by one delta.
    std::vector<double> previous_level{0.0};
    for (size_t m = 1; m <= max_variable_mods_per_peptide_; ++m)
    {
      std::vector<double> next_level;
      next_level.reserve(previous_level.size() * eligible_deltas.size());
      for (double prev : previous_level)
      {
        for (double d : eligible_deltas)
        {
          next_level.push_back(prev + d);
        }
      }
      // Dedup within next_level and against result.
      std::sort(next_level.begin(), next_level.end());
      next_level.erase(
          std::unique(next_level.begin(), next_level.end(),
                      [](double a, double b) { return std::abs(a - b) < 1e-6; }),
          next_level.end());
      for (double v : next_level)
      {
        // Insert into result if not already present (within tolerance).
        auto it = std::lower_bound(result.begin(), result.end(), v - 1e-6);
        if (it == result.end() || std::abs(*it - v) >= 1e-6)
        {
          result.insert(it, v);
        }
      }
      previous_level = std::move(next_level);
    }

    return result;
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
    // Thin wrapper: forward class flags to the series-explicit implementation.
    generateFragmentsForSeries_(fragments, sequence, seq_len, peptide_idx,
                                n_term_mod_mass, c_term_mod_mass, residue_mod_masses,
                                add_b_ions_, add_a_ions_, add_c_ions_,
                                add_y_ions_, add_x_ions_, add_z_ions_);
  }

  void FragmentIndex::generateFragmentsForSeries_(
    std::vector<Fragment>& fragments,
    const char* sequence,
    size_t seq_len,
    UInt32 peptide_idx,
    double n_term_mod_mass,
    double c_term_mod_mass,
    const double* residue_mod_masses,
    bool add_b,
    bool add_a,
    bool add_c,
    bool add_y,
    bool add_x,
    bool add_z) const
  {
    const double proton = Constants::PROTON_MASS_U;
    const auto& table = residue_mass_table_;

    // Generate prefix ions (b, a, c) - left to right cumulative sum
    // Fragment charge is always 1 for the index (matching original TSG call)
    // Ion index is 1-based: i=0 produces ion 1 (b1/a1/c1), so skip when (i+1) <= min_ion_index_
    if (add_b || add_a || add_c)
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

          if (add_b)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.b_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_a)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.a_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_c)
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
    if (add_y || add_x || add_z)
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

          if (add_y)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.y_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_x)
          {
            float mz = static_cast<float>((cumulative + ion_offsets_.x_offset) / z);
            if (mz >= fragment_min_mz_ && mz <= fragment_max_mz_)
              fragments.emplace_back(peptide_idx, mz);
          }
          if (add_z)
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
    // swap-to-empty ensures the underlying heap capacity is actually released.
    // std::vector::clear() alone only resets size; capacity stays resident, which
    // defeats the M1 optimization of freeing the fragment index before downstream
    // PeptideIndexing / Aho-Corasick peaks (fi_fragments_ alone is hundreds of MB
    // on human-proteome builds). The swap idiom is the only portable way to force
    // deallocation across libstdc++/libc++/MSVC.
    std::vector<Fragment>().swap(fi_fragments_);
    std::vector<Peptide>().swap(fi_peptides_);
    std::vector<float>().swap(bucket_min_mz_);
    std::vector<uint32_t>().swap(protein_lengths_);
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

    // Apply variable modifications from bitmask.
    // SNES defensive masking: bit 31 (SNES_KIND_BIT_MASK) encodes Single-C-ness,
    // not a slot. Mask it off before iterating slot bits so the loop remains
    // correct when a SNES mother reaches this function. In non-SNES mode
    // mod_bitmask_ is never written with bit 31, so masking is a zero-cost
    // no-op there.
    const uint32_t slot_bits = peptide.mod_bitmask_ & SNES_SLOT_MASK;
    if (slot_bits != 0)
    {
      const char* seq_ptr = protein_seq.c_str() + peptide.sequence_.first;
      size_t seq_len = peptide.sequence_.second;
      bool is_prot_nterm = (peptide.sequence_.first == 0);
      bool is_prot_cterm = (peptide.sequence_.first + seq_len == protein_seq.size());
      ModSlot slots[MAX_MOD_SLOTS];
      size_t n_slots = buildModSlots_(seq_ptr, seq_len, slots, is_prot_nterm, is_prot_cterm);

      for (size_t s = 0; s < n_slots; ++s)
      {
        if (!(slot_bits & (1u << s))) continue;

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

  int FragmentIndex::realizeSNESLength(const Peptide& mother,
                                       const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                       double target_mh_plus,
                                       double tolerance_lower_magnitude,
                                       double tolerance_upper_magnitude,
                                       bool tolerance_ppm) const
  {
    if (!is_snes_mode_) return -1;

    const std::string& protein_seq = fasta_entries[mother.protein_idx].sequence;
    const bool is_single_c = isSingleCMother(mother.mod_bitmask_);
    const size_t mother_start = mother.sequence_.first;
    const size_t mother_length = mother.sequence_.second;

    static const double water = Residue::getInternalToFull().getMonoWeight();
    const double base = water + Constants::PROTON_MASS_U
                        + fixed_nterm_delta_ + fixed_cterm_delta_;
    const double tol_lo_da = tolerance_ppm
      ? (tolerance_lower_magnitude * std::max(target_mh_plus, 0.0) * 1e-6)
      : tolerance_lower_magnitude;
    const double tol_hi_da = tolerance_ppm
      ? (tolerance_upper_magnitude * std::max(target_mh_plus, 0.0) * 1e-6)
      : tolerance_upper_magnitude;

    // Scan residues from the anchored terminus outward, accumulating the
    // residue-mass sum incrementally. Single-N scans left-to-right (the first
    // residue is mother_start); Single-C scans right-to-left (the first residue
    // is mother_start + mother_length - 1). At each prefix length k in
    // [peptide_min_length_, mother_length], the corresponding realized sub-peptide
    // has mass (base + cumulative). Pick the length whose signed delta
    // (realized_mass - target) falls in [-tol_lo_da, +tol_hi_da] and is closest
    // to the target in magnitude. Asymmetric bounds preserve calibrated
    // windows (e.g. [100 ppm, 5 ppm] after bias correction).
    double cumulative = 0.0;
    double best_abs_delta = std::numeric_limits<double>::max();
    int best_length = -1;
    for (size_t k = 0; k < mother_length; ++k)
    {
      const size_t res_idx = is_single_c
        ? mother_start + (mother_length - 1 - k)
        : mother_start + k;
      const unsigned char aa = static_cast<unsigned char>(protein_seq[res_idx]);
      cumulative += residue_mass_table_[aa] + fixed_mod_deltas_[aa];

      const size_t length = k + 1;
      if (length < peptide_min_length_) continue;

      const double realized_mass = base + cumulative;
      const double delta = realized_mass - target_mh_plus;        // signed
      if (delta < -tol_lo_da || delta > tol_hi_da) continue;
      const double abs_delta = std::abs(delta);
      if (abs_delta < best_abs_delta)
      {
        best_abs_delta = abs_delta;
        best_length = static_cast<int>(length);
      }
    }
    return best_length;
  }

  AASequence FragmentIndex::reconstructRealizedSubSequence(
      const Peptide& mother,
      const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
      size_t realized_length,
      uint32_t subset_bitmask) const
  {
    const std::string& protein_seq = fasta_entries[mother.protein_idx].sequence;
    const bool is_single_c = isSingleCMother(mother.mod_bitmask_);
    const size_t mother_start = mother.sequence_.first;
    const size_t mother_length = mother.sequence_.second;

    // Single-N: realized = [mother_start, mother_start + realized_length)
    // Single-C: realized = [mother_start + mother_length - realized_length, mother_start + mother_length)
    const size_t realized_start = is_single_c
      ? mother_start + mother_length - realized_length
      : mother_start;

    AASequence seq = AASequence::fromString(protein_seq.substr(realized_start, realized_length));

    const bool has_mods = !(modifications_fixed_.empty() && modifications_variable_.empty());
    if (!has_mods && subset_bitmask == 0) return seq;

    // Fixed residue mods — applied to every residue of the realized sub-peptide.
    for (size_t i = 0; i < seq.size(); ++i)
    {
      const unsigned char aa = static_cast<unsigned char>(protein_seq[realized_start + i]);
      if (fixed_mod_ptrs_[aa] != nullptr)
      {
        seq.setModification(i, fixed_mod_ptrs_[aa]);
      }
    }
    if (fixed_nterm_mod_ptr_ != nullptr) seq.setNTerminalModification(fixed_nterm_mod_ptr_);
    if (fixed_cterm_mod_ptr_ != nullptr) seq.setCTerminalModification(fixed_cterm_mod_ptr_);

    // SNES v1.1: apply variable mods from subset_bitmask.
    if (subset_bitmask != 0)
    {
      const char* seq_ptr = protein_seq.c_str() + realized_start;
      const bool is_prot_nterm = (realized_start == 0);
      const bool is_prot_cterm = (realized_start + realized_length == protein_seq.size());
      ModSlot slots[MAX_MOD_SLOTS];
      size_t n_slots = buildModSlots_(seq_ptr, realized_length, slots, is_prot_nterm, is_prot_cterm);

      for (size_t s = 0; s < n_slots; ++s)
      {
        if (!(subset_bitmask & (1u << s))) continue;
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

  void FragmentIndex::generateSNESMothers_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    // Residue-mass table already initialized by the caller (generatePeptides).
    // Still need the modification tables for fixed-mod deltas.
    const bool has_modifications = !(modifications_fixed_.empty() && modifications_variable_.empty());
    const bool has_variable_mods = !modifications_variable_.empty();

    if (has_modifications)
    {
      initModificationTables_();
    }

    if (has_variable_mods)
    {
      // SNES v1 restriction: variable modifications are not enumerated on mother
      // peptides. Rationale: which slots survive realization depends on the C-terminus
      // (for Single-N) or N-terminus (for Single-C), which is only resolved in the
      // realization step; tracking per-length slot masks would negate the SNES size
      // win. Terminal + internal variable-mod support is earmarked for v1.2+.
      std::string ignored = ListUtils::concatenate(modifications_variable_, ", ");
      OPENMS_LOG_INFO << "[FragmentIndex] SNES v1.1: variable modifications are not enumerated on mother peptides "
                      << "(mother index stores unmodified mothers); query-time subset enumeration "
                      << "applies them to realized sub-peptides using configured mods: " << ignored << std::endl;
    }

    std::atomic<size_t> skipped_peptides{0};

    OPENMS_LOG_INFO << "Generating SNES mother peptides..." << std::endl;

#ifdef _OPENMP
    const int num_threads = omp_get_max_threads();
#else
    const int num_threads = 1;
#endif
    vector<vector<Peptide>> thread_peptides(num_threads);
    // Heuristic reserve: ~2 * avg_protein_length mothers per protein (Single-N + Single-C).
    const size_t est_per_thread =
        (fasta_entries.size() * 2 * std::max<size_t>(peptide_max_length_, 20)) / num_threads + 1;
    for (int t = 0; t < num_threads; ++t) thread_peptides[t].reserve(est_per_thread);

    // Mother mass = water + proton + sum(residue_masses + fixed residue-mod deltas)
    //             + fixed_nterm_delta + fixed_cterm_delta
    //
    // This is always an upper bound on the realized sub-peptide mass (trimming
    // residues off one end only removes mass; fixed terminal mods either apply to
    // the realized peptide too or are absent — they never add mass to the realization
    // that wasn't in the mother). The one-sided precursor filter
    // `mother_mass >= observed_P - tol` therefore admits every realizable candidate;
    // the ProSEAlgorithm realization step exactly rematches the observed precursor
    // and removes the false positives that survived the filter.
    static const double water = Residue::getInternalToFull().getMonoWeight();
    const double base_sum_constants = water + Constants::PROTON_MASS_U
                                      + fixed_nterm_delta_ + fixed_cterm_delta_;

    #pragma omp parallel for
    for (SignedSize protein_idx = 0; protein_idx < (SignedSize)fasta_entries.size(); ++protein_idx)
    {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      const FASTAFile::FASTAEntry& protein = fasta_entries[protein_idx];
      const std::string& seq = protein.sequence;
      const size_t L = seq.size();
      if (L < peptide_min_length_) continue;

      // Fast path: for the overwhelmingly common case of a protein with no
      // ambiguous residues, a single linear scan tells us we can skip the per-
      // mother X/B/Z check entirely. For the rare protein with X/B/Z, we fall
      // back to checking only the next-bad-position via std::string::find_first_of
      // starting from the mother's start index — still O(length) in the worst
      // case, but only if the mother actually overlaps a bad region.
      //
      // Follow-up: a protein-wide skip here is overly conservative. A max-length
      // mother spanning an `X` is discarded, but shorter realizations with the
      // same anchor *not* covering the `X` could be valid. CodeRabbit #1.
      // Correct fix: split the protein into contiguous unambiguous spans and
      // generate mothers per-span. Deferred to v1.1.
      const bool protein_has_ambiguous = seq.find_first_of("XBZ") != std::string::npos;

      // Honor peptide:max_size=0 as "no maximum" (the documented semantics of
      // the non-SNES path). Using raw peptide_max_length_ in std::min would give
      // length 0 and an empty SNES index.
      const size_t effective_max_length = (peptide_max_length_ == 0) ? L : peptide_max_length_;

      auto emitMother = [&](size_t start, size_t length, bool is_single_c)
      {
        if (length < peptide_min_length_) return;
        const char* seq_ptr = seq.c_str() + start;

        // Reject mothers containing ambiguous codes — any realized sub-peptide
        // spanning an X/B/Z would fail AASequence::fromString downstream.
        if (protein_has_ambiguous)
        {
          const size_t bad = seq.find_first_of("XBZ", start);
          if (bad != std::string::npos && bad < start + length)
          {
            skipped_peptides.fetch_add(1);
            return;
          }
        }

        double mass = base_sum_constants;
        for (size_t k = 0; k < length; ++k)
        {
          const unsigned char aa = static_cast<unsigned char>(seq_ptr[k]);
          mass += residue_mass_table_[aa] + fixed_mod_deltas_[aa];
        }
        const float mz = static_cast<float>(mass);
        // Only the lower bound is safe at mother-generation time: shorter
        // realizations of a mother whose total mass exceeds peptide_max_mass_
        // can still fall within the user's configured mass range. Enforce the
        // upper bound at realization time via the precursor-tolerance window
        // (which is always <= peptide_max_mass_ for observed spectra). CodeRabbit #5.
        if (mz < peptide_min_mass_) return;

        const uint32_t kind_bits = is_single_c ? SNES_KIND_BIT_MASK : 0u;
        thread_peptides[tid].emplace_back(
            static_cast<UInt32>(protein_idx),
            kind_bits,
            std::make_pair(static_cast<uint16_t>(start), static_cast<uint16_t>(length)),
            mz);
      };

      // Single-N mothers: anchored at position i, span the longest possible peptide
      // starting there (capped at effective_max_length). i sweeps [0, L - min_length].
      for (size_t i = 0; i + peptide_min_length_ <= L; ++i)
      {
        const size_t length = std::min<size_t>(effective_max_length, L - i);
        emitMother(i, length, /*is_single_c=*/false);
      }

      // Single-C mothers: anchored at position j (last residue), span the longest
      // possible peptide ending there. j sweeps [min_length - 1, L - 1].
      // When j + 1 <= effective_max_length the mother happens to coincide with a
      // Single-N mother at position 0 — that's harmless redundancy: both emit a
      // different ion series into the index, so there's no duplicate fragment.
      // Guard against min_length=0: j would wrap to SIZE_MAX. Clamp to 1 locally;
      // this is the only code path sensitive to the min_length=0 edge case.
      const size_t snes_min_length = std::max<size_t>(1, peptide_min_length_);
      for (size_t j = snes_min_length - 1; j < L; ++j)
      {
        const size_t length = std::min<size_t>(effective_max_length, j + 1);
        const size_t start = j + 1 - length;
        emitMother(start, length, /*is_single_c=*/true);
      }
    }

    // Merge per-thread buckets (same shape as generatePeptides).
    size_t total = 0;
    for (int t = 0; t < num_threads; ++t) total += thread_peptides[t].size();
    fi_peptides_.reserve(total);
    for (int t = 0; t < num_threads; ++t)
    {
      fi_peptides_.insert(fi_peptides_.end(), thread_peptides[t].begin(), thread_peptides[t].end());
      vector<Peptide>().swap(thread_peptides[t]);
    }

    OPENMS_LOG_INFO << "Sorting SNES mother peptides..." << std::endl;
    sort(fi_peptides_.begin(), fi_peptides_.end(), [](const Peptide& a, const Peptide& b)
    {
      return std::tie(a.precursor_mz_, a.protein_idx) < std::tie(b.precursor_mz_, b.protein_idx);
    });

    OPENMS_LOG_INFO << "Generated " << fi_peptides_.size() << " SNES mothers ("
                    << skipped_peptides.load() << " skipped due to ambiguous residues or mass filter)." << std::endl;
  }

  void FragmentIndex::generatePeptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
      initResidueMassTable_();

      // SNES mode dispatch: for non-specific searches with snes_enabled, switch to
      // mother-peptide indexing instead of the O(L^2) sub-peptide enumeration below.
      // The SNES path has its own mod-table init; everything else (fragment emission,
      // query layer, build-level bucketing) reads the populated fi_peptides_ uniformly.
      if (is_snes_mode_)
      {
        generateSNESMothers_(fasta_entries);
        return;
      }

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
      protein_lengths_.clear();
      protein_lengths_.reserve(fasta_entries.size());
      for (const auto& e : fasta_entries)
      {
        protein_lengths_.push_back(static_cast<uint32_t>(e.sequence.size()));
      }

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

        if (is_snes_mode_)
        {
          // SNES mother: emit a single ion series determined by bit 31 of mod_bitmask_.
          // Variable modifications are disabled in SNES v1 (see generateSNESMothers_);
          // only fixed residue/terminal modifications are applied. The "free" terminus
          // (C-term for Single-N, N-term for Single-C) does not receive its fixed
          // terminal modification at index time — the fragment series we emit doesn't
          // reach that terminus, and the sub-peptide mass is recomputed at realization.
          const bool is_single_c = isSingleCMother(pep.mod_bitmask_);

          vector<double> mod_masses(seq_len, 0.0);
          bool has_residue_mods = false;
          for (size_t i = 0; i < seq_len; ++i)
          {
            const double delta = fixed_mod_deltas_[static_cast<unsigned char>(seq_ptr[i])];
            if (delta != 0.0)
            {
              mod_masses[i] = delta;
              has_residue_mods = true;
            }
          }

          const double n_term_mod = is_single_c ? 0.0 : fixed_nterm_delta_;
          const double c_term_mod = is_single_c ? fixed_cterm_delta_ : 0.0;

          // SNES candidate lookup in querySpectrumSNES_ only targets b-ions (for
          // Single-N mothers) and y-ions (for Single-C mothers) — the other ion
          // series (a/c/x/z) are not targeted and indexing them would be wasted
          // storage at best and source of silent data loss at worst if a user
          // disabled the primary series via ion toggles. Force b-only/y-only
          // regardless of the class add_*_ions_ flags. CodeRabbit #6.
          generateFragmentsForSeries_(
            thread_fragments[tid], seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
            n_term_mod, c_term_mod,
            has_residue_mods ? mod_masses.data() : nullptr,
            /*add_b=*/!is_single_c,
            /*add_a=*/false,
            /*add_c=*/false,
            /*add_y=*/ is_single_c,
            /*add_x=*/false,
            /*add_z=*/false);
        }
        else if (!has_modifications || pep.mod_bitmask_ == 0)
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

          // Rebuild mod slots and apply variable mods from bitmask.
          // SNES defensive masking: bit 31 (SNES_KIND_BIT_MASK) encodes
          // Single-C-ness, not a slot. SNES mothers go through a different
          // branch above, so masking is a no-op in current code, but applying
          // it here makes the invariant self-documenting and prevents a latent
          // bug if a future refactor routes an SNES-kind-marked peptide through
          // this path. In non-SNES mode mod_bitmask_ is never written with bit
          // 31, so masking is zero-cost.
          const uint32_t slot_bits = pep.mod_bitmask_ & SNES_SLOT_MASK;
          const string& prot_seq = fasta_entries[pep.protein_idx].sequence;
          bool is_prot_nterm = (pep.sequence_.first == 0);
          bool is_prot_cterm = (pep.sequence_.first + seq_len == prot_seq.size());
          ModSlot slots[MAX_MOD_SLOTS];
          size_t n_slots = buildModSlots_(seq_ptr, seq_len, slots, is_prot_nterm, is_prot_cterm);
          for (size_t s = 0; s < n_slots; ++s)
          {
            if (!(slot_bits & (1u << s))) continue;

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

    // SNES mode uses querySpectrumSNES_ directly (dispatched in querySpectrum);
    // this function is only reached for non-SNES searches.
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

  void FragmentIndex::querySpectrumSNES_(const MSSpectrum& spectrum,
                                          const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                          SpectrumMatchesTopN& sms)
  {
    // Preconditions checked by the public entry point querySpectrum.

    const auto& precursor = spectrum.getPrecursors()[0];
    vector<uint16_t> charges;
    // Precursor::getCharge() returns signed Int; treat non-positive (0 = unset,
    // rare negative-mode encodings) as "unknown" and fall back to the
    // configured min..max range. Without this guard, static_cast<uint16_t>(-1)
    // would wrap to 65535 and be used as an actual charge downstream.
    if (precursor.getCharge() > 0)
    {
      charges.push_back(static_cast<uint16_t>(precursor.getCharge()));
    }
    else
    {
      for (uint16_t z = min_precursor_charge_; z <= max_precursor_charge_; ++z)
      {
        charges.push_back(z);
      }
    }

    // Phase 1 — byte-count scoring across ALL mothers.
    //
    // Thread-local buffer reused across calls; sized to fi_peptides_.size() and
    // zeroed via .assign. Avoids per-spectrum allocation and keeps the table hot
    // in cache for large indices. Saturates at UINT16_MAX (far above any realistic
    // matched-peak count) to protect against pathological inputs without branch-on-overflow.
    thread_local std::vector<uint16_t> score_table;
    score_table.assign(fi_peptides_.size(), 0);

    // Fragment-charge upper bound for the byte scan. Use the max charge in the
    // `charges` list (the spectrum's known charge, or max_precursor_charge_ when
    // unknown), capped by max_fragment_charge_. This matches the per-iteration
    // `std::min(precursor_charge, max_fragment_charge_)` policy of non-SNES
    // queryPeaks — using the static max_precursor_charge_ would inflate byte
    // counts with matches at fragment charges above the actual precursor charge.
    uint16_t byte_scan_max_frag_charge = 0;
    for (uint16_t c : charges) byte_scan_max_frag_charge = std::max(byte_scan_max_frag_charge, c);
    byte_scan_max_frag_charge = std::min(byte_scan_max_frag_charge, max_fragment_charge_);

    for (const Peak1D& peak : spectrum)
    {
      for (uint16_t frag_charge = 1; frag_charge <= byte_scan_max_frag_charge; ++frag_charge)
      {
        const float adjusted_mass = static_cast<float>(peak.getMZ()) * frag_charge
          - (frag_charge - 1) * static_cast<float>(Constants::PROTON_MASS_U);
        const float frag_tol = fragment_mz_tolerance_unit_ppm_
          ? Math::ppmToMass<float>(static_cast<float>(fragment_mz_tolerance_), adjusted_mass)
          : static_cast<float>(fragment_mz_tolerance_);

        // Bucket-range lookup mirrors query(): bucket_min_mz_ holds the smallest
        // fragment_mz of each bucket, so a peak with mz ∈ [bucket_min, next_bucket_min)
        // falls inside the bucket starting at bucket_min.
        auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(),
                                        adjusted_mass - frag_tol);
        auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(),
                                         adjusted_mass + frag_tol);
        if (left_it != bucket_min_mz_.begin()) --left_it;

        const size_t bucket_begin = std::distance(bucket_min_mz_.begin(), left_it);
        const size_t bucket_end = std::distance(bucket_min_mz_.begin(), right_it);

        for (size_t j = bucket_begin; j < bucket_end; ++j)
        {
          const auto slice_begin = fi_fragments_.begin() + (j * bucketsize_);
          const auto slice_end = ((j + 1) * bucketsize_) >= fi_fragments_.size()
            ? fi_fragments_.end()
            : (fi_fragments_.begin() + ((j + 1) * bucketsize_));

          // No peptide_idx pre-filter: every peptide in the bucket is a potential
          // match regardless of its mother mass. The precursor filter is applied
          // downstream via the fragment-bin-as-precursor trick.
          for (auto it = slice_begin; it != slice_end; ++it)
          {
            if (adjusted_mass >= it->fragment_mz_ - frag_tol
                && adjusted_mass <= it->fragment_mz_ + frag_tol)
            {
              auto& cell = score_table[it->peptide_idx_];
              if (cell < std::numeric_limits<uint16_t>::max()) ++cell;
            }
          }
        }
      }
    }

    // Phase 2 — candidate collection via the fragment-index-as-precursor-filter trick.
    //
    // A Single-N mother of any realized length k has b_k m/z = (M_sub)+H+ − water.
    // A Single-C mother of any realized length k has y_k m/z = (M_sub)+H+ directly.
    // So mothers that could realize the observed precursor mass are exactly the
    // ones with an indexed fragment in a narrow m/z window around those targets.
    //
    // We walk the fragment buckets at each target once per (charge, iso_err) and
    // emit the dedup'd set of matching mother ids whose phase-1 byte score meets
    // the minimum-matched-peaks threshold.
    static const double water = Residue::getInternalToFull().getMonoWeight();

    // Dedup guard (thread_local, sized once per query). The actual per-
    // (charge, iso_err, sigma) reset happens inside the Σ loops below via
    // std::fill — this assign() is only for size-safe initialization of
    // the thread-local buffer when fi_peptides_.size() changes between calls.
    thread_local std::vector<uint8_t> emitted;
    emitted.assign(fi_peptides_.size(), 0);

    // Helper: compute the iso-shifted observed (M+H)+ for a given (charge, iso_err).
    // Used by the subset-enumeration post-pass to reconstruct the realization target.
    auto shifted_mh_for = [&](uint16_t charge, int16_t iso_err) -> float {
      const float mh_plus = static_cast<float>(precursor.getMZ()) * charge
        - (charge - 1) * static_cast<float>(Constants::PROTON_MASS_U);
      return mh_plus + static_cast<float>(iso_err) * static_cast<float>(Constants::C13C12_MASSDIFF_U);
    };

    // Asymmetric tolerance: tol_lo <= 0 (low-side magnitude, sign-flipped),
    // tol_hi >= 0. A match requires (fragment_mz - target_mz) ∈ [tol_lo, tol_hi].
    // Preserves calibrated windows like [100 ppm, 5 ppm] where the symmetric
    // max-collapse over-admitted ~20× on the tighter side.
    auto collect_candidates =
      [&](float target_mz, float tol_lo, float tol_hi, bool expect_single_c,
          int16_t iso_err, uint16_t charge,
          SnesAnchor require_anchor, float sigma_tag)
    {
      auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(),
                                      target_mz + tol_lo);
      auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(),
                                       target_mz + tol_hi);
      if (left_it != bucket_min_mz_.begin()) --left_it;

      const size_t bucket_begin = std::distance(bucket_min_mz_.begin(), left_it);
      const size_t bucket_end = std::distance(bucket_min_mz_.begin(), right_it);

      for (size_t j = bucket_begin; j < bucket_end; ++j)
      {
        const auto slice_begin = fi_fragments_.begin() + (j * bucketsize_);
        const auto slice_end = ((j + 1) * bucketsize_) >= fi_fragments_.size()
          ? fi_fragments_.end()
          : (fi_fragments_.begin() + ((j + 1) * bucketsize_));

        for (auto it = slice_begin; it != slice_end; ++it)
        {
          const float delta = it->fragment_mz_ - target_mz;
          if (delta < tol_lo || delta > tol_hi) continue;

          const UInt32 id = it->peptide_idx_;
          if (emitted[id]) continue;

          const auto& mother = fi_peptides_[id];
          if (isSingleCMother(mother.mod_bitmask_) != expect_single_c) continue;

          // SNES v1.1: anchor-specific filter for PROTEIN_N/C_TERM mod walks.
          if (require_anchor == SnesAnchor::PROT_NTERM && mother.sequence_.first != 0) continue;
          if (require_anchor == SnesAnchor::PROT_CTERM)
          {
            const uint32_t prot_len = protein_lengths_[mother.protein_idx];
            if (static_cast<uint32_t>(mother.sequence_.first) + mother.sequence_.second != prot_len) continue;
          }

          if (score_table[id] < min_matched_peaks_) continue;

          emitted[id] = 1;
          SpectrumMatch sm;
          sm.peptide_idx_ = id;
          sm.num_matched_ = score_table[id];
          sm.isotope_error_ = iso_err;
          sm.precursor_charge_ = charge;
          sm.sigma_delta_ = sigma_tag;
          sms.hits_.push_back(sm);
        }
      }
    };

    // SNES v1.1: precompute the set-difference Σ values that only appear in
    // the protein-anchored-only sets. These drive extra bin walks gated by
    // the SnesAnchor filter in collect_candidates. When no protein-term
    // variable mods are configured, these vectors are empty and the extra
    // walks are skipped entirely.
    auto set_difference_tol = [](const std::vector<double>& A, const std::vector<double>& B) {
      std::vector<double> out;
      out.reserve(A.size());
      for (double a : A)
      {
        bool found = false;
        for (double b : B)
        {
          if (std::abs(a - b) < 1e-6) { found = true; break; }
        }
        if (!found) out.push_back(a);
      }
      return out;
    };
    const std::vector<double> prot_nterm_extra = set_difference_tol(
        snes_sigma_delta_set_with_prot_nterm_, snes_sigma_delta_set_);
    const std::vector<double> prot_cterm_extra = set_difference_tol(
        snes_sigma_delta_set_with_prot_cterm_, snes_sigma_delta_set_);

    for (uint16_t charge : charges)
    {
      const float mh_plus = static_cast<float>(precursor.getMZ()) * charge
        - (charge - 1) * static_cast<float>(Constants::PROTON_MASS_U);

      // Open-search mode (very wide precursor tolerance auto-detected in
      // isOpenSearchMode_()) collapses the isotope-error iteration to a
      // single iso_err == 0 pass — at open-search windows, adding multiples
      // of C13C12_MASSDIFF_U to the target is a no-op on candidate admission
      // (the window already spans many isotope peaks) and just inflates the
      // hit list with duplicate-labelled candidates. Mirrors the non-SNES
      // `queryPeaks` path (FragmentIndex.cpp:1478-1480).
      const bool open_mode = isOpenSearchMode_();
      const int16_t iso_lo = open_mode ? 0 : min_isotope_error_;
      const int16_t iso_hi = open_mode ? 0 : max_isotope_error_;
      for (int16_t iso_err = iso_lo; iso_err <= iso_hi; ++iso_err)
      {
        const float shifted_mh = mh_plus
          + static_cast<float>(iso_err) * static_cast<float>(Constants::C13C12_MASSDIFF_U);

        // Asymmetric precursor tolerance for the bin walk — use computeMassWindow_
        // to get signed (lo <= 0, hi >= 0) Da bounds that respect calibrated
        // asymmetric windows. Previously collapsed to max(lower, upper), which
        // over-admitted by ~20× on the tighter side for calibrated configs
        // like [100 ppm, 5 ppm]. Same reference m/z (shifted_mh) for all
        // targets within this (charge, iso_err) iteration.
        const auto prec_window = computeMassWindow_(shifted_mh);
        const float prec_tol_lo = prec_window.first;   // <= 0
        const float prec_tol_hi = prec_window.second;  // >= 0

        // Baseline Σ loop: walks every mother regardless of protein anchor.
        // Each (charge, iso_err, sigma) triple is independent — reset the dedup
        // guard per iteration so the same mother can emit distinct matches at
        // different Σ values (each represents a distinct variable-mod assignment).
        for (double sigma : snes_sigma_delta_set_)
        {
          // Reset dedup per (charge, iso_err, sigma) combo so the same mother
          // can re-emit at distinct sigma values (each is a distinct match).
          std::fill(emitted.begin(), emitted.end(), 0);

          const float s = static_cast<float>(sigma);

          // Target m/z derivation, accounting for the fact that SNES fragment
          // generation (build()) emits Single-N b-ions with c_term_mod = 0 and
          // Single-C y-ions with n_term_mod = 0:
          //   Realized (M+H)+ = water + proton + fixed_nterm + fixed_cterm + Σ_internal
          //   b_k (Single-N)  = proton + fixed_nterm + Σ_internal
          //                   => (M+H)+ − water − fixed_cterm − Σ
          //   y_k (Single-C)  = water + proton + fixed_cterm + Σ_internal
          //                   => (M+H)+ − fixed_nterm − Σ
          // Missing the fixed-term offsets would shift the lookup target by the
          // corresponding delta and silently miss all candidates when a user
          // configures Acetyl (N-term) / Amidated (C-term) / similar.
          collect_candidates(shifted_mh - static_cast<float>(water)
                                        - static_cast<float>(fixed_cterm_delta_) - s,
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/false, iso_err, charge,
                             SnesAnchor::NONE, s);
          collect_candidates(shifted_mh - static_cast<float>(fixed_nterm_delta_) - s,
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/true, iso_err, charge,
                             SnesAnchor::NONE, s);

          // Supplementary full-length realization at this Σ — recover b_L and y_L
          // cases that aren't in the fragment index (v1 indexes b_1..b_{L-1}
          // and y_1..y_{L-1} only). Match mothers whose precursor_mz_ equals
          // shifted_mh - s within the asymmetric prec tolerance.
          //
          // This matters especially for (a) sub-peptides of length == max_length,
          // which have no "longer-mother" partial-realization alternative, and
          // (b) proteins shorter than max_length, where every mother is at full
          // protein length.
          {
            const float target = shifted_mh - s;
            auto lb = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                        target + prec_tol_lo,
                                        [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
            auto ub = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                        target + prec_tol_hi,
                                        [](float b, const Peptide& a) { return b < a.precursor_mz_; });
            for (auto it = lb; it != ub; ++it)
            {
              const UInt32 id = static_cast<UInt32>(std::distance(fi_peptides_.begin(), it));
              if (emitted[id]) continue;
              if (score_table[id] < min_matched_peaks_) continue;
              emitted[id] = 1;
              SpectrumMatch sm;
              sm.peptide_idx_ = id;
              sm.num_matched_ = score_table[id];
              sm.isotope_error_ = iso_err;
              sm.precursor_charge_ = charge;
              sm.sigma_delta_ = s;
              sms.hits_.push_back(sm);
            }
          }
        }

        // Extra walks for PROTEIN_N_TERM-only Σ values (Single-N mothers at
        // protein position 0 only). Empty when no PROTEIN_N_TERM variable mods
        // are configured.
        for (double sigma : prot_nterm_extra)
        {
          std::fill(emitted.begin(), emitted.end(), 0);
          const float s = static_cast<float>(sigma);
          collect_candidates(shifted_mh - static_cast<float>(water)
                                        - static_cast<float>(fixed_cterm_delta_) - s,
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/false, iso_err, charge,
                             SnesAnchor::PROT_NTERM, s);
          // Supplementary full-length at this Σ, PROT_NTERM-gated.
          {
            const float target = shifted_mh - s;
            auto lb = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                        target + prec_tol_lo,
                                        [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
            auto ub = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                        target + prec_tol_hi,
                                        [](float b, const Peptide& a) { return b < a.precursor_mz_; });
            for (auto it = lb; it != ub; ++it)
            {
              const UInt32 id = static_cast<UInt32>(std::distance(fi_peptides_.begin(), it));
              if (emitted[id]) continue;
              if (fi_peptides_[id].sequence_.first != 0) continue; // PROT_NTERM anchor
              if (isSingleCMother(fi_peptides_[id].mod_bitmask_)) continue; // Single-N only
              if (score_table[id] < min_matched_peaks_) continue;
              emitted[id] = 1;
              SpectrumMatch sm;
              sm.peptide_idx_ = id;
              sm.num_matched_ = score_table[id];
              sm.isotope_error_ = iso_err;
              sm.precursor_charge_ = charge;
              sm.sigma_delta_ = s;
              sms.hits_.push_back(sm);
            }
          }
        }

        // Extra walks for PROTEIN_C_TERM-only Σ values.
        for (double sigma : prot_cterm_extra)
        {
          std::fill(emitted.begin(), emitted.end(), 0);
          const float s = static_cast<float>(sigma);
          collect_candidates(shifted_mh - static_cast<float>(fixed_nterm_delta_) - s,
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/true, iso_err, charge,
                             SnesAnchor::PROT_CTERM, s);
          // Supplementary full-length at this Σ, PROT_CTERM-gated.
          {
            const float target = shifted_mh - s;
            auto lb = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                        target + prec_tol_lo,
                                        [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
            auto ub = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                        target + prec_tol_hi,
                                        [](float b, const Peptide& a) { return b < a.precursor_mz_; });
            for (auto it = lb; it != ub; ++it)
            {
              const UInt32 id = static_cast<UInt32>(std::distance(fi_peptides_.begin(), it));
              if (emitted[id]) continue;
              if (!isSingleCMother(fi_peptides_[id].mod_bitmask_)) continue; // Single-C only
              const uint32_t prot_len = protein_lengths_[fi_peptides_[id].protein_idx];
              if (static_cast<uint32_t>(fi_peptides_[id].sequence_.first)
                  + fi_peptides_[id].sequence_.second != prot_len) continue;
              if (score_table[id] < min_matched_peaks_) continue;
              emitted[id] = 1;
              SpectrumMatch sm;
              sm.peptide_idx_ = id;
              sm.num_matched_ = score_table[id];
              sm.isotope_error_ = iso_err;
              sm.precursor_charge_ = charge;
              sm.sigma_delta_ = s;
              sms.hits_.push_back(sm);
            }
          }
        }
      }
    }

    // SNES v1.1 subset enumeration: expand each (mother, Σ) hit in sms.hits_
    // into one SpectrumMatch per valid variable-mod subset on the realized
    // sub-peptide. Σ=0 hits pass through unchanged (bitmask=0). Per-mother
    // cap of 16 subsets across all (k, Σ) tuples prevents degenerate blowup.
    {
      std::vector<SpectrumMatch> expanded;
      expanded.reserve(sms.hits_.size());
      std::unordered_map<size_t, size_t> subsets_per_mother;

      const auto& fasta_entries_ref = fasta_entries;

      for (const SpectrumMatch& sm_raw : sms.hits_)
      {
        if (sm_raw.sigma_delta_ == 0.0f)
        {
          // No variable mods: pass through unchanged, bitmask already 0.
          expanded.push_back(sm_raw);
          continue;
        }

        const Peptide& mother = fi_peptides_[sm_raw.peptide_idx_];
        const double iso_shifted_target =
            static_cast<double>(shifted_mh_for(sm_raw.precursor_charge_, sm_raw.isotope_error_))
            - static_cast<double>(sm_raw.sigma_delta_);
        const int realized_len = realizeSNESLength(
            mother, fasta_entries_ref, iso_shifted_target,
            precursor_mass_tolerance_lower_,
            precursor_mass_tolerance_upper_,
            precursor_mass_tolerance_unit_ppm_);
        if (realized_len < 0) continue;

        const std::string& protein_seq = fasta_entries_ref[mother.protein_idx].sequence;
        const bool is_single_c = isSingleCMother(mother.mod_bitmask_);
        const size_t sub_start = is_single_c
            ? mother.sequence_.first + mother.sequence_.second - static_cast<size_t>(realized_len)
            : mother.sequence_.first;
        const size_t sub_len = static_cast<size_t>(realized_len);
        const char* seq_ptr = protein_seq.c_str() + sub_start;
        const bool is_prot_nterm = (sub_start == 0);
        const bool is_prot_cterm = (sub_start + sub_len == protein_seq.size());

        ModSlot slots[MAX_MOD_SLOTS];
        const size_t n_slots = buildModSlots_(seq_ptr, sub_len, slots, is_prot_nterm, is_prot_cterm);

        // Enumerate bitmask subsets 1..(2^n_slots - 1) with constraints:
        //   - popcount ≤ max_variable_mods_per_peptide_
        //   - no two active bits share a residue position
        //   - Σ_subset ≈ sigma_delta_ within 1e-6 Da
        // Cap: ≤ 16 subsets per mother (across all k, Σ tuples in this query).
        if (n_slots == 0) continue;
        // Enumerate all non-empty bitmasks in [1, 2^n_slots - 1]. Use uint64_t for
        // the upper bound to avoid (1u << 32) UB at n_slots=32, and to include
        // bitmask 0xFFFFFFFF at that boundary. The iteration variable is still
        // uint32_t since n_slots ≤ 32 (bounded by MAX_MOD_SLOTS).
        const uint64_t max_bitmask64 = (n_slots >= 32)
            ? (uint64_t{1} << 32)
            : (uint64_t{1} << n_slots);
        // bm is uint64_t so the terminating increment past UINT32_MAX doesn't
        // wrap to 0 and re-enter the loop (n_slots == 32 isn't reachable in
        // SNES mode — bit 31 is reserved for the kind flag — but the
        // defensive width keeps the loop terminating cleanly on any
        // future widening of MAX_MOD_SLOTS).
        for (uint64_t bm = 1; bm < max_bitmask64; ++bm)
        {
          if (static_cast<size_t>(std::popcount(bm)) > max_variable_mods_per_peptide_) continue;

          // Position-conflict check. Use 1ULL for the shift so n_slots up to
          // 63 remain well-defined after the bm → uint64_t widening above.
          bool conflict = false;
          for (size_t a = 0; a < n_slots && !conflict; ++a)
          {
            if (!(bm & (uint64_t{1} << a))) continue;
            for (size_t b = a + 1; b < n_slots; ++b)
            {
              if (!(bm & (uint64_t{1} << b))) continue;
              if (slots[a].position == slots[b].position) { conflict = true; break; }
            }
          }
          if (conflict) continue;

          // Σ match check.
          // sigma_delta_ is stored as float; tolerate float→double rounding
          // by using 1e-4 Da (~0.1 mDa) rather than 1e-6. Minimum modification
          // delta separation in Unimod is ≥1 mDa, so 0.1 mDa is safe.
          double subset_sigma = 0.0;
          for (size_t s = 0; s < n_slots; ++s)
          {
            if (bm & (uint64_t{1} << s)) subset_sigma += slots[s].delta_mass;
          }
          if (std::abs(subset_sigma - static_cast<double>(sm_raw.sigma_delta_)) >= 1e-4) continue;

          // Per-mother cap.
          size_t& count = subsets_per_mother[sm_raw.peptide_idx_];
          if (count >= 16)
          {
            OPENMS_LOG_DEBUG << "[FragmentIndex] SNES per-mother subset cap "
                             << "hit for mother_idx=" << sm_raw.peptide_idx_
                             << " at sigma_delta=" << sm_raw.sigma_delta_ << std::endl;
            break;
          }

          SpectrumMatch sm_variant = sm_raw;
          // subset_bitmask_ is uint32_t; bm is uint64_t for termination safety
          // but bm's value is always < 2^n_slots ≤ 2^32, so this narrowing is
          // value-preserving.
          sm_variant.subset_bitmask_ = static_cast<uint32_t>(bm);
          expanded.push_back(sm_variant);
          ++count;
        }
      }
      sms.hits_ = std::move(expanded);
    }

    // Cap the candidate set at `max_processed_hits_` (same policy as the
    // non-SNES path via queryPeaks→trimHits). Without this cap every candidate
    // pays TSG + AASequence + HyperScore cost in scoreSpectraAgainstIndex_;
    // ProSE's coarser fragment bucketing (sqrt(N) per bucket vs MetaMorpheus's
    // 1-mDa fixed bins) admits more candidates per spectrum than MM's tight
    // bin lookup, so the cap is necessary to keep per-spectrum work bounded.
    // Top-K by num_matched is safe here because the fragment-index-as-precursor
    // filter has already tightly constrained candidates to mothers compatible
    // with the observed precursor — no length bias like in the v1 design.
    trimHits(sms);
  }

  void FragmentIndex::querySpectrum(const OpenMS::MSSpectrum& spectrum,
                                    OpenMS::FragmentIndex::SpectrumMatchesTopN& sms)
  {
    // Backward-compatible 2-arg overload. Delegates to the 3-arg overload
    // with an empty FASTA. Safe for non-SNES and SNES-without-var-mods
    // callers — the subset-enumeration block in querySpectrumSNES_ only
    // dereferences fasta_entries when sm.sigma_delta_ != 0, which cannot
    // occur when modifications_variable_ is empty. SNES + var-mods callers
    // must use the 3-arg overload; this guard rejects them explicitly
    // rather than producing undefined behavior.
    if (is_snes_mode_ && !modifications_variable_.empty())
    {
      OPENMS_LOG_ERROR << "[FragmentIndex] querySpectrum called without FASTA in SNES mode "
                          "with variable modifications — results would be undefined. "
                          "Use querySpectrum(spectrum, fasta_entries, sms) instead.\n";
      return;
    }
    static const std::vector<FASTAFile::FASTAEntry> empty_fasta;
    querySpectrum(spectrum, empty_fasta, sms);
  }

  void FragmentIndex::querySpectrum(const OpenMS::MSSpectrum& spectrum,
                                    const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
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

      if (is_snes_mode_)
      {
        querySpectrumSNES_(spectrum, fasta_entries, sms);
        return;
      }

      // Non-SNES path: fasta_entries not needed.
      // two posible modes. Precursor has a charge or we test all possible charges
      vector<size_t> charges;
      if (precursor[0].getCharge())
      {
        charges.push_back(precursor[0].getCharge());
      }
      else
      {
        for (uint16_t i = min_precursor_charge_; i <= max_precursor_charge_; i++)
        {
          charges.push_back(i);
        }
      }

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
    // Default iso range [0, +2]: Orbitrap/QExactive/tims monoisotopic peak picking
    // fails predominantly *upward* (picks the +1 or +2 isotope instead of the true
    // monoisotopic). Symmetric ranges like [-1, +1] waste a query slot on the
    // rare downward mispick. Matches MetaMorpheus/MSFragger defaults.
    defaults_.setValue("precursor:isotope_error_min", 0, "Minimum allowed precursor isotope error");
    defaults_.setValue("precursor:isotope_error_max", 2, "Maximum allowed precursor isotope error");

    // SNES (Speedy Non-specific Enzyme Search): only takes effect when
    // peptide:enzyme_specificity is "none". For full/semi tryptic searches this flag
    // has no effect — the standard digestion + precursor-window lookup is always used.
    // v1 ships opt-in (default false); will flip to default true once external
    // workloads are validated end-to-end.
    defaults_.setValue("snes_enabled", "false",
      "[experimental, v1 opt-in] When peptide:enzyme_specificity=none, use mother-"
      "peptide indexing (Single-N + Single-C) instead of naïve O(L^2) sub-peptide "
      "enumeration. Orders-of-magnitude smaller index and faster search on "
      "non-specific workloads (immunopeptidomics). Ignored for specific/semi-"
      "specific enzymes. Variable modifications are applied via query-time subset enumeration (v1.1).");
    defaults_.setValidStrings("snes_enabled", {"true", "false"});
    
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

    // Derive SNES mode: snes_enabled switch AND enzyme_specificity == SPEC_NONE.
    // snes_enabled is a no-op for specific/semi-specific searches — its only purpose
    // is to provide a debug escape hatch for regression-testing the legacy
    // non-specific enumeration path. Keep snes_enabled_ stored raw so callers can
    // distinguish "SNES turned off" from "SNES not applicable for this enzyme".
    snes_enabled_ = param_.getValue("snes_enabled").toString() == "true";
    is_snes_mode_ = snes_enabled_ && (enzyme_specificity_ == EnzymaticDigestion::SPEC_NONE);

    // SNES v1 indexes b-ions for Single-N mothers and y-ions for Single-C mothers
    // regardless of the user's ions:add_b_ions / ions:add_y_ions toggles (the
    // candidate lookup in querySpectrumSNES_ hard-codes b/y precursor-equivalent
    // targets). Downstream scoring (ProSEAlgorithm) does honor those toggles
    // when building theoretical spectra, so turning b or y off leaves the SNES
    // filter admitting candidates that cannot be scored well — silent quality
    // degradation. Reject the configuration explicitly in v1.
    if (is_snes_mode_ && (!add_b_ions_ || !add_y_ions_))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "SNES mode (snes_enabled=true with enzyme_specificity=none) requires both "
        "ions:add_b_ions=true and ions:add_y_ions=true in v1. Additional ion "
        "series (a/c/x/z) may be enabled freely for downstream scoring.");
    }

    if (isOpenSearchMode_())
    {
      OPENMS_LOG_WARN << "[FragmentIndex] Open-search mode auto-triggered: window [-"
                      << precursor_mass_tolerance_lower_ << ", +"
                      << precursor_mass_tolerance_upper_ << "] "
                      << (precursor_mass_tolerance_unit_ppm_ ? "ppm" : "Da")
                      << " exceeds threshold. Isotope-error iteration collapses to [0, 0]."
                      << std::endl;
    }

    // Re-initialize modification tables to reflect the current
    // modifications_fixed_ / modifications_variable_ values.
    // Reset the guard so initModificationTables_() re-runs unconditionally.
    mod_tables_initialized_ = false;
    initModificationTables_();

    // SNES v1.1: precompute Σ_delta enumeration for the query path.
    // Three sets support anchor-dependent bin walks:
    //   baseline:            ANYWHERE + N_TERM + C_TERM variable mods
    //   with_prot_nterm:     baseline + PROTEIN_N_TERM variable mods
    //   with_prot_cterm:     baseline + PROTEIN_C_TERM variable mods
    // Non-SNES queries never consult these; populated unconditionally (cheap)
    // so that toggling snes_enabled at runtime does not require a rebuild.
    snes_sigma_delta_set_ = computeSnesSigmaDeltaSet_(false, false);
    snes_sigma_delta_set_with_prot_nterm_ = computeSnesSigmaDeltaSet_(true, false);
    snes_sigma_delta_set_with_prot_cterm_ = computeSnesSigmaDeltaSet_(false, true);

    const size_t largest_set = std::max({snes_sigma_delta_set_.size(),
                                          snes_sigma_delta_set_with_prot_nterm_.size(),
                                          snes_sigma_delta_set_with_prot_cterm_.size()});
    if (is_snes_mode_ && largest_set > 64)
    {
      OPENMS_LOG_WARN << "[FragmentIndex] SNES Σ_delta set has "
                      << largest_set << " entries — query performance will "
                      << "scale linearly with this. Consider reducing "
                      << "modifications:variable or variable_max_per_peptide.\n";
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
