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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <functional>
#include <map>
#include <mutex>
#include <set>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <boost/sort/sort.hpp>
#if defined(__linux__)
  #include <cstdint>
  #include <type_traits>
  #include <sys/mman.h> // madvise (releasePagesInParallel, adviseHugePages)
  #include <unistd.h>
#endif
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
  #include <intrin.h> // _mm_prefetch (queryPeaks)
#endif

using namespace std;



namespace OpenMS
{

  // Static member definitions
  std::array<double, 128> FragmentIndex::residue_mass_table_{};
  std::once_flag FragmentIndex::mass_table_once_flag_;
  FragmentIndex::IonOffsets FragmentIndex::ion_offsets_{};

  namespace
  {
    // Characters a peptide in the index may contain: one-letter codes of residues with a known
    // elemental formula. Excludes the ambiguous codes B, X and Z, stop codons ('*') and any other
    // symbol. AASequence parses '*' as a weightless X, so a peptide containing one could be
    // indexed but never scored.
    const std::array<bool, 256>& indexableResidues()
    {
      static const std::array<bool, 256> table = []
      {
        std::array<bool, 256> indexable{};
        const ResidueDB* rdb = ResidueDB::getInstance();
        for (char c = 'A'; c <= 'Z'; ++c)
        {
          const Residue* r = rdb->getResidue(static_cast<unsigned char>(c));
          indexable[static_cast<unsigned char>(c)] = (r != nullptr && !r->getFormula().isEmpty());
        }
        return indexable;
      }();
      return table;
    }
  }

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
      ion_offsets_.zp1_offset = Residue::getInternalToZp1Ion().getMonoWeight();
    });
  }

  void FragmentIndex::checkFixedModifications(const StringList& fixed_modifications)
  {
    if (fixed_modifications.empty()) return;
    // The index has one N- and one C-terminal fixed mass for all peptides (fixed_nterm_delta_ / fixed_cterm_delta_,
    // used for the precursor mass, the fragments and the reconstructed sequence alike).
    const ResidueModification* fixed_terminal[2] = {nullptr, nullptr}; // N-, C-terminus
    for (const auto& [mod_ptr, residue_ptr] : ModifiedPeptideGenerator::getModifications(fixed_modifications).val)
    {
      const ResidueModification::TermSpecificity term_spec = mod_ptr->getTermSpecificity();
      if (term_spec == ResidueModification::ANYWHERE) continue;
      const bool n_term = term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::PROTEIN_N_TERM;
      const std::string terminus = n_term ? "N-terminus" : "C-terminus";
      const std::string use_variable = ", but a fixed terminal modification is applied to every peptide. "
                                       "Specify it as a variable modification (modifications:variable) instead.";
      if (term_spec == ResidueModification::PROTEIN_N_TERM || term_spec == ResidueModification::PROTEIN_C_TERM)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Fixed modification '" + mod_ptr->getFullId() + "' applies to the protein " + terminus + " only" + use_variable);
      }
      const char origin = mod_ptr->getOrigin();
      if (origin != 'X' && origin != '.')
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Fixed modification '" + mod_ptr->getFullId() + "' applies only to peptides with " + std::string(1, origin)
          + " at the " + terminus + use_variable);
      }
      const ResidueModification*& previous = fixed_terminal[n_term ? 0 : 1];
      if (previous != nullptr)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Fixed modifications '" + previous->getFullId() + "' and '" + mod_ptr->getFullId() + "' both modify the peptide "
          + terminus + ", which carries one modification. Keep one of them as a fixed modification.");
      }
      previous = mod_ptr;
    }
  }

  namespace
  {
    /// 0 (N) or 1 (C) for a modification of a whole peptide or protein terminus (no residue preference), else -1
    int wholeTerminus(const ResidueModification& mod)
    {
      if (mod.getOrigin() != 'X' && mod.getOrigin() != '.') return -1;
      switch (mod.getTermSpecificity())
      {
        case ResidueModification::N_TERM:
        case ResidueModification::PROTEIN_N_TERM:
          return 0;
        case ResidueModification::C_TERM:
        case ResidueModification::PROTEIN_C_TERM:
          return 1;
        default:
          return -1;
      }
    }
  } // namespace

  StringList FragmentIndex::shadowedVariableTerminalModifications(const StringList& fixed_modifications,
                                                                  const StringList& variable_modifications)
  {
    StringList shadowed;
    if (fixed_modifications.empty() || variable_modifications.empty()) return shadowed;
    bool fixed_terminus[2] = {false, false}; // N-, C-terminus
    for (const auto& [mod_ptr, residue_ptr] : ModifiedPeptideGenerator::getModifications(fixed_modifications).val)
    {
      if (const int t = wholeTerminus(*mod_ptr); t >= 0) fixed_terminus[t] = true;
    }
    for (const std::string& name : variable_modifications) // in the given order (getModifications() returns a hash map)
    {
      for (const auto& [mod_ptr, residue_ptr] : ModifiedPeptideGenerator::getModifications({name}).val)
      {
        if (const int t = wholeTerminus(*mod_ptr); t >= 0 && fixed_terminus[t]) shadowed.push_back(name);
      }
    }
    return shadowed;
  }

  void FragmentIndex::initModificationTables_()
  {
    if (mod_tables_initialized_) return;
    checkFixedModifications(modifications_fixed_);

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
          // Residue-specific fixed mod (e.g., Carbamidomethyl on C): applies at all matching positions.
          // Residue-specific terminal ones are rejected by checkFixedModifications() above (they would need a
          // per-peptide terminal mass); the branches below are kept for completeness.
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
          // Pure terminal mod (no specific AA). A terminus carries one modification, so it is not applied where a
          // fixed terminal modification sits (as in ModifiedPeptideGenerator; shadowedVariableTerminalModifications()).
          // Applying both would give the index the sum of both masses, while the reconstructed sequence keeps one.
          if (term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::PROTEIN_N_TERM)
          {
            if (fixed_nterm_mod_ptr_ == nullptr) variable_nterm_mods_.push_back(entry);
          }
          else if (term_spec == ResidueModification::C_TERM || term_spec == ResidueModification::PROTEIN_C_TERM)
          {
            if (fixed_cterm_mod_ptr_ == nullptr) variable_cterm_mods_.push_back(entry);
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

  template <typename FragmentSink>
  void FragmentIndex::generateFragmentsLightweight_(
    FragmentSink& fragments,
    FragmentSink& electron_fragments,
    const char* sequence,
    size_t seq_len,
    UInt32 peptide_idx,
    double n_term_mod_mass,
    double c_term_mod_mass,
    const double* residue_masses,
    const double* residue_mod_masses) const
  {
    // Thin wrapper: forward class flags to the series-explicit implementation.
    generateFragmentsForSeries_(fragments, sequence, seq_len, peptide_idx,
                                n_term_mod_mass, c_term_mod_mass, residue_masses, residue_mod_masses,
                                add_b_ions_, add_a_ions_, add_c_ions_,
                                add_y_ions_, add_x_ions_, add_z_ions_, add_zp1_ions_);
    // ions:electron_ions: the c and z+1 ions that the series above lack, in a set of their own
    if (electron_ions_ && !(add_c_ions_ && add_zp1_ions_))
    {
      generateFragmentsForSeries_(electron_fragments, sequence, seq_len, peptide_idx,
                                  n_term_mod_mass, c_term_mod_mass, residue_masses, residue_mod_masses,
                                  /*add_b=*/false, /*add_a=*/false, /*add_c=*/!add_c_ions_,
                                  /*add_y=*/false, /*add_x=*/false, /*add_z=*/false,
                                  /*add_zp1=*/!add_zp1_ions_);
    }
  }

  template <typename FragmentSink>
  void FragmentIndex::generateFragmentsForSeries_(
    FragmentSink& fragments,
    const char* sequence,
    size_t seq_len,
    UInt32 peptide_idx,
    double n_term_mod_mass,
    double c_term_mod_mass,
    const double* residue_masses,
    const double* residue_mod_masses,
    bool add_b,
    bool add_a,
    bool add_c,
    bool add_y,
    bool add_x,
    bool add_z,
    bool add_zp1) const
  {
    const double proton = Constants::PROTON_MASS_U;
    const float min_mz = fragment_min_mz_;
    const float max_mz = fragment_max_mz_;
    const size_t min_ion_index = min_ion_index_;
    const IonOffsets offsets = ion_offsets_;
    const auto add = [&](double mass)
    {
      const float mz = static_cast<float>(mass);
      if (mz >= min_mz && mz <= max_mz) fragments.emplace_back(peptide_idx, mz);
    };

    // Prefix ions (b, a, c) sum up the residues from the left, suffix ions (y, x, z, z+1) from the right:
    // two independent sums, advanced together. Step i yields ion number i + 1 of both (i = 0: b1 / y1).
    // Fragment charge is always 1 for the index (matching original TSG call)
    double prefix = proton + n_term_mod_mass;
    double suffix = proton + c_term_mod_mass;
    for (size_t i = 0; i + 1 < seq_len; ++i)
    {
      const size_t j = seq_len - 1 - i;
      double prefix_residue = residue_masses[static_cast<unsigned char>(sequence[i])];
      double suffix_residue = residue_masses[static_cast<unsigned char>(sequence[j])];
      if (residue_mod_masses)
      {
        prefix_residue += residue_mod_masses[i];
        suffix_residue += residue_mod_masses[j];
      }
      prefix += prefix_residue;
      suffix += suffix_residue;

      if (i + 1 <= min_ion_index) continue; // skip ions below min_ion_index

      if (add_b) add(prefix + offsets.b_offset);
      if (add_a) add(prefix + offsets.a_offset);
      if (add_c) add(prefix + offsets.c_offset);
      if (add_y) add(suffix + offsets.y_offset);
      if (add_x) add(suffix + offsets.x_offset);
      if (add_z) add(suffix + offsets.z_offset);
      if (add_zp1) add(suffix + offsets.zp1_offset);
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

  namespace
  {
    // Gives the memory pages of a large buffer back to the operating system with the OpenMP
    // threads, so that the deallocation that follows has no pages left to unmap. Unmapping the
    // ~2 GB of a proteome-wide index is otherwise a serial step of 0.2 - 0.4 s per search.
    // The content of the buffer is lost: the caller deallocates it right away and does not read
    // it before, so no result depends on this. Failing calls leave the pages to the deallocation.
    template <typename T>
    void releasePagesInParallel(std::vector<T>& buffer)
    {
#if defined(__linux__) && defined(_OPENMP)
      static_assert(std::is_trivially_destructible<T>::value, "destructors would read the released pages");
      const std::uintptr_t page = static_cast<std::uintptr_t>(sysconf(_SC_PAGESIZE));
      const std::uintptr_t slice = (std::uintptr_t(32) << 20) / page * page;
      // whole pages inside the buffer only: its first and last page may hold other data
      const std::uintptr_t begin = (reinterpret_cast<std::uintptr_t>(buffer.data()) + page - 1) / page * page;
      const std::uintptr_t end = (reinterpret_cast<std::uintptr_t>(buffer.data()) + buffer.capacity() * sizeof(T)) / page * page;
      if (omp_get_max_threads() < 2 || omp_in_parallel() || end < begin + 2 * slice) { return; }
      const SignedSize slices = static_cast<SignedSize>((end - begin + slice - 1) / slice);
      #pragma omp parallel for schedule(dynamic)
      for (SignedSize i = 0; i < slices; ++i)
      {
        const std::uintptr_t from = begin + static_cast<std::uintptr_t>(i) * slice;
        madvise(reinterpret_cast<void*>(from), std::min(slice, end - from), MADV_DONTNEED);
      }
#else
      (void) buffer;
#endif
    }

    // Asks Linux to back a large buffer with transparent huge pages (with THP in "madvise" mode, only buffers advised so
    // get them): writing it then takes a 512th of the page faults, and releasing it a fraction of the time. Only the
    // 2 MB pages inside the buffer are advised, and the caller writes all of it right after, so this backs no memory the
    // buffer would not take anyway.
    template <typename T>
    void adviseHugePages(std::vector<T>& buffer)
    {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
      constexpr std::uintptr_t huge_page = std::uintptr_t(2) << 20;
      const std::uintptr_t begin = (reinterpret_cast<std::uintptr_t>(buffer.data()) + huge_page - 1) / huge_page * huge_page;
      const std::uintptr_t end = (reinterpret_cast<std::uintptr_t>(buffer.data()) + buffer.size() * sizeof(T)) / huge_page * huge_page;
      if (end > begin) { madvise(reinterpret_cast<void*>(begin), end - begin, MADV_HUGEPAGE); }
#else
      (void) buffer;
#endif
    }
  }

  void FragmentIndex::clear()
  {
    releasePagesInParallel(fi_fragments_);
    releasePagesInParallel(electron_fragments_);
    releasePagesInParallel(fi_peptides_);
    // swap-to-empty ensures the underlying heap capacity is actually released.
    // std::vector::clear() alone only resets size; capacity stays resident, which
    // defeats the M1 optimization of freeing the fragment index before downstream
    // PeptideIndexing / Aho-Corasick peaks (fi_fragments_ alone is hundreds of MB
    // on human-proteome builds). The swap idiom is the only portable way to force
    // deallocation across libstdc++/libc++/MSVC.
    std::vector<Fragment>().swap(fi_fragments_);
    std::vector<Fragment>().swap(electron_fragments_);
    std::vector<Peptide>().swap(fi_peptides_);
    std::vector<RemovedOccurrence>().swap(removed_occurrences_);
    std::vector<float>().swap(bucket_min_mz_);
    std::vector<float>().swap(electron_bucket_min_mz_);
    std::vector<UInt32>().swap(bucket_skip_);
    std::vector<UInt32>().swap(electron_bucket_skip_);
    std::vector<uint32_t>().swap(protein_lengths_);
    is_build_ = false;
    mod_tables_initialized_ = false;
  }

  bool FragmentIndex::isProteinNTerminal_(const std::string& protein, Size start) const
  {
    return start == 0 || (clip_nterm_methionine_ && start == 1 && ! protein.empty() && protein[0] == 'M');
  }

  AASequence FragmentIndex::reconstructModifiedSequence(
    const Peptide& peptide,
    const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const
  {
    const string& protein_seq = fasta_entries[peptide.protein_idx].sequence;
    AASequence seq = AASequence::fromString(StringUtils::substr(protein_seq, peptide.sequence_.first, peptide.sequence_.second));

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
      bool is_prot_nterm = isProteinNTerminal_(protein_seq, peptide.sequence_.first);
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

    AASequence seq = AASequence::fromString(StringUtils::substr(protein_seq, realized_start, realized_length));

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
      const bool is_prot_nterm = isProteinNTerminal_(protein_seq, realized_start);
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

#if defined(__GLIBCXX__) && defined(_OPENMP)
  namespace
  {
    // libstdc++'s std::sort (bits/stl_algo.h, std::__sort) is an introsort: __introsort_loop partitions a range
    // around the median of three elements, recurses into the right part, loops on the left part and stops at
    // ranges of at most _S_threshold = 16 elements (heap sort once 2 * floor(log2(n)) partitionings are nested);
    // one insertion sort over the whole range (__final_insertion_sort) finishes the job. The functions below make
    // the same comparisons and swaps, so they leave the same permutation - also among elements that compare
    // equal - but hand the right parts to OpenMP tasks:
    //  - a partitioning step reads and writes only its own range (the scans of __unguarded_partition are stopped
    //    by elements of the range), and the two parts are disjoint, so the order in which they are processed
    //    does not matter;
    //  - after the partitioning the ranges of <= 16 elements ("leaves") are in ascending order relative to each
    //    other, so the final insertion sort never moves an element across a leaf boundary: it is an insertion
    //    sort of each leaf, which can be done as soon as the leaf is known. (__insertion_sort, used for the first
    //    16 elements, and __unguarded_linear_insert differ in their comparisons, not in where an element ends up.)
    // Both points need what std::sort needs anyway: a comparison that is a strict weak ordering and has no state.
    constexpr std::ptrdiff_t STD_SORT_LEAF_SIZE = 16; // _S_threshold

    template <typename T, typename Less>
    void introsortLoopLikeStdSort(T* first, T* last, int depth_limit, Less less, std::ptrdiff_t min_task_size)
    {
      while (last - first > STD_SORT_LEAF_SIZE)
      {
        if (depth_limit == 0)
        {
          // std::__partial_sort(first, last, last, comp) == __make_heap + __sort_heap; the range is sorted afterwards
          std::make_heap(first, last, less);
          std::sort_heap(first, last, less);
          return;
        }
        --depth_limit;
        // std::__unguarded_partition_pivot: __move_median_to_first(first, first + 1, mid, last - 1) ...
        T* const a = first + 1;
        T* const b = first + (last - first) / 2;
        T* const c = last - 1;
        if (less(*a, *b))
        {
          if (less(*b, *c)) std::iter_swap(first, b);
          else if (less(*a, *c)) std::iter_swap(first, c);
          else std::iter_swap(first, a);
        }
        else if (less(*a, *c)) std::iter_swap(first, a);
        else if (less(*b, *c)) std::iter_swap(first, c);
        else std::iter_swap(first, b);
        // ... and __unguarded_partition(first + 1, last, pivot = first)
        T* cut = first + 1;
        T* hi = last;
        while (true)
        {
          while (less(*cut, *first)) ++cut;
          --hi;
          while (less(*first, *hi)) --hi;
          if (!(cut < hi)) break;
          std::iter_swap(cut, hi);
          ++cut;
        }
        // right part [cut, last): recursion in std::sort, a task here if it is large enough to be worth one
        if (last - cut > min_task_size)
        {
          #pragma omp task default(none) firstprivate(cut, last, depth_limit, less, min_task_size)
          introsortLoopLikeStdSort(cut, last, depth_limit, less, min_task_size);
        }
        else
        {
          introsortLoopLikeStdSort(cut, last, depth_limit, less, min_task_size);
        }
        last = cut; // left part [first, cut): next iteration
      }
      // leaf: the part of __final_insertion_sort that concerns [first, last)
      for (T* i = first + 1; i < last; ++i)
      {
        T value = std::move(*i);
        T* pos = i;
        for (; pos != first && less(value, *(pos - 1)); --pos) *pos = std::move(*(pos - 1));
        *pos = std::move(value);
      }
    }
  } // namespace
#endif

  void FragmentIndex::sortPeptides_(std::vector<Peptide>& peptides, size_t min_task_size)
  {
    const auto by_mz_then_protein = [](const Peptide& a, const Peptide& b)
    {
      return std::tie(a.precursor_mz_, a.protein_idx) < std::tie(b.precursor_mz_, b.protein_idx);
    };
#if defined(__GLIBCXX__) && defined(_OPENMP)
    // Relies on libstdc++'s std::sort algorithm (see introsortLoopLikeStdSort); FragmentIndex_test compares the two.
    // Any other standard library, and a small input, gets the plain std::sort call.
    if (peptides.size() > min_task_size)
    {
      Peptide* const first = peptides.data();
      Peptide* const last = first + peptides.size();
      const int depth_limit = 2 * (static_cast<int>(std::bit_width(peptides.size())) - 1); // std::__lg(n) * 2
      const std::ptrdiff_t min_size = static_cast<std::ptrdiff_t>(min_task_size);
      #pragma omp parallel
      #pragma omp single nowait
      introsortLoopLikeStdSort(first, last, depth_limit, by_mz_then_protein, min_size);
      return;
    }
#else
    (void)min_task_size;
#endif
    std::sort(peptides.begin(), peptides.end(), by_mz_then_protein);
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
    const std::array<bool, 256>& indexable = indexableResidues();

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

      // Position of the first residue at or after `from` that cannot be indexed
      // (X/B/Z, a stop codon or any other symbol), or npos.
      auto findUnindexable = [&seq, &indexable](size_t from)
      {
        for (size_t i = from; i < seq.size(); ++i)
        {
          if (!indexable[static_cast<unsigned char>(seq[i])]) return i;
        }
        return std::string::npos;
      };

      // Honor peptide:max_size=0 as "no maximum" (the documented semantics of
      // the non-SNES path). Using raw peptide_max_length_ in std::min would give
      // length 0 and an empty SNES index.
      const size_t effective_max_length = (peptide_max_length_ == 0) ? L : peptide_max_length_;

      // Mass-compute + filter + emit. No residue check here: the dispatch below
      // either calls sweepSpan(0, L) on a protein that can be indexed as a whole
      // or splits at X/B/Z, stop codons and other symbols, so span boundaries
      // structurally prevent any such residue from reaching this lambda.
      auto emitMother = [&](size_t start, size_t length, bool is_single_c)
      {
        if (length < peptide_min_length_) return;
        const char* seq_ptr = seq.c_str() + start;

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

      // Single-N mothers anchored at every position in [s, e - min_length], length
      // capped at effective_max_length and at the span end. Single-C mothers
      // anchored at every position j in [s + snes_min_length - 1, e - 1] with the
      // same length cap. snes_min_length guards the peptide_min_length_=0 corner
      // case (j would wrap to SIZE_MAX otherwise).
      const size_t snes_min_length = std::max<size_t>(1, peptide_min_length_);
      auto sweepSpan = [&](size_t s, size_t e)
      {
        if (e <= s || e - s < peptide_min_length_)
        {
          if (e > s) skipped_peptides.fetch_add(1);
          return;
        }
        for (size_t i = s; i + peptide_min_length_ <= e; ++i)
        {
          const size_t length = std::min<size_t>(effective_max_length, e - i);
          emitMother(i, length, /*is_single_c=*/false);
        }
        for (size_t j = s + snes_min_length - 1; j < e; ++j)
        {
          const size_t length = std::min<size_t>(effective_max_length, j + 1 - s);
          const size_t start = j + 1 - length;
          emitMother(start, length, /*is_single_c=*/true);
        }
      };

      // No X/B/Z (or stop codon, or other symbol) anywhere: sweep the whole
      // protein as a single span. Otherwise: split into contiguous unambiguous
      // spans and sweep each.
      // Issue #9192 item 2: previously the whole mother was dropped on any
      // X/B/Z overlap; truncating to the unambiguous prefix/suffix at the same
      // anchor preserves valid shorter realizations.
      const size_t first_bad = findUnindexable(0);
      if (first_bad == std::string::npos)
      {
        sweepSpan(0, L);
      }
      else
      {
        size_t p = 0;
        size_t bad = first_bad;
        while (true)
        {
          sweepSpan(p, bad);
          p = bad + 1;
          if (p >= L) break;  // protein ended with X/B/Z — no tail span
          bad = findUnindexable(p);
          if (bad == std::string::npos) { sweepSpan(p, L); break; }  // last span — no more X/B/Z
        }
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
    sortPeptides_(fi_peptides_);

    OPENMS_LOG_INFO << "Generated " << fi_peptides_.size() << " SNES mothers ("
                    << skipped_peptides.load() << " spans skipped — shorter than peptide:min_size)." << std::endl;
  }

  namespace
  {
    // Cleavage rule of the form "(?<=[KRX])" or "(?<=[KRX])(?!P)" (Trypsin/P, Trypsin, Lys-C, Arg-C, ...):
    // cleave after a residue of the first set unless the next residue is the blocking one.
    struct SimpleCleavageRule
    {
      std::array<bool, 256> after{};      ///< residues N-terminal of a cleavage site
      std::array<bool, 256> not_before{}; ///< residues C-terminal of a site that prevent the cleavage
    };

    // True (and @p rule filled) if the enzyme regex is exactly "(?<=[" + A-Z letters + "])", optionally followed
    // by "(?!" + one A-Z letter + ")". Everything else is left to the regex-based library digest.
    bool parseSimpleCleavageRule(const std::string& regex, SimpleCleavageRule& rule)
    {
      const auto is_letter = [](char c) { return c >= 'A' && c <= 'Z'; };
      if (regex.compare(0, 5, "(?<=[") != 0) return false;
      size_t i = 5;
      for (; i < regex.size() && is_letter(regex[i]); ++i) rule.after[static_cast<unsigned char>(regex[i])] = true;
      if (i == 5 || regex.compare(i, 2, "])") != 0) return false;
      i += 2;
      if (i == regex.size()) return true;
      if (regex.size() != i + 5 || regex.compare(i, 3, "(?!") != 0 || !is_letter(regex[i + 3]) || regex[i + 4] != ')') return false;
      rule.not_before[static_cast<unsigned char>(regex[i + 3])] = true;
      return true;
    }

    // Fully specific digest of a non-empty @p seq with a SimpleCleavageRule, without regex, substr copy or std::set.
    // Appends to @p out exactly the (start, length) spans, in the same order, that generatePeptides() otherwise gets
    // from EnzymaticDigestion::digestUnmodified() followed by its initial-Met-loss block:
    //  - tokenize_() yields 0 and every position 0 < p < size matched by the (zero-width) regex,
    //  - digestAfterTokenize_() emits the products with 0, then 1, 2, ... missed cleavages, each from N- to C-terminus,
    //  - the Met-loss block appends the N-terminal products of seq.substr(1) not yet present as (1, length).
    // @p sites is scratch space (kept by the caller to avoid an allocation per protein).
    void digestSimpleCleavage(const SimpleCleavageRule& rule, const std::string& seq, size_t missed_cleavages,
                              size_t min_length, size_t max_length, bool clip_nterm_methionine,
                              std::vector<size_t>& sites, std::vector<std::pair<size_t, size_t>>& out)
    {
      const size_t n = seq.size();
      sites.clear();
      sites.push_back(0);
      for (size_t p = 1; p < n; ++p)
      {
        if (rule.after[static_cast<unsigned char>(seq[p - 1])] && !rule.not_before[static_cast<unsigned char>(seq[p])]) sites.push_back(p);
      }
      const size_t count = sites.size();
      sites.push_back(n); // sentinel: the last product of every missed-cleavage level ends at the sequence end

      // digestUnmodified(): a maximum of 0 or beyond the sequence length means "no upper limit"
      const auto upper_limit = [max_length](size_t length) { return (max_length == 0 || max_length > length) ? length : max_length; };
      if (n >= min_length)
      {
        const size_t max_len = upper_limit(n);
        for (size_t mc = 0; mc <= missed_cleavages && mc < count; ++mc)
        {
          for (size_t j = 1; j + mc <= count; ++j)
          {
            const size_t l = sites[j + mc] - sites[j - 1];
            if (l >= min_length && l <= max_len) out.emplace_back(sites[j - 1], l);
          }
        }
      }

      // Initial-Met loss. Removing the first residue moves every cleavage site p >= 2 to p - 1 (the rule only looks
      // at the residues next to a site) and turns a site at 1 into the new start, so the N-terminal products of
      // seq.substr(1) end at the sites behind position 1 (or at the sequence end).
      if (clip_nterm_methionine && n > 1 && seq[0] == 'M' && n - 1 >= min_length)
      {
        const size_t n_spans = out.size();
        const size_t first = (sites[1] == 1) ? 2 : 1;     // index of the first site (or the end) behind position 1
        const size_t clipped_count = count - (first - 1); // number of sites of seq.substr(1), including its start
        const size_t max_len = upper_limit(n - 1);
        for (size_t mc = 0; mc <= missed_cleavages && mc < clipped_count; ++mc)
        {
          const size_t l = sites[first + mc] - 1;
          if (l < min_length || l > max_len) continue;
          // spans starting at 1 exist already only if 1 is a cleavage site (enzyme cleaving after M)
          const auto spans_end = out.begin() + n_spans;
          if (first == 2 && std::find(out.begin(), spans_end, std::make_pair(size_t(1), l)) != spans_end) continue;
          out.emplace_back(1, l);
        }
      }
    }

    // Variable-modification slots a non-SNES peptide may use: bit 31 of mod_bitmask_ marks Single-C mothers in SNES
    // mode, and the reconstruction masks it off (SNES_SLOT_MASK) in either mode.
    constexpr size_t MAX_ENUMERATED_SLOTS = 31;

    // The smallest x' >= x with at most max_set_bits bits set, or a value >= end if there is none below end
    // (end <= 2^32). Every number in [x, x + lowest set bit of x) keeps all bits of x and has at least as many set;
    // so while x has too many, the next candidate is x plus its lowest set bit.
    uint64_t nextSubsetWithin(uint64_t x, size_t max_set_bits, uint64_t end)
    {
      while (x < end && static_cast<size_t>(std::popcount(x)) > max_set_bits) x += x & (~x + 1);
      return x;
    }
  } // namespace

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
      size_t capped_peptides = 0; // peptides with more variable-modification slots than MAX_ENUMERATED_SLOTS

      ProteaseDigestion digestor;
      digestor.setEnzyme(digestion_enzyme_);
      digestor.setMissedCleavages(missed_cleavages_);
      digestor.setSpecificity(enzyme_specificity_);

      // Regex-free digest for fully specific searches with a simple cleavage rule (same spans, same order);
      // every other enzyme / specificity keeps using the library digest.
      SimpleCleavageRule cleavage_rule;
      const bool simple_digest = enzyme_specificity_ == EnzymaticDigestion::SPEC_FULL
                                 && digestor.getEnzymeName() != EnzymaticDigestion::UnspecificCleavage
                                 && parseSimpleCleavageRule(ProteaseDB::getInstance()->getEnzyme(digestion_enzyme_)->getRegEx(), cleavage_rule);

      OPENMS_LOG_INFO << "Generating peptides..." << std::endl;

      // Per-thread peptide vectors to avoid omp critical
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
#else
      const int num_threads = 1;
#endif
      // One cache line per thread: emplace_back() rewrites the vector's end pointer for every peptide, and adjacent
      // std::vector objects (24 bytes each) would make the threads bounce the same line back and forth.
      struct alignas(64) PaddedPeptides : vector<Peptide> {};
      vector<PaddedPeptides> thread_peptides(num_threads);
      // A tryptic digest with two missed cleavages yields about one peptide per three residues. Reserving one per two
      // spares the per-thread vectors their reallocation copies in the common case; untouched pages cost nothing.
      size_t total_residues = 0;
      for (const FASTAFile::FASTAEntry& entry : fasta_entries) total_residues += entry.sequence.size();
      const size_t est_per_thread = std::max(fasta_entries.size() * 5, total_residues / 2) / num_threads + 1;
      for (int t = 0; t < num_threads; ++t)
        thread_peptides[t].reserve(est_per_thread);

      const std::array<bool, 256>& indexable = indexableResidues();
      const auto is_unindexable = [&indexable](char c) { return !indexable[static_cast<unsigned char>(c)]; };

      vector<pair<size_t, size_t>> digested_peptides;
      vector<size_t> cleavage_sites;
      #pragma omp parallel for private(digested_peptides, cleavage_sites)
      for (SignedSize protein_idx = 0; protein_idx < (SignedSize)fasta_entries.size(); ++protein_idx)
      {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        digested_peptides.clear();
        const FASTAFile::FASTAEntry& protein = fasta_entries[protein_idx];
        if (simple_digest && !protein.sequence.empty())
        {
          digestSimpleCleavage(cleavage_rule, protein.sequence, missed_cleavages_, peptide_min_length_, peptide_max_length_,
                               clip_nterm_methionine_, cleavage_sites, digested_peptides);
        }
        else
        {
          digestor.digestUnmodified(protein.sequence, digested_peptides, peptide_min_length_, peptide_max_length_);
          if (clip_nterm_methionine_ && protein.sequence.size() > 1 && protein.sequence[0] == 'M'
              && enzyme_specificity_ != EnzymaticDigestion::SPEC_NONE)
          {
            // Digest the mature sequence separately so length and missed-cleavage limits
            // apply AFTER loss of the initial Met. Keep only its N-terminal spans:
            // internal peptides already exist in the ordinary digest.
            vector<pair<size_t, size_t>> clipped_peptides;
            digestor.digestUnmodified(protein.sequence.substr(1), clipped_peptides, peptide_min_length_, peptide_max_length_);
            std::set<size_t> existing_lengths;
            for (const auto& span : digested_peptides)
            {
              if (span.first == 1) { existing_lengths.insert(span.second); }
            }
            for (const auto& span : clipped_peptides)
            {
              if (span.first == 0 && existing_lengths.insert(span.second).second) { digested_peptides.emplace_back(1, span.second); }
            }
          }
        }

        // getProteinOccurrences() relies on this loop treating spans with equal residues alike: whether a span is
        // skipped and which entries it gets depend on its residues only (and, for protein-terminal modifications,
        // which hasProteinOccurrences() excludes, on its position in the protein).
        for (const pair<size_t, size_t>& digested_peptide : digested_peptides)
        {
          // skip peptides containing unknown or ambiguous AA codes (X, B, Z), stop codons ('*')
          // or other symbols
          {
            const std::string_view sub(protein.sequence.data() + digested_peptide.first, digested_peptide.second);
            if (std::any_of(sub.begin(), sub.end(), is_unindexable))
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
            bool is_prot_nterm = isProteinNTerminal_(protein.sequence, digested_peptide.first);
            bool is_prot_cterm = (digested_peptide.first + seq_len == protein.sequence.size());
            ModSlot slots[MAX_MOD_SLOTS];
            size_t n_slots = buildModSlots_(seq_ptr, seq_len, slots, is_prot_nterm, is_prot_cterm);
            if (n_slots > MAX_ENUMERATED_SLOTS)
            {
              // buildModSlots_() stops at MAX_MOD_SLOTS; the first 31 slots keep their bits
              n_slots = MAX_ENUMERATED_SLOTS;
              #pragma omp atomic
              capped_peptides++;
            }

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

              // Enumerate the slot subsets with at most max_variable_mods_per_peptide_ slots, in increasing bitmask
              // order (the emission order fixes the order of equal-mass variants in the index). nextSubsetWithin()
              // skips the subsets with more slots instead of visiting all 2^n_slots of them.
              const uint64_t end_bitmask = uint64_t{1} << n_slots;
              for (uint64_t subset = 0; subset < end_bitmask; subset = nextSubsetWithin(subset + 1, max_variable_mods_per_peptide_, end_bitmask))
              {
                const uint32_t bitmask = static_cast<uint32_t>(subset);

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
        OPENMS_LOG_WARN << skipped_peptides << " peptides skipped due to unknown or ambiguous AA (X/B/Z), stop codons or other symbols\n";
      }
      if (capped_peptides > 0)
      {
        OPENMS_LOG_WARN << capped_peptides << " peptide(s) have more than " << MAX_ENUMERATED_SLOTS
                        << " sites for variable modifications; only the first " << MAX_ENUMERATED_SLOTS
                        << " sites are considered for them (see modifications:variable)." << std::endl;
      }

      // Merge per-thread peptide vectors, in thread order: the sort below keys only on
      // (precursor_mz_, protein_idx) — which does NOT cover mod_bitmask_ / sequence_ — so
      // equal-key peptides are distinguishable and their relative order (hence every downstream
      // peptide index) depends on this concatenation. Peptide() writes nothing, so the resize only
      // allocates, and each thread copies (and first touches) the part that it generated.
      std::vector<size_t> merge_offsets(num_threads + 1, fi_peptides_.size());
      for (int t = 0; t < num_threads; ++t) merge_offsets[t + 1] = merge_offsets[t] + thread_peptides[t].size();
      fi_peptides_.resize(merge_offsets[num_threads]);
      #pragma omp parallel for schedule(static, 1)
      for (int t = 0; t < num_threads; ++t)
      {
        std::copy(thread_peptides[t].begin(), thread_peptides[t].end(), fi_peptides_.begin() + merge_offsets[t]);
        vector<Peptide>().swap(thread_peptides[t]);
      }

      OPENMS_LOG_INFO << "Sorting peptides..." << std::endl;
      sortPeptides_(fi_peptides_);
      OPENMS_LOG_INFO << "done." << std::endl;
  }

  FragmentIndex::PeptidoformRendering_ FragmentIndex::peptidoformRendering_() const
  {
    // AASequence::toString() writes the N-terminal modification, one token per residue (its one-letter code, or its
    // modification's toString()) and the C-terminal modification. Modifications are compared by that rendering:
    // different ResidueModification objects may render alike (Acetyl (N-term) and Acetyl (Protein N-term) are both
    // ".(Acetyl)"). Classes start at 256, above the residue bytes; 0 is "no terminal modification".
    PeptidoformRendering_ rendering;
    if (modifications_fixed_.empty() && modifications_variable_.empty()) return rendering;
    struct Member
    {
      const ResidueModification* mod;
      bool variable;
    };
    std::map<std::string, std::vector<Member>> members; // rendering -> modifications
    const auto add = [&](const ResidueModification* mod, const bool variable)
    {
      if (mod == nullptr) return;
      auto& list = members[mod->toString()];
      if (std::none_of(list.begin(), list.end(), [&](const Member& m) { return m.mod == mod && m.variable == variable; }))
      {
        list.push_back({mod, variable});
      }
    };
    for (const ResidueModification* mod : fixed_mod_ptrs_) add(mod, false);
    add(fixed_nterm_mod_ptr_, false);
    add(fixed_cterm_mod_ptr_, false);
    const auto add_variable = [&](const VarModEntry& entry)
    {
      add(entry.mod_ptr, true);
      rendering.slots_depend_on_context |= (entry.term_spec == ResidueModification::PROTEIN_N_TERM
                                            || entry.term_spec == ResidueModification::PROTEIN_C_TERM);
    };
    for (const auto& entries : variable_mod_table_)
    {
      for (const VarModEntry& entry : entries) add_variable(entry);
    }
    for (const VarModEntry& entry : variable_nterm_mods_) add_variable(entry);
    for (const VarModEntry& entry : variable_cterm_mods_) add_variable(entry);

    uint32_t cls = 256;
    for (const auto& [text, list] : members)
    {
      const auto& first = list.front();
      size_t variable_count = 0;
      for (const Member& m : list)
      {
        if (std::none_of(rendering.mod_class.begin(), rendering.mod_class.end(), [&](const auto& c) { return c.first == m.mod; }))
        {
          rendering.mod_class.emplace_back(m.mod, cls);
        }
        variable_count += m.variable ? 1 : 0;
        // Alike rendered modifications that differ in mass or residue, or a fixed one that a variable one renders
        // like (its mass is added on top), can give equal peptidoforms different precursor m/z or residues.
        rendering.runs_hold_peptidoforms &= (m.mod->getDiffMonoMass() == first.mod->getDiffMonoMass()
                                             && m.mod->getOrigin() == first.mod->getOrigin()
                                             && m.variable == first.variable);
      }
      rendering.unique_variable_renderings &= (variable_count <= 1);
      ++cls;
    }
    return rendering;
  }

  void FragmentIndex::renderPeptidoform_(const Peptide& peptide,
                                         const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                         const PeptidoformRendering_& rendering,
                                         std::vector<uint32_t>& tokens) const
  {
    // [N-term, residue 0, ..., residue len-1, C-term], with the modifications reconstructModifiedSequence() sets
    const std::string& protein = fasta_entries[peptide.protein_idx].sequence;
    const size_t start = peptide.sequence_.first;
    const size_t len = peptide.sequence_.second;
    const char* seq = protein.data() + start;
    tokens.assign(len + 2, 0);
    for (size_t i = 0; i < len; ++i) tokens[i + 1] = static_cast<unsigned char>(seq[i]);
    if (modifications_fixed_.empty() && modifications_variable_.empty()) return;
    tokens[0] = rendering.classOf(fixed_nterm_mod_ptr_);
    tokens[len + 1] = rendering.classOf(fixed_cterm_mod_ptr_);
    for (size_t i = 0; i < len; ++i)
    {
      const ResidueModification* fixed = fixed_mod_ptrs_[static_cast<unsigned char>(seq[i])];
      if (fixed != nullptr) tokens[i + 1] = rendering.classOf(fixed);
    }
    const uint32_t slot_bits = peptide.mod_bitmask_ & SNES_SLOT_MASK;
    if (slot_bits == 0) return;
    ModSlot slots[MAX_MOD_SLOTS];
    const size_t n_slots = buildModSlots_(seq, len, slots, isProteinNTerminal_(protein, start), start + len == protein.size());
    for (size_t s = 0; s < n_slots; ++s)
    {
      if (!(slot_bits & (1u << s))) continue;
      const uint32_t cls = rendering.classOf(slots[s].mod_ptr); // replaces a fixed modification at the same place
      if (slots[s].position == ModSlot::NTERM_SLOT) tokens[0] = cls;
      else if (slots[s].position == ModSlot::CTERM_SLOT) tokens[len + 1] = cls;
      else tokens[slots[s].position + 1] = cls;
    }
  }

  bool FragmentIndex::samePeptidoform_(const Peptide& a, const Peptide& b,
                                       const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                       const PeptidoformRendering_& rendering,
                                       std::vector<uint32_t>& tokens_a, std::vector<uint32_t>& tokens_b) const
  {
    // With rendering.runs_hold_peptidoforms, equal renderings imply equal residues (every rendered modification
    // names one residue).
    const size_t len = a.sequence_.second;
    if (b.sequence_.second != len) return false;
    const std::string& protein_a = fasta_entries[a.protein_idx].sequence;
    const std::string& protein_b = fasta_entries[b.protein_idx].sequence;
    const size_t start_a = a.sequence_.first;
    const size_t start_b = b.sequence_.first;
    if (std::memcmp(protein_a.data() + start_a, protein_b.data() + start_b, len) != 0) return false;
    // No variable modification: the fixed modifications depend on the residues only.
    if (a.mod_bitmask_ == 0 && b.mod_bitmask_ == 0) return true;
    // The same slots (buildModSlots_() sees the same residues and protein-terminal context): the same active slots
    // set the same modifications; different ones set differently rendered ones, unless two variable modifications
    // render alike.
    if (!rendering.slots_depend_on_context
        || (isProteinNTerminal_(protein_a, start_a) == isProteinNTerminal_(protein_b, start_b)
            && (start_a + len == protein_a.size()) == (start_b + len == protein_b.size())))
    {
      if (a.mod_bitmask_ == b.mod_bitmask_) return true;
      if (rendering.unique_variable_renderings) return false;
    }
    renderPeptidoform_(a, fasta_entries, rendering, tokens_a);
    renderPeptidoform_(b, fasta_entries, rendering, tokens_b);
    return tokens_a == tokens_b;
  }

  Size FragmentIndex::deduplicateByString_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    // The definition of peptide:deduplicate, literally (OpenMS #10394): group by a hash of the rendered sequence,
    // compare the strings within a group, keep the first entry of every string. build() uses it only for
    // configurations in which equal renderings need not have equal precursor m/z (see PeptidoformRendering_).
    std::vector<std::pair<size_t, Size>> fingerprints(fi_peptides_.size());
#pragma omp parallel for default(none) shared(fingerprints, fasta_entries)
    for (SignedSize i = 0; i < static_cast<SignedSize>(fi_peptides_.size()); ++i)
    {
      fingerprints[i] = {std::hash<std::string> {}(reconstructModifiedSequence(fi_peptides_[i], fasta_entries).toString()), static_cast<Size>(i)};
    }
    // Original index breaks hash ties so the first representative is stable.
    std::sort(fingerprints.begin(), fingerprints.end());
    std::vector<uint8_t> duplicate(fi_peptides_.size(), 0);
    // Removed entries as (index of their kept entry, index of the removed entry), both before compaction
    std::vector<std::pair<Size, Size>> removed_entries;
    for (Size begin = 0; begin < fingerprints.size();)
    {
      Size end = begin + 1;
      while (end < fingerprints.size() && fingerprints[end].first == fingerprints[begin].first)
      {
        ++end;
      }
      if (end - begin > 1)
      {
        std::unordered_map<std::string, Size> first_of; // rendered peptidoform -> its first (kept) entry
        for (Size i = begin; i < end; ++i)
        {
          const Size index = fingerprints[i].second;
          const auto [it, inserted] = first_of.emplace(reconstructModifiedSequence(fi_peptides_[index], fasta_entries).toString(), index);
          duplicate[index] = ! inserted;
          if (! inserted) removed_entries.emplace_back(it->second, index);
        }
      }
      begin = end;
    }
    // Ordered by kept entry, then by the position of the removed entry: getRemovedOccurrences()
    std::sort(removed_entries.begin(), removed_entries.end());
    removed_occurrences_.clear();
    removed_occurrences_.reserve(removed_entries.size());
    auto next_removed = removed_entries.begin();
    Size retained = 0;
    for (Size i = 0; i < fi_peptides_.size(); ++i)
    {
      if (duplicate[i]) continue;
      // removed entries come after their kept entry, which is therefore still in place here
      for (; next_removed != removed_entries.end() && next_removed->first == i; ++next_removed)
      {
        const Peptide& occurrence = fi_peptides_[next_removed->second];
        removed_occurrences_.push_back({static_cast<UInt32>(retained), occurrence.protein_idx, occurrence.sequence_.first});
      }
      fi_peptides_[retained++] = fi_peptides_[i];
    }
    const Size removed = fi_peptides_.size() - retained;
    fi_peptides_.erase(fi_peptides_.begin() + retained, fi_peptides_.end());
    return removed;
  }

  namespace
  {
    // build() generates the fragments partitioned into m/z bins and sortAndBucketFragments_() finishes the
    // bins independently of each other. A bin holds the fragments whose m/z bit patterns agree
    // in all but the lowest MZ_BIN_SHIFT bits (1 Th wide at m/z 512-1024, 2 Th at 1024-2048, ...). For floats
    // that are neither negative nor NaN the IEEE-754 bit pattern orders like the value, so the bins ascend
    // in m/z and the low bits ("cell") order the m/z values inside a bin.
    constexpr int MZ_BIN_SHIFT = 14;
    constexpr uint32_t MZ_BIN_END = 0x7F800000u >> MZ_BIN_SHIFT; // first bin of infinity, NaN and negative m/z
    constexpr size_t MAX_MZ_BINS = 8192;

    template <typename FragmentT>
    inline uint32_t mzBin(const FragmentT& fragment)
    {
      return std::bit_cast<uint32_t>(fragment.fragment_mz_) >> MZ_BIN_SHIFT;
    }

    // Bin of an m/z, counted from first_bin; whatever lies outside of [first_bin, first_bin + last] goes
    // to bin last.
    inline uint32_t mzBinFrom(float mz, uint32_t first_bin, uint32_t last)
    {
      return std::min((std::bit_cast<uint32_t>(mz) >> MZ_BIN_SHIFT) - first_bin, last);
    }

    inline void prefetchForRead(const void* address, std::ptrdiff_t byte_offset); // defined below, with queryPeaks()
  }

  void FragmentIndex::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    build(fasta_entries, {});
  }

  void FragmentIndex::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                            const std::function<const MSExperiment*(Size)>& searched_spectra)
  {
      // A rebuild replaces the previous database. generatePeptides() appends, so stale
      // peptides would otherwise be kept and their coordinates interpreted against the
      // new FASTA. Also leaves isBuild() false if this build throws.
      clear();
      protein_lengths_.reserve(fasta_entries.size());
      for (const auto& e : fasta_entries)
      {
        // Peptide coordinates (start offset, length) are stored as 16-bit values in
        // Peptide::sequence_, so a FASTA entry must not exceed 65535 residues. Beyond that,
        // the start-offset cast to uint16_t in generatePeptides()/generateSNESMothers_ would
        // wrap modulo 65536 and silently index fragments from the wrong subsequence. Fail loud
        // instead: long metaproteomic contigs / six-frame-translated frames must be split first.
        if (e.sequence.size() > std::numeric_limits<uint16_t>::max())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "FragmentIndex: FASTA entry '" + e.identifier + "' has " + std::to_string(e.sequence.size())
            + " residues, exceeding the supported maximum of 65535 (peptide offsets are stored as 16-bit). "
            "Split long contigs / six-frame-translated frames into windows of at most 65535 residues "
            "(with overlap >= peptide:max_size so no peptide is lost across a split) before building the index.");
        }
        protein_lengths_.push_back(static_cast<uint32_t>(e.sequence.size()));
      }

      /// generate all Peptides (also initializes residue mass table and mod tables)
      generatePeptides(fasta_entries);

      const bool has_modifications = !(modifications_fixed_.empty() && modifications_variable_.empty());

      OPENMS_LOG_INFO << "Generating fragments..." << std::endl;

#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
#else
      const int num_threads = 1;
#endif

      // The fragments are generated twice: the first pass counts them per m/z bin, the second writes each
      // one straight to its place in fi_fragments_ / electron_fragments_, which end up partitioned into
      // ascending m/z bins as sortAndBucketFragments_() wants them - no intermediate copy of the fragments.
      // Both passes work on the same portions of the peptides (fixed, independent of the number of
      // threads); in every bin a portion writes behind the portions before it, so that a bin is filled
      // in ascending peptide order.
      // The bins start at fragment_min_mz_, but not below 1 Th (no singly charged ion is lighter), and end
      // at fragment_max_mz_ or after MAX_MZ_BINS bins. A fragment outside of them would be put into the last
      // bin and leave the order to the general way of sortAndBucketFragments_(), as does the single bin
      // used if the limits allow for no fragment in any bin.
      const float lowest_mz = std::max(fragment_min_mz_, 1.0f);
      const bool binnable = (lowest_mz <= fragment_max_mz_);
      const uint32_t first_bin = binnable ? mzBinFrom(lowest_mz, 0, MZ_BIN_END) : 0;
      const uint32_t last = binnable ? mzBinFrom(fragment_max_mz_, first_bin, MAX_MZ_BINS - 1) : 0; // last bin, counted from first_bin
      const size_t num_bins = size_t(last) + 1;
      const size_t portion_size = std::max<size_t>(4 * num_bins, 1024); // peptides; the tables below take 2 bytes per peptide

      // peptide:deduplicate - protein occurrences are not distinct peptide hypotheses: of every peptidoform
      // (reconstructModifiedSequence(...).toString()) only the first entry is kept. The protein occurrences of the others
      // are recorded (getRemovedOccurrences()), so that getProteinOccurrences() still lists every protein occurrence,
      // including target/decoy shared sequences. SNES entries are mother peptides with different anchors, not scored
      // forms.
      // The entries of a peptidoform normally have equal residues and bitwise equal precursor_mz_ (generatePeptides()
      // adds the same masses in the same order), so they lie in one run of equal precursor_mz_. The first pass meets
      // them one after the other, their sequences in the cache, and keeps the first entry of every peptidoform
      // (samePeptidoform_()): its portions start at runs, and each moves its kept entries to its front, in place
      // (fresh memory would cost more page faults than the comparisons cost time). The portions are then closed up,
      // and the second pass works on what is left of each. Otherwise (a modification configured fixed and variable,
      // or alike rendered ones of different mass or residue) the rendered strings are compared beforehand.
      const bool deduplicate_requested = param_.getValue("peptide:deduplicate").toBool() && !is_snes_mode_;
      const PeptidoformRendering_ rendering = deduplicate_requested ? peptidoformRendering_() : PeptidoformRendering_{};
      const bool deduplicate_by_string = deduplicate_requested && !rendering.runs_hold_peptidoforms;
      if (deduplicate_by_string)
      {
        const Size removed = deduplicateByString_(fasta_entries);
        OPENMS_LOG_INFO << "Collapsed " << removed << " repeated peptidoform occurrences." << std::endl;
      }

      // Peptides that no spectrum can reach need no fragments (the open search windows reach nearly all peptides).
      // After the deduplication by strings: equal renderings may have different precursor m/z there, so the kept
      // entry of a peptidoform has to be chosen among all of its entries, not only among those in a window. The
      // deduplication in runs of equal precursor m/z below keeps or drops a run as a whole and may follow.
      if (searched_spectra && !is_snes_mode_ && !isOpenSearchMode_())
      {
        if (const MSExperiment* spectra = searched_spectra(fi_peptides_.size())) { keepPeptidesInPrecursorWindows_(*spectra); }
      }
      const bool deduplicate = deduplicate_requested && !deduplicate_by_string;
      const size_t num_peptides = fi_peptides_.size();
      const SignedSize num_portions = static_cast<SignedSize>((num_peptides + portion_size - 1) / portion_size);

      // Per portion and bin: the number of fragments (first pass), then the position of the next one (second pass)
      vector<size_t> positions(num_portions * num_bins, 0);
      vector<size_t> electron_positions(electron_ions_ ? num_portions * num_bins : 0, 0); // ions:electron_ions

      // What generateFragments...() write to in the first and in the second pass
      struct BinCounter
      {
        size_t* count;
        uint32_t first_bin, last;
        void emplace_back(UInt32 /*peptide_idx*/, float mz) { ++count[mzBinFrom(mz, first_bin, last)]; }
        void flush() {}
      };
      struct BinWriter
      {
        Fragment* destination;
        size_t* next;
        uint32_t first_bin, last;
        // The fragments are collected here and placed a buffer at a time, in the order of their arrival
        // (the scattered writes are faster in a loop of their own than mixed into the generation)
        enum : size_t { capacity = 1024 };
        size_t size = 0;
        Fragment buffer[capacity];
        void emplace_back(UInt32 peptide_idx, float mz)
        {
          buffer[size++] = Fragment(peptide_idx, mz);
          if (size == capacity) flush();
        }
        void flush()
        {
          for (size_t k = 0; k < size; ++k) destination[next[mzBinFrom(buffer[k].fragment_mz_, first_bin, last)]++] = buffer[k];
          size = 0;
        }
      };

      // Residue masses with the fixed modifications included, for peptides without variable modifications.
      // They yield the same fragment masses as the array of per-residue deltas used otherwise: with that,
      // the generator adds the delta to the residue mass (the sum stored here) or, for a residue without
      // fixed modification, +0.0, which leaves the mass as it is (no mass in the table is -0.0).
      std::array<double, 128> fixed_residue_masses = residue_mass_table_;
      for (size_t aa = 0; aa < fixed_residue_masses.size(); ++aa)
      {
        if (fixed_mod_deltas_[aa] != 0.0) fixed_residue_masses[aa] += fixed_mod_deltas_[aa];
      }

      // Unified fragment generation path for all cases.
      // For modified peptides: reconstruct per-residue deltas from bitmask + mod tables.
      // No AASequence construction, no ModifiedPeptideGenerator.
      // mod_masses: buffer for the per-residue mass deltas, kept by the caller (no allocation per peptide)
      const auto generate_fragments_of = [&](const size_t peptide_idx, vector<double>& mod_masses, auto& fragment_sink, auto& electron_sink)
      {
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

          const double n_term_mod = is_single_c ? 0.0 : fixed_nterm_delta_;
          const double c_term_mod = is_single_c ? fixed_cterm_delta_ : 0.0;

          // SNES candidate lookup in querySpectrumSNES_ only targets b-ions (for
          // Single-N mothers) and y-ions (for Single-C mothers) — the other ion
          // series (a/c/x/z/z+1) are not targeted and indexing them would be wasted
          // storage at best and source of silent data loss at worst if a user
          // disabled the primary series via ion toggles. Force b-only/y-only
          // regardless of the class add_*_ions_ flags. CodeRabbit #6.
          generateFragmentsForSeries_(
            fragment_sink, seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
            n_term_mod, c_term_mod,
            fixed_residue_masses.data(), nullptr,
            /*add_b=*/!is_single_c,
            /*add_a=*/false,
            /*add_c=*/false,
            /*add_y=*/ is_single_c,
            /*add_x=*/false,
            /*add_z=*/false,
            /*add_zp1=*/false);
        }
        else if (!has_modifications || pep.mod_bitmask_ == 0)
        {
          // No modifications or unmodified variant: just fixed mod deltas (if any)
          if (!has_modifications)
          {
            generateFragmentsLightweight_(
              fragment_sink, electron_sink, seq_ptr, seq_len,
              static_cast<UInt32>(peptide_idx), 0.0, 0.0, residue_mass_table_.data(), nullptr);
          }
          else
          {
            // Fixed mods only — their deltas are part of fixed_residue_masses
            generateFragmentsLightweight_(
              fragment_sink, electron_sink, seq_ptr, seq_len,
              static_cast<UInt32>(peptide_idx), fixed_nterm_delta_, fixed_cterm_delta_,
              fixed_residue_masses.data(), nullptr);
          }
        }
        else
        {
          // Variable modifications active: reconstruct delta array from bitmask
          mod_masses.assign(seq_len, 0.0);
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
          bool is_prot_nterm = isProteinNTerminal_(prot_seq, pep.sequence_.first);
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
            fragment_sink, electron_sink, seq_ptr, seq_len,
            static_cast<UInt32>(peptide_idx), n_term_mod, c_term_mod,
            residue_mass_table_.data(), has_residue_mods ? mod_masses.data() : nullptr);
        }
      };

      // Portion p takes the peptides [portion_start[p], portion_start[p + 1]). With deduplication, the first pass's
      // portions start at the first run of equal precursor_mz_ that begins at or after their nominal start.
      vector<size_t> portion_start(num_portions + 1, num_peptides);
      for (SignedSize portion = 0; portion < num_portions; ++portion)
      {
        size_t start = static_cast<size_t>(portion) * portion_size;
        while (deduplicate && start > 0 && start < num_peptides && fi_peptides_[start].precursor_mz_ == fi_peptides_[start - 1].precursor_mz_) ++start;
        portion_start[portion] = start;
      }
      vector<size_t> kept_count(deduplicate ? num_portions : 0); // first pass: entries kept of each portion
      // First pass: the entries each portion removed, with the position of their kept entry within the portion,
      // in the order met (getRemovedOccurrences())
      vector<vector<RemovedOccurrence>> removed_of_portion(deduplicate ? num_portions : 0);
      // Key of a peptide within its run: a hash of its length and of up to 8 residues after the first and before the
      // last one, and of its active modification slots where these name the peptidoform (the same residues give the
      // same slots, and no two variable modifications render alike). Equal peptidoforms agree in them, the peptides
      // of a run mostly differ in them (distinct ones of one composition, or the sites of a modification). 32 bits,
      // so that the comparison with all keys of a run vectorises; equal keys are only candidates.
      const bool key_has_slots = !rendering.slots_depend_on_context && rendering.unique_variable_renderings;
      const auto runKey = [&fasta_entries, key_has_slots](const Peptide& peptide)
      {
        uint64_t residues = 0;
        const size_t length = peptide.sequence_.second;
        const char* sequence = fasta_entries[peptide.protein_idx].sequence.data() + peptide.sequence_.first;
        if (length >= sizeof(residues) + 2) std::memcpy(&residues, sequence + 1, sizeof(residues));
        else if (length > 2) std::memcpy(&residues, sequence + 1, length - 2);
        const uint64_t slots = key_has_slots ? peptide.mod_bitmask_ : 0;
        const uint64_t mixed = (residues ^ (static_cast<uint64_t>(length) << 48)) * 0x9E3779B97F4A7C15ull + slots * 0xC2B2AE3D27D4EB4Full;
        return static_cast<uint32_t>(mixed >> 32);
      };

      // One pass over all peptides (first_pass: the counting one); make_sink(fragments, table row of the portion)
      // returns what a set of fragments is written to. The peptides are sorted by mass and heavier ones have more
      // fragments: the portions are handed out one by one so that all threads stay busy.
      const auto generate_fragments = [&](const auto make_sink, const bool first_pass)
      {
        #pragma omp parallel
        {
          vector<double> mod_masses;
          vector<uint32_t> run_keys;           // keys of the kept entries of the current run
          vector<Peptide> run_entries;         // the kept entries of the current run
          vector<UInt32> run_positions;        // their positions within the portion's kept entries
          vector<uint32_t> tokens_a, tokens_b; // samePeptidoform_()
          const size_t size = fi_peptides_.size();
          #pragma omp for schedule(dynamic)
          for (SignedSize portion = 0; portion < num_portions; ++portion)
          {
            auto fragment_sink = make_sink(fi_fragments_, positions.data() + portion * num_bins);
            auto electron_sink = make_sink(electron_fragments_, electron_positions.empty() ? nullptr : electron_positions.data() + portion * num_bins);
            const size_t begin = portion_start[portion];
            const size_t end = portion_start[portion + 1];
            const bool deduplicating = deduplicate && first_pass;
            // Deduplicating portions rewrite their entries: look ahead only within the own portion then.
            const size_t ahead_end = deduplicating ? end : size;
            size_t kept = 0;
            for (size_t peptide_idx = begin; peptide_idx < end; ++peptide_idx)
            {
              // Sorted by mass, the peptides come from the proteins in no order: fetch the sequences of the
              // next ones (first the string object, then its characters) while this one is worked on
              if (peptide_idx + 16 < ahead_end)
              {
                prefetchForRead(&fasta_entries[fi_peptides_[peptide_idx + 16].protein_idx].sequence, 0);
                const Peptide& ahead = fi_peptides_[peptide_idx + 8];
                prefetchForRead(fasta_entries[ahead.protein_idx].sequence.data(), ahead.sequence_.first);
              }
              if (!deduplicating)
              {
                generate_fragments_of(peptide_idx, mod_masses, fragment_sink, electron_sink);
                continue;
              }
              const Peptide peptide = fi_peptides_[peptide_idx];
              if (peptide_idx == begin || peptide.precursor_mz_ != run_entries.front().precursor_mz_)
              {
                run_keys.clear();
                run_entries.clear();
                run_positions.clear();
              }
              const uint32_t key = runKey(peptide);
              // Long runs are mostly distinct peptides of one composition: test all keys at once, without branches,
              // and compare entries only where a key matches.
              bool candidate = false;
              for (const uint32_t earlier : run_keys) candidate |= (earlier == key);
              size_t repeated = run_keys.size(); // the kept entry this one repeats, if any
              for (size_t k = 0; candidate && k < run_keys.size(); ++k)
              {
                if (run_keys[k] == key && samePeptidoform_(run_entries[k], peptide, fasta_entries, rendering, tokens_a, tokens_b))
                {
                  repeated = k;
                  break;
                }
              }
              if (repeated < run_keys.size())
              {
                removed_of_portion[portion].push_back({run_positions[repeated], peptide.protein_idx, peptide.sequence_.first});
                continue;
              }
              run_keys.push_back(key);
              run_entries.push_back(peptide);
              run_positions.push_back(static_cast<UInt32>(kept));
              fi_peptides_[begin + kept] = peptide; // at or before peptide_idx
              generate_fragments_of(begin + kept, mod_masses, fragment_sink, electron_sink); // the counter ignores the index
              ++kept;
            }
            if (deduplicating) kept_count[portion] = kept;
            fragment_sink.flush();
            electron_sink.flush();
          }
        }
      };

      // First pass: count
      generate_fragments([&](vector<Fragment>& /*fragments*/, size_t* count) { return BinCounter{count, first_bin, last}; }, true);

      if (deduplicate)
      {
        // Close up the portions: the kept entries of each are numbered on from those of the portions before. Each
        // portion moves to the left, onto kept entries of earlier portions only, after these have moved: in waves, a
        // portion one wave after the latest of the portions whose entries it covers.
        const vector<size_t> source = portion_start;
        for (SignedSize portion = 0; portion < num_portions; ++portion)
        {
          portion_start[portion + 1] = portion_start[portion] + kept_count[portion];
        }
        const size_t kept_total = portion_start[num_portions];
        vector<int> wave(num_portions, -1); // -1: stays
        int num_waves = 0;
        for (SignedSize portion = 0; portion < num_portions; ++portion)
        {
          const size_t to = portion_start[portion];
          if (kept_count[portion] == 0 || to == source[portion]) continue;
          int w = 0;
          for (SignedSize q = portion - 1; q >= 0 && source[q] + kept_count[q] > to; --q)
          {
            if (wave[q] >= 0 && source[q] < to + kept_count[portion]) w = std::max(w, wave[q] + 1);
          }
          wave[portion] = w;
          num_waves = std::max(num_waves, w + 1);
        }
        for (int w = 0; w < num_waves; ++w)
        {
          #pragma omp parallel for schedule(dynamic)
          for (SignedSize portion = 0; portion < num_portions; ++portion)
          {
            if (wave[portion] != w) continue;
            const auto from = fi_peptides_.begin() + source[portion];
            std::copy(from, from + kept_count[portion], fi_peptides_.begin() + portion_start[portion]);
          }
        }
        OPENMS_LOG_INFO << "Collapsed " << (num_peptides - kept_total) << " repeated peptidoform occurrences." << std::endl;
        fi_peptides_.erase(fi_peptides_.begin() + kept_total, fi_peptides_.end());

        // The removed occurrences under the final index of their kept entry; stable, so that the occurrences of an
        // entry stay in the order met (a run, and so every kept entry with its repeats, lies in one portion)
        removed_occurrences_.clear();
        removed_occurrences_.reserve(num_peptides - kept_total);
        for (SignedSize portion = 0; portion < num_portions; ++portion)
        {
          for (RemovedOccurrence occurrence : removed_of_portion[portion])
          {
            occurrence.peptide_idx += static_cast<UInt32>(portion_start[portion]);
            removed_occurrences_.push_back(occurrence);
          }
          vector<RemovedOccurrence>().swap(removed_of_portion[portion]);
        }
        std::stable_sort(removed_occurrences_.begin(), removed_occurrences_.end(),
                         [](const RemovedOccurrence& a, const RemovedOccurrence& b) { return a.peptide_idx < b.peptide_idx; });
      }

      // Turn the counts into the position at which each portion fills each bin; returns the number of fragments
      const auto counts_to_positions = [num_bins](vector<size_t>& table)
      {
        vector<size_t> bin_position(num_bins, 0);
        for (size_t row = 0; row < table.size(); row += num_bins)
        {
          for (size_t bin = 0; bin < num_bins; ++bin) bin_position[bin] += table[row + bin];
        }
        size_t total = 0;
        for (size_t& position : bin_position)
        {
          const size_t count = position;
          position = total;
          total += count;
        }
        for (size_t row = 0; row < table.size(); row += num_bins)
        {
          for (size_t bin = 0; bin < num_bins; ++bin)
          {
            const size_t count = table[row + bin];
            table[row + bin] = bin_position[bin];
            bin_position[bin] += count;
          }
        }
        return total;
      };
      // Fragment's default constructor leaves the new elements uninitialised: the second pass writes them first
      fi_fragments_.resize(counts_to_positions(positions));
      electron_fragments_.resize(counts_to_positions(electron_positions));
      adviseHugePages(fi_fragments_);
      adviseHugePages(electron_fragments_);

      // Second pass: write
      generate_fragments([&](vector<Fragment>& fragments, size_t* next) { return BinWriter{fragments.data(), next, first_bin, last}; }, false);

      OPENMS_LOG_INFO << "Sorting fragments..." << std::endl;

      // Empty database (no peptide passed length / mass / motif filters): nothing to bucket.
      // Mark as built and return — guards against the OMP loop below dividing by zero
      // when bucketsize_ becomes 0. This is a real risk for immunopeptidomics FASTAs that
      // contain entries shorter than peptide:min_size.
      if (fi_fragments_.empty() && electron_fragments_.empty())
      {
        bucketsize_ = 1; // keep non-zero to preserve bucket-walking loop invariants
        OPENMS_LOG_INFO << "[FragmentIndex] No fragments generated — index is empty." << std::endl;
        is_build_ = true;
        return;
      }

      /// Bucket size chosen to approximate MSFragger's fixed ~0.02 Da fragment-bin density:
      /// in the dense 500-1500 Da region a typical tryptic/immunopeptidomics index holds
      /// ~4-8k fragments per 0.02 Da window, so a bucket covers roughly one query tolerance
      /// window instead of the much wider sqrt(N) span.
      bucketsize_ = 4096;
      OPENMS_LOG_INFO << "Creating DB with bucket_size " << bucketsize_ << endl;

      sortAndBucketFragments_(fi_fragments_, bucket_min_mz_, num_threads);
      // ions:electron_ions: the c and z+1 ions get buckets of their own, walked only on request
      sortAndBucketFragments_(electron_fragments_, electron_bucket_min_mz_, num_threads);
      buildSkipTables_();

      is_build_ = true;
      OPENMS_LOG_INFO << "Fragment index built!" << endl;
  }

  namespace
  {
    // Does the work of sortAndBucketFragments_() for fragments that arrive partitioned into ascending m/z
    // bins (see build()): identical result, linear time, each bin handled in the cache.
    // Returns false if the fragments are not partitioned like that; it may have reordered them by then.
    //
    // The result is fixed by the fragments alone: bucket k holds those of rank [k * bucketsize,
    // (k + 1) * bucketsize) in the order (m/z, peptide index), sorted by (peptide index, m/z), and
    // bucket_min_mz[k] is the m/z of rank k * bucketsize. Fragments that agree in both are bitwise
    // identical, so it does not matter how these are computed. Here, a bin is first brought into peptide
    // order; walking it in that order, the number of fragments with a smaller m/z plus the number of those
    // with this m/z seen so far is the rank of a fragment, which names its bucket, and the fragments of a
    // bucket are met in peptide order.
    template <typename FragmentT>
    bool sortAndBucketBinnedFragments(std::vector<FragmentT>& fragments, std::vector<float>& bucket_min_mz, const size_t bucketsize)
    {
      const size_t size = fragments.size();
      FragmentT* const data = fragments.data();
      const auto by_mz_then_peptide = [](const FragmentT& a, const FragmentT& b)
      {
        return std::tie(a.fragment_mz_, a.peptide_idx_) < std::tie(b.fragment_mz_, b.peptide_idx_);
      };
      const auto by_peptide_then_mz = [](const FragmentT& a, const FragmentT& b)
      {
        return std::tie(a.peptide_idx_, a.fragment_mz_) < std::tie(b.peptide_idx_, b.fragment_mz_);
      };
      const int bucket_shift = std::has_single_bit(bucketsize) ? std::countr_zero(bucketsize) : -1;
      const auto bucket_of = [bucketsize, bucket_shift](size_t rank) { return bucket_shift >= 0 ? rank >> bucket_shift : rank / bucketsize; };

      // Locate the bins. Whether every fragment lies in its bin is checked when the bin is processed.
      struct Bin
      {
        uint32_t id;
        size_t begin;
      };
      std::vector<Bin> bins;
      for (size_t i = 0; i < size;)
      {
        const uint32_t id = mzBin(data[i]);
        if (id >= MZ_BIN_END || bins.size() == MAX_MZ_BINS || (!bins.empty() && id <= bins.back().id)) return false;
        bins.push_back({id, i});
        i = std::partition_point(data + i + 1, data + size, [id](const FragmentT& f) { return mzBin(f) <= id; }) - data;
      }
      const SignedSize num_bins = static_cast<SignedSize>(bins.size());
      bins.push_back({MZ_BIN_END, size});

      bucket_min_mz.resize((size + bucketsize - 1) / bucketsize);

      constexpr size_t num_cells = size_t(1) << MZ_BIN_SHIFT; // distinct m/z values of a bin
      constexpr int radix_bits = 12;
      constexpr size_t radix_size = size_t(1) << radix_bits;
      constexpr int radix_passes = 3; // 12 + 12 + 8 bits of the peptide index
      bool binned = true;
      #pragma omp parallel
      {
        std::vector<FragmentT> scratch;
        std::vector<uint32_t> cell_rank(num_cells);
        std::vector<uint32_t> digit_count(radix_passes * radix_size);
        std::vector<size_t> bucket_next;

        #pragma omp for schedule(dynamic)
        for (SignedSize b = 0; b < num_bins; ++b)
        {
          const uint32_t id = bins[b].id;
          const size_t begin = bins[b].begin;
          const size_t n = bins[b + 1].begin - begin;
          FragmentT* const bin = data + begin;
          // The bin starts in bucket first_bucket, lead fragments after the start of that bucket. Ranks are
          // counted from the start of that bucket.
          const size_t first_bucket = bucket_of(begin);
          const size_t lead = begin - first_bucket * bucketsize;

          // (a bin too large for 32-bit ranks is treated like fragments outside their bin: the general way)
          bool in_bin = (n < std::numeric_limits<uint32_t>::max() - bucketsize);
          if (in_bin && n < num_cells / 4)
          {
            // Few fragments: not worth the tables below, sort by comparison
            for (size_t k = 0; k < n; ++k) in_bin &= (mzBin(bin[k]) == id);
            if (in_bin)
            {
              std::sort(bin, bin + n, by_mz_then_peptide);
              for (size_t piece = 0; piece < n;)
              {
                const size_t piece_end = std::min(n, (bucket_of(lead + piece) + 1) * bucketsize - lead);
                if (bucket_of(lead + piece) * bucketsize == lead + piece) bucket_min_mz[bucket_of(lead + piece) + first_bucket] = bin[piece].fragment_mz_;
                std::sort(bin + piece, bin + piece_end, by_peptide_then_mz);
                piece = piece_end;
              }
              continue;
            }
          }

          // Fragments per m/z value; is the bin in peptide order already?
          bool peptide_order = true;
          if (in_bin)
          {
            std::fill(cell_rank.begin(), cell_rank.end(), 0u);
            uint32_t previous = 0;
            for (size_t k = 0; k < n; ++k)
            {
              const uint32_t bits = std::bit_cast<uint32_t>(bin[k].fragment_mz_);
              in_bin &= ((bits >> MZ_BIN_SHIFT) == id);
              ++cell_rank[bits & (num_cells - 1)];
              peptide_order &= (bin[k].peptide_idx_ >= previous);
              previous = bin[k].peptide_idx_;
            }
          }
          if (!in_bin)
          {
            #pragma omp critical (FragmentIndex_sortAndBucketBinnedFragments)
            binned = false;
            continue;
          }

          if (scratch.size() < n) scratch.resize(n);
          FragmentT* from = bin;
          FragmentT* to = scratch.data();
          if (!peptide_order)
          {
            // LSD radix sort by peptide index; a digit in which all indices agree needs no pass
            std::fill(digit_count.begin(), digit_count.end(), 0u);
            for (size_t k = 0; k < n; ++k)
            {
              const uint32_t peptide = from[k].peptide_idx_;
              for (int pass = 0; pass < radix_passes; ++pass) ++digit_count[pass * radix_size + ((peptide >> (pass * radix_bits)) & (radix_size - 1))];
            }
            for (int pass = 0; pass < radix_passes; ++pass)
            {
              uint32_t* next = &digit_count[pass * radix_size];
              const int shift = pass * radix_bits;
              if (next[(from[0].peptide_idx_ >> shift) & (radix_size - 1)] == n) continue;
              uint32_t sum = 0;
              for (size_t digit = 0; digit < radix_size; ++digit)
              {
                const uint32_t count = next[digit];
                next[digit] = sum;
                sum += count;
              }
              for (size_t k = 0; k < n; ++k) to[next[(from[k].peptide_idx_ >> shift) & (radix_size - 1)]++] = from[k];
              std::swap(from, to);
            }
          }

          // Turn the counts into the rank of the first fragment of each m/z value. The m/z value that covers
          // the first rank of a bucket is the smallest of that bucket.
          uint32_t rank = static_cast<uint32_t>(lead);
          size_t bucket = (lead == 0) ? 0 : 1; // next bucket to start, relative to first_bucket
          for (size_t cell = 0; cell < num_cells; ++cell)
          {
            const uint32_t count = cell_rank[cell];
            cell_rank[cell] = rank;
            rank += count;
            for (; bucket * bucketsize < rank; ++bucket)
            {
              bucket_min_mz[first_bucket + bucket] = std::bit_cast<float>(static_cast<uint32_t>((id << MZ_BIN_SHIFT) | cell));
            }
          }

          // Place the fragments, in peptide order, into the next free slot of their bucket
          bucket_next.resize(bucket_of(lead + n - 1) + 1);
          bucket_next[0] = 0;
          for (size_t k = 1; k < bucket_next.size(); ++k) bucket_next[k] = k * bucketsize - lead;
          for (size_t k = 0; k < n; ++k)
          {
            const uint32_t fragment_rank = cell_rank[std::bit_cast<uint32_t>(from[k].fragment_mz_) & (num_cells - 1)]++;
            to[bucket_next[bucket_of(fragment_rank)]++] = from[k];
          }
          if (to != bin) std::copy(to, to + n, bin);

          // Several fragments of one peptide in a bucket are next to each other now: order them by m/z
          for (size_t k = 1; k < n; ++k)
          {
            for (size_t i = k; i > 0 && bin[i - 1].peptide_idx_ == bin[i].peptide_idx_ && bin[i].fragment_mz_ < bin[i - 1].fragment_mz_; --i)
            {
              std::swap(bin[i - 1], bin[i]);
            }
          }
        }
      }
      if (!binned) return false;

      // A bucket that extends over several bins consists of one sorted piece per bin: merge them
      #pragma omp parallel for schedule(dynamic)
      for (SignedSize b = 1; b < num_bins; ++b)
      {
        const size_t bucket_begin = bucket_of(bins[b].begin) * bucketsize;
        if (bins[b].begin == bucket_begin || bins[b - 1].begin > bucket_begin) continue; // not cut, or done with an earlier bin
        const size_t bucket_end = std::min(bucket_begin + bucketsize, size);
        for (SignedSize piece = b; piece < num_bins && bins[piece].begin < bucket_end; ++piece)
        {
          std::inplace_merge(data + bucket_begin, data + bins[piece].begin, data + std::min(bins[piece + 1].begin, bucket_end), by_peptide_then_mz);
        }
      }
      return true;
    }
  }

  void FragmentIndex::sortAndBucketFragments_(std::vector<Fragment>& fragments,
                                              std::vector<float>& bucket_min_mz,
                                              int num_threads)
  {
      if (fragments.empty()) return;

      // build() delivers the fragments partitioned into m/z bins, which are finished without a global sort.
      // The general way below remains for fragments in any other order and defines the result.
      if (sortAndBucketBinnedFragments(fragments, bucket_min_mz, bucketsize_)) return;

      /// 1.) First all Fragments are sorted by their own mass (parallel via Boost.Sort).
      /// Boost defaults to std::thread::hardware_concurrency() threads, which ignores both
      /// the tool's -threads option and the process CPU affinity: on a 384-core node a
      /// single-threaded run spawned 384 sort threads. Use the same budget as the OpenMP
      /// regions around it.
      boost::sort::block_indirect_sort(fragments.begin(), fragments.end(), [](const Fragment& a, const Fragment& b)
      {
        return std::tie(a.fragment_mz_, a.peptide_idx_) < std::tie(b.fragment_mz_, b.peptide_idx_);
      }, static_cast<uint32_t>(num_threads));

      /// 2.) Within each fragment-m/z bucket, re-sort the fragments by their originating peptide
      /// index so that query() can binary-search a candidate peptide range inside a bucket.
      ///
      /// bucket_min_mz[k] is the smallest fragment m/z in bucket k. Because fragments is
      /// already globally sorted by fragment_mz_, these per-bucket minima are monotonically
      /// non-decreasing, so we write them directly by bucket index — no omp critical and no
      /// trailing re-sort of bucket_min_mz is required.
      const size_t num_buckets = (fragments.size() + bucketsize_ - 1) / bucketsize_;
      bucket_min_mz.resize(num_buckets);

      // Per-thread scratch buffers reused as the LSD-radix ping-pong destination across buckets.
      vector<vector<Fragment>> radix_scratch(num_threads);
      for (auto& s : radix_scratch) s.reserve(bucketsize_);

      #pragma omp parallel for
      for (SignedSize b = 0; b < (SignedSize)num_buckets; ++b)
      {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        const size_t i = static_cast<size_t>(b) * bucketsize_;
        bucket_min_mz[b] = fragments[i].fragment_mz_;

        Fragment* base = fragments.data() + i;
        const size_t n = std::min<size_t>(bucketsize_, fragments.size() - i);

        // LSD radix sort of the bucket by peptide_idx_ (uint32). peptide_idx_ spans the full
        // peptide range, so a value-range counting sort is not applicable; instead we do a few
        // 8-bit passes — only as many bytes as the largest index in the bucket needs (3 for a
        // ~2M-peptide database). Stable, branch-free, and ~3-4x faster than std::sort on these
        // dense 4096-element buckets. (Replaces the per-bucket std::sort.)
        vector<Fragment>& scratch = radix_scratch[tid];
        scratch.resize(n);

        uint32_t max_idx = 0;
        for (size_t k = 0; k < n; ++k) max_idx = std::max(max_idx, base[k].peptide_idx_);
        int num_passes = 1;
        for (uint32_t m = max_idx; m >>= 8; ) ++num_passes;

        Fragment* src = base;
        Fragment* dst = scratch.data();
        for (int p = 0; p < num_passes; ++p)
        {
          const int shift = p * 8;
          uint32_t count[256] = {0};
          for (size_t k = 0; k < n; ++k) ++count[(src[k].peptide_idx_ >> shift) & 0xFFu];
          uint32_t sum = 0;
          for (int c = 0; c < 256; ++c) { uint32_t t = count[c]; count[c] = sum; sum += t; }
          for (size_t k = 0; k < n; ++k)
          {
            const uint32_t radix = (src[k].peptide_idx_ >> shift) & 0xFFu;
            dst[count[radix]++] = src[k];
          }
          std::swap(src, dst);
        }
        // After an odd number of passes the sorted data lives in scratch — copy it back in place.
        if (src != base) std::copy(src, src + n, base);
      }
  }

  void FragmentIndex::buildSkipTables_()
  {
    // Sample the peptide_idx_ of every SKIP_STRIDE_-th fragment of each (peptide-sorted) bucket
    // and of its last fragment: sample k is the fragment at min(k * SKIP_STRIDE_, size - 1).
    // queryPeaks() finds the start of the candidate ranges in these few cache-resident entries
    // instead of binary-searching the 32 kB bucket itself. Costs 4 bytes per SKIP_STRIDE_
    // fragments (< 1% of the index).
    skip_per_bucket_ = (bucketsize_ + SKIP_STRIDE_ - 1) / SKIP_STRIDE_ + 1;
    auto sample = [this](const std::vector<Fragment>& fragments, std::vector<UInt32>& skip)
    {
      const size_t num_buckets = (fragments.size() + bucketsize_ - 1) / bucketsize_;
      skip.resize(num_buckets * skip_per_bucket_);
      #pragma omp parallel for
      for (SignedSize b = 0; b < (SignedSize)num_buckets; ++b)
      {
        const size_t begin = static_cast<size_t>(b) * bucketsize_;
        const size_t last = std::min<size_t>(bucketsize_, fragments.size() - begin) - 1;
        UInt32* samples = skip.data() + static_cast<size_t>(b) * skip_per_bucket_;
        for (size_t k = 0; k < skip_per_bucket_; ++k)
        {
          samples[k] = fragments[begin + std::min(k * SKIP_STRIDE_, last)].peptide_idx_;
        }
      }
    };
    sample(fi_fragments_, bucket_skip_);
    sample(electron_fragments_, electron_bucket_skip_);
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

  void FragmentIndex::keepPeptidesInPrecursorWindows_(const MSExperiment& spectra)
  {
    // The windows of querySpectrum(), with its charges and isotope errors and the float arithmetic of
    // searchDifferentPrecursorRanges() and getPeptidesInMassWindow(). The margin only covers a different
    // rounding of the same expressions (e.g. fused multiply-adds in one place only); it is far smaller
    // than a window.
    std::vector<std::pair<float, float>> windows;
    std::vector<uint16_t> charges;
    for (const MSSpectrum& spectrum : spectra)
    {
      if (spectrum.empty() || spectrum.getMSLevel() != 2 || spectrum.getPrecursors().size() != 1) { continue; } // not searched
      const Precursor& precursor = spectrum.getPrecursors()[0];
      charges.clear();
      if (precursor.getCharge()) { charges.push_back(static_cast<uint16_t>(precursor.getCharge())); }
      else
      {
        for (uint16_t charge = min_precursor_charge_; charge <= max_precursor_charge_; ++charge) { charges.push_back(charge); }
      }
      for (const uint16_t charge : charges)
      {
        const float precursor_mass = (float)precursor.getMZ() * charge - ((charge - 1) * Constants::PROTON_MASS_U);
        for (int16_t isotope_error = min_isotope_error_; isotope_error <= max_isotope_error_; ++isotope_error)
        {
          const float shifted_mass = precursor_mass
            + static_cast<float>(isotope_error) * static_cast<float>(Constants::C13C12_MASSDIFF_U);
          const auto window = computeMassWindow_(shifted_mass);
          const float margin = 1e-3f + 1e-6f * std::fabs(shifted_mass);
          const float lo = shifted_mass + window.first - margin;
          const float hi = shifted_mass + window.second + margin;
          // a NaN bound gives the query an arbitrary peptide range: keep them all
          if (std::isnan(lo) || std::isnan(hi)) { return; }
          windows.emplace_back(lo, hi);
        }
      }
    }
    if (windows.empty()) { return; } // nothing will be searched
    std::sort(windows.begin(), windows.end());

    // Moves the peptides of the merged windows to the front, in their order (the windows ascend).
    const size_t num_peptides = fi_peptides_.size();
    struct KeptRange { size_t begin, end, new_begin; }; // [begin, end) moved to new_begin
    std::vector<KeptRange> kept_ranges;
    auto kept_end = fi_peptides_.begin();
    auto from = fi_peptides_.begin();
    for (size_t w = 0; w < windows.size();)
    {
      const float lo = windows[w].first;
      float hi = windows[w].second;
      for (++w; w < windows.size() && windows[w].first <= hi; ++w) { hi = std::max(hi, windows[w].second); }
      const auto first = std::lower_bound(from, fi_peptides_.end(), lo, [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
      from = std::upper_bound(first, fi_peptides_.end(), hi, [](float b, const Peptide& a) { return b < a.precursor_mz_; });
      if (first == from) { continue; }
      kept_ranges.push_back({static_cast<size_t>(first - fi_peptides_.begin()), static_cast<size_t>(from - fi_peptides_.begin()),
                             static_cast<size_t>(kept_end - fi_peptides_.begin())});
      kept_end = (kept_end == first) ? from : std::copy(first, from, kept_end);
    }
    fi_peptides_.erase(kept_end, fi_peptides_.end());

    // Occurrences removed by a deduplication before (ordered by their kept entry): those of the kept entries follow
    // them to their new index, the others go with them.
    if (!removed_occurrences_.empty())
    {
      size_t num_kept = 0;
      auto range = kept_ranges.cbegin();
      for (const RemovedOccurrence& occurrence : removed_occurrences_)
      {
        while (range != kept_ranges.cend() && range->end <= occurrence.peptide_idx) { ++range; }
        if (range == kept_ranges.cend()) { break; }
        if (occurrence.peptide_idx < range->begin) { continue; }
        RemovedOccurrence moved = occurrence;
        moved.peptide_idx = static_cast<UInt32>(range->new_begin + (occurrence.peptide_idx - range->begin));
        removed_occurrences_[num_kept++] = moved;
      }
      removed_occurrences_.resize(num_kept);
    }
    OPENMS_LOG_INFO << "The precursor windows of the spectra reach " << fi_peptides_.size() << " of " << num_peptides << " peptides." << std::endl;
  }

  vector<FragmentIndex::Hit> FragmentIndex::query(const OpenMS::Peak1D& peak,
                                                  const pair<size_t, size_t>& peptide_idx_range,
                                                  uint16_t peak_charge)
  {
      float adjusted_mass = peak.getMZ() * (float)peak_charge -((peak_charge-1) * Constants::PROTON_MASS_U);

      float frag_tol = fragment_mz_tolerance_unit_ppm_ ? Math::ppmToMass(fragment_mz_tolerance_, adjusted_mass) : fragment_mz_tolerance_;

      // The fragments are matched against the same bounds that select the buckets: a test written
      // differently (adjusted_mass - frag_tol <= fragment_mz_ is not the same in float arithmetic as
      // adjusted_mass >= fragment_mz_ - frag_tol) would let a fragment within an ulp of a bound count
      // only if its bucket happens to be visited, i.e. depend on the bucket layout.
      const float mz_lo = adjusted_mass - frag_tol;
      const float mz_hi = adjusted_mass + frag_tol;
      auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_lo);
      auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_hi);

      if (left_it != bucket_min_mz_.begin()) --left_it;

      auto in_range_buckets = make_pair(std::distance(bucket_min_mz_.begin(), left_it), std::distance(bucket_min_mz_.begin(), right_it));

      // Public API entry point; the internal search path (queryPeaks) counts matches
      // inline instead of materializing one Hit vector per (peak, fragment charge).
      // Sizing the reservation from the candidate range is not viable: the range is
      // O(index) in open search (millions of peptides) while the number of fragments
      // inside one tolerance window is at most a few thousand.
      vector<FragmentIndex::Hit> hits;


      for (UInt32 j = in_range_buckets.first; j < in_range_buckets.second; j++)
      {
        auto slice_begin = fi_fragments_.begin() + (j*bucketsize_);
        auto slice_end = ((j+1) * bucketsize_) >= fi_fragments_.size() ? fi_fragments_.end() : (fi_fragments_.begin() + ((j+1) * bucketsize_)) ;

        auto left_iter = std::lower_bound(slice_begin, slice_end, peptide_idx_range.first, [](Fragment a, UInt32 b) { return a.peptide_idx_ < b;} );

        while (left_iter != slice_end) // sequential scan
        {
          // peptide_idx_range is half-open [first, second) — stop BEFORE index second.
          if (left_iter->peptide_idx_ >= peptide_idx_range.second) break;

          if (left_iter->fragment_mz_ >= mz_lo && left_iter->fragment_mz_ <= mz_hi)
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

  namespace
  {
    /// Hint that the cache line @p byte_offset bytes from @p address is read soon (build, queryPeaks).
    /// The line need not exist: a prefetch never faults.
    inline void prefetchForRead(const void* address, std::ptrdiff_t byte_offset)
    {
      const std::uintptr_t line = reinterpret_cast<std::uintptr_t>(address) + byte_offset;
#if defined(__GNUC__) || defined(__clang__)
      __builtin_prefetch(reinterpret_cast<const void*>(line));
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
      _mm_prefetch(reinterpret_cast<const char*>(line), _MM_HINT_T0);
#else
      (void)line;
#endif
    }
  }

  void FragmentIndex::queryPeaks(SpectrumMatchesTopN& candidates, const MSSpectrum& spectrum,
                                const std::vector<CandidateBlock_>& blocks,
                                const uint16_t precursor_charge,
                                const bool with_electron_ions)
  {
      // One call == all (isotope error) blocks of one precursor charge: count matched fragments
      // per candidate peptide and APPEND, block by block, only the candidates that clear the
      // emit threshold below. Materializing the dense [first, second) ranges instead would cost
      // one 24-byte zero entry per peptide in the precursor window (millions of them per spectrum
      // in open search) that trimHits drops again right afterwards.
      //
      // The blocks share ONE walk over the fragment buckets: which buckets a peak reaches and
      // which fragments it matches depends on the peak and the fragment charge only, not on the
      // isotope error. A block-by-block walk would visit every bucket once per block.
      if (blocks.empty()) return;

      // blocks are disjoint and ascending (see searchDifferentPrecursorRanges): their cells lie
      // back to back in the count table and [lo, hi) spans all of them.
      const UInt32 lo = static_cast<UInt32>(blocks.front().first);
      const size_t hi = blocks.back().second;
      const size_t window = blocks.back().cell_offset + (blocks.back().second - blocks.back().first);

      // Persistent thread-local buffers indexed BLOCK-RELATIVE (cell = cell_offset + peptide_idx
      // - first). Their capacity is the high-water candidate window seen on this thread and only
      // ever grows, so a closed search costs a few kB per thread no matter how large the index
      // is. Nothing here is derived from fi_peptides_.size(), so a rebuilt, cleared or chunked
      // index cannot invalidate them.
      //
      // Invariant: every nonzero cell is listed in touched_ids. The reset below
      // therefore restores the all-zero state in O(touched) rather than an O(window) memset,
      // and it does so for ANY subsequent window size.
      thread_local std::vector<uint32_t> match_counts;   // matched-peak count per cell
      thread_local std::vector<UInt32> touched_ids;      // cells written since the last reset
      thread_local std::vector<UInt32> emit_ids;         // subset of touched_ids that is emitted

      // Cells left over from a wider previous window still index within the table and
      // are still nonzero, so they must be cleared here regardless of the current window.
      for (UInt32 cell : touched_ids) match_counts[cell] = 0;
      touched_ids.clear();

      // Release the high-water table once this thread enters genuinely-small-window
      // territory (a closed search following an open search): the open-search table is
      // O(index) per thread and would otherwise stay resident for the thread's lifetime.
      // Gated on closed-search mode AND the current window being small: high-mass
      // precursors in the thin tail of the index can produce sub-threshold windows
      // MID-open-search, and releasing there would make the next ordinary spectrum
      // re-allocate the table. A fresh vector is all-zero, preserving the reset invariant.
      constexpr size_t small_window = 32768;       // closed-search windows are ~1e3-1e4
      constexpr size_t release_bytes = 1u << 20;   // keep tables below ~1 MB regardless
      if (!isOpenSearchMode_() && window < small_window
          && match_counts.size() * sizeof(uint32_t) > release_bytes)
      {
        std::vector<uint32_t>(window).swap(match_counts);
        touched_ids.shrink_to_fit();
        emit_ids.clear();
        emit_ids.shrink_to_fit();
      }
      // Grow-only otherwise: resize value-initializes just the new tail, and the reset
      // above already restored every older cell to zero.
      else if (window > match_counts.size()) match_counts.resize(window);

      // Loop-invariant: same cap as the previous per-peak std::min(). A precursor charge of
      // 0 yields an empty fragment-charge loop, hence no candidates — as before.
      const uint16_t actual_max = std::min(precursor_charge, max_fragment_charge_);

      // the thread-local buffers, resolved once for the loops below
      uint32_t* const counts = match_counts.data();
      std::vector<UInt32>& touched = touched_ids;

      // A bucket visit: the bucket, where in it the candidate ranges are expected to start,
      // and the m/z window of the peak its fragments are matched against (the bounds that
      // selected the bucket).
      struct Visit
      {
        const Fragment* begin;
        const Fragment* guess;
        const Fragment* end;
        float mz_lo;
        float mz_hi;
      };

      // Tolerance window and half-open peptide-range test are identical to
      // FragmentIndex::query(). A bucket is sorted by peptide_idx_ and the blocks are disjoint
      // and ascending, so one forward pass serves every block: each fragment whose peptide lies
      // in a block and whose m/z matches is counted once into that block — exactly the
      // (peak, fragment) pairs a walk per block counts.
      auto scan = [&](const Visit& v)
      {
          // it = std::lower_bound(v.begin, v.end, lo), reached from the guess: exact for any guess
          const Fragment* it = v.guess;
          if (it != v.end && it->peptide_idx_ < lo)
          {
            do { ++it; } while (it != v.end && it->peptide_idx_ < lo);
          }
          else
          {
            while (it != v.begin && (it - 1)->peptide_idx_ >= lo) --it;
          }
          for (const CandidateBlock_& block : blocks)
          {
            uint32_t* const block_counts = counts + block.cell_offset;
            while (it != v.end && it->peptide_idx_ < block.first) ++it;   // between two blocks

            // candidate ranges are half-open [first, second) — stop BEFORE index second.
            for (; it != v.end && it->peptide_idx_ < block.second; ++it)
            {
              if (it->fragment_mz_ >= v.mz_lo && it->fragment_mz_ <= v.mz_hi)
              {
                uint32_t& count = block_counts[it->peptide_idx_ - block.first];
                if (count == 0) touched.push_back(static_cast<UInt32>(&count - counts));
                ++count;   // uint32_t, same type as SpectrumMatch::num_matched_ — no saturation
              }
            }
            if (it == v.end) return;
          }
      };

      // The counts are integers, so the order in which the bucket visits are scanned does not
      // matter. A visit is queued (and the fragments around its guess are prefetched) and scanned
      // only after the next pipeline_depth visits were located: the cache misses of that many
      // buckets overlap instead of being taken one after the other.
      constexpr size_t pipeline_depth = 16;
      Visit pipeline[pipeline_depth];
      size_t num_queued = 0;
      size_t num_scanned = 0;

      // Bucket range of a peak as in FragmentIndex::query() — same buckets visited.
      auto count_matches = [&](const std::vector<Fragment>& fragments, const std::vector<float>& bucket_min_mz,
                               const std::vector<UInt32>& skip, float adjusted_mass, float frag_tol)
      {
          const float mz_lo = adjusted_mass - frag_tol;
          const float mz_hi = adjusted_mass + frag_tol;
          auto left_it = std::lower_bound(bucket_min_mz.begin(), bucket_min_mz.end(), mz_lo);
          auto right_it = std::upper_bound(bucket_min_mz.begin(), bucket_min_mz.end(), mz_hi);

          if (left_it != bucket_min_mz.begin()) --left_it;

          const size_t bucket_begin = std::distance(bucket_min_mz.begin(), left_it);
          const size_t bucket_end = std::distance(bucket_min_mz.begin(), right_it);

          for (size_t j = bucket_begin; j < bucket_end; j++)
          {
            const Fragment* slice_begin = fragments.data() + (j*bucketsize_);
            const Fragment* slice_end = fragments.data() + std::min((j+1) * bucketsize_, fragments.size());

            // Where the candidate ranges start in the bucket is looked up in its skip table (a
            // few cache-resident entries) instead of by a binary search over the 32 kB bucket.
            // The samples only steer: skipped buckets hold no fragment of a candidate, and scan
            // corrects the guess to the exact lower bound.
            const UInt32* samples = skip.data() + j * skip_per_bucket_;
            UInt32 num_below = 0;
            for (size_t k = 0; k < skip_per_bucket_; ++k) num_below += (samples[k] < lo);
            if (num_below == skip_per_bucket_) continue;   // last fragment of the bucket < lo
            const Fragment* guess = slice_begin;
            if (num_below == 0)
            {
              if (samples[0] >= hi) continue;              // first fragment of the bucket >= hi
            }
            else
            {
              // the lower bound lies between two sampled fragments: interpolate its position
              const size_t last = static_cast<size_t>(slice_end - slice_begin) - 1;
              const size_t pos_below = std::min((num_below - 1) * SKIP_STRIDE_, last);
              const size_t pos_above = std::min(num_below * SKIP_STRIDE_, last);
              const float fraction = static_cast<float>(lo - samples[num_below - 1])
                                   / static_cast<float>(samples[num_below] - samples[num_below - 1]);
              guess += pos_below + 1 + static_cast<size_t>(fraction * static_cast<float>(pos_above - pos_below - 1));
            }
            // the lower bound is a few fragments off the guess, the candidates follow it
            prefetchForRead(guess, -32);
            prefetchForRead(guess, 32);
            prefetchForRead(guess, 96);

            if (num_queued - num_scanned == pipeline_depth) scan(pipeline[num_scanned++ % pipeline_depth]);
            pipeline[num_queued++ % pipeline_depth] = Visit{slice_begin, guess, slice_end, mz_lo, mz_hi};
          }
      };

      for (const Peak1D& peak : spectrum)
      {
        for (uint16_t fragment_charge = 1; fragment_charge <= actual_max; fragment_charge++)
        {
          float adjusted_mass = peak.getMZ() * (float)fragment_charge -((fragment_charge-1) * Constants::PROTON_MASS_U);

          float frag_tol = fragment_mz_tolerance_unit_ppm_ ? Math::ppmToMass(fragment_mz_tolerance_, adjusted_mass) : fragment_mz_tolerance_;

          count_matches(fi_fragments_, bucket_min_mz_, bucket_skip_, adjusted_mass, frag_tol);
          // ions:electron_ions: the c and z+1 ions count only when the caller asks for them
          if (with_electron_ions) count_matches(electron_fragments_, electron_bucket_min_mz_, electron_bucket_skip_, adjusted_mass, frag_tol);
        }
      }
      while (num_scanned < num_queued) scan(pipeline[num_scanned++ % pipeline_depth]);

      // Prefilter totals of all candidates with at least one matched fragment (SpectrumMatchesTopN): touched lists
      // exactly these cells, and each cell's count is its number of matched (peak, fragment) pairs.
      uint64_t matched_pairs = 0;
      for (UInt32 cell : touched) matched_pairs += counts[cell];
      candidates.scored_candidates_ += touched.size();
      candidates.matched_peaks_ += matched_pairs;

      // trimHits sorts by num_matched_ descending first and then drops everything below
      // min_matched_peaks_, so a below-threshold candidate can neither displace an
      // above-threshold one from the top-N nor survive the trim: thresholding here leaves
      // the surviving set unchanged. The threshold is clamped to 1 because a candidate that
      // matched no peak at all carries no information — with min_matched_peaks_ == 0 the
      // unclamped filter would emit the entire precursor window.
      const uint32_t emit_min = std::max<uint32_t>(min_matched_peaks_, 1u);

      auto emit = [&](const CandidateBlock_& block, size_t cell)
      {
        SpectrumMatch& sm = candidates.hits_.emplace_back();
        sm.num_matched_ = counts[cell];
        sm.precursor_charge_ = precursor_charge;
        sm.isotope_error_ = block.isotope_error;
        sm.peptide_idx_ = block.first + (cell - block.cell_offset);
      };

      // Blocks are emitted in the order given (ascending isotope error) and ascending peptide
      // index within a block: the order the dense per-block array had. Emitting in that order
      // keeps the fi_peptides_ / fasta_entries accesses of the downstream scoring pass
      // sequential; correctness no longer rides on it, because trimHits' comparator now ends
      // in peptide_idx_ and therefore admits no equal-key candidates at all on this path — the
      // top-N cut is the same whatever order the blocks appended in. Cells ascend with the
      // block and, inside a block, with the peptide index, so ascending cells are that order.
      constexpr size_t sweep_factor = 8;
      if (window <= sweep_factor * touched.size())
      {
        // A good part of the cells was written (always so in a closed search, with its few-kB
        // table): reading the cells in order is cheaper than gathering and sorting the survivors.
        for (const CandidateBlock_& block : blocks)
        {
          const size_t cell_end = block.cell_offset + (block.second - block.first);
          for (size_t cell = block.cell_offset; cell < cell_end; ++cell)
          {
            if (counts[cell] >= emit_min) emit(block, cell);
          }
        }
      }
      else
      {
        // Threshold BEFORE ordering: touched_ids holds every candidate with at least one matched
        // fragment (up to millions in open search) while the survivors are orders of magnitude
        // fewer, and touched_ids itself must stay complete for the next call's reset.
        emit_ids.clear();
        for (UInt32 cell : touched)
        {
          if (counts[cell] >= emit_min) emit_ids.push_back(cell);
        }
        std::sort(emit_ids.begin(), emit_ids.end());

        const CandidateBlock_* block = blocks.data();
        for (UInt32 cell : emit_ids)
        {
          while (cell >= block->cell_offset + (block->second - block->first)) ++block;
          emit(*block, cell);
        }
      }
  }

  void FragmentIndex::trimHits(OpenMS::FragmentIndex::SpectrumMatchesTopN& init_hits) const
  {
      // Single ranking predicate for both branches below (they used to carry two copies of it,
      // which is exactly how they would drift apart). Keys, most significant first:
      //   1. num_matched_      descending — more matched fragments wins
      //   2. |isotope_error_|  ascending  — prefer the assignment closest to monoisotopic
      //   3. isotope_error_    ascending  — resolve -k vs +k by sign instead of by luck
      //   4. precursor_charge_ ascending
      //   5. peptide_idx_      ascending  — final key, added so that neither std::sort nor
      //      std::partial_sort (both unstable) can let the arrangement of the input decide
      //      which of several equally scored candidates survives the top-N cut. On the
      //      non-SNES path this makes the comparator a strict total order: two hits sharing
      //      all five keys would have to be the same peptide at the same isotope error and
      //      charge, i.e. the same (charge, iso) block, and queryPeaks emits each peptide at
      //      most once per block. querySpectrumSNES_ can still emit one mother twice within a
      //      block at two different Σ (sigma_delta_) values; those remain tied here, but their
      //      emission order is itself fixed, so the outcome is reproducible either way.
      auto by_rank = [](const SpectrumMatch& a, const SpectrumMatch& b)
      {
        if (a.num_matched_ != b.num_matched_)
        {
          return a.num_matched_ > b.num_matched_;
        }
        // Prefer isotope_error close to 0: abs(isotope_error), then isotope_error, then precursor_charge
        const auto abs_iso_a = a.isotope_error_ < 0 ? -a.isotope_error_ : a.isotope_error_;
        const auto abs_iso_b = b.isotope_error_ < 0 ? -b.isotope_error_ : b.isotope_error_;
        if (abs_iso_a != abs_iso_b) return abs_iso_a < abs_iso_b;
        if (a.isotope_error_ != b.isotope_error_) return a.isotope_error_ < b.isotope_error_;
        if (a.precursor_charge_ != b.precursor_charge_) return a.precursor_charge_ < b.precursor_charge_;
        return a.peptide_idx_ < b.peptide_idx_;
      };

      if (init_hits.hits_.size() > max_processed_hits_)
      {
        std::partial_sort(init_hits.hits_.begin(), init_hits.hits_.begin() + max_processed_hits_,
                          init_hits.hits_.end(), by_rank);

        init_hits.hits_.resize(max_processed_hits_);
      }
      else
      {
        std::sort(init_hits.hits_.begin(), init_hits.hits_.end(), by_rank);
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
                                                     uint16_t charge,
                                                     bool with_electron_ions)
  {
    // Open mode absorbs isotope shifts into the wide window — no per-isotope iteration.
    const bool open_mode = isOpenSearchMode_();
    const int16_t iso_lo = open_mode ? 0 : min_isotope_error_;
    const int16_t iso_hi = open_mode ? 0 : max_isotope_error_;

    // SNES mode uses querySpectrumSNES_ directly (dispatched in querySpectrum);
    // this function is only reached for non-SNES searches.
    //
    // The isotope-error blocks of this charge are collected and handed to queryPeaks together,
    // which walks the fragment buckets once for all of them and appends its (already compacted
    // and threshold-filtered) matches block by block, in ascending isotope error, directly to
    // the caller's accumulator.
    thread_local std::vector<CandidateBlock_> blocks;
    blocks.clear();
    for (int16_t isotope_error = iso_lo; isotope_error <= iso_hi; ++isotope_error)
    {
      const float shifted_mass = precursor_mass
        + static_cast<float>(isotope_error) * static_cast<float>(Constants::C13C12_MASSDIFF_U);

      const auto window = computeMassWindow_(shifted_mass);

      // candidates_range is half-open [first, second) — the scan in queryPeaks stops
      // strictly before peptide_idx == second.
      auto candidates_range = getPeptidesInMassWindow(shifted_mass, window);
      if (candidates_range.first >= candidates_range.second) continue;   // empty block: no candidates

      // queryPeaks serves all blocks in one forward pass per bucket, which needs them disjoint
      // and ascending. The windows of consecutive isotope errors are (they are ~1 Da apart);
      // with a precursor tolerance wide enough to make them overlap, the blocks collected so
      // far are searched first — the appended blocks stay in ascending isotope error.
      if (!blocks.empty() && candidates_range.first < blocks.back().second)
      {
        queryPeaks(sms, spectrum, blocks, charge, with_electron_ions);
        blocks.clear();
      }
      const size_t cell_offset = blocks.empty() ? 0 : blocks.back().cell_offset + (blocks.back().second - blocks.back().first);
      blocks.push_back({candidates_range.first, candidates_range.second, cell_offset, isotope_error});
    }
    queryPeaks(sms, spectrum, blocks, charge, with_electron_ions);
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
    const size_t n_mothers = fi_peptides_.size();
    const size_t n_words = (n_mothers + 63) / 64;

    // Persistent thread-local buffers, sized once. Between spectra they are
    // restored to all-zero by O(touched) resets (touched_ids / emitted_touched),
    // avoiding a full-index memset per spectrum (the score_table.assign + the
    // per-(charge,iso,sigma) std::fill(emitted) were the dominant non-scan cost
    // at proteome scale).
    thread_local std::vector<uint16_t> score_table;     // matched-peak count per viable mother
    thread_local std::vector<uint64_t> viable_words;    // precursor-viability bitset, 1 bit / mother
    thread_local std::vector<uint8_t> emitted;          // Phase-2 per-(charge,iso,sigma) dedup guard
    thread_local std::vector<UInt32> touched_ids;       // mothers marked viable (== only ids ever written)
    thread_local std::vector<UInt32> emitted_touched;   // ids set in emitted since its last reset

    // Restore the buffers to all-zero for this spectrum. They persist across queries
    // AND across different / rebuilt FragmentIndex instances on the same thread, so the
    // index size can change between calls (e.g. a later, smaller chunk in a chunked
    // search). On a size change we must NOT walk the stale touched lists — their ids
    // index the previous (possibly larger) size, so they would read/write out of bounds.
    if (score_table.size() != n_mothers || viable_words.size() != n_words || emitted.size() != n_mothers)
    {
      score_table.assign(n_mothers, 0);   // full zero-fill makes the touched lists redundant
      viable_words.assign(n_words, 0);
      emitted.assign(n_mothers, 0);
      touched_ids.clear();
      emitted_touched.clear();
    }
    else
    {
      // Same index as the previous query: restore only the touched entries (O(touched)
      // instead of a full-index memset — the whole point of the optimization).
      for (UInt32 id : touched_ids) { score_table[id] = 0; viable_words[id >> 6] &= ~(uint64_t{1} << (id & 63)); }
      touched_ids.clear();
      for (UInt32 id : emitted_touched) emitted[id] = 0;
      emitted_touched.clear();
    }

    auto viable_test = [&](UInt32 id) -> bool { return (viable_words[id >> 6] >> (id & 63)) & uint64_t{1}; };
    auto viable_set  = [&](UInt32 id) { uint64_t& w = viable_words[id >> 6]; const uint64_t b = uint64_t{1} << (id & 63);
                                        if (!(w & b)) { w |= b; touched_ids.push_back(id); } };
    auto emit_mark   = [&](UInt32 id) { emitted[id] = 1; emitted_touched.push_back(id); };
    auto reset_emitted_touched = [&]() { for (UInt32 id : emitted_touched) emitted[id] = 0; emitted_touched.clear(); };

    // Single source of truth for the SNES precursor-derived bin-walk targets, shared by
    // the viability pre-pass AND Phase-2 so the two cannot drift (the superset guarantee
    // depends on byte-identical target m/z). SNES build() emits Single-N b-ions with
    // c_term_mod=0 and Single-C y-ions with n_term_mod=0, so for a realized (M+H)+:
    //   b_k (Single-N) = (M+H)+ - water - fixed_cterm - Sigma ; y_k (Single-C) = (M+H)+ - fixed_nterm - Sigma
    const float snes_water = static_cast<float>(Residue::getInternalToFull().getMonoWeight());
    auto target_single_n = [&](float shifted_mh, float s) { return shifted_mh - snes_water - static_cast<float>(fixed_cterm_delta_) - s; };
    auto target_single_c = [&](float shifted_mh, float s) { return shifted_mh - static_cast<float>(fixed_nterm_delta_) - s; };
    auto target_full     = [&](float shifted_mh, float s) { return shifted_mh - s; };

    // Fragment-charge upper bound for the byte scan. Use the max charge in the
    // `charges` list (the spectrum's known charge, or max_precursor_charge_ when
    // unknown), capped by max_fragment_charge_. This matches the per-iteration
    // `std::min(precursor_charge, max_fragment_charge_)` policy of non-SNES
    // queryPeaks — using the static max_precursor_charge_ would inflate byte
    // counts with matches at fragment charges above the actual precursor charge.
    uint16_t byte_scan_max_frag_charge = 0;
    for (uint16_t c : charges) byte_scan_max_frag_charge = std::max(byte_scan_max_frag_charge, c);
    byte_scan_max_frag_charge = std::min(byte_scan_max_frag_charge, max_fragment_charge_);

    // ===================== Pre-pass: mark precursor-viable mothers =====================
    // The Phase-1 byte scan below only writes a count for a mother that this pre-pass
    // marked viable. The marked set is a strict SUPERSET of the mothers Phase-2 can
    // emit: we walk the SAME precursor-derived fragment targets (Single-N b_k,
    // Single-C y_k) and the SAME full-length precursor_mz_ ranges, but drop the
    // single_c / protein-anchor / score-threshold filters (over-approximation). So no
    // mother that Phase-2 emits is ever missing its count, while the non-viable mothers
    // no longer take a (cache-missing) random write per matched fragment. The bin-walk
    // targets come from the shared target_single_n / target_single_c / target_full
    // helpers used by Phase-2 too, so the two paths cannot drift out of sync.
    {
      const bool open_mode_pp = isOpenSearchMode_();
      const int16_t iso_lo_pp = open_mode_pp ? 0 : min_isotope_error_;
      const int16_t iso_hi_pp = open_mode_pp ? 0 : max_isotope_error_;

      std::vector<double> sigma_union = snes_sigma_delta_set_;
      sigma_union.insert(sigma_union.end(), snes_sigma_delta_set_with_prot_nterm_.begin(), snes_sigma_delta_set_with_prot_nterm_.end());
      sigma_union.insert(sigma_union.end(), snes_sigma_delta_set_with_prot_cterm_.begin(), snes_sigma_delta_set_with_prot_cterm_.end());

      auto mark_bucket_range = [&](float target, float tol_lo, float tol_hi) {
        // same bounds for the bucket selection and the fragment test (see query())
        const float mz_lo = target + tol_lo;
        const float mz_hi = target + tol_hi;
        auto lb = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_lo);
        auto rb = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_hi);
        if (lb != bucket_min_mz_.begin()) --lb;
        const size_t jb = std::distance(bucket_min_mz_.begin(), lb);
        const size_t je = std::distance(bucket_min_mz_.begin(), rb);
        for (size_t j = jb; j < je; ++j)
        {
          const auto sb = fi_fragments_.begin() + (j * bucketsize_);
          const auto se = ((j + 1) * bucketsize_) >= fi_fragments_.size()
            ? fi_fragments_.end() : (fi_fragments_.begin() + ((j + 1) * bucketsize_));
          for (auto it = sb; it != se; ++it)
          {
            if (it->fragment_mz_ >= mz_lo && it->fragment_mz_ <= mz_hi) viable_set(it->peptide_idx_);
          }
        }
      };
      auto mark_precursor_range = [&](float target, float tol_lo, float tol_hi) {
        auto lb = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(), target + tol_lo,
                                   [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
        auto ub = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(), target + tol_hi,
                                   [](float b, const Peptide& a) { return b < a.precursor_mz_; });
        for (auto it = lb; it != ub; ++it)
          viable_set(static_cast<UInt32>(std::distance(fi_peptides_.begin(), it)));
      };

      for (uint16_t charge : charges)
      {
        const float mh_plus = static_cast<float>(precursor.getMZ()) * charge
          - (charge - 1) * static_cast<float>(Constants::PROTON_MASS_U);
        for (int16_t iso_err = iso_lo_pp; iso_err <= iso_hi_pp; ++iso_err)
        {
          const float shifted_mh = mh_plus
            + static_cast<float>(iso_err) * static_cast<float>(Constants::C13C12_MASSDIFF_U);
          const auto prec_window = computeMassWindow_(shifted_mh);
          const float tlo = prec_window.first;   // <= 0
          const float thi = prec_window.second;  // >= 0
          for (double sigma : sigma_union)
          {
            const float s = static_cast<float>(sigma);
            mark_bucket_range(target_single_n(shifted_mh, s), tlo, thi); // Single-N b_k
            mark_bucket_range(target_single_c(shifted_mh, s), tlo, thi); // Single-C y_k
            mark_precursor_range(target_full(shifted_mh, s), tlo, thi);  // full-length
          }
        }
      }
    }

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
        // falls inside the bucket starting at bucket_min. The fragments are matched
        // against the same bounds.
        const float mz_lo = adjusted_mass - frag_tol;
        const float mz_hi = adjusted_mass + frag_tol;
        auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_lo);
        auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_hi);
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
            if (it->fragment_mz_ >= mz_lo && it->fragment_mz_ <= mz_hi)
            {
              const UInt32 id = it->peptide_idx_;
              if (!viable_test(id)) continue;   // precursor-prefilter: skip the non-viable mothers
              auto& cell = score_table[id];
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
    // the minimum-matched-peaks threshold. The target m/z are computed via the
    // shared target_single_n / target_single_c / target_full helpers (setup block).

    // Dedup guard `emitted` is declared and reset (O(touched), via emitted_touched)
    // in the setup block above. The per-(charge, iso_err, sigma) reset below is also
    // O(touched) rather than a full-index std::fill.

    // Helper: compute the iso-shifted observed (M+H)+ for a given (charge, iso_err).
    // Used by the subset-enumeration post-pass to reconstruct the realization target.
    auto shifted_mh_for = [&](uint16_t charge, int16_t iso_err) -> float {
      const float mh_plus = static_cast<float>(precursor.getMZ()) * charge
        - (charge - 1) * static_cast<float>(Constants::PROTON_MASS_U);
      return mh_plus + static_cast<float>(iso_err) * static_cast<float>(Constants::C13C12_MASSDIFF_U);
    };

    // Asymmetric tolerance: tol_lo <= 0 (low-side magnitude, sign-flipped),
    // tol_hi >= 0. A match requires fragment_mz ∈ [target_mz + tol_lo, target_mz + tol_hi],
    // the bounds that also select the buckets.
    // Preserves calibrated windows like [100 ppm, 5 ppm] where the symmetric
    // max-collapse over-admitted ~20× on the tighter side.
    auto collect_candidates =
      [&](float target_mz, float tol_lo, float tol_hi, bool expect_single_c,
          int16_t iso_err, uint16_t charge,
          SnesAnchor require_anchor, float sigma_tag)
    {
      const float mz_lo = target_mz + tol_lo;
      const float mz_hi = target_mz + tol_hi;
      auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_lo);
      auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), mz_hi);
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
          if (it->fragment_mz_ < mz_lo || it->fragment_mz_ > mz_hi) continue;

          const UInt32 id = it->peptide_idx_;
          if (emitted[id]) continue;

          const auto& mother = fi_peptides_[id];
          if (isSingleCMother(mother.mod_bitmask_) != expect_single_c) continue;

          // SNES v1.1: anchor-specific filter for PROTEIN_N/C_TERM mod walks.
          if (require_anchor == SnesAnchor::PROT_NTERM && !isProteinNTerminal_(fasta_entries[mother.protein_idx].sequence, mother.sequence_.first))
          {
            continue;
          }
          if (require_anchor == SnesAnchor::PROT_CTERM)
          {
            const uint32_t prot_len = protein_lengths_[mother.protein_idx];
            if (static_cast<uint32_t>(mother.sequence_.first) + mother.sequence_.second != prot_len) continue;
          }

          if (score_table[id] < min_matched_peaks_) continue;

          emit_mark(id);
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
          reset_emitted_touched();

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
          collect_candidates(target_single_n(shifted_mh, s),
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/false, iso_err, charge,
                             SnesAnchor::NONE, s);
          collect_candidates(target_single_c(shifted_mh, s),
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
            const float target = target_full(shifted_mh, s);
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
              emit_mark(id);
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
        // protein position 0, or 1 after enabled Met clipping). Empty when no PROTEIN_N_TERM variable mods
        // are configured.
        for (double sigma : prot_nterm_extra)
        {
          reset_emitted_touched();
          const float s = static_cast<float>(sigma);
          collect_candidates(target_single_n(shifted_mh, s),
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/false, iso_err, charge,
                             SnesAnchor::PROT_NTERM, s);
          // Supplementary full-length at this Σ, PROT_NTERM-gated.
          {
            const float target = target_full(shifted_mh, s);
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
              if (!isProteinNTerminal_(fasta_entries[fi_peptides_[id].protein_idx].sequence, fi_peptides_[id].sequence_.first))
              {
                continue; // PROT_NTERM anchor
              }
              if (isSingleCMother(fi_peptides_[id].mod_bitmask_)) continue; // Single-N only
              if (score_table[id] < min_matched_peaks_) continue;
              emit_mark(id);
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
          reset_emitted_touched();
          const float s = static_cast<float>(sigma);
          collect_candidates(target_single_c(shifted_mh, s),
                             prec_tol_lo, prec_tol_hi, /*expect_single_c=*/true, iso_err, charge,
                             SnesAnchor::PROT_CTERM, s);
          // Supplementary full-length at this Σ, PROT_CTERM-gated.
          {
            const float target = target_full(shifted_mh, s);
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
              emit_mark(id);
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
        const bool is_prot_nterm = isProteinNTerminal_(protein_seq, sub_start);
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
    querySpectrum(spectrum, fasta_entries, sms, false);
  }

  void FragmentIndex::querySpectrum(const OpenMS::MSSpectrum& spectrum,
                                    const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                    OpenMS::FragmentIndex::SpectrumMatchesTopN& sms,
                                    bool with_electron_ions)
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

      // Charge outer, isotope error inner (inside searchDifferentPrecursorRanges); each
      // block appends to sms in ascending peptide index. A per-charge staging container
      // would only add a second full copy of every candidate before the same tail insert.
      for (uint16_t charge : charges)
      {
        float mz;
        mz = (float)precursor[0].getMZ() * charge - ((charge-1) * Constants::PROTON_MASS_U);
        searchDifferentPrecursorRanges(spectrum, mz, sms, charge, with_electron_ions);
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
    
    defaults_.setValue("ions:add_z_ions", "false", "Add peaks of z-ions (y - NH3) to the spectrum. For ETD, EThcD or ETciD spectra, use ions:add_zp1_ions instead.");
    defaults_.setValidStrings("ions:add_z_ions", {"true","false"});

    defaults_.setValue("ions:add_zp1_ions", "false", "Add peaks of z+1 ions (z-dot, y - NH2) to the spectrum, the main C-terminal fragments of ETD, EThcD and ETciD spectra.");
    defaults_.setValidStrings("ions:add_zp1_ions", {"true","false"});

    defaults_.setValue("ions:electron_ions", "false", "Also index c and z+1 ions, the main fragments of ETD, ECD, EThcD and ETciD spectra, apart from the ion series above: querySpectrum() matches a spectrum against them only when asked to, e.g. for an electron-activated spectrum. Series enabled above are not indexed twice. Ignored in SNES mode.");
    defaults_.setValidStrings("ions:electron_ions", {"true","false"});
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

    vector<std::string> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'. A terminus carries one modification: a variable modification of the whole terminus (e.g. 'Acetyl (Protein N-term)') is not searched where a fixed one sits on it (e.g. 'TMT6plex (N-term)').");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");

    vector<std::string> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    defaults_.setValue("enzyme", "Trypsin", "Enzyme for digestion");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));


    defaults_.setValue("peptide:missed_cleavages", 1, "Missed cleavages for digestion");
    defaults_.setValue(
      "peptide:clip_nterm_methionine", "false",
      "Also consider loss of the initial M of a protein. Length and missed-cleavage limits apply to the clipped peptide, "
      "which remains eligible for protein N-terminal variable modifications. Non-specific searches already include these sequences.");
    defaults_.setValidStrings("peptide:clip_nterm_methionine", {"true", "false"});
    defaults_.setValue("peptide:deduplicate", "false",
                       "Index each exact modified peptide once, retaining one representative protein coordinate. "
                       "Callers must recover complete protein mappings separately. Does not apply to SNES mother indices.",
                       {"advanced"});
    defaults_.setValidStrings("peptide:deduplicate", {"true", "false"});
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
    // Default iso range [-2, 0]: Orbitrap/QExactive/tims monoisotopic peak picking
    // fails predominantly *upward* (picks the +1 or +2 isotope instead of the true
    // monoisotopic). The query matches observed mass + isotope_error * C13C12, so an
    // upward mispick needs a negative isotope error. Symmetric ranges like [-1, +1]
    // waste a query slot on the rare downward mispick. Matches the MetaMorpheus/MSFragger
    // defaults (0, +1, +2 there, counted as observed minus theoretical).
    defaults_.setValue("precursor:isotope_error_min", -2,
                       "Minimum precursor isotope error searched, in 13C spacings (1.00336 Da) added to the observed "
                       "precursor mass: -1 finds a peptide whose first 13C isotope peak was selected as the precursor.");
    defaults_.setValue("precursor:isotope_error_max", 0,
                       "Maximum precursor isotope error searched, with the sign of precursor:isotope_error_min: +1 finds "
                       "a peptide whose precursor was selected one 13C spacing below its monoisotopic peak.");

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
    // Kept in sync with ProSEAlgorithm's own "annotate:PSM" valid-strings list
    // (ProSEAlgorithm.cpp) — both declare this parameter since the merged Search:
    // subtree is passed down to FragmentIndex as well.
    defaults_.setValidStrings("annotate:PSM",
                              std::vector<std::string>{
                                "ALL",
                                Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM,
                                Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
                                Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
                                Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION,
                                Constants::UserParam::NUM_MATCHED_PEAKS,
                                Constants::UserParam::MATCHED_PREFIX_IONS,
                                Constants::UserParam::MATCHED_SUFFIX_IONS,
                                Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE,
                                Constants::UserParam::MATCHED_ION_CURRENT,
                                Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM,
                                Constants::UserParam::HYPERSCORE_ZSCORE,
                                Constants::UserParam::LN_NUM_CANDIDATES,
                                Constants::UserParam::MATCHED_ION_CURRENT_FRACTION,
                                Constants::UserParam::COMPLEMENTARY_IONS_FRACTION}
    );
    defaults_.setValue("report:top_hits", 1, "Maximum number of top scoring hits per spectrum that are reported.");
    defaults_.setSectionDescription("report", "Reporting Options");
    defaults_.setValue("peptide:motif", "", "If set, only peptides that contain this motif (provided as RegEx) will be considered.");
    defaults_.setSectionDescription("peptide", "Peptide Options");

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
    add_zp1_ions_ = param_.getValue("ions:add_zp1_ions").toBool();
    electron_ions_ = param_.getValue("ions:electron_ions").toBool();
    digestion_enzyme_ = param_.getValue("enzyme").toString();
    clip_nterm_methionine_ = param_.getValue("peptide:clip_nterm_methionine").toBool();
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
        "series (a/c/x/z/z+1) may be enabled freely for downstream scoring.");
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

  bool FragmentIndex::hasProteinOccurrences(const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const
  {
    if (!is_build_ || is_snes_mode_ || protein_lengths_.size() != fasta_entries.size()) return false;
    for (Size i = 0; i < fasta_entries.size(); ++i)
    {
      if (protein_lengths_[i] != fasta_entries[i].sequence.size()) return false;
    }
    // The configured modifications themselves (not the tables built from them): a protein-terminal one, fixed or
    // variable, may make the entries of a span depend on its position in the protein.
    for (const StringList* mods : {&modifications_fixed_, &modifications_variable_})
    {
      for (const auto& mod_residue : ModifiedPeptideGenerator::getModifications(*mods).val)
      {
        const ResidueModification::TermSpecificity term = mod_residue.first->getTermSpecificity();
        if (term == ResidueModification::PROTEIN_N_TERM || term == ResidueModification::PROTEIN_C_TERM) return false;
      }
    }
    return true;
  }

  void FragmentIndex::getProteinOccurrences(Size peptide_idx, const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                            std::vector<std::pair<UInt32, UInt32>>& occurrences) const
  {
    const Size first_new = occurrences.size();
    const Peptide& peptide = fi_peptides_[peptide_idx];
    const float mz = peptide.precursor_mz_;
    const size_t length = peptide.sequence_.second;
    const char* residues = fasta_entries[peptide.protein_idx].sequence.data() + peptide.sequence_.first;
    // The entries of the spans with these residues and slots: same precursor_mz_ (bitwise), length, mod_bitmask_ and
    // residues, in the run of equal precursor_mz_ around peptide_idx (fi_peptides_ is sorted by precursor_mz_).
    // With peptide:deduplicate, that is the entry itself.
    size_t first = peptide_idx;
    while (first > 0 && fi_peptides_[first - 1].precursor_mz_ == mz) --first;
    for (size_t i = first; i < fi_peptides_.size() && fi_peptides_[i].precursor_mz_ == mz; ++i)
    {
      const Peptide& entry = fi_peptides_[i];
      if (entry.sequence_.second == length && entry.mod_bitmask_ == peptide.mod_bitmask_
          && std::memcmp(fasta_entries[entry.protein_idx].sequence.data() + entry.sequence_.first, residues, length) == 0)
      {
        occurrences.emplace_back(entry.protein_idx, entry.sequence_.first);
      }
    }
    // The occurrences of its peptidoform that peptide:deduplicate removed (removed_occurrences_ is ordered by the kept
    // entry). A span with several entries of the peptidoform (a modification configured fixed and variable renders
    // alike) is listed once.
    const UInt32 kept = static_cast<UInt32>(peptide_idx);
    const auto removed_begin = std::lower_bound(removed_occurrences_.begin(), removed_occurrences_.end(), kept,
      [](const RemovedOccurrence& occurrence, const UInt32 index) { return occurrence.peptide_idx < index; });
    if (removed_begin == removed_occurrences_.end() || removed_begin->peptide_idx != kept) return;
    for (auto it = removed_begin; it != removed_occurrences_.end() && it->peptide_idx == kept; ++it)
    {
      occurrences.emplace_back(it->protein_idx, it->start);
    }
    std::sort(occurrences.begin() + first_new, occurrences.end());
    occurrences.erase(std::unique(occurrences.begin() + first_new, occurrences.end()), occurrences.end());
  }

}
