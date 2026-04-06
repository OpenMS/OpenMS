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
#include <functional>
#include <mutex>


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
    if (add_y_ions_ || add_x_ions_ || add_z_ions_)
    {
      {
        constexpr int z = 1;
        double base_mass = proton * z + c_term_mod_mass;
        double cumulative = base_mass;

        for (size_t j = seq_len; j > 1; --j)
        {
          double res_mass = table[static_cast<unsigned char>(sequence[j - 1])];
          if (residue_mod_masses) res_mass += residue_mod_masses[j - 1];
          cumulative += res_mass;

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

      ModifiedPeptideGenerator::MapToResidueType fixed_modifications;
      ModifiedPeptideGenerator::MapToResidueType variable_modifications;
      if (has_modifications)
      {
        fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
        variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
      }

      size_t skipped_peptides = 0;

      ProteaseDigestion digestor;
      digestor.setEnzyme(digestion_enzyme_);
      digestor.setMissedCleavages(missed_cleavages_);

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

          if (has_modifications)
          {
            // Modified path: still needs AASequence for modification application
            AASequence unmod_peptide = AASequence::fromString(protein.sequence.substr(digested_peptide.first, seq_len));
            vector<AASequence> modified_peptides;
            AASequence mod_peptide = AASequence(unmod_peptide);

            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, max_variable_mods_per_peptide_, modified_peptides);

            UInt32 modification_idx = 0;
            for (const AASequence& modified_peptide : modified_peptides)
            {
              float modified_mz = modified_peptide.getMZ(1);
              if (modified_mz < peptide_min_mass_ || modified_mz > peptide_max_mass_)
              {
                continue;
              }
              thread_peptides[tid].emplace_back(static_cast<UInt32>(protein_idx),
                                      modification_idx,
                                      std::make_pair(static_cast<uint16_t>(digested_peptide.first),
                                                     static_cast<uint16_t>(seq_len)),
                                      modified_mz);
              ++modification_idx;
            }
          }
          else
          {
            // Unmodified path: compute precursor m/z directly from lookup table
            float unmodified_mz = computePrecursorMzFromChars(seq_ptr, seq_len, residue_mass_table_);
            if (peptide_min_mass_ <= unmodified_mz && unmodified_mz <= peptide_max_mass_)
            {
              thread_peptides[tid].emplace_back(static_cast<UInt32>(protein_idx), 0,
                                      std::make_pair(static_cast<uint16_t>(digested_peptide.first),
                                                     static_cast<uint16_t>(seq_len)),
                                      unmodified_mz);
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
      /// generate all Peptides (also initializes residue mass table)
      generatePeptides(fasta_entries);

      const bool has_modifications = !(modifications_fixed_.empty() && modifications_variable_.empty());

      // Only needed for the modified path
      ModifiedPeptideGenerator::MapToResidueType fixed_modifications;
      ModifiedPeptideGenerator::MapToResidueType variable_modifications;
      if (has_modifications)
      {
        fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
        variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
      }

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

      if (has_modifications)
      {
        // Modified path: reconstruct AASequence to apply modifications,
        // then extract per-residue mass deltas for lightweight generation
        vector<AASequence> mod_peptides;

        #pragma omp parallel for private(mod_peptides)
        for (SignedSize peptide_idx = 0; peptide_idx < (SignedSize)fi_peptides_.size(); peptide_idx++)
        {
#ifdef _OPENMP
          const int tid = omp_get_thread_num();
#else
          const int tid = 0;
#endif
          const Peptide& pep = fi_peptides_[peptide_idx];
          mod_peptides.clear();

          const string& protein_seq = fasta_entries[pep.protein_idx].sequence;
          const char* seq_ptr = protein_seq.c_str() + pep.sequence_.first;
          size_t seq_len = pep.sequence_.second;

          AASequence unmod_peptide = AASequence::fromString(protein_seq.substr(pep.sequence_.first, seq_len));
          AASequence mod_peptide = AASequence(unmod_peptide);
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, max_variable_mods_per_peptide_, mod_peptides);

          const AASequence& final_peptide = mod_peptides[pep.modification_idx_];

          // Extract per-residue modification mass deltas
          vector<double> mod_masses(seq_len, 0.0);
          bool has_residue_mods = false;
          for (size_t i = 0; i < seq_len; ++i)
          {
            double modified_mass = final_peptide[i].getMonoWeight(Residue::Internal);
            double unmodified_mass = residue_mass_table_[static_cast<unsigned char>(seq_ptr[i])];
            double delta = modified_mass - unmodified_mass;
            if (delta != 0.0)
            {
              mod_masses[i] = delta;
              has_residue_mods = true;
            }
          }

          double n_term_mod = 0.0;
          if (final_peptide.hasNTerminalModification())
            n_term_mod = final_peptide.getNTerminalModification()->getDiffMonoMass();

          double c_term_mod = 0.0;
          if (final_peptide.hasCTerminalModification())
            c_term_mod = final_peptide.getCTerminalModification()->getDiffMonoMass();

          generateFragmentsLightweight_(
            thread_fragments[tid], seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
            n_term_mod, c_term_mod, has_residue_mods ? mod_masses.data() : nullptr);
        }
      }
      else
      {
        // Unmodified path: compute fragments directly from protein sequence chars.
        // No AASequence parsing, no TheoreticalSpectrumGenerator needed.
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

          generateFragmentsLightweight_(
            thread_fragments[tid], seq_ptr, seq_len, static_cast<UInt32>(peptide_idx),
            0.0, 0.0, nullptr);
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

      /// 1.) First all Fragments are sorted by their own mass!
      sort(fi_fragments_.begin(), fi_fragments_.end(), [](const Fragment& a, const Fragment& b)
      {
        return std::tie(a.fragment_mz_, a.peptide_idx_) < std::tie(b.fragment_mz_, b.peptide_idx_);
      });

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

  std::pair<size_t, size_t > FragmentIndex::getPeptidesInPrecursorRange(float precursor_mass,
                                                                       const std::pair<float, float>& window)
  {
      float prec_tol = precursor_mz_tolerance_unit_ppm_ ? Math::ppmToMass(precursor_mz_tolerance_, precursor_mass) : precursor_mz_tolerance_ ;

      auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass - prec_tol + window.first, [](const Peptide& a, float b) { return a.precursor_mz_ < b;});
      auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass + prec_tol + window.second, [](float b, const Peptide& a) { return b < a.precursor_mz_;});
      return make_pair(std::distance(fi_peptides_.begin(), left_it), std::distance(fi_peptides_.begin(), right_it));
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
          if(left_iter->peptide_idx_ > peptide_idx_range.second) break;

          if ((adjusted_mass >= left_iter->fragment_mz_ - frag_tol ) && adjusted_mass <= (left_iter->fragment_mz_+ frag_tol))
          {

            hits.emplace_back(left_iter->peptide_idx_, left_iter->fragment_mz_);
            #ifdef DEBUG_FRAGMENT_INDEX
            if (left_iter->peptide_idx_ < peptide_idx_range.first || left_iter->peptide_idx_ > peptide_idx_range.second)
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
      int16_t  min_isotope_error_applied;
      int16_t  max_isotope_error_applied;
      float precursor_window_upper_applied;
      float precursor_window_lower_applied;
      if (isOpenSearchMode_())
      {
        min_isotope_error_applied = 0;
        max_isotope_error_applied = 0;
        precursor_window_upper_applied = open_precursor_window_upper_;
        precursor_window_lower_applied = open_precursor_window_lower_;
      }
      else
      {
        min_isotope_error_applied = min_isotope_error_;
        max_isotope_error_applied = max_isotope_error_;
        precursor_window_upper_applied = 0;
        precursor_window_lower_applied = 0;
      }

      for (int16_t isotope_error = min_isotope_error_applied; isotope_error <= max_isotope_error_applied; isotope_error++)
      {
        SpectrumMatchesTopN candidates_iso_error;
        float precursor_mass_isotope_error = precursor_mass + ((float)isotope_error * (float)Constants::C13C12_MASSDIFF_U);
        auto candidates_range = getPeptidesInPrecursorRange(precursor_mass_isotope_error, {precursor_window_lower_applied, precursor_window_upper_applied}); // for the simple search we do not apply any modification window!!
        candidates_iso_error.hits_.resize(candidates_range.second - candidates_range.first + 1);

        queryPeaks(candidates_iso_error, spectrum, candidates_range, isotope_error, charge);

        // take only top 50 hits
        //trimHits(candidates_iso_error);

        sms += candidates_iso_error;
      }
      //trimHits(sms);
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


    defaults_.setValue("precursor:mass_tolerance", 10.0, "Tolerance for precursor-m/z in search");
    std::vector<std::string> precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

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
    // Open search window bounds (used when tolerance > 1 Da or > 1000 ppm)
    defaults_.setValue("precursor:open_window_lower", -100.0, "lower bound of the open precursor window");
    defaults_.setValue("precursor:open_window_upper", 200.0, "upper bound of the open precursor window");

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
    missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_min_mass_ = param_.getValue("peptide:min_mass");
    peptide_max_mass_ = param_.getValue("peptide:max_mass");
    peptide_min_length_ = param_.getValue("peptide:min_size");
    peptide_max_length_ = param_.getValue("peptide:max_size");
    fragment_min_mz_ = param_.getValue("fragment:min_mz");
    fragment_max_mz_ = param_.getValue("fragment:max_mz");
    
    precursor_mz_tolerance_ = param_.getValue("precursor:mass_tolerance");
    fragment_mz_tolerance_ = param_.getValue("fragment:mass_tolerance");
    precursor_mz_tolerance_unit_ppm_ = param_.getValue("precursor:mass_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment:mass_tolerance_unit").toString() == "ppm";
    
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
    // Open search mode is automatically determined in isOpenSearchMode_()
    if (isOpenSearchMode_())
    {
      OPENMS_LOG_INFO << "[FragmentIndex] Open-search mode enabled because precursor mass tolerance ("
                      << precursor_mz_tolerance_ << " "
                      << (precursor_mz_tolerance_unit_ppm_ ? "ppm" : "Da")
                      << ") exceeds threshold (1000 ppm or 1 Da)." << std::endl;
    }
    open_precursor_window_lower_ = param_.getValue("precursor:open_window_lower");
    open_precursor_window_upper_ = param_.getValue("precursor:open_window_upper");
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
