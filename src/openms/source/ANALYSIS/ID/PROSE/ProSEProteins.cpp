// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>
#include "ProSEImpl.h"

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <functional>
#include <limits>
#include <map>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  namespace
  {
    // ProteaseDigestion::isValidProduct(protein, pos, length, true, true, false) of a fully specific digestion with the
    // enzyme of @p rule, without the regex: the start is a site (or the protein start, or within the first three residues
    // of a protein that starts with M), and the end is a site or the protein end. (isValidProduct_() tokenizes
    // [pos - 2, end + 2): every site p inside it with pos - 2 < p < end + 2 sees the residues p - 1 and p.)
    // Requires 0 <= pos, 0 < length, pos + length <= protein.size().
    bool isValidFullySpecific(const SegmentRule& rule, const std::string& protein, Size pos, Size length)
    {
      const auto site = [&rule, &protein](Size p)
      {
        return rule.after[static_cast<unsigned char>(protein[p - 1])] && !rule.not_before[static_cast<unsigned char>(protein[p])];
      };
      const Size end = pos + length;
      const bool spec_n = pos == 0 || (pos <= 2 && protein[0] == 'M') || site(pos);
      const bool spec_c = end == protein.size() || site(end);
      return spec_n && spec_c;
    }

  // =====================================================================
  // Protein mapping from the fragment index (peptide:protein_mapping = index):
  // what PeptideIndexing::run() writes, without PeptideIndexing (the spans of the hits from the index, or with a
  // protein-terminal modification from a pass over the database for the hits' sequences)
  // =====================================================================
    bool isAmbiguousResidue(const char c)
    {
      return c == 'B' || c == 'J' || c == 'Z' || c == 'X';
    }

    // The residues AhoCorasickAmbiguous matches an ambiguous protein residue with: B = D or N, J = I or L,
    // Z = E or Q, X = any of its unambiguous amino acids (every other letter, O and U included)
    const std::string& residuesMatchedBy(const char ambiguous)
    {
      static const std::string b = "DN", j = "IL", z = "EQ";
      static const std::string x = []
      {
        std::string all;
        for (char c = 'A'; c <= 'Z'; ++c)
        {
          if (!AA(c).isAmbiguous()) all += c;
        }
        return all;
      }();
      switch (ambiguous)
      {
        case 'B': return b;
        case 'J': return j;
        case 'Z': return z;
        default: return x;
      }
    }

    // The first @p length (<= 8) residues of @p residues as a number that sorts like them
    uint64_t residueKey(const std::string_view residues, const Size length)
    {
      uint64_t key = 0;
      for (Size i = 0; i < length; ++i) key = (key << 8) | static_cast<unsigned char>(residues[i]);
      return key;
    }

    using Occurrence = std::pair<UInt32, UInt32>; // {protein index, position}

    /**
      @brief The spans that step 3 of buildProteinMapping_() looks up for one ambiguous residue of a protein: those
      with no ambiguous residue before it and 1 to aaa_max ambiguous residues from it on (so every span over ambiguous
      residues belongs to its first one).

      With r_0 = position < r_1 < ... the ambiguous residues from the position on, a span [start, end) with start <=
      r_0 contains exactly k of them if r_{k-1} < end <= r_k (r_k = the protein end if there is none). Each span is
      looked up once per combination of the residues its ambiguous ones stand for. Used for both the count (before
      any lookup, to decide on the fallback) and the lookups, so the two agree.
    */
    class AmbiguousSpans
    {
    public:
      /// PeptideIndexing's maximum of aaa_max
      static constexpr Size MAX_AMBIGUOUS = 10;

      AmbiguousSpans(const std::string& protein, const Size position, const Size max_length, const Size aaa_max) :
        protein_(protein), first_start_(position), aaa_max_(std::min(aaa_max, MAX_AMBIGUOUS))
      {
        // the first start: after the previous ambiguous residue, and close enough to reach the position
        while (first_start_ > 0 && position - first_start_ + 1 < max_length && !isAmbiguousResidue(protein[first_start_ - 1]))
        {
          --first_start_;
        }
        // the ambiguous residues from the position on that a span can reach: up to aaa_max + 1 (the last one only
        // bounds the spans with aaa_max)
        const Size reach = std::min(protein.size(), position + max_length);
        ambiguous_[0] = position;
        for (Size i = position + 1; i < reach && count_ <= aaa_max_; ++i)
        {
          if (isAmbiguousResidue(protein[i])) ambiguous_[count_++] = i;
        }
      }

      /// The largest number of ambiguous residues in a span
      Size maxAmbiguous() const { return std::min(aaa_max_, count_); }

      /// The k-th ambiguous residue from the position on (k = 0: the position)
      Size ambiguous(const Size k) const { return ambiguous_[k]; }

      /// The starts [first, last] of the spans of @p length with exactly @p k (1..maxAmbiguous()) ambiguous residues
      /// (first > last if there are none)
      std::pair<SignedSize, SignedSize> starts(const Size length, const Size k) const
      {
        const SignedSize end_after = static_cast<SignedSize>(ambiguous_[k - 1]);
        const SignedSize end_at_most = static_cast<SignedSize>(k < count_ ? ambiguous_[k] : protein_.size());
        const SignedSize len = static_cast<SignedSize>(length);
        return {std::max(static_cast<SignedSize>(first_start_), end_after + 1 - len),
                std::min(static_cast<SignedSize>(ambiguous_[0]), end_at_most - len)};
      }

      /// The sequences looked up for the spans of @p length (combinations saturate at 2^32, far above any budget)
      Size lookups(const Size length) const
      {
        Size total = 0, combinations = 1;
        for (Size k = 1; k <= maxAmbiguous(); ++k)
        {
          combinations = std::min(combinations * residuesMatchedBy(protein_[ambiguous_[k - 1]]).size(), Size(1) << 32);
          const auto [first, last] = starts(length, k);
          if (first <= last) total += combinations * static_cast<Size>(last - first + 1);
        }
        return total;
      }

    private:
      const std::string& protein_;
      Size first_start_;
      Size aaa_max_;
      std::array<Size, MAX_AMBIGUOUS + 1> ambiguous_{};
      Size count_ = 1; ///< entries of ambiguous_
    };
  }

  ProSEAlgorithm::ProteinMapping_ ProSEAlgorithm::buildProteinMapping_(const FragmentIndex& index,
                                                                      const std::vector<FASTAFile::FASTAEntry>& db,
                                                                      std::vector<UInt32> candidates,
                                                                      const Param& indexer_parameters)
  {
    ProteinMapping_ mapping;

    // PeptideIndexing settings the mapping reproduces (see ProteinMapping_)
    if (static_cast<Int>(indexer_parameters.getValue("mismatches_max")) != 0
        || indexer_parameters.getValue("IL_equivalent").toBool()
        || !indexer_parameters.getValue("allow_nterm_protein_cleavage").toBool()
        || indexer_parameters.getValue("write_protein_sequence").toBool()
        || indexer_parameters.getValue("write_protein_description").toBool()
        || indexer_parameters.getValue("keep_unreferenced_proteins").toBool()
        || indexer_parameters.getValue("decoy_string").toString().empty())
    {
      mapping.fallback_reason = "PeptideIndexing settings other than no mismatches, no I/L equivalence, N-terminal methionine "
                                "cleavage, no protein sequences or descriptions, no unreferenced proteins and a given decoy string";
      return mapping;
    }
    const std::string enzyme_name = indexer_parameters.getValue("enzyme:name").toString();
    const std::string specificity = indexer_parameters.getValue("enzyme:specificity").toString();
    const Param index_parameters = index.getParameters();
    if (specificity != EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL]
        || enzyme_name == EnzymaticDigestion::UnspecificCleavage || enzyme_name == EnzymaticDigestion::NoCleavage
        || index_parameters.getValue("enzyme").toString() != enzyme_name
        || index_parameters.getValue("peptide:enzyme_specificity").toString() != specificity)
    {
      mapping.fallback_reason = "not a fully specific digestion with the enzyme of the fragment index";
      return mapping;
    }
    // The spans with the residues of a hit come from the index entries if they list them all; with a protein-terminal
    // modification they do not (the entries of a span depend on its position in its protein), and the spans are found
    // by a pass over the database instead (1. below).
    const bool index_occurrences = index.hasProteinOccurrences(db);
    if (!index_occurrences && !index.isBuiltFrom(db))
    {
      mapping.fallback_reason = "the fragment index does not list every protein occurrence of its peptides "
                                "(low_memory, or not built from this database)";
      return mapping;
    }
    const Size aaa_max = static_cast<Size>(static_cast<Int>(indexer_parameters.getValue("aaa_max")));
    ProteaseDigestion enzyme;
    enzyme.setEnzyme(enzyme_name);
    enzyme.setSpecificity(EnzymaticDigestion::SPEC_FULL);

    // The candidates, once each, and the lengths of their sequences. A hit with an ambiguous residue leaves the
    // mapping to PeptideIndexing; checked first, since it needs no pass over the database.
#ifdef PROSE_GNU_PARALLEL
    if (candidates.size() >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::sort(candidates.begin(), candidates.end()); }
    else
#endif
    { std::sort(candidates.begin(), candidates.end()); }
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.empty()) return mapping; // nothing to map (applyProteinMapping_() decides)
    const std::vector<FragmentIndex::Peptide>& peptides = index.getPeptides();
    if (candidates.back() >= peptides.size())
    {
      mapping.fallback_reason = "a candidate the fragment index does not hold";
      return mapping;
    }
    const auto residues_of = [&](const UInt32 candidate)
    {
      const FragmentIndex::Peptide& peptide = peptides[candidate];
      return std::string_view(db[peptide.protein_idx].sequence).substr(peptide.sequence_.first, peptide.sequence_.second);
    };
    bool unambiguous_candidates = true;
#pragma omp parallel for schedule(dynamic, 256) reduction(&& : unambiguous_candidates)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      const std::string_view residues = residues_of(candidates[i]);
      unambiguous_candidates = unambiguous_candidates && std::none_of(residues.begin(), residues.end(), isAmbiguousResidue);
    }
    if (!unambiguous_candidates)
    {
      mapping.fallback_reason = "a hit with an ambiguous residue";
      return mapping;
    }
    Size min_length = std::numeric_limits<Size>::max(), max_length = 0;
    for (const UInt32 candidate : candidates)
    {
      min_length = std::min<Size>(min_length, peptides[candidate].sequence_.second);
      max_length = std::max<Size>(max_length, peptides[candidate].sequence_.second);
    }
    std::vector<uint8_t> length_used(max_length + 1, 0);
    for (const UInt32 candidate : candidates) length_used[peptides[candidate].sequence_.second] = 1;

    // The database: letters A-Z only (PeptideIndexing removes '*' and skips other symbols, which shifts positions),
    // no stretch of more than aaa_max X (PeptideIndexing splits proteins at such stretches), and at most a budget of
    // sequences that 3. looks up for the spans over ambiguous residues (beyond it, PeptideIndexing does this better).
    // The lookups are counted here, without looking anything up, so that a fallback costs little more than this pass.
    constexpr Size lookup_budget = Size(1) << 22;
    std::array<uint8_t, 256> residue_class{}; // 0: unambiguous letter, 1: ambiguous letter, 2: anything else
    residue_class.fill(2);
    for (char c = 'A'; c <= 'Z'; ++c) residue_class[static_cast<unsigned char>(c)] = isAmbiguousResidue(c) ? 1 : 0;
    std::vector<uint8_t> with_spans(db.size(), 0); // has spans for 3.
    bool letters_only = true, short_stretches = true;
    std::atomic<Size> lookups{0}; // exact while at most lookup_budget (then counting stops)
#pragma omp parallel for schedule(dynamic, 256) reduction(&& : letters_only, short_stretches)
    for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
    {
      const std::string& protein = db[p].sequence;
      uint8_t classes = 0; // OR of its residue classes
      for (const char c : protein) classes |= residue_class[static_cast<unsigned char>(c)];
      letters_only = letters_only && classes < 2;
      if (classes != 1) continue;
      Size stretch = 0;
      for (const char c : protein)
      {
        stretch = (c == 'X') ? stretch + 1 : 0;
        short_stretches = short_stretches && stretch <= aaa_max;
      }
      Size protein_lookups = 0;
      for (Size position = 0; position < protein.size(); ++position)
      {
        if (!isAmbiguousResidue(protein[position])) continue;
        if (protein_lookups > lookup_budget || lookups.load(std::memory_order_relaxed) > lookup_budget) break;
        const AmbiguousSpans spans(protein, position, max_length, aaa_max);
        for (Size length = min_length; length <= max_length; ++length)
        {
          if (length_used[length]) protein_lookups += spans.lookups(length);
        }
      }
      with_spans[p] = protein_lookups > 0;
      lookups.fetch_add(std::min(protein_lookups, lookup_budget + 1), std::memory_order_relaxed);
    }
    if (!letters_only || !short_stretches)
    {
      mapping.fallback_reason = letters_only ? "the database has stretches of more than aaa_max X"
                                             : "the database has symbols other than the letters A-Z";
      return mapping;
    }
    if (lookups.load() > lookup_budget)
    {
      mapping.fallback_reason = "too many spans over ambiguous residues in the database";
      return mapping;
    }

    // 1. the spans of the digest: per candidate, then per sequence (from the index; else after the sequences, below)
    std::vector<std::vector<Occurrence>> candidate_occurrences(candidates.size());
    if (index_occurrences)
    {
#pragma omp parallel for schedule(dynamic, 64)
      for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
      {
        index.getProteinOccurrences(candidates[i], db, candidate_occurrences[i]);
      }
    }
    // One entry per sequence (the candidates of a sequence differ in their modifications only: they have the same
    // spans). The sequences go to ProteinMapping_::SHARDS hash tables by their hash, filled in parallel; the entries
    // are numbered table by table, in the order of the candidates.
    constexpr Size shards = ProteinMapping_::SHARDS;
    std::vector<Size> shard_of(candidates.size());
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      shard_of[i] = std::hash<std::string_view>{}(residues_of(candidates[i])) % shards;
    }
    std::vector<std::vector<Size>> shard_candidates(shards);
    for (Size i = 0; i < candidates.size(); ++i) shard_candidates[shard_of[i]].push_back(i);
    mapping.entries.resize(shards);
    std::vector<std::vector<Size>> shard_entries(shards); // per table: the candidate of each of its entries
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize shard = 0; shard < static_cast<SignedSize>(shards); ++shard)
    {
      auto& table = mapping.entries[shard];
      table.reserve(shard_candidates[shard].size());
      for (const Size i : shard_candidates[shard])
      {
        const auto [entry, added] = table.try_emplace(residues_of(candidates[i]), shard_entries[shard].size());
        if (added)
        {
          shard_entries[shard].push_back(i);
        }
        else
        { // another modified form of the same residues: the same spans (merged and made unique below)
          std::vector<Occurrence>& first = candidate_occurrences[shard_entries[shard][entry->second]];
          first.insert(first.end(), candidate_occurrences[i].begin(), candidate_occurrences[i].end());
        }
      }
    }
    std::vector<Size> shard_start(shards + 1, 0);
    for (Size shard = 0; shard < shards; ++shard) shard_start[shard + 1] = shard_start[shard] + shard_entries[shard].size();
    mapping.sequences.resize(shard_start[shards]);
    mapping.occurrences.resize(shard_start[shards]);
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize shard = 0; shard < static_cast<SignedSize>(shards); ++shard)
    {
      for (auto& [sequence, entry] : mapping.entries[shard]) entry += shard_start[shard];
      for (Size k = 0; k < shard_entries[shard].size(); ++k)
      {
        const Size i = shard_entries[shard][k];
        mapping.sequences[shard_start[shard] + k] = residues_of(candidates[i]);
        mapping.occurrences[shard_start[shard] + k] = std::move(candidate_occurrences[i]);
      }
    }

    std::vector<std::tuple<Size, UInt32, UInt32>> extra; // spans of 1. (database pass), 2. and 3.: {sequence entry, protein, position}

    // 1. without the occurrences of the index: every span of the database with the residues of a sequence that
    //    PeptideIndexing accepts (isValidProduct() as in its search, which matches the residues literally here: the
    //    spans over ambiguous residues are those of 3.). Positions are prefiltered by their first residues (a bit set over
    //    the sequences' first prefix_length residues as a base-26 number), then looked up by the key of their first
    //    key_length residues. The database holds the letters A-Z only (checked above).
    if (!index_occurrences)
    {
      // isValidProduct() tokenizes with the regex per call (~20 us); a simple rule is evaluated directly
      SegmentRule rule;
      const bool simple_rule = parseSegmentRule(ProteaseDB::getInstance()->getEnzyme(enzyme_name)->getRegEx(), rule);
      const Size key_length = std::min<Size>(min_length, 8);
      std::vector<std::pair<uint64_t, Size>> keyed(mapping.sequences.size());
#pragma omp parallel for schedule(static)
      for (SignedSize entry = 0; entry < static_cast<SignedSize>(keyed.size()); ++entry)
      {
        keyed[entry] = {residueKey(mapping.sequences[entry], key_length), static_cast<Size>(entry)};
      }
#ifdef PROSE_GNU_PARALLEL
      if (keyed.size() >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::sort(keyed.begin(), keyed.end()); } // (distinct pairs)
      else
#endif
      { std::sort(keyed.begin(), keyed.end()); }
      const Size prefix_length = std::min<Size>(key_length, 5);
      Size prefix_count = 1;
      for (Size i = 0; i < prefix_length; ++i) prefix_count *= 26;
      const Size prefix_shift = prefix_count / 26; // 26^(prefix_length - 1)
      const auto prefix_of = [prefix_length](const char* residues)
      {
        Size code = 0;
        for (Size i = 0; i < prefix_length; ++i) code = code * 26 + static_cast<Size>(residues[i] - 'A');
        return code;
      };
      std::vector<uint64_t> prefix_bits((prefix_count + 63) / 64, 0);
      for (const std::string_view sequence : mapping.sequences)
      {
        const Size code = prefix_of(sequence.data());
        prefix_bits[code >> 6] |= uint64_t(1) << (code & 63);
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme; // as PeptideIndexing: one per thread
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
#pragma omp for schedule(dynamic, 16) nowait
        for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
        {
          const std::string& protein = db[p].sequence;
          if (protein.size() < min_length) continue;
          Size code = prefix_of(protein.data()) / 26; // the first prefix_length - 1 residues
          for (Size position = 0; position + min_length <= protein.size(); ++position)
          {
            code = (code % prefix_shift) * 26 + static_cast<Size>(protein[position + prefix_length - 1] - 'A');
            if (((prefix_bits[code >> 6] >> (code & 63)) & 1) == 0) continue;
            const uint64_t key = residueKey(std::string_view(protein).substr(position), key_length);
            for (auto k = std::lower_bound(keyed.begin(), keyed.end(), std::make_pair(key, Size(0)));
                 k != keyed.end() && k->first == key; ++k)
            {
              const std::string_view sequence = mapping.sequences[k->second];
              if (position + sequence.size() <= protein.size() && protein.compare(position, sequence.size(), sequence) == 0
                  && (simple_rule ? isValidFullySpecific(rule, protein, position, sequence.size())
                                  : thread_enzyme.isValidProduct(protein, static_cast<int>(position), static_cast<int>(sequence.size()), true, true, false)))
              {
                found.emplace_back(k->second, static_cast<UInt32>(p), static_cast<UInt32>(position));
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // 2. positions 1 and 2 of proteins that start with M: found through their first residues
    {
      // {key of the residues from the position on, protein, position}, sorted. Filled in parallel by chunks of
      // proteins, grouped by the first residue (the highest byte of the key), then each group sorted on its own.
      const Size key_length = std::min<Size>(min_length, 8);
      using Start = std::tuple<uint64_t, UInt32, UInt32>;
      const auto is_start = [&](const std::string& protein, Size position)
      {
        return position + min_length <= protein.size() && protein[0] == 'M';
      };
      constexpr Size letters = 26;
      const Size num_chunks = std::max<Size>(1, std::min<Size>(256, db.size() / 1024));
      std::vector<Size> offset(num_chunks * letters + 1, 0); // [letter][chunk]: count, then where the chunk fills it
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize chunk = 0; chunk < static_cast<SignedSize>(num_chunks); ++chunk)
      {
        for (Size p = db.size() * chunk / num_chunks; p < db.size() * (chunk + 1) / num_chunks; ++p)
        {
          const std::string& protein = db[p].sequence;
          for (Size position = 1; position <= 2; ++position)
          {
            if (is_start(protein, position)) ++offset[(protein[position] - 'A') * num_chunks + chunk + 1];
          }
        }
      }
      for (Size i = 1; i < offset.size(); ++i) offset[i] += offset[i - 1];
      std::vector<Start> starts(offset.back());
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize chunk = 0; chunk < static_cast<SignedSize>(num_chunks); ++chunk)
      {
        for (Size p = db.size() * chunk / num_chunks; p < db.size() * (chunk + 1) / num_chunks; ++p)
        {
          const std::string& protein = db[p].sequence;
          for (Size position = 1; position <= 2; ++position)
          {
            if (!is_start(protein, position)) continue;
            starts[offset[(protein[position] - 'A') * num_chunks + chunk]++] =
              Start(residueKey(std::string_view(protein).substr(position), key_length), static_cast<UInt32>(p), static_cast<UInt32>(position));
          }
        }
      }
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize letter = 0; letter < static_cast<SignedSize>(letters); ++letter)
      {
        // after the fill, offset[letter * num_chunks + num_chunks - 1] is where the letter's group ends
        const Size begin = letter == 0 ? 0 : offset[letter * num_chunks - 1];
        const Size end = offset[(letter + 1) * num_chunks - 1];
        std::sort(starts.begin() + begin, starts.begin() + end);
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme; // as PeptideIndexing: one per thread
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
#pragma omp for schedule(dynamic, 256) nowait
        for (SignedSize entry = 0; entry < static_cast<SignedSize>(mapping.sequences.size()); ++entry)
        {
          const std::string_view sequence = mapping.sequences[entry];
          const uint64_t key = residueKey(sequence, key_length);
          for (auto s = std::lower_bound(starts.begin(), starts.end(), std::make_tuple(key, UInt32(0), UInt32(0)));
               s != starts.end() && std::get<0>(*s) == key; ++s)
          {
            const UInt32 protein = std::get<1>(*s), position = std::get<2>(*s);
            const std::string& protein_sequence = db[protein].sequence;
            if (position + sequence.size() <= protein_sequence.size() && protein_sequence.compare(position, sequence.size(), sequence) == 0
                && thread_enzyme.isValidProduct(protein_sequence, static_cast<int>(position), static_cast<int>(sequence.size()), true, true, false))
            {
              found.emplace_back(static_cast<Size>(entry), protein, position);
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // 3. spans over ambiguous residues (AmbiguousSpans; within the budget checked above): each looked up with every
    //    combination of the residues its ambiguous ones stand for
    {
      // work items: an ambiguous residue and a length with spans (a few residues with several ambiguous ones around
      // them would otherwise be the critical path); at most one per two lookups, so bounded by the budget
      struct Item { UInt32 protein, position, length; };
      std::vector<Item> items;
#pragma omp parallel
      {
        std::vector<Item> thread_items;
#pragma omp for schedule(dynamic, 16) nowait
        for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
        {
          if (!with_spans[p]) continue;
          const std::string& protein = db[p].sequence;
          for (Size position = 0; position < protein.size(); ++position)
          {
            if (!isAmbiguousResidue(protein[position])) continue;
            const AmbiguousSpans spans(protein, position, max_length, aaa_max);
            for (Size length = min_length; length <= max_length; ++length)
            {
              if (length_used[length] && spans.lookups(length) > 0)
              {
                thread_items.push_back({static_cast<UInt32>(p), static_cast<UInt32>(position), static_cast<UInt32>(length)});
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        items.insert(items.end(), thread_items.begin(), thread_items.end());
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme;
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
        std::string window;
        std::array<Size, AmbiguousSpans::MAX_AMBIGUOUS> choice{};
#pragma omp for schedule(dynamic, 1) nowait
        for (SignedSize i = 0; i < static_cast<SignedSize>(items.size()); ++i)
        {
          const Item& item = items[i];
          const std::string& protein_sequence = db[item.protein].sequence;
          const AmbiguousSpans spans(protein_sequence, item.position, max_length, aaa_max);
          for (Size k = 1; k <= spans.maxAmbiguous(); ++k)
          {
            const auto [first, last] = spans.starts(item.length, k);
            for (SignedSize start = first; start <= last; ++start)
            {
              window.assign(protein_sequence, static_cast<Size>(start), item.length);
              std::fill_n(choice.begin(), k, 0);
              while (true)
              {
                for (Size w = 0; w < k; ++w)
                {
                  window[spans.ambiguous(w) - static_cast<Size>(start)] = residuesMatchedBy(protein_sequence[spans.ambiguous(w)])[choice[w]];
                }
                const Size entry = mapping.find(window);
                if (entry < mapping.sequences.size()
                    && thread_enzyme.isValidProduct(protein_sequence, static_cast<int>(start), static_cast<int>(item.length), true, true, false))
                {
                  found.emplace_back(entry, item.protein, static_cast<UInt32>(start));
                }
                Size w = 0; // the next combination
                for (; w < k; ++w)
                {
                  if (++choice[w] < residuesMatchedBy(protein_sequence[spans.ambiguous(w)]).size()) break;
                  choice[w] = 0;
                }
                if (w == k) break;
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // all spans of every sequence, ascending and once (2. and 3. may find spans of 1. again)
    // (grouped by entry through counts, then appended per entry in parallel: in the order of `extra` within an entry,
    // which the sort below makes irrelevant)
    if (!extra.empty())
    {
      std::vector<Size> first(mapping.occurrences.size() + 1, 0);
      for (const auto& found : extra) ++first[std::get<0>(found) + 1];
      for (Size k = 1; k < first.size(); ++k) first[k] += first[k - 1];
      std::vector<Occurrence> grouped(extra.size());
      {
        std::vector<Size> next(first.begin(), first.end() - 1);
        for (const auto& [entry, protein, position] : extra) grouped[next[entry]++] = Occurrence(protein, position);
      }
      std::vector<std::tuple<Size, UInt32, UInt32>>().swap(extra);
      FirstError error; // (an allocation failure)
#pragma omp parallel for schedule(dynamic, 256)
      for (SignedSize k = 0; k < static_cast<SignedSize>(mapping.occurrences.size()); ++k)
      {
        if (first[k] == first[k + 1]) continue;
        try
        {
          mapping.occurrences[k].insert(mapping.occurrences[k].end(), grouped.begin() + first[k], grouped.begin() + first[k + 1]);
        }
        catch (...)
        {
          error.set(k);
        }
      }
      error.rethrow();
    }
#pragma omp parallel for schedule(dynamic, 256)
    for (SignedSize k = 0; k < static_cast<SignedSize>(mapping.occurrences.size()); ++k)
    {
      std::vector<Occurrence>& occurrences = mapping.occurrences[k];
      std::sort(occurrences.begin(), occurrences.end());
      occurrences.erase(std::unique(occurrences.begin(), occurrences.end()), occurrences.end());
    }
    return mapping;
  }

  bool ProSEAlgorithm::applyProteinMapping_(const ProteinMapping_& mapping,
                                            const std::vector<FASTAFile::FASTAEntry>& db,
                                            const Param& indexer_parameters,
                                            std::vector<ProteinIdentification>& protein_ids,
                                            PeptideIdentificationList& peptide_ids,
                                            std::string& fallback_reason)
  {
    if (!mapping.fallback_reason.empty())
    {
      fallback_reason = mapping.fallback_reason;
      return false;
    }
    if (protein_ids.size() != 1)
    {
      fallback_reason = "not exactly one identification run";
      return false;
    }
    {
      // PeptideIndexing changes its cleavage rules for results of X! Tandem and MS-GF+
      std::string engine = protein_ids[0].getOriginalSearchEngineName();
      StringUtils::toUpper(engine);
      const ProteinIdentification::SearchParameters& search_parameters = protein_ids[0].getSearchParameters();
      if (engine == "XTANDEM" || engine == "MS-GF+" || engine == "MSGFPLUS"
          || search_parameters.metaValueExists("SE:XTandem") || search_parameters.metaValueExists("SE:MS-GF+"))
      {
        fallback_reason = "results of X! Tandem or MS-GF+";
        return false;
      }
    }

    // the sequence entry of every hit (PeptideIndexing maps the unmodified sequence)
    std::vector<Size> first_hit(peptide_ids.size() + 1, 0);
    for (Size i = 0; i < peptide_ids.size(); ++i) first_hit[i + 1] = first_hit[i] + peptide_ids[i].getHits().size();
    if (first_hit.back() == 0)
    {
      fallback_reason = "no peptide hits"; // PeptideIndexing's own handling (PEPTIDE_IDS_EMPTY)
      return false;
    }
    std::vector<Size> hit_entry(first_hit.back());
    bool all_mapped = true;
#pragma omp parallel for schedule(dynamic, 64) reduction(&& : all_mapped)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      const std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        const Size entry = mapping.find(hits[h].getSequence().toUnmodifiedString());
        const bool mapped = entry < mapping.sequences.size() && !mapping.occurrences[entry].empty();
        all_mapped = all_mapped && mapped;
        if (mapped) hit_entry[first_hit[i] + h] = entry;
      }
    }
    if (!all_mapped)
    {
      fallback_reason = "a hit whose sequence was not mapped from the fragment index";
      return false;
    }

    const std::string decoy_string = indexer_parameters.getValue("decoy_string").toString();
    const bool decoy_prefix = indexer_parameters.getValue("decoy_string_position").toString() == "prefix";
    std::vector<uint8_t> is_decoy(db.size());
#pragma omp parallel for schedule(static)
    for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
    {
      is_decoy[p] = decoy_prefix ? StringUtils::hasPrefix(db[p].identifier, decoy_string) : StringUtils::hasSuffix(db[p].identifier, decoy_string);
    }
    // the proteins of the hits
    std::vector<uint8_t> entry_used(mapping.occurrences.size(), 0);
    for (const Size entry : hit_entry) entry_used[entry] = 1;
    std::vector<uint8_t> protein_used(db.size(), 0);
    for (Size entry = 0; entry < entry_used.size(); ++entry)
    {
      if (!entry_used[entry]) continue;
      for (const Occurrence& occurrence : mapping.occurrences[entry]) protein_used[occurrence.first] = 1;
    }
    std::vector<UInt32> proteins; // in database order (PeptideIndexing's std::set of protein indices)
    bool any_decoy = false;
    for (Size p = 0; p < db.size(); ++p)
    {
      if (!protein_used[p]) continue;
      proteins.push_back(static_cast<UInt32>(p));
      any_decoy = any_decoy || is_decoy[p];
    }
    if (!any_decoy && indexer_parameters.getValue("missing_decoy_action").toString() != "silent")
    {
      fallback_reason = "no hit maps to a decoy"; // PeptideIndexing's warning or error
      return false;
    }

    // The meta value names in the order in which PeptideIndexing registers them (its first hit is mapped): MetaInfo
    // keeps, and idXML writes, the values of a hit ordered by registry index
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const UInt target_decoy = registry.registerName("target_decoy");
    const UInt protein_references = registry.registerName("protein_references");
    const DataValue target("target"), decoy("decoy"), target_and_decoy("target+decoy"), unique("unique"), non_unique("non-unique");

    // the evidences of every hit, in PeptideIndexing's order (protein, position)
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        PeptideHit& hit = hits[h];
        const std::vector<Occurrence>& occurrences = mapping.occurrences[hit_entry[first_hit[i] + h]];
        const Size length = hit.getSequence().size();
        std::vector<PeptideEvidence> evidences;
        evidences.reserve(occurrences.size());
        bool in_target = false, in_decoy = false;
        Size protein_count = 0;
        UInt32 last_protein = std::numeric_limits<UInt32>::max();
        for (const auto& [protein, position] : occurrences)
        {
          const std::string& sequence = db[protein].sequence;
          const char before = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : sequence[position - 1];
          const char after = (position + length >= sequence.size()) ? PeptideEvidence::C_TERMINAL_AA : sequence[position + length];
          evidences.emplace_back(db[protein].identifier, static_cast<Int>(position), static_cast<Int>(position + length) - 1, before, after);
          if (protein != last_protein)
          {
            last_protein = protein;
            ++protein_count;
          }
          (is_decoy[protein] ? in_decoy : in_target) = true;
        }
        hit.setPeptideEvidences(std::move(evidences));
        hit.setMetaValue(target_decoy, (in_target && in_decoy) ? target_and_decoy : (in_target ? target : decoy));
        hit.setMetaValue(protein_references, protein_count == 1 ? unique : non_unique);
      }
    }

    // the protein hits: the proteins of the hits (keep_unreferenced_proteins = false), in database order
    std::vector<ProteinHit> protein_hits(proteins.size());
#pragma omp parallel for schedule(static)
    for (SignedSize k = 0; k < static_cast<SignedSize>(proteins.size()); ++k)
    {
      protein_hits[k].setAccession(db[proteins[k]].identifier);
      protein_hits[k].setMetaValue(target_decoy, is_decoy[proteins[k]] ? decoy : target);
    }
    protein_ids[0].getHits() = std::move(protein_hits);

    // PeptideIndexing's settings, as it records them
    ProteaseDigestion enzyme;
    enzyme.setEnzyme(indexer_parameters.getValue("enzyme:name").toString());
    enzyme.setSpecificity(ProteaseDigestion::getSpecificityByName(indexer_parameters.getValue("enzyme:specificity").toString()));
    ProteinIdentification::SearchParameters search_parameters = protein_ids[0].getSearchParameters();
    search_parameters.setMetaValue("PeptideIndexer:decoy_string", decoy_string);
    search_parameters.setMetaValue("PeptideIndexer:decoy_string_position", decoy_prefix ? "prefix" : "suffix");
    search_parameters.setMetaValue("PeptideIndexer:enzyme", enzyme.getEnzymeName());
    search_parameters.setMetaValue("PeptideIndexer:enzyme_specificity", EnzymaticDigestion::NamesOfSpecificity[enzyme.getSpecificity()]);
    search_parameters.setMetaValue("PeptideIndexer:aaa_max", static_cast<Int>(indexer_parameters.getValue("aaa_max")));
    search_parameters.setMetaValue("PeptideIndexer:mismatches_max", static_cast<Int>(indexer_parameters.getValue("mismatches_max")));
    search_parameters.setMetaValue("PeptideIndexer:IL_equivalent", indexer_parameters.getValue("IL_equivalent").toBool() ? "true" : "false");
    search_parameters.setMetaValue("PeptideIndexer:allow_nterm_protein_cleavage",
                                   indexer_parameters.getValue("allow_nterm_protein_cleavage").toBool() ? "true" : "false");
    search_parameters.setMetaValue("PeptideIndexer:unmatched_action", indexer_parameters.getValue("unmatched_action").toString());
    search_parameters.setMetaValue("PeptideIndexer:missing_decoy_action", indexer_parameters.getValue("missing_decoy_action").toString());
    protein_ids[0].setSearchParameters(std::move(search_parameters));
    return true;
  }

  // =====================================================================
  // Shared protein-FDR finalization for a COMPLETE protein set. Single source
  // of truth for both the file-based search() below and the ProSE TOPP tool's
  // single-file path (the two previously held byte-identical copies).
  // =====================================================================
  // static
  std::string ProSEAlgorithm::inferProteins(std::vector<ProteinIdentification>& protein_ids,
                                            const PeptideIdentificationList& peptide_ids,
                                            const PeptideIdentificationList& inference_only_psms)
  {
    if (protein_ids.empty()) { return ""; }
    ProteinIdentification& run = protein_ids[0];

    // Copies of the best hits with what the inference reads (sequence, charge, peptide evidences, target/decoy and
    // the scores), in the run: the inference switches their scores and filters them, while the PSMs stay as they are.
    // PSMs without hits are left out, as the score switch reads the first PSM.
    PeptideIdentificationList psms;
    psms.reserve(peptide_ids.size() + inference_only_psms.size());
    std::vector<std::string> keys;
    const auto add = [&](const PeptideIdentification& id)
    {
      if (id.getHits().empty()) { return; }
      const PeptideHit* best = &id.getHits()[0];
      for (const PeptideHit& hit : id.getHits())
      {
        if (id.isHigherScoreBetter() ? hit.getScore() > best->getScore() : hit.getScore() < best->getScore()) { best = &hit; }
      }
      PeptideHit hit(best->getScore(), best->getRank(), best->getCharge(), best->getSequence());
      hit.setPeptideEvidences(best->getPeptideEvidences());
      best->getKeys(keys);
      for (const std::string& key : keys)
      {
        if (key == Constants::UserParam::TARGET_DECOY || Scores::isKnownScoreType(key)) { hit.setMetaValue(key, best->getMetaValue(key)); }
      }
      PeptideIdentification& copy = psms.emplace_back();
      copy.setIdentifier(run.getIdentifier());
      copy.setScoreType(id.getScoreType());
      copy.setHigherScoreBetter(id.isHigherScoreBetter());
      copy.setHits({std::move(hit)});
    };
    for (const PeptideIdentification& id : peptide_ids)
    {
      if (id.getIdentifier() == run.getIdentifier()) { add(id); } // the PSMs of other runs are not the run's evidence
    }
    for (const PeptideIdentification& id : inference_only_psms) { add(id); }

    // The protein score: the best PEP of its PSMs, if they have one (Percolator's, also when its q-value is the main
    // score), as ProteomicsLFQ scores proteins; otherwise the best search engine score (e.g. 'ln(hyperscore)', also
    // after FDR:PSM made the q-value the main score). The q-values themselves rank the proteins coarsely: many PSMs
    // share the best q-value.
    // Every PSM must have the score: when Percolator rescored only some files of a pooled set, only their PSMs have a
    // PEP, and all of them are then scored by ProSE's own score, which all have ('ln(hyperscore)', the main score, or
    // a meta value after rescoring or FDR:PSM).
    std::string score_type;
    if (!psms.empty())
    {
      const IDScoreSwitcherAlgorithm switcher;
      const auto all_have = [&](IDScoreSwitcherAlgorithm::ScoreType type)
      {
        return std::all_of(psms.begin(), psms.end(), [&](const PeptideIdentification& id)
                           { return !switcher.findScoreType(id, type).score_name.empty(); });
      };
      if (all_have(IDScoreSwitcherAlgorithm::ScoreType::PEP)) { score_type = "PEP"; }
      // (after useSearchEngineScore_(), ProSE's score is the main score of the copies, which "RAW" finds)
      else if (useSearchEngineScore_(psms) || all_have(IDScoreSwitcherAlgorithm::ScoreType::RAW)) { score_type = "RAW"; }
    }
    BasicProteinInferenceAlgorithm bpia;
    Param bpia_param = bpia.getParameters();
    bpia_param.setValue("score_type", score_type);
    bpia.setParameters(bpia_param);
    bpia.run(psms, run);
    return score_type;
  }

  // static
  void ProSEAlgorithm::removeProteinsReferencedByNeither_(std::vector<ProteinIdentification>& protein_ids,
                                                         const PeptideIdentificationList& peptide_ids,
                                                         const PeptideIdentificationList& also_referencing)
  {
    // the accessions each run's PSMs of either list reference (as IDFilter::removeUnreferencedProteins(), by identifier)
    std::map<std::string, std::unordered_set<std::string>> referenced;
    for (const PeptideIdentificationList* list : {&peptide_ids, &also_referencing})
    {
      for (const PeptideIdentification& id : *list)
      {
        std::unordered_set<std::string>& accessions = referenced[id.getIdentifier()];
        for (const PeptideHit& hit : id.getHits())
        {
          for (const PeptideEvidence& evidence : hit.getPeptideEvidences()) { accessions.insert(evidence.getProteinAccession()); }
        }
      }
    }
    for (ProteinIdentification& run : protein_ids)
    {
      const std::unordered_set<std::string>& accessions = referenced[run.getIdentifier()];
      std::vector<ProteinHit>& hits = run.getHits();
      hits.erase(std::remove_if(hits.begin(), hits.end(), [&](const ProteinHit& hit) { return !accessions.contains(hit.getAccession()); }),
                 hits.end());
    }
  }

  // static
  void ProSEAlgorithm::removeUnreportedProteins(std::vector<ProteinIdentification>& protein_ids,
                                                const PeptideIdentificationList& peptide_ids,
                                                const PeptideIdentificationList& also_referencing)
  {
    if (also_referencing.empty()) { IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids); }
    else { removeProteinsReferencedByNeither_(protein_ids, peptide_ids, also_referencing); }
    for (ProteinIdentification& run : protein_ids)
    {
      IDFilter::updateProteinGroups(run.getIndistinguishableProteins(), run.getHits());
      IDFilter::updateProteinGroups(run.getProteinGroups(), run.getHits());
    }
  }

  // static
  bool ProSEAlgorithm::applyCompleteSetProteinFDR(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      const std::string& decoy_string,
      bool decoy_is_prefix,
      double protein_fdr,
      bool remove_decoys,
      const PeptideIdentificationList& inference_only_psms)
  {
    // Defensive: an empty decoy marker would make every protein match the prefix/suffix guard
    // below (hasPrefix(x, "") is always true) and corrupt picked-protein FDR. Callers gate on
    // have_decoys (so decoy_string is non-empty in practice), but this is a public static helper.
    if (decoy_string.empty())
    {
      OPENMS_LOG_WARN << "[ProSE] applyCompleteSetProteinFDR called with an empty decoy marker; "
                      << "skipping protein FDR (cannot identify decoys)." << std::endl;
      return false;
    }

    // Aggregate the best PSM score (PEP, or the search engine score) per peptide per protein over the reported PSMs and
    // the inference-only ones, then apply picked-protein FDR (Savitski et al. 2015) over the full target+decoy protein
    // set. The inference uses the best hit of each PSM, and so do the reported PSMs from here on.
    IDFilter::keepNBestHits(peptide_ids, 1);
    const std::string score_type = inferProteins(protein_ids, peptide_ids, inference_only_psms);

    const auto is_decoy = [&](const ProteinHit& ph)
    {
      return decoy_is_prefix ? StringUtils::hasPrefix(ph.getAccession(), decoy_string)
                             : StringUtils::hasSuffix(ph.getAccession(), decoy_string);
    };
    const auto& inferred = protein_ids[0].getHits();
    const Size decoy_proteins = static_cast<Size>(std::count_if(inferred.begin(), inferred.end(), is_decoy));
    const auto with_hits = [](const PeptideIdentification& id) { return !id.getHits().empty(); };
    const Size reported = static_cast<Size>(std::count_if(peptide_ids.begin(), peptide_ids.end(), with_hits));
    const Size inference_only = static_cast<Size>(std::count_if(inference_only_psms.begin(), inference_only_psms.end(), with_hits));
    OPENMS_LOG_INFO << "[ProSE] Protein inference on " << reported + inference_only << " PSMs (" << inference_only
                    << " of them up to FDR:protein_inference_PSM, not reported), protein score "
                    << (score_type.empty() ? "the main score" : score_type) << ": " << inferred.size() - decoy_proteins
                    << " target and " << decoy_proteins << " decoy proteins." << std::endl;

    // Picked-protein FDR needs identified decoy PROTEINS, not merely a decoy database: if PSM-level
    // filtering removed all decoy evidence, q-values would be estimated from targets only. The
    // callers gate on a DB-level / pre-inference "decoys present" check, which does not catch this
    // post-inference target-only edge -- guard it here.
    if (decoy_proteins == 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Protein FDR requested but no decoy proteins remain after inference "
                      << "(marker '" << decoy_string << "'); skipping picked-protein FDR to avoid "
                      << "target-only q-values." << std::endl;
      removeUnreportedProteins(protein_ids, peptide_ids); // the proteins only inference-only PSMs support
      return false;
    }

    FalseDiscoveryRate fdr;
    if (!remove_decoys)
    {
      // By default the picked-protein FDR drops the decoy proteins from the run; keep them, with q-values, so that the
      // filter below applies to them as to the targets (as mergeFileResults() does). The targets' q-values are the same.
      Param fdr_param = fdr.getParameters();
      fdr_param.setValue("add_decoy_proteins", "true");
      fdr.setParameters(fdr_param);
    }
    fdr.applyPickedProteinFDR(protein_ids[0], decoy_string, decoy_is_prefix);
    IDFilter::filterHitsByScore(protein_ids, protein_fdr);
    // The indistinguishable groups have q-values of their own (the picked FDR of the groups, which counts a group
    // once, however many targets it holds): a group whose proteins pass can still be above the cutoff. Filter the
    // groups at the cutoff too, and keep only the proteins of the groups that pass, so that every reported protein
    // and group is within it.
    std::vector<ProteinIdentification::ProteinGroup>& groups = protein_ids[0].getIndistinguishableProteins();
    if (!groups.empty())
    {
      IDFilter::filterGroupsByScore(groups, protein_fdr, protein_ids[0].isHigherScoreBetter());
      IDFilter::removeUngroupedProteins(groups, protein_ids[0].getHits());
    }

    // Decoy removal is the finalization step. removeDecoyHits strips decoy PSMs; decoy
    // proteins then fall out as unreferenced. Repair indistinguishable-protein and
    // protein-group references and drop peptide evidence pointing at removed proteins,
    // else idXML storage fails on dangling references.
    if (remove_decoys)
    {
      IDFilter::removeDecoyHits(peptide_ids);
      IDFilter::removeEmptyIdentifications(peptide_ids);
    }
    // Proteins (and groups) without a reported PSM go: the decoys once their PSMs are removed, and the proteins that
    // only inference-only PSMs support.
    removeUnreportedProteins(protein_ids, peptide_ids);
    // A PSM all of whose proteins failed the protein FDR goes too (as in ProteomicsLFQ): it would be reported without
    // any protein (no protein accession, start, end or flanking residues), so every reported PSM maps to a reported
    // protein. The proteins of the PSMs that stay are not affected.
    IDFilter::removeDanglingProteinReferences(peptide_ids, protein_ids, /*remove_peptides_without_reference=*/ true);
    IDFilter::removeEmptyIdentifications(peptide_ids);

    OPENMS_LOG_INFO << "[ProSE] Protein inference + picked-protein FDR: "
                    << protein_ids[0].getHits().size() << " proteins at "
                    << protein_fdr * 100 << "% FDR." << std::endl;
    return true;
  }
} // namespace OpenMS
