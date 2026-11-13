// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>
#include "ProSEImpl.h"

#include <OpenMS/ANALYSIS/ID/CoIsolationDetector.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  // =====================================================================
  // Build a SearchContext (decoys + FragmentIndex) from a FASTA database.
  // Hoisted out of the in-memory search() body so callers can build the
  // index once and reuse it across many spectrum files.
  // =====================================================================
  Param ProSEAlgorithm::fragmentIndexParameters_(bool electron_ions) const
  {
    Param p = getParameters();
    // FragmentIndex has its own boolean 'decoys' flag {true,false}; ProSE's enum
    // value (auto/generate/ignore) would fail its validation. ProSE builds the
    // decoy-augmented database itself, so the FragmentIndex must never generate.
    p.remove("decoys");
    p.setValue("decoys", "false");
    p.remove("decoy_method"); // ProSE's own (the decoys of the database it builds)
    p.remove("decoy_parity_report");
    // ions:by_activation is resolved by the caller, which knows the spectra. The c and z+1 ions
    // for electron-activated spectra go into a set of their own, which the other spectra are not
    // matched against (see scoreSpectraAgainstIndex_).
    p.remove("ions:by_activation");
    p.setValue("ions:electron_ions", electron_ions ? "true" : "false");
    // annotation-only parameters of ProSE (see annotateIonPriors_)
    p.remove("annotate:self_trained_ion_priors");
    p.removeAll("annotate:ion_prior_");
    p.remove("peptide:protein_mapping"); // ProSE's own (FragmentIndex::getProteinOccurrences() serves it)
    p.remove("database:peff_variants"); // ProSE's own (FragmentIndex::setPeptideRegions() serves it)
    return p;
  }

  namespace
  {
    /**
      @brief DecoyHelper::findDecoyString(@p db, quiet = true), with the OpenMP threads for large databases.

      DecoyHelper::countDecoys() sums, over the proteins, the matches of the prefix and the suffix pattern and keeps,
      per match, the spelling of its last occurrence (within a protein the suffix after the prefix). Chunks of
      consecutive proteins counted on their own and joined in their order give the same. The decision is the one of
      findDecoyString() for these statistics; it does not depend on the order in which it looks at the matches,
      since at most one prefix and one suffix can reach 80% of all prefixes or suffixes.

      At most 8 threads: the step (about 60 ms serially for 50,000 entries) gains little beyond 4 threads, and a
      larger team costs CPU time (measured at 64 threads: 2-4 CPU-s per run for no shorter wall time).

      Most accessions cannot match a pattern, and searching for the suffix pattern tries every position: an accession
      is searched only if mayMatchDecoyPattern() says it can match.
    */
    DecoyHelper::Result findDecoyStringInParallel(const std::vector<FASTAFile::FASTAEntry>& db)
    {
      // Whether DecoyHelper's prefix (or suffix) pattern can match the lower-case accession: "^(decoy_*|dec_*|...)"
      // matches exactly the accessions that start with an affix, "(_decoy*|_dec*|...)$" exactly those that end with
      // '_', an affix without its last character and this character any number of times. Boost's ^ and $ also match
      // at the line breaks inside a text: such an accession is always searched.
      const auto mayMatchDecoyPattern = [](const std::string& lower, bool suffix)
      {
        if (lower.find_first_of("\n\r\f") != std::string::npos) { return true; }
        for (const std::string& affix : DecoyHelper::affixes)
        {
          if (!suffix)
          {
            if (StringUtils::hasPrefix(lower, affix)) { return true; }
            continue;
          }
          const std::string_view stem(affix.data(), affix.size() - 1);
          for (std::string_view head(lower); ; head.remove_suffix(1))
          {
            if (head.size() > stem.size() && head.substr(head.size() - stem.size()) == stem && head[head.size() - stem.size() - 1] == '_')
            {
              return true;
            }
            if (head.empty() || head.back() != affix.back()) { break; }
          }
        }
        return false;
      };
#ifdef _OPENMP
      constexpr int max_threads = 8;
      const int threads = omp_in_parallel() ? 1 : std::min(omp_get_max_threads(), max_threads);
#else
      const int threads = 1;
#endif
      if (threads < 2 || db.size() < 8192)
      {
        FASTAContainer<TFI_Vector> container(db);
        return DecoyHelper::findDecoyString(container, /*quiet=*/true);
      }
      struct Counts
      {
        std::map<std::string, std::pair<Size, Size>> count; ///< match -> occurrences as prefix, as suffix
        std::map<std::string, std::string> spelling;        ///< match -> spelling of its last occurrence
        Size prefixes = 0, suffixes = 0;
      };
      const SignedSize n = static_cast<SignedSize>(db.size());
      const SignedSize num_chunks = 4 * static_cast<SignedSize>(threads);
      std::vector<Counts> chunks(num_chunks);
#pragma omp parallel num_threads(threads)
      {
        const RegularExpression prefix_pattern(DecoyHelper::regexstr_prefix);
        const RegularExpression suffix_pattern(DecoyHelper::regexstr_suffix);
        std::string lower, match;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize c = 0; c < num_chunks; ++c)
        {
          Counts& chunk = chunks[c];
          for (SignedSize i = n * c / num_chunks; i < n * (c + 1) / num_chunks; ++i)
          {
            const std::string& identifier = db[i].identifier;
            lower = identifier;
            StringUtils::toLower(lower);
            if (mayMatchDecoyPattern(lower, false) && prefix_pattern.search(lower, &match))
            {
              ++chunk.prefixes;
              ++chunk.count[match].first;
              chunk.spelling[match] = StringUtils::prefix(identifier, match.length());
            }
            if (mayMatchDecoyPattern(lower, true) && suffix_pattern.search(lower, &match))
            {
              ++chunk.suffixes;
              ++chunk.count[match].second;
              chunk.spelling[match] = StringUtils::suffix(identifier, match.length());
            }
          }
        }
      }
      Counts all;
      for (const Counts& chunk : chunks) // in order: later chunks hold later proteins
      {
        all.prefixes += chunk.prefixes;
        all.suffixes += chunk.suffixes;
        for (const auto& [match, count] : chunk.count)
        {
          all.count[match].first += count.first;
          all.count[match].second += count.second;
        }
        for (const auto& [match, spelling] : chunk.spelling) all.spelling[match] = spelling;
      }
      // DecoyHelper::findDecoyString()'s decision
      const double proteins = static_cast<double>(db.size());
      if (static_cast<double>(all.prefixes + all.suffixes) < 0.4 * proteins || all.prefixes == all.suffixes)
      {
        return {false, "?", true};
      }
      for (const auto& [match, count] : all.count)
      {
        if (static_cast<double>(count.first) / static_cast<double>(all.prefixes) >= 0.8 && static_cast<double>(count.first) / proteins >= 0.4)
        {
          return {true, all.spelling[match], true};
        }
      }
      for (const auto& [match, count] : all.count)
      {
        if (static_cast<double>(count.second) / static_cast<double>(all.suffixes) >= 0.8 && static_cast<double>(count.second) / proteins >= 0.4)
        {
          return {true, all.spelling[match], false};
        }
      }
      return {false, "?", true};
    }
  }

  // static
  Size ProSEAlgorithm::removeDuplicateProteins(std::vector<FASTAFile::FASTAEntry>& db, const std::string& source)
  {
    // accession -> position of its first entry; the entries are compared by their sequences
    std::unordered_map<std::string_view, Size> first;
    first.reserve(db.size());
    std::vector<bool> repeated(db.size(), false);
    Size n_repeated = 0;
    for (Size i = 0; i < db.size(); ++i)
    {
      const auto [it, inserted] = first.emplace(db[i].identifier, i);
      if (inserted) { continue; }
      if (db[it->second].sequence != db[i].sequence)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, source,
          "The protein database '" + source + "' has two entries of the accession '" + db[i].identifier + "' with "
          "different sequences (entries " + StringUtils::toStr(it->second + 1) + " and " + StringUtils::toStr(i + 1)
          + "). An accession must name one protein: rename or remove one of them.");
      }
      repeated[i] = true;
      ++n_repeated;
    }
    if (n_repeated == 0) { return 0; }
    first.clear(); // (views of the identifiers, which are moved below)
    Size kept = 0;
    for (Size i = 0; i < db.size(); ++i)
    {
      if (repeated[i]) { continue; }
      if (kept != i) { db[kept] = std::move(db[i]); }
      ++kept;
    }
    db.resize(kept);
    OPENMS_LOG_WARN << "[ProSE] The protein database '" << source << "' repeats " << n_repeated << " entr"
                    << (n_repeated == 1 ? "y" : "ies") << " (same accession and sequence as an entry before, e.g. the "
                    << "file or the contaminants were added twice); searching each protein once." << std::endl;
    return n_repeated;
  }

  ProSEAlgorithm::DecoyStrategy_
  ProSEAlgorithm::resolveDecoyStrategy_(const std::vector<FASTAFile::FASTAEntry>& db) const
  {
    // Detect pre-existing decoys. First the common-marker heuristic
    // (DecoyHelper: "decoy", "rev", "xxx", ... as prefix or suffix), then a
    // literal fall-back to the configured decoy_prefix so custom markers
    // outside the vocabulary are still recognised.
    // quiet=true: a target-only database is a normal case here (auto/generate
    // then synthesise decoys), so suppress DecoyHelper's "unable to determine
    // decoy string" ERROR/WARN noise — we handle the negative result ourselves.
    const DecoyHelper::Result det = findDecoyStringInParallel(db);

    bool existing = det.success;
    std::string ext_string = det.success ? det.name : decoy_prefix_;
    bool ext_is_prefix = det.success ? det.is_prefix : true;

    Size n_decoys = 0; // the entries with the marker (as PeptideIndexing finds them: this spelling exactly)
    if (!existing && !decoy_prefix_.empty())
    {
      Size n_prefix = 0, n_suffix = 0;
      for (const FASTAFile::FASTAEntry& e : db)
      {
        n_prefix += StringUtils::hasPrefix(e.identifier, decoy_prefix_);
        n_suffix += StringUtils::hasSuffix(e.identifier, decoy_prefix_);
      }
      if (n_prefix + n_suffix > 0)
      {
        existing = true;
        ext_string = decoy_prefix_;
        ext_is_prefix = n_prefix > 0; // prefer prefix when both orientations occur
        n_decoys = ext_is_prefix ? n_prefix : n_suffix;
      }
    }
    else if (existing)
    {
      n_decoys = countDecoyEntries_(db, ext_string, ext_is_prefix);
    }

    DecoyStrategy_ s;
    switch (decoy_mode_)
    {
      case DecoyMode_::AUTO:
        if (existing)
        {
          // Reuse the decoys already in the database; never double-generate. Target-decoy FDR estimates the false
          // targets by the decoys, one decoy per target: with far fewer (e.g. a few entries of the decoy_prefix found
          // literally, which DecoyHelper's 40 % rule does not detect) or far more decoys, every q-value is off by
          // about that factor. Refuse such a database rather than reporting wrong FDRs.
          const Size n_targets = db.size() - n_decoys;
          if (static_cast<double>(n_decoys) < MIN_DECOYS_PER_TARGET * static_cast<double>(n_targets)
              || static_cast<double>(n_decoys) > MAX_DECOYS_PER_TARGET * static_cast<double>(n_targets))
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              paramName_("decoys") + "=auto: the protein database has " + StringUtils::toStr(n_decoys) + " decoy protein(s) (marker '"
              + ext_string + "', " + (ext_is_prefix ? "prefix" : "suffix") + ") for " + StringUtils::toStr(n_targets)
              + " target protein(s). Target-decoy FDR needs about one decoy per target (ProSE reuses decoys at "
              + StringUtils::toStr(MIN_DECOYS_PER_TARGET) + " to " + StringUtils::toStr(MAX_DECOYS_PER_TARGET)
              + " per target); with these its q-values would be wrong. Set '" + paramName_("decoys") + "' to 'generate' to replace them "
              "with decoys generated from the targets, or use a database with one decoy per target.");
          }
          s.generate = false; s.strip_existing = false; s.have_decoys = true;
          s.decoy_string = ext_string; s.is_prefix = ext_is_prefix;
        }
        else
        {
          // No decoys present: synthesise them so target-decoy FDR is possible.
          s.generate = true; s.strip_existing = false; s.have_decoys = true;
          s.decoy_string = decoy_prefix_; s.is_prefix = true;
        }
        break;

      case DecoyMode_::GENERATE:
        if (existing)
        {
          OPENMS_LOG_WARN << "[ProSE] decoys=generate: removing pre-existing decoys (marker '"
                          << ext_string << "') and regenerating from the target proteins.\n";
        }
        else
        {
          // Nothing detected to strip. If the database still carries some decoy-like accessions
          // with a non-standard or too-sparse marker (below DecoyHelper's threshold, and not the
          // configured decoy_prefix), they are treated as targets and reversed into decoys-of-decoys.
          // Warn so the user can set -Search:decoy_prefix to the real marker or pre-clean the FASTA.
          FASTAContainer<TFI_Vector> stats_container(db);
          const DecoyHelper::DecoyStatistics stats = DecoyHelper::countDecoys(stats_container);
          if (stats.all_prefix_occur + stats.all_suffix_occur > 0)
          {
            OPENMS_LOG_WARN << "[ProSE] decoys=generate: " << (stats.all_prefix_occur + stats.all_suffix_occur)
                            << " accession(s) carry a decoy-like marker but too few to auto-detect, so "
                            << "they are not stripped and will be treated as targets (risking "
                            << "decoys-of-decoys). Set -Search:decoy_prefix to the actual marker or "
                            << "pre-clean the database.\n";
          }
        }
        s.generate = true; s.strip_existing = existing; s.have_decoys = true;
        s.decoy_string = decoy_prefix_; s.is_prefix = true;
        s.strip_string = ext_string; s.strip_is_prefix = ext_is_prefix;
        break;

      case DecoyMode_::IGNORE:
        if (existing)
        {
          OPENMS_LOG_WARN << "[ProSE] decoys=ignore: removing pre-existing decoys (marker '"
                          << ext_string << "'); searching the target proteins only.\n";
        }
        s.generate = false; s.strip_existing = existing; s.have_decoys = false;
        s.decoy_string.clear(); s.is_prefix = true;
        s.strip_string = ext_string; s.strip_is_prefix = ext_is_prefix;
        break;
    }
    return s;
  }

  namespace
  {
    // The segments {start, length} of @p sequence: digestor.digestUnmodified(sequence, spans, 1, 0) (fully specific, no
    // missed cleavages), with @p rule without the regex (~1 ms per protein with the regex): its tokenize_() yields 0 and
    // every site 0 < p < size, and every segment between consecutive ones (and the last one up to the end) is a product.
    void digestSegments(const ProteaseDigestion& digestor, const SegmentRule* rule, const std::string& sequence,
                        std::vector<std::pair<Size, Size>>& spans)
    {
      if (rule == nullptr)
      {
        digestor.digestUnmodified(sequence, spans, 1, 0);
        return;
      }
      spans.clear();
      if (sequence.empty()) return;
      Size start = 0;
      for (Size p = 1; p < sequence.size(); ++p)
      {
        if (rule->after[static_cast<unsigned char>(sequence[p - 1])] && !rule->not_before[static_cast<unsigned char>(sequence[p])])
        {
          spans.emplace_back(start, p - start);
          start = p;
        }
      }
      spans.emplace_back(start, sequence.size() - start);
    }

    // decoy_method 'keep_termini': the residues that stay in place inside a segment besides its first and last one,
    // derived from the cleavage rule @p regex of the search enzyme: every residue letter it names (the rules of the
    // ProteaseDB name them in lookarounds and character classes, ambiguity codes included) and P. A rule that tells
    // residues apart only by the letters it names then finds the target's sites in the decoy: inside a segment every
    // residue stays or trades places with one the rule treats alike. Trypsin and Trypsin/P: K, R, X and P, so the
    // decoys of DecoyDatabase's -Decoy:non_shuffle_pattern K,P,R except where an X precedes a P (Trypsin only; an X ends
    // a segment otherwise). Chymotrypsin: F, Y, W, L, J, X and P (K, R and P for every enzyme, as before, moved an F,
    // Y, W or L before P, which mostly made a new site).
    std::string keepTerminiFixedResidues(const std::string& regex)
    {
      std::string residues = "P";
      for (const char c : regex)
      {
        if (c >= 'A' && c <= 'Z' && residues.find(c) == std::string::npos) residues += c;
      }
      return residues;
    }

    // decoy_method 'keep_termini': OpenMS DecoyDatabase's construction (method reverse, Decoy:keepPeptideNTerm and
    // keepPeptideCTerm true; Decoy:non_shuffle_pattern @p fixed_residues of keepTerminiFixedResidues(), its default is
    // empty) with @p digestor (the search enzyme, no missed cleavages, full specificity): every segment is
    // MRMDecoy::reversePeptideLight(segment, keepN, keepC, fixed_residues). Residues not in any segment (none for a
    // digest that tiles the protein) stay in place.
    std::string keepTerminiDecoy(const std::string& protein, const ProteaseDigestion& digestor, const SegmentRule* rule,
                                 const std::string& fixed_residues, std::vector<std::pair<Size, Size>>& spans)
    {
      digestSegments(digestor, rule, protein, spans);
      std::string decoy = protein;
      for (const auto& [start, length] : spans)
      {
        decoy.replace(start, length, MRMDecoy::reversePeptideLight(protein.substr(start, length), {}, true, true, fixed_residues).first);
      }
      return decoy;
    }

    // FNV-1a 64 of the generated decoy sequences, each followed by '\n', in target order (decoy_parity_report)
    void fnv1a(UInt64& h, const std::string& s)
    {
      for (const unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
      h ^= static_cast<unsigned char>('\n'); h *= 1099511628211ULL;
    }
  }

  // Build the searched database according to @p strategy: optionally strip
  // pre-existing decoys, optionally generate fresh decoys from the targets.
  // Produces FASTA entries only — does not build a FragmentIndex.
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const DecoyStrategy_& strategy) const
  {
    return buildDecoyAugmentedDB_(std::vector<FASTAFile::FASTAEntry>(fasta_db), strategy);
  }

  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      std::vector<FASTAFile::FASTAEntry>&& fasta_db,
      const DecoyStrategy_& strategy) const
  {
    // The entries are filtered in place instead of being copied one by one: the same entries
    // in the same order.
    std::vector<FASTAFile::FASTAEntry> db = std::move(fasta_db);
    // For in-memory databases; one read from a file had its repeats removed by loadFASTA() already (and its conflicts
    // reported with the file name), so this finds none.
    removeDuplicateProteins(db, "the protein database");

    // decoys=auto reusing pre-existing decoys logs nothing otherwise; surface the auto-detected
    // marker so a rare DecoyHelper misdetection (a target DB whose accessions start with a decoy
    // affix) is diagnosable. Single emission: buildDecoyAugmentedDB_ runs once per search.
    const bool reuses_decoys = decoy_mode_ == DecoyMode_::AUTO && !strategy.generate && strategy.have_decoys;
    if (reuses_decoys)
    {
      OPENMS_LOG_INFO << "[ProSE] decoys=auto: reusing existing decoys detected in the database "
                      << "(marker '" << strategy.decoy_string << "', "
                      << (strategy.is_prefix ? "prefix" : "suffix") << ")." << std::endl;
    }
    // Initial-Met clipping adds N-terminal peptides of every protein that starts with M. Generated decoys keep the
    // initial Met (below); supplied decoys often do not (a reversed protein ends with it), and then the clipped
    // peptides enlarge the target space only. Count both kinds while the entries are filtered.
    const bool check_met_symmetry = reuses_decoys && param_.getValue("peptide:clip_nterm_methionine").toBool();
    Size targets = 0, targets_with_met = 0, decoys = 0, decoys_with_met = 0;

    // 1. Keep targets, dropping pre-existing decoys when requested. A stop codon ('*') that ends
    //    a sequence, as in databases translated from genomes (e.g. SGD), is not a residue: remove
    //    it so the C-terminal peptide stays searchable and decoys are built from the protein
    //    alone. FragmentIndex skips peptides that contain a stop codon inside the sequence.
    //    An entry left without residues has nothing to search, and decoy generation needs
    //    residues: drop it. Lowercase residues are read as uppercase ones, as the FASTA format
    //    defines them (NCBI); they are neither indexed nor reversed otherwise.
    const auto existing_decoy = [&strategy](const FASTAFile::FASTAEntry& e)
    {
      return strategy.strip_existing &&
        (strategy.strip_is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.strip_string)
                                  : StringUtils::hasSuffix(e.identifier, strategy.strip_string));
    };
    // (the lowercase check reads every residue: in parallel, 10-20 ms on one thread for a human target-decoy database)
    Size lowercase_entries = 0;
#pragma omp parallel for schedule(dynamic, 256) reduction(+: lowercase_entries) if (db.size() >= 1024)
    for (SignedSize i = 0; i < static_cast<SignedSize>(db.size()); ++i)
    {
      std::string& sequence = db[i].sequence;
      if (!existing_decoy(db[i]) && std::any_of(sequence.begin(), sequence.end(), [](char c) { return c >= 'a' && c <= 'z'; }))
      {
        StringUtils::toUpper(sequence);
        ++lowercase_entries;
      }
    }
    auto kept = db.begin();
    for (FASTAFile::FASTAEntry& e : db)
    {
      const bool is_existing_decoy = existing_decoy(e);
      if (!is_existing_decoy)
      {
        std::string& sequence = e.sequence;
        while (!sequence.empty() && sequence.back() == '*') { sequence.pop_back(); }
        if (!sequence.empty())
        {
          if (check_met_symmetry)
          {
            const bool is_decoy = strategy.is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.decoy_string)
                                                     : StringUtils::hasSuffix(e.identifier, strategy.decoy_string);
            (is_decoy ? decoys : targets) += 1;
            (is_decoy ? decoys_with_met : targets_with_met) += (sequence.size() > 1 && sequence[0] == 'M') ? 1 : 0;
          }
          if (&*kept != &e) { *kept = std::move(e); }
          ++kept;
        }
      }
    }
    db.erase(kept, db.end());
    if (lowercase_entries > 0)
    {
      OPENMS_LOG_INFO << "[ProSE] " << lowercase_entries << " database entries hold lowercase residues; they are searched "
                      << "as uppercase ones." << std::endl;
    }
    // Warn when the share of decoys that start with M is below half the share of targets that do.
    if (check_met_symmetry && targets_with_met > 0 && 2 * decoys_with_met * targets < targets_with_met * decoys)
    {
      OPENMS_LOG_WARN << "[ProSE] peptide:clip_nterm_methionine: " << targets_with_met << " of " << targets
                      << " target proteins but only " << decoys_with_met << " of " << decoys << " decoy proteins "
                      << "in the database start with M. Removing the initial Met adds N-terminal peptides to the "
                      << "targets that have no decoy counterpart, which makes target-decoy FDR estimates slightly "
                      << "optimistic. Use '-Search:decoys generate' (generated decoys keep the initial Met) or "
                      << "'-Search:peptide:clip_nterm_methionine false' for a symmetric search space." << std::endl;
    }
    // The options that rely on symmetric decoys, with decoys taken from the database: say so (non-default options only)
    if (reuses_decoys && decoy_keep_termini_)
    {
      OPENMS_LOG_WARN << "[ProSE] decoy_method keep_termini has no effect: the decoys of the database are used "
                      << "(decoys 'auto'). Use '-Search:decoys generate' to replace them." << std::endl;
    }
    // The effective N-terminal rescoring options of the search whose database this builds: those requested, unless
    // the decoys cannot support them (below). Re-derived at every database build, so a later search of the same object
    // starts from the requested values again. (One search at a time per object, as for the other search state.)
    terminal_features_ = param_.getValue("rescoring:terminal_features").toString() == "true";
    stratum_q_values_ = param_.getValue("rescoring:stratum_q_values").toString() == "true";
    nterm_options_disabled_ = false;
    if (check_met_symmetry && (terminal_features_ || stratum_q_values_)
        && targets_with_met > 0 && 2 * decoys_with_met * targets < targets_with_met * decoys)
    {
      // With these decoys the features would tell targets from decoys (the mechanism of OpenMS issue #10726), so the
      // requested options are switched off for this search rather than only warned about. The check is the share of
      // decoys that start with M against the targets' (the asymmetry of whole-protein reversal); it does not prove
      // that other decoys populate every N-terminal class like the targets.
      OPENMS_LOG_WARN << "[ProSE] rescoring:terminal_features / rescoring:stratum_q_values with these decoys: the "
                      << "Met-clipped (and protein N-terminal) peptides have (almost) no decoys, so the features would tell "
                      << "targets from decoys and the stratum q-values would be too low. The requested option(s) are switched "
                      << "off for this search (recorded in the output's search parameters); use '-Search:decoys generate "
                      << "-Search:decoy_method keep_termini' for decoys that support them." << std::endl;
      terminal_features_ = false;
      stratum_q_values_ = false;
      nterm_options_disabled_ = true;
    }

    // 1b. The windows of the PEFF sequence variants (database:peff_variants), as further targets: the decoys below are
    //     generated for them as for the proteins
    std::shared_ptr<ProSEPEFF> peff;
    if (peff_variants_)
    {
      ProSEPEFF::Settings peff_settings;
      peff_settings.enzyme = enzyme_;
      peff_settings.specificity = peptide_enzyme_specificity_;
      peff_settings.peptide_max_length = peptide_max_size_;
      // per cleavage segment (keep_termini, or reversePeptides() with a specific digestion) or the whole entry reversed
      peff_settings.segment_decoys = peptide_enzyme_specificity_ != EnzymaticDigestion::SPEC_NONE;
      peff_settings.clip_nterm_methionine = param_.getValue("peptide:clip_nterm_methionine").toBool();
      peff = std::make_shared<ProSEPEFF>();
      const Size n_windows = peff->addVariantWindows(db, peff_settings);
      if (n_windows > 0 && reuses_decoys)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          paramName_("database:peff_variants") + " needs generated decoys: the decoys of the database have none for the "
          "sequence variants. Use '" + paramName_("decoys") + " generate'.");
      }
      OPENMS_LOG_INFO << "[ProSE] database:peff_variants: " << n_windows << " sequence variants (\\VariantSimple) added "
                      << "to the database." << std::endl;
      if (n_windows == 0) { peff.reset(); }
    }

    // 2. Generate decoys by reversing the (remaining) target proteins.
    if (strategy.generate)
    {
      // decoy_string is the prefix prepended to generated accessions; an empty prefix would
      // produce decoys with the same accession as their targets, silently breaking FDR.
      if (strategy.decoy_string.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          paramName_("decoy_prefix") + " must be non-empty to generate decoys (" + paramName_("decoys") + "=auto without existing "
          "decoys, or " + paramName_("decoys") + "=generate).");
      }
      DecoyGenerator decoy_generator;
      const size_t old_size = db.size();
      db.reserve(old_size * 2);
      // decoy_method 'keep_termini' (fully specific digestion only: the candidates are unions of whole segments; otherwise
      // 'reverse')
      const bool keep_termini = decoy_keep_termini_ && peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_FULL
                                && enzyme_ != EnzymaticDigestion::UnspecificCleavage && enzyme_ != EnzymaticDigestion::NoCleavage;
      if (decoy_keep_termini_ && !keep_termini)
      {
        OPENMS_LOG_WARN << "[ProSE] decoy_method keep_termini needs a fully specific digestion with a cleaving enzyme; "
                        << "the decoys are generated with 'reverse'." << std::endl;
      }
      ProteaseDigestion segments;
      Size site_mismatches = 0;
      UInt64 checksum = 14695981039346656037ULL;
      SegmentRule segment_rule;
      bool simple_segments = false;
      std::string fixed_residues;
      if (keep_termini)
      {
        segments.setEnzyme(enzyme_);
        segments.setMissedCleavages(0);
        segments.setSpecificity(EnzymaticDigestion::SPEC_FULL);
        const std::string& regex = ProteaseDB::getInstance()->getEnzyme(enzyme_)->getRegEx();
        simple_segments = parseSegmentRule(regex, segment_rule);
        fixed_residues = keepTerminiFixedResidues(regex);
      }
      const SegmentRule* const segment_rule_ptr = simple_segments ? &segment_rule : nullptr;
      // The decoy sequences in parallel (each from its target alone; serially 1.2-1.6 s per run for a human database, as
      // AASequence parsing and the digestion per protein dominate), then appended in target order as before.
      const bool clip_nterm_methionine = param_.getValue("peptide:clip_nterm_methionine").toBool();
      std::vector<std::string> decoy_sequences(old_size);
      std::vector<char> site_mismatch(old_size, 0);
      std::vector<std::exception_ptr> decoy_errors(old_size); // an exception must not leave the parallel region
      // Only sequences of plain residue letters in parallel: AASequence::fromString() of a sequence with modification
      // syntax may register modifications (shared state) and warn; then serially, stopping at the first error, as before.
      bool plain_letters = true;
#pragma omp parallel for schedule(dynamic, 256) reduction(&&: plain_letters) if (old_size >= 256)
      for (SignedSize si = 0; si < static_cast<SignedSize>(old_size); ++si)
      {
        const std::string& sequence = db[static_cast<size_t>(si)].sequence;
        plain_letters = plain_letters && std::all_of(sequence.begin(), sequence.end(), [](char c) { return c >= 'A' && c <= 'Z'; });
      }
      const bool decoys_in_parallel = plain_letters && old_size >= 256;
      bool serial_stop = false; // (written only without the parallel region)
#pragma omp parallel if (decoys_in_parallel)
      {
        std::vector<std::pair<Size, Size>> target_spans, decoy_spans;
        // monotonic: without the parallel region (serial fallback) the proteins are processed in order
#if _OPENMP >= 201511
#pragma omp for schedule(monotonic: dynamic, 64)
#else
#pragma omp for schedule(dynamic, 64) // OpenMP 2.0 (MSVC): no modifier, and dynamic is monotonic
#endif
        for (SignedSize si = 0; si < static_cast<SignedSize>(old_size); ++si)
        {
          const size_t i = static_cast<size_t>(si);
          if (serial_stop) { continue; }
          try
          {
            const std::string& sequence = db[i].sequence;
            if (keep_termini)
            {
              decoy_sequences[i] = keepTerminiDecoy(sequence, segments, segment_rule_ptr, fixed_residues, target_spans);
              // the guarantee holds only if the decoy has the target's cleavage sites (every rule that tells residues
              // apart only by the letters it names, as those of the ProteaseDB; checked)
              digestSegments(segments, segment_rule_ptr, decoy_sequences[i], decoy_spans);
              site_mismatch[i] = decoy_spans != target_spans;
              continue;
            }
            // Keep an initial Met on generated decoys as well: otherwise clipping
            // expands the target N-terminal search space without its decoy counterpart.
            const bool preserve_met = clip_nterm_methionine && sequence.size() > 1 && sequence[0] == 'M';
            const std::string target = preserve_met ? sequence.substr(1) : sequence;
            std::string decoy = peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_NONE
              ? decoy_generator.reverseProtein(AASequence::fromString(target)).toString()
              : decoy_generator.reversePeptides(AASequence::fromString(target), enzyme_).toString();
            if (preserve_met) { decoy.insert(decoy.begin(), 'M'); }
            decoy_sequences[i] = std::move(decoy);
          }
          catch (...)
          {
            decoy_errors[i] = std::current_exception();
            if (!decoys_in_parallel) { serial_stop = true; } // (serially: the first error stops, as before)
          }
        }
      }
      for (const std::exception_ptr& error : decoy_errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = db[i];
        e.sequence = std::move(decoy_sequences[i]);
        site_mismatches += site_mismatch[i];
        if (keep_termini || decoy_parity_report_) fnv1a(checksum, e.sequence);
        e.identifier = strategy.decoy_string + e.identifier;  // decoy_string is the prefix to add
        db.push_back(std::move(e));
      }
      std::vector<std::string>().swap(decoy_sequences);
      if (keep_termini)
      {
        if (site_mismatches > 0)
        {
          OPENMS_LOG_WARN << "[ProSE] decoy_method keep_termini: " << site_mismatches << " of " << old_size << " decoy "
                          << "proteins do not have the cleavage sites of their target (enzyme '" << enzyme_ << "'); their "
                          << "peptides are not paired with the target's." << std::endl;
        }
      }
      if (decoy_parity_report_)
      {
        OPENMS_LOG_INFO << "[ProSE] decoy parity: " << old_size << " decoys generated (" << (keep_termini ? "keep_termini" : "reverse")
                        << "), cleavage-site mismatches " << (keep_termini ? StringUtils::toStr(site_mismatches) : std::string("not checked"))
                        << ", FNV-1a 64 of the decoy sequences in target order " << std::hex << checksum << std::dec << std::endl;
      }
      Math::RandomShuffler shuffler(42);  // fixed seed for reproducible decoy ordering across runs/files
      shuffler.portable_random_shuffle(db.begin(), db.end());

      // Under the default decoys=auto, decoys are generated whenever the database lacks them —
      // including no-ProSE-FDR runs, because the decoys are typically consumed by a downstream or
      // external FDR step (FalseDiscoveryRate / Percolator) or a merged run. This roughly doubles
      // index build + search; for a genuinely target-only search, use decoys=ignore.
      if (decoy_mode_ == DecoyMode_::AUTO && fdr_psm_ == 0.0 && fdr_protein_ == 0.0)
      {
        OPENMS_LOG_INFO << "[ProSE] decoys=auto generated " << old_size << " decoys (the database had "
                        << "none) for downstream target-decoy FDR. Use '-Search:decoys ignore' for a "
                        << "target-only search if you do not need FDR." << std::endl;
      }
      if (peff) { peff->addDecoys(db, strategy.decoy_string); }
    }
    peff_ = std::move(peff);
    return db;
  }

  void ProSEAlgorithm::setPeptideRegions_(FragmentIndex& index, const std::shared_ptr<const ProSEPEFF>& peff)
  {
    index.setPeptideRegions(peff ? peff->peptideRegions() : nullptr);
  }

  // =====================================================================
  // Strided calibration sample. Used by the chunked search paths so that
  // calibration always sees a representative, suitably-sized subset of the
  // full DB — independent of database_chunk_size_. Prevents the first-chunk
  // starvation mode where small chunk_size values (e.g. 500) produce a
  // calibration pool below calibration_min_psms_ and silently disable
  // calibration. See issue #9182.
  // =====================================================================
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildCalibrationSample_(
      const std::vector<FASTAFile::FASTAEntry>& full_db) const
  {
    // Calibration FI size is bounded by the user's memory budget — which is
    // exactly what they declared via database_chunk_size_. This keeps peak RSS
    // during calibration identical to peak RSS during the main search
    // (one chunk FI at a time), regardless of digestion mode.
    //
    // Critical for immunopeptidomics (non-specific digestion, MHC-I/II variant
    // databases): a naive fixed-size 5000-protein cal sample can generate
    // 50+ GB of fragment index because each protein produces thousands of
    // candidate peptides — defeating the memory savings chunking is designed
    // to deliver. Tying cal_target to chunk_size keeps the promise that
    // "user's memory budget = one chunk's FI."
    //
    // Strided (not contiguous / random) to give the calibration pool even
    // coverage of the full DB.
    //
    // Residual gap vs non-chunked calibration:
    //   - Trypsin / moderate chunk_size (e.g. 5000): cal sample is adequate;
    //     yield is within ~2% of non-chunked.
    //   - Small chunk_size (e.g. 100 proteins for MHC-II variants): cal pool
    //     can fall below calibration_min_psms_ → runCalibrationPass_ returns
    //     success=false, calibration silently no-ops. The main search uses
    //     user-configured tolerances.
    // Follow-up issue: pool calibration PSMs across the first K chunks of the
    // main search loop, so calibration quality becomes independent of
    // chunk_size.
    const Size cal_target = std::min(
        std::max<Size>(database_chunk_size_, Size(1)),
        full_db.size());
    if (cal_target >= full_db.size()) return full_db;

    const Size stride = std::max<Size>(1, full_db.size() / cal_target);
    std::vector<FASTAFile::FASTAEntry> cal_db;
    cal_db.reserve(cal_target);
    for (Size i = 0; i < full_db.size() && cal_db.size() < cal_target; i += stride)
    {
      cal_db.push_back(full_db[i]);
    }
    return cal_db;
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db) const
  {
    return prepareContext(fasta_db, false);
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db, bool electron_ions) const
  {
    return prepareContext_(std::vector<FASTAFile::FASTAEntry>(fasta_db), electron_ions);
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext_(
      std::vector<FASTAFile::FASTAEntry>&& fasta_db, bool electron_ions,
      const std::function<const PeakMap*(Size)>& searched_spectra) const
  {
    SearchContext ctx;

    startProgress(0, 1, "Generate decoys...");
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    ctx.db = buildDecoyAugmentedDB_(std::move(fasta_db), strategy);
    ctx.decoy_string = strategy.decoy_string;
    ctx.decoy_is_prefix = strategy.is_prefix;
    ctx.have_decoys = strategy.have_decoys;
    ctx.peff = peff_;
    endProgress();

    // build fragment index
    startProgress(0, 1, "Building fragment index...");
    buildFragmentIndex_(ctx, fragmentIndexParameters_(electron_ions), searched_spectra);
    ctx.electron_ions = electron_ions;
    endProgress();
    reportDecoyParity_(ctx);

    return ctx;
  }

  double ProSEAlgorithm::measuredError_(double offset, double sigma)
  {
    return std::fabs(offset) + 4.0 * sigma;
  }

  double ProSEAlgorithm::calibrationPassWindow_(double configured, bool ppm)
  {
    // At least 20 ppm, which measures an error of up to about 20 ppm without censoring it. In Da at least twice the
    // configured window, up to 0.5 Da: a fixed 0.5 Da would search high-resolution MS2 configured at 0.02 Da 25 times
    // wider, with more chance matches among the PSMs the error is measured from; ion-trap windows (0.25 Da or more)
    // get 0.5 Da as before. The configured window if wider.
    return std::max(configured, ppm ? 20.0 : std::min(0.5, 2.0 * configured));
  }

  void ProSEAlgorithm::buildFragmentIndex_(SearchContext& ctx, const Param& params,
                                           const std::function<const PeakMap*(Size)>& searched_spectra) const
  {
    setPeptideRegions_(ctx.fragment_index, ctx.peff);
    if (!searched_spectra || !calibration_enabled_)
    {
      ctx.fragment_index.setParameters(params);
      ctx.fragment_index.build(ctx.db, searched_spectra);
      return;
    }
    // With calibration, the peptides the calibration pass can reach: it searches with the widest precursor windows, and
    // the calibrated windows of the main search never exceed the configured ones. With calibration:recalibrate offset
    // the main search's windows lie around the corrected precursor m/z, moved by a median error of the pass (within its
    // window on each side): the configured window is added. The tolerances are query-time parameters of the index
    // (FragmentIndex::updateMembers_()), so the configured ones are set after the build.
    const bool ppm = precursor_mass_tolerance_unit_ == "ppm";
    const bool offset = calibration_recalibrate_ == "offset";
    Param pass = params;
    pass.setValue("precursor:mass_tolerance_lower", calibrationPassWindow_(precursor_mass_tolerance_lower_, ppm)
                                                    + (offset ? precursor_mass_tolerance_lower_ : 0.0));
    pass.setValue("precursor:mass_tolerance_upper", calibrationPassWindow_(precursor_mass_tolerance_upper_, ppm)
                                                    + (offset ? precursor_mass_tolerance_upper_ : 0.0));
    ctx.fragment_index.setParameters(pass);
    ctx.fragment_index.build(ctx.db, searched_spectra);
    ctx.fragment_index.setParameters(params);
  }

  void ProSEAlgorithm::reportDecoyParity_(const SearchContext& ctx) const
  {
    if (!decoy_parity_report_) return;
    const std::vector<FASTAFile::FASTAEntry>& db = ctx.db;
    const bool clip = param_.getValue("peptide:clip_nterm_methionine").toBool();
    std::vector<char> decoy(db.size(), 0);
    if (!ctx.decoy_string.empty())
    {
      for (Size i = 0; i < db.size(); ++i)
      {
        decoy[i] = ctx.decoy_is_prefix ? StringUtils::hasPrefix(db[i].identifier, ctx.decoy_string)
                                       : StringUtils::hasSuffix(db[i].identifier, ctx.decoy_string);
      }
    }
    // subclass of an occurrence: 0 internal, 1 protein N-term starting with M, 2 protein N-term otherwise, 3 Met-clipped
    const auto subclass = [&](UInt32 p, Size start) -> int
    {
      const std::string& q = db[p].sequence;
      if (start == 0) return (!q.empty() && q[0] == 'M') ? 1 : 2;
      if (clip && start == 1 && q[0] == 'M') return 3;
      return 0;
    };
    Size entries[2][4] = {}, occurrences[2][4] = {}, acetyl[2][4] = {};
    // protein N-terminal acetylation (a variable modification such as Acetyl (Protein N-term)) of a kept entry
    std::vector<signed char> acetylated;
    const auto is_acetylated = [&](Size idx) -> bool
    {
      const FragmentIndex::Peptide& pe = ctx.fragment_index.getPeptides()[idx];
      if (pe.mod_bitmask_ == 0) return false;
      if (acetylated.empty()) acetylated.assign(ctx.fragment_index.getPeptides().size(), -1);
      if (acetylated[idx] < 0)
      {
        const AASequence seq = ctx.fragment_index.reconstructModifiedSequence(pe, db);
        acetylated[idx] = (seq.hasNTerminalModification()
                           && seq.getNTerminalModificationName().find("Acetyl") != std::string::npos) ? 1 : 0;
      }
      return acetylated[idx] == 1;
    };
    std::set<std::string_view> distinct[2][4];
    std::set<std::string_view> any_label[2];
    const std::vector<FragmentIndex::Peptide>& peptides = ctx.fragment_index.getPeptides();
    const auto count = [&](UInt32 p, Size start, Size length, bool unmodified)
    {
      const int l = decoy[p] ? 1 : 0, c = subclass(p, start);
      ++occurrences[l][c];
      if (unmodified)
      {
        const std::string_view seq(db[p].sequence.data() + start, length);
        distinct[l][c].insert(seq);
        any_label[l].insert(seq);
      }
    };
    for (Size idx = 0; idx < peptides.size(); ++idx)
    {
      const FragmentIndex::Peptide& pe = peptides[idx];
      const int l = decoy[pe.protein_idx] ? 1 : 0, c = subclass(pe.protein_idx, pe.sequence_.first);
      ++entries[l][c];
      count(pe.protein_idx, pe.sequence_.first, pe.sequence_.second, pe.mod_bitmask_ == 0);
      if (c != 0 && is_acetylated(idx)) ++acetyl[l][c];
    }
    for (const FragmentIndex::RemovedOccurrence& ro : ctx.fragment_index.getRemovedOccurrences())
    {
      const FragmentIndex::Peptide& kept = peptides[ro.peptide_idx];
      count(ro.protein_idx, ro.start, kept.sequence_.second, kept.mod_bitmask_ == 0);
      const int c = subclass(ro.protein_idx, ro.start);
      if (c != 0 && is_acetylated(ro.peptide_idx)) ++acetyl[decoy[ro.protein_idx] ? 1 : 0][c];
    }
    Size shared = 0;
    for (const std::string_view& seq : any_label[0]) shared += any_label[1].count(seq);
    static const char* names[4] = {"internal", "nterm_M", "nterm_other", "metclip"};
    OPENMS_LOG_INFO << "[ProSE] decoy parity (fragment index; target / decoy): subclass, peptidoform entries kept, "
                    << "protein occurrences of peptidoforms, distinct unmodified sequences" << std::endl;
    for (int c = 0; c < 4; ++c)
    {
      OPENMS_LOG_INFO << "[ProSE] decoy parity " << names[c] << ": entries " << entries[0][c] << " / " << entries[1][c]
                      << ", occurrences " << occurrences[0][c] << " / " << occurrences[1][c]
                      << ", sequences " << distinct[0][c].size() << " / " << distinct[1][c].size()
                      << ", N-terminally acetylated occurrences " << acetyl[0][c] << " / " << acetyl[1][c]
                      << (occurrences[0][c] == occurrences[1][c] && acetyl[0][c] == acetyl[1][c] ? " (occurrences equal)" : " (occurrences DIFFER)") << std::endl;
    }
    OPENMS_LOG_INFO << "[ProSE] decoy parity: distinct unmodified sequences target " << any_label[0].size() << ", decoy "
                    << any_label[1].size() << ", in both " << shared << std::endl;
  }

  void ProSEAlgorithm::scoreSpectraAgainstIndex_(
      const PeakMap& spectra,
      FragmentIndex& fi,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const SpectrumGenerators_& generators,
      double effective_fragment_tol,
      bool fragment_mass_tolerance_unit_ppm,
      bool open_search_mode,
      std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
      std::vector<CandidatePoolStats_>& pool_stats,
      const std::string& progress_label,
      const PeakMap* query_spectra) const
  {
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }
    if (query_spectra != nullptr && query_spectra->size() != spectra.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, query_spectra->size(),
        "Candidate-retrieval spectra must hold one entry per scoring spectrum.");
    }
    startProgress(0, spectra.size(), progress_label);
    // incremented by every thread and read by the master thread for the progress display: atomic for both
    std::atomic<Size> count_spectra{0};
    std::atomic<int> scoring_threads{0}; // the team of the first region, for the log (see the end)
    const double proton_mass_u = Constants::PROTON_MASS_U;
    const Size keep = std::max(searchTopHits_(), Size(2)); // keep ≥2 for delta score
    const bool deduplicate_chunks = param_.getValue("peptide:deduplicate").toBool() && database_chunk_size_ > 0;

    // An exception on a worker thread (a modification AASequence rejects, a scorer's parameter check, ...) must not
    // leave the parallel region, which would terminate the process. The one of the SMALLEST failing scan index is
    // rethrown after the loop: what the scan-order loop on one thread reports. After a failure only the spectra before
    // it are still scored (to find an earlier failure); every other spectrum is skipped. (Spectra after the failing one
    // that came earlier in the scoring order may have been scored: the results are discarded with the exception.)
    std::exception_ptr scoring_error;
    const SignedSize n_spectra = static_cast<SignedSize>(spectra.size());
    std::atomic<SignedSize> first_failed{n_spectra};
    // The spectra are scored in the order of their precursor mass (scoringOrder_()): neighbouring iterations then probe
    // the same precursor-window slices of the fragment index, which stay in the caches. Every spectrum is scored on its
    // own and its results go to its scan index: the same results in any order.
    const std::vector<SignedSize> order = scoringOrder_(spectra, query_spectra);
    // While a background read lends threads of this thread, in slices, between which this thread takes them back once
    // the read has ended (see BackgroundSpectra). The slices run over positions in the scoring order.
    const SignedSize slice = readLendsThreads() ? std::max<SignedSize>(1, (n_spectra + 15) / 16) : std::max<SignedSize>(1, n_spectra);
    for (SignedSize slice_begin = 0; slice_begin < n_spectra; slice_begin += slice)
    {
    takeBackThreadsOfEndedRead();
    const SignedSize slice_end = std::min(n_spectra, slice_begin + slice);
    // Runs of neighbouring masses per core (dynamic chunks of up to 64 spectra), but at least 8 chunks per thread of the
    // slice so that small inputs and slices stay balanced. In peptide-mass mode (low_memory) a spectrum costs milliseconds
    // (every candidate in its precursor windows is realised and its fragments generated), not microseconds, and the cost
    // grows with the mass, so the last chunks of the ascending-mass order left the other threads idle (at 64 threads the
    // last 2 of 10 s ran on 37 -> 12 cores): one spectrum per chunk then, the dispatch costs nothing against such spectra.
#ifdef _OPENMP
    const int chunk = fi.isPeptideMassMode() ? 1 : static_cast<int>(std::clamp<SignedSize>(
      (slice_end - slice_begin) / (SignedSize(8) * std::max(1, omp_get_max_threads())), 1, 64));
#else
    [[maybe_unused]] const int chunk = 1; // (only the OpenMP directive uses it)
#endif
#pragma omp parallel for schedule(dynamic, chunk) default(none) shared(annotated_hits, pool_stats, query_spectra, count_spectra, scoring_threads, fi, generators, db, fragment_mass_tolerance_unit_ppm, spectra, open_search_mode, proton_mass_u, effective_fragment_tol, keep, deduplicate_chunks, scoring_error, first_failed, order, slice_begin, slice_end, chunk)
    for (SignedSize position = slice_begin; position < slice_end; ++position)
    {
      const SignedSize scan_index = order[position];
      if (scan_index >= first_failed.load(std::memory_order_relaxed)) continue;
      try
      {
      const Size spectra_done = count_spectra.fetch_add(1, std::memory_order_relaxed) + 1;
      IF_MASTERTHREAD
      {
        setProgress(spectra_done);
#ifdef _OPENMP
        if (scoring_threads.load(std::memory_order_relaxed) == 0) { scoring_threads.store(omp_get_num_threads(), std::memory_order_relaxed); }
#endif
      }

      const MSSpectrum& exp_spectrum = spectra[scan_index];
      const TheoreticalSpectrumGenerator& spectrum_generator = generators.forSpectrum(exp_spectrum);
      ScoringScratch_& scratch = scoringScratch_();
      // The candidates of this spectrum: the query appends to the hits, so the hits of the thread's previous spectrum
      // are cleared first (clear() keeps the capacity)
      FragmentIndex::SpectrumMatchesTopN& top_sms = scratch.top_sms;
      top_sms.clear();
      // ions:by_activation: only electron-activated spectra are matched against the c and z+1 ions,
      // so that these ions do not change which candidates the other spectra keep
      const MSSpectrum& query = query_spectra != nullptr ? (*query_spectra)[scan_index] : exp_spectrum;
      fi.querySpectrum(query, db, top_sms, generators.electronIons(exp_spectrum));

      const bool peptide_mass_mode = fi.isPeptideMassMode();

      // Reused across the candidates of this spectrum and, per thread, across spectra. Avoids per-candidate heap
      // churn of a fresh PeakSpectrum + its DataArrays (TSG's add_metainfo fills
      // StringDataArrays with ion names — these are a notable allocation hot spot
      // when the candidate count per spectrum is in the tens/hundreds). Both generator calls below start by
      // clearing the spectrum, so what the previous candidate or spectrum left in it does not matter.
      PeakSpectrum& theo_spectrum = scratch.theo_spectrum;
      std::vector<uint8_t>& theo_ion_series = scratch.theo_ion_series;

      for (const auto& sms : top_sms.hits_)
      {
        const FragmentIndex::Peptide& sms_pep = fi.getPeptides()[sms.peptide_idx_];

        // Peptide-mass mode: the candidate is the prefix of the mother, with the variable modifications, that the index
        // matched to the precursor
        AASequence mod_candidate = peptide_mass_mode
          ? fi.reconstructRealizedSubSequence(sms_pep, db, sms.realized_length_, sms.subset_bitmask_)
          : fi.reconstructModifiedSequence(sms_pep, db);

        // Index construction already removes repeated occurrences within a
        // chunk. Across chunks, skip the same hypothesis before updating either
        // scores or pool statistics. Charge and isotope hypotheses stay distinct.
        if (deduplicate_chunks)
        {
          std::string key = mod_candidate.toString();
          key += "\t" + std::to_string(sms.precursor_charge_) + "\t" + std::to_string(sms.isotope_error_);
          if (! pool_stats[scan_index].seen_candidates.insert(std::move(key)).second) { continue; }
        }

        HyperScore::PSMDetail detail;
        double score;
        // The peaks of getSpectrum() with the ion series per peak instead of the ion names (no strings, no data arrays):
        // the same scores
        if (spectrum_generator.getIonSeriesSpectrum(theo_spectrum, theo_ion_series, mod_candidate, 1, scoringMaxCharge_(sms.precursor_charge_)))
        {
          score = mass_accuracy_score_
            ? HyperScore::computeMassAccuracy(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                              theo_ion_series, mass_error_sd_ppm_, detail)
            : HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                            theo_ion_series, detail);
        }
        else
        {
          // Clear peaks + data arrays (ion names / charges) before refilling for the
          // next candidate; getSpectrum appends to whatever is there. clear(false) keeps the capacity of the
          // peak vector (clear(true) would shrink it); the data arrays are emptied either way.
          theo_spectrum.clear(false);
          theo_spectrum.getPrecursors().clear(); // (clear(false) keeps the meta data; getSpectrum() appends a precursor per call)
          spectrum_generator.getSpectrum(theo_spectrum, mod_candidate, 1, scoringMaxCharge_(sms.precursor_charge_));
          // Note: TSG emits sorted output when add_metainfo=true (see the
          // sortByPositionPresorted() call at the tail of getSpectrum_); the extra
          // sortByPosition() pass here was a redundant O(N) scan per candidate.
          score = mass_accuracy_score_
            ? HyperScore::computeMassAccuracy(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                              mass_error_sd_ppm_, detail)
            : HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum, detail);
        }

        // Summarise the candidate before it can be dropped below or pruned at the
        // end of the loop: the pool-derived PSM features describe the whole search
        // space, so a candidate that scored 0 (no fragment matched) still counts
        // and still belongs in the null distribution the z-score is measured against.
        // Each scan_index is owned by exactly one thread, so this needs no guard.
        pool_stats[scan_index].add(score);

        if (score == 0) continue;

        AnnotatedHit_ ah;
        ah.sequence = std::move(mod_candidate);
        ah.score = score;
        // Account for the additional charge hypotheses in the ion-count fractions.
        double seq_length = static_cast<double>(ah.sequence.size()) * scoringMaxCharge_(sms.precursor_charge_);
        ah.prefix_fraction = static_cast<float>(detail.matched_prefix_ions / seq_length);
        ah.suffix_fraction = static_cast<float>(detail.matched_suffix_ions / seq_length);
        ah.mean_error = static_cast<float>(detail.mean_error);
        ah.matched_prefix_ions = static_cast<uint16_t>(detail.matched_prefix_ions);
        ah.matched_suffix_ions = static_cast<uint16_t>(detail.matched_suffix_ions);
        ah.isotope_error = sms.isotope_error_;
        ah.applied_charge = sms.precursor_charge_;
        ah.index_peptide = static_cast<UInt32>(sms.peptide_idx_);
        ah.delta_mass = 0.0;
        if (open_search_mode)
        {
          double theo_mh_plus = ah.sequence.getMZ(1);
          double exp_mz = exp_spectrum.getPrecursors()[0].getMZ();
          double exp_mh_plus = exp_mz * sms.precursor_charge_ - ((sms.precursor_charge_ - 1) * proton_mass_u);
          ah.delta_mass = exp_mh_plus - theo_mh_plus;
        }

        annotated_hits[scan_index].push_back(std::move(ah));
      }

      // Prune to top-N per spectrum to bound memory: up to
      // scoring:max_candidates_per_spectrum hits (each owning a heap AASequence)
      // would otherwise stay resident until postProcessHits_ truncates them.
      // Correct because all of this spectrum's candidates are appended above and
      // scores are independent of the chunk a candidate came from: a hit outside
      // the top-N here can never re-enter it.
      auto& hits = annotated_hits[scan_index];
#ifdef _OPENMP
      pool_stats[scan_index].scoring_thread = omp_get_thread_num();
#endif
      if (hits.size() > keep)
      {
        std::partial_sort(hits.begin(), hits.begin() + keep, hits.end(), AnnotatedHit_::hasBetterScore);
        hits.resize(keep);
        hits.shrink_to_fit();
      }
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_scoring_error)
        {
          if (scan_index < first_failed.load(std::memory_order_relaxed))
          {
            scoring_error = std::current_exception();
            first_failed.store(scan_index, std::memory_order_relaxed);
          }
        }
      }
    }
    }
    OPENMS_LOG_DEBUG << "[ProSE] Scored " << n_spectra << " spectra with " << std::max(1, scoring_threads.load(std::memory_order_relaxed)) << " OpenMP threads" << std::endl;
    endProgress();
    if (scoring_error) std::rethrow_exception(scoring_error);
  }

  std::vector<SignedSize> ProSEAlgorithm::scoringOrder_(const PeakMap& spectra, const PeakMap* query_spectra)
  {
    // Key: the neutral precursor mass of the first precursor of the candidate-retrieval spectrum (an unknown or
    // non-positive charge counts as 1; only the order depends on it). Spectra without a finite key (no precursor, a
    // non-finite m/z) go to the end in scan order. (key, scan index) is a total order: any sort gives the same order.
    const SignedSize n = static_cast<SignedSize>(spectra.size());
    std::vector<SignedSize> order(n);
    if (n < 2)
    {
      for (SignedSize i = 0; i < n; ++i) order[i] = i;
      return order;
    }
    std::vector<std::pair<double, SignedSize>> keyed(n);
    std::vector<unsigned char> valid(n);
#pragma omp parallel for schedule(static) if (n >= 16384)
    for (SignedSize i = 0; i < n; ++i)
    {
      const MSSpectrum& s = query_spectra != nullptr ? (*query_spectra)[i] : spectra[i];
      double key = std::numeric_limits<double>::quiet_NaN();
      if (!s.getPrecursors().empty())
      {
        const Precursor& precursor = s.getPrecursors()[0];
        const int charge = precursor.getCharge() > 0 ? precursor.getCharge() : 1;
        key = (precursor.getMZ() - Constants::PROTON_MASS_U) * charge;
      }
      valid[i] = std::isfinite(key) ? 1 : 0;
      keyed[i] = {valid[i] ? key : 0.0, i};
    }
    // the spectra with a key first (in scan order), then the others (in scan order)
    SignedSize n_valid = 0;
    for (SignedSize i = 0; i < n; ++i) { if (valid[i]) keyed[n_valid++] = keyed[i]; }
    SignedSize next_invalid = n_valid;
    for (SignedSize i = 0; i < n; ++i) { if (!valid[i]) order[next_invalid++] = i; }
    const auto by_key = [](const std::pair<double, SignedSize>& a, const std::pair<double, SignedSize>& b)
    {
      return a.first < b.first || (!(b.first < a.first) && a.second < b.second);
    };
#ifdef PROSE_GNU_PARALLEL
    if (n_valid >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::sort(keyed.begin(), keyed.begin() + n_valid, by_key); }
    else
#endif
    { std::sort(keyed.begin(), keyed.begin() + n_valid, by_key); }
    for (SignedSize p = 0; p < n_valid; ++p) order[p] = keyed[p].second;
    return order;
  }

  namespace
  {
    // chimeric:enabled: the peaks of @p spectrum that the theoretical ions of @p sequence (fragment charges up to
    // @p max_charge) claim: every peak within the fragment tolerance of an ion. Merged into the ascending peak indices
    // @p claimed.
    void claimPeaks(const MSSpectrum& spectrum, const AASequence& sequence, const TheoreticalSpectrumGenerator& generator,
                    int max_charge, double fragment_tol, bool fragment_tol_ppm, PeakSpectrum& theo, std::vector<UInt32>& claimed)
    {
      theo.clear(true);
      generator.getSpectrum(theo, sequence, 1, max_charge);
      for (const Peak1D& ion : theo)
      {
        const double tol = fragment_tol_ppm ? ion.getMZ() * fragment_tol * 1e-6 : fragment_tol;
        for (auto it = spectrum.MZBegin(ion.getMZ() - tol); it != spectrum.end() && it->getMZ() <= ion.getMZ() + tol; ++it)
        {
          claimed.push_back(static_cast<UInt32>(it - spectrum.begin()));
        }
      }
      std::sort(claimed.begin(), claimed.end());
      claimed.erase(std::unique(claimed.begin(), claimed.end()), claimed.end());
    }
  }

  MSSpectrum ProSEAlgorithm::residualSpectrum_(const MSSpectrum& spectrum, const std::vector<UInt32>& removed, const Precursor& precursor)
  {
    std::vector<Size> kept;
    kept.reserve(spectrum.size());
    auto next_removed = removed.begin();
    for (Size p = 0; p < spectrum.size(); ++p)
    {
      if (next_removed != removed.end() && *next_removed == p)
      {
        ++next_removed;
        continue;
      }
      kept.push_back(p);
    }
    MSSpectrum residual = spectrum;
    residual.select(kept);
    residual.removeMetaValue(CoIsolationDetector::META_MZ);
    residual.removeMetaValue(CoIsolationDetector::META_CHARGE);
    residual.removeMetaValue(OBSERVED_COISOLATED_MZ);
    residual.setPrecursors({precursor});
    return residual;
  }

  void ProSEAlgorithm::searchCoIsolated_(const PeakMap& spectra,
                                         const std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
                                         FragmentIndex& fi,
                                         const std::vector<FASTAFile::FASTAEntry>& db,
                                         const SpectrumGenerators_& generators,
                                         double fragment_tol,
                                         bool fragment_tol_ppm,
                                         bool open_search_mode,
                                         std::vector<CoIsolatedPSM_>& psms,
                                         std::vector<std::vector<AnnotatedHit_>>& hits,
                                         std::vector<CandidatePoolStats_>& pool_stats) const
  {
    psms.clear();
    hits.clear();
    pool_stats.clear();
    const Size max_n = static_cast<Size>(static_cast<int>(param_.getValue("chimeric:max_coisolated")));
    const Size max_hits = std::min(report_top_hits_, static_cast<Size>(static_cast<int>(param_.getValue("chimeric:top_hits"))));

    // The spectra with co-isolated precursors and at least one hit: their best hit, the peaks claimed so far and
    // the winners found so far (indices into the round results below)
    struct Scan
    {
      Size spectrum = 0;
      DoubleList mzs;
      DoubleList observed_mzs; // the mzs as observed (OBSERVED_COISOLATED_MZ), if calibration:recalibrate corrected them
      IntList charges;
      const AnnotatedHit_* primary = nullptr;
      std::vector<UInt32> claimed;
      std::vector<std::pair<Size, Size>> winners; // (round, position in that round's results)
    };
    std::vector<Scan> scans;
    for (Size i = 0; i < spectra.size() && i < annotated_hits.size(); ++i)
    {
      const MSSpectrum& spectrum = spectra[i];
      if (annotated_hits[i].empty() || spectrum.empty() || spectrum.getPrecursors().empty()
          || !spectrum.metaValueExists(CoIsolationDetector::META_MZ) || !spectrum.metaValueExists(CoIsolationDetector::META_CHARGE))
      {
        continue;
      }
      Scan scan;
      scan.spectrum = i;
      scan.mzs = spectrum.getMetaValue(CoIsolationDetector::META_MZ);
      if (spectrum.metaValueExists(OBSERVED_COISOLATED_MZ)) { scan.observed_mzs = spectrum.getMetaValue(OBSERVED_COISOLATED_MZ); }
      scan.charges = spectrum.getMetaValue(CoIsolationDetector::META_CHARGE);
      const Size n = std::min({scan.mzs.size(), scan.charges.size(), max_n});
      scan.mzs.resize(n);
      scan.charges.resize(n);
      if (n == 0) continue;
      scan.primary = &*std::min_element(annotated_hits[i].begin(), annotated_hits[i].end(), AnnotatedHit_::hasBetterScore);
      scans.push_back(std::move(scan));
    }
    if (scans.empty()) return;

#pragma omp parallel for schedule(dynamic, 64) default(none) shared(scans, spectra, generators, fragment_tol, fragment_tol_ppm)
    for (SignedSize s = 0; s < static_cast<SignedSize>(scans.size()); ++s)
    {
      Scan& scan = scans[s];
      const MSSpectrum& spectrum = spectra[scan.spectrum];
      PeakSpectrum theo;
      claimPeaks(spectrum, scan.primary->sequence, generators.forSpectrum(spectrum), scoringMaxCharge_(scan.primary->applied_charge),
                 fragment_tol, fragment_tol_ppm, theo, scan.claimed);
    }

    // One round per co-isolated precursor rank; a round's residual spectra are scored in batches
    constexpr Size batch_size = 8192;
    struct Result
    {
      Size scan = 0;
      CoIsolatedPSM_ psm;
      std::vector<AnnotatedHit_> hits; // the winner first
      CandidatePoolStats_ stats;
      bool emitted = false;
    };
    std::vector<std::vector<Result>> rounds;
    Size n_searched = 0, n_dropped = 0;
    PeakSpectrum theo;
    for (Size k = 0; k < max_n; ++k)
    {
      std::vector<Size> in_round;
      for (Size s = 0; s < scans.size(); ++s)
      {
        if (k < scans[s].mzs.size() && scans[s].charges[k] > 0) in_round.push_back(s);
      }
      if (in_round.empty()) break;
      std::vector<Result>& results = rounds.emplace_back();
      for (Size begin = 0; begin < in_round.size(); begin += batch_size)
      {
        const Size end = std::min(begin + batch_size, in_round.size());
        PeakMap residuals;
        std::vector<CoIsolatedPSM_> batch_psms;
        for (Size b = begin; b < end; ++b)
        {
          const Scan& scan = scans[in_round[b]];
          const MSSpectrum& spectrum = spectra[scan.spectrum];
          CoIsolatedPSM_ psm;
          psm.parent = scan.spectrum;
          psm.precursor = spectrum.getPrecursors()[0];
          psm.precursor.setMZ(scan.mzs[k]);
          // the m/z as observed, which postProcessHits_() reports (not the parent's; none if it was not corrected)
          if (k < scan.observed_mzs.size()) { psm.precursor.setMetaValue(OBSERVED_PRECURSOR_MZ, scan.observed_mzs[k]); }
          else { psm.precursor.removeMetaValue(OBSERVED_PRECURSOR_MZ); }
          psm.precursor.setCharge(scan.charges[k]);
          psm.removed = scan.claimed;
          residuals.addSpectrum(residualSpectrum_(spectrum, psm.removed, psm.precursor));
          batch_psms.push_back(std::move(psm));
        }
        std::vector<std::vector<AnnotatedHit_>> batch_hits(residuals.size());
        std::vector<CandidatePoolStats_> batch_stats(residuals.size());
        scoreSpectraAgainstIndex_(residuals, fi, db, generators, fragment_tol, fragment_tol_ppm, open_search_mode,
                                  batch_hits, batch_stats, "Scoring co-isolated precursors...");
        n_searched += residuals.size();
        residuals.clear(true);
        for (Size b = begin; b < end; ++b)
        {
          Scan& scan = scans[in_round[b]];
          // the hits that are not the spectrum's best peptide, the best (the winner) first
          std::vector<AnnotatedHit_>& candidates = batch_hits[b - begin];
          candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                          [&scan](const AnnotatedHit_& h) { return h.sequence == scan.primary->sequence; }),
                           candidates.end());
          if (candidates.empty()) continue;
          std::sort(candidates.begin(), candidates.end(), AnnotatedHit_::hasBetterScore);
          if (candidates.size() > max_hits) candidates.resize(max_hits);
          const auto winner = candidates.begin();
          const MSSpectrum& spectrum = spectra[scan.spectrum];
          claimPeaks(spectrum, winner->sequence, generators.forSpectrum(spectrum), scoringMaxCharge_(winner->applied_charge),
                     fragment_tol, fragment_tol_ppm, theo, scan.claimed);
          // dropped if the spectrum already reports this peptide at this charge: its distinct precursor m/z would
          // make it a precursor of its own in Percolator's target-decoy competition and count it twice
          auto same = [&winner](const AnnotatedHit_& h) { return h.applied_charge == winner->applied_charge && h.sequence == winner->sequence; };
          bool repeated = std::any_of(annotated_hits[scan.spectrum].begin(), annotated_hits[scan.spectrum].end(), same);
          for (const auto& [round, position] : scan.winners)
          {
            repeated = repeated || same(rounds[round][position].hits.front());
          }
          Result result;
          result.scan = in_round[b];
          result.psm = std::move(batch_psms[b - begin]);
          result.hits = std::move(candidates);
          result.stats = std::move(batch_stats[b - begin]);
          result.emitted = !repeated;
          if (repeated) ++n_dropped;
          scan.winners.emplace_back(k, results.size());
          results.push_back(std::move(result));
        }
      }
    }

    // In spectrum order, within a spectrum in detection order
    for (const Scan& scan : scans)
    {
      for (const auto& [round, position] : scan.winners)
      {
        Result& result = rounds[round][position];
        if (!result.emitted) continue;
        psms.push_back(std::move(result.psm));
        hits.push_back(std::move(result.hits));
        pool_stats.push_back(std::move(result.stats));
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Chimeric search: " << n_searched << " co-isolated precursors searched, " << psms.size()
                    << " with a PSM (" << n_dropped << " more repeat a PSM of their spectrum and are dropped)." << std::endl;
  }

  // =====================================================================
  // Chunked search implementation. Takes a pre-built decoy-augmented DB
  // (from buildDecoyAugmentedDB_) and splits it into chunks for scoring.
  // Called from the single-file search() wrapper and the multi-file wrapper.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::searchChunked_(
      PeakMap& spectra,
      std::vector<FASTAFile::FASTAEntry>& full_db,
      const DecoyStrategy_& strategy,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    const Size n_chunks = (full_db.size() + database_chunk_size_ - 1) / database_chunk_size_;
    if (chimeric_enabled_)
    {
      OPENMS_LOG_WARN << "[ProSE] chimeric:enabled is not supported with database:chunk_size > 0 and is ignored." << std::endl;
    }
    OPENMS_LOG_INFO << "[ProSE] Database chunking enabled: " << full_db.size()
                    << " proteins (incl. decoys), chunk_size=" << database_chunk_size_
                    << " → " << n_chunks << " chunks." << std::endl;

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
    bool open_search_mode = isOpenSearchMode_();
    FragmentIonLikelihoodModel::PeakLists ion_evidence; // annotate:self_trained_ion_priors
    PeakMap evidence_spectra, query_spectra;
    PeakMap* query_ptr = query_raw_spectrum_ ? &query_spectra : nullptr;
    PeakMap* evidence_ptr = param_.getValue("annotate:local_fragment_evidence").toBool() ? &evidence_spectra : nullptr;
    preprocessSearchedSpectra_(spectra, "", ion_evidence, evidence_ptr, query_ptr);

    // ions:by_activation: the chunk indices hold c and z+1 ions if electron-activated spectra are searched
    const Size n_electron_activated = countElectronActivated_(spectra);
    const bool electron_ions = n_electron_activated > 0;
    if (electron_ions)
    {
      OPENMS_LOG_INFO << "[ProSE] " << n_electron_activated << " of " << spectra.size()
                      << " spectra are electron-activated (ETD, ECD, EThcD or ETciD) and are also scored"
                      << " with c and z+1 ions." << std::endl;
    }

    // Effective tolerances — may be narrowed by calibration below. Kept asymmetric
    // (lower / upper separately) so the FragmentIndex query can exploit the full
    // calibrated window rather than the looser max() collapse.
    double effective_precursor_tol_lower = precursor_mass_tolerance_lower_;
    double effective_precursor_tol_upper = precursor_mass_tolerance_upper_;
    double effective_fragment_tol = fragment_mass_tolerance_;

    // Save originals so we can restore algo-level members after search (they may be
    // mutated below for the OSMA mod-match tolerance computation).
    const double orig_prec_tol_lower = precursor_mass_tolerance_lower_;
    const double orig_prec_tol_upper = precursor_mass_tolerance_upper_;
    bool calibration_applied = false;
    // Restore algo-level tolerance members on every return path, and when an exception leaves this function (e.g. a
    // worker exception that the scoring loop rethrows). Calibration may mutate them below; leaving them mutated would
    // poison subsequent searches that reuse this ProSEAlgorithm instance (e.g. the multi-file wrapper, or a caller
    // that catches the exception).
    auto restore_tolerances = [&]()
    {
      if (calibration_applied)
      {
        calibration_applied = false;
        precursor_mass_tolerance_lower_ = orig_prec_tol_lower;
        precursor_mass_tolerance_upper_ = orig_prec_tol_upper;
      }
    };
    const OnScopeExit restore_on_exit(restore_tolerances);

    // Optional calibration on a strided sample of the full DB. The sample size is
    // bounded (see buildCalibrationSample_) so calibration memory stays O(chunk)
    // regardless of database_chunk_size_. Strided rather than first-N so small
    // chunk_size values don't starve the calibration pool (#9182).
    CalibrationResult_ cal; // (its recalibration, if any, is recorded in the search parameters below)
    if (calibration_enabled_ && !open_search_mode)
    {
      std::vector<FASTAFile::FASTAEntry> cal_db = buildCalibrationSample_(full_db);
      FragmentIndex cal_fi;
      cal_fi.setParameters(fragmentIndexParameters_(electron_ions));
      setPeptideRegions_(cal_fi, peff_);
      cal_fi.build(cal_db);

      cal = runCalibrationPass_(spectra, cal_fi, cal_db, query_ptr, strategy.have_decoys ? strategy.decoy_string : std::string(),
                                strategy.is_prefix);
      if (cal.success)
      {
        applyRecalibration_(cal, spectra, query_ptr, evidence_ptr);
        effective_fragment_tol = cal.fragment_tolerance;
        if (!cal.extreme_bias)
        {
          effective_precursor_tol_lower = cal.cal_lower;
          effective_precursor_tol_upper = cal.cal_upper;
          precursor_mass_tolerance_lower_ = cal.cal_lower;
          precursor_mass_tolerance_upper_ = cal.cal_upper;
          calibration_applied = true;
          OPENMS_LOG_INFO << "[ProSE] Calibration (chunked, strided sample): shift=" << cal.precursor_shift
                          << " " << precursor_mass_tolerance_unit_
                          << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                          << " fragment=" << cal.fragment_tolerance << std::endl;
        }
        else
        {
          OPENMS_LOG_WARN << "[ProSE] Calibration: extreme bias, precursor calibration discarded. "
                          << "Fragment calibration applied." << std::endl;
        }
      }
    }

    // 3. Prepare spectrum generators (once).
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // 4. Allocate per-spectrum hit accumulator (persists across chunks).
    // not reserved: most spectra keep few hits or none, and report:top_hits per spectrum ahead of time is a lot
    vector<vector<AnnotatedHit_>> annotated_hits(spectra.size());
    // Accumulates across chunks: every chunk contributes its own candidates for the
    // same spectrum, and the pool features must see all of them.
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    // 5. Chunk loop.
    const Size chunk_size = database_chunk_size_;
    Size chunk_idx = 0;
    for (Size start = 0; start < full_db.size(); start += chunk_size)
    {
      ++chunk_idx;
      const Size end = std::min(start + chunk_size, full_db.size());

      OPENMS_LOG_INFO << "[ProSE] Chunk " << chunk_idx << ": proteins "
                      << start << "–" << (end - 1) << " (" << (end - start) << " proteins)" << std::endl;

      // Build FragmentIndex for this chunk only.
      std::vector<FASTAFile::FASTAEntry> chunk_db(full_db.begin() + start, full_db.begin() + end);
      FragmentIndex chunk_fi;
      {
        Param fi_params = fragmentIndexParameters_(electron_ions);
        // Apply calibrated tolerances (if calibration succeeded above). Asymmetric
        // lower/upper preserved — collapsing to max() would re-open the tight side
        // of the calibrated window and admit spurious decoy candidates.
        fi_params.setValue("fragment:mass_tolerance", effective_fragment_tol);
        fi_params.setValue("precursor:mass_tolerance_lower", effective_precursor_tol_lower);
        fi_params.setValue("precursor:mass_tolerance_upper", effective_precursor_tol_upper);
        chunk_fi.setParameters(fi_params);
      }
      setPeptideRegions_(chunk_fi, peff_);
      chunk_fi.build(chunk_db);

      // Score all spectra against this chunk's index.
      scoreSpectraAgainstIndex_(spectra, chunk_fi, chunk_db, generators,
                                effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                                open_search_mode, annotated_hits, pool_stats,
                                "Scoring chunk " + StringUtils::toStr(chunk_idx) + "...", query_ptr);

      // Prune to top-N per spectrum after each chunk to bound memory growth.
      // Without this, K chunks × T candidates per spectrum could accumulate
      // K*T hits before postProcessHits_ truncates — problematic for 100M+
      // peptide databases with many chunks. The pruning is correct because
      // scores are independent across chunks: a hit that fails to place in the
      // current top-N cannot improve when more chunks are added.
      const Size keep = std::max(searchTopHits_(), Size(2)); // keep ≥2 for delta score
#pragma omp parallel for default(none) shared(annotated_hits, keep)
      for (SignedSize si = 0; si < (SignedSize)annotated_hits.size(); ++si)
      {
        if (annotated_hits[si].size() > keep)
        {
          std::partial_sort(annotated_hits[si].begin(),
                            annotated_hits[si].begin() + keep,
                            annotated_hits[si].end(),
                            AnnotatedHit_::hasBetterScore);
          annotated_hits[si].resize(keep);
        }
      }
    } // end chunk loop
    query_spectra.clear(true);
    for (auto& stats : pool_stats)
    {
      std::unordered_set<std::string>().swap(stats.seen_candidates);
    }

    // 6. Post-process merged hits (sort, annotate, PeptideIndexing, FDR).
    //    This runs once on ALL hits accumulated across all chunks.
    //    Set mod-match tolerance for open-search modification analysis downstream.
    last_mod_match_tolerance_used_ = computeModMatchTolerance_();

    startProgress(0, 1, "Post-processing PSMs...");
    postProcessHits_(spectra,
      annotated_hits,
      pool_stats,
      protein_ids,
      peptide_ids,
      searchTopHits_(),
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      std::max(precursor_mass_tolerance_lower_, precursor_mass_tolerance_upper_),
      effective_fragment_tol,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      "", // no database filename for in-memory search
      evidence_ptr,
      nullptr,
      false,
      cal.appliedFragmentRecalibration()
      );
    recordRecalibration_(cal, protein_ids);
    // Release the hits and the additional peak list before PeptideIndexing.
    freeOnOwnerThreads(annotated_hits, [&pool_stats](SignedSize i) { return pool_stats[i].scoring_thread; }); // (where scored)
    freeInParallel(pool_stats);
    freeInParallel(evidence_spectra);
    endProgress();

    annotatePeptDeepFeatures_(spectra, protein_ids, peptide_ids);

    // 7. PeptideIndexing against the FULL database (not per-chunk).
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    // Use the effective decoy marker/position that was actually searched. A
    // non-empty string (decoy_prefix_ when target-only) avoids PeptideIndexing's
    // own auto-detection; missing_decoy_action=silent keeps it quiet when none match.
    param_pi.setValue("decoy_string", strategy.decoy_string.empty() ? decoy_prefix_ : strategy.decoy_string);
    param_pi.setValue("decoy_string_position", strategy.is_prefix ? "prefix" : "suffix");
    param_pi.setValue("enzyme:name", enzyme_);
    param_pi.setValue("enzyme:specificity",
                      EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
    param_pi.setValue("missing_decoy_action", "silent");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(full_db, protein_ids, peptide_ids);
    if (peff_ && indexer_exit == PeptideIndexing::ExitCodes::EXECUTION_OK) { peff_->mapToParents(protein_ids, peptide_ids); }

    if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
    {
      restore_tolerances();
      if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
        return ExitCodes::INPUT_FILE_EMPTY;
      else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
        return ExitCodes::UNEXPECTED_RESULT;
      else
        return ExitCodes::UNKNOWN_ERROR;
    }

    // 8. PSM-level FDR only. Protein-level FDR is the caller's responsibility
    //    (matching the non-chunked search(spectra, ctx, ...) semantics — that
    //    method also does PSM FDR only; the file-based search() and the single-file
    //    searchWithModificationAnalysis overloads apply protein FDR post-call).
    const bool has_decoys = strategy.have_decoys;

    // Optional per-run ion priors need the target/decoy labels and the native scores: after
    // PeptideIndexing, before FDR overwrites the scores.
    if (self_trained_ion_priors_)
    {
      annotateIonPriors_(spectra, ion_evidence, protein_ids, peptide_ids);
      ion_evidence.clear();
    }

    // Pre-FDR stats (target/decoy counts + HyperScore distribution).
    capturePreFdrStats_(peptide_ids, last_run_stats_);

    // PSM-level FDR, after Percolator rescoring with 'rescore' (see applyPsmFDR_)
    last_percolator_ = applyPsmFDR_(protein_ids, peptide_ids, has_decoys, last_run_stats_, nullptr, &inference_only_peptide_ids_);

    restore_tolerances();

    collectRunStatistics_(spectra, protein_ids, peptide_ids, last_run_stats_);

    return ExitCodes::EXECUTION_OK;
  }
} // namespace OpenMS
