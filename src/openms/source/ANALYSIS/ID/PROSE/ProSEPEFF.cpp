// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEPEFF.h>

#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/PEFFFile.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>

namespace OpenMS
{
  const std::string& ProSEPEFF::variantMetaValue()
  {
    static const std::string name = "PEFF_variants";
    return name;
  }

  std::string ProSEPEFF::kindName(const VariantKind kind)
  {
    switch (kind)
    {
      case VariantKind::SIMPLE: return "simple";
      case VariantKind::STOP_GAIN: return "stop_gain";
      case VariantKind::COMPLEX: return "complex";
    }
    return "simple";
  }

  std::vector<ProSEPEFF::EquivalentSubstitution> ProSEPEFF::referenceEquivalentSubstitutions(const StringList& modifications)
  {
    std::vector<EquivalentSubstitution> out;
    if (modifications.empty()) { return out; }
    static const std::string residues = "ACDEFGHIKLMNPQRSTVWY";
    const ResidueDB* rdb = ResidueDB::getInstance();
    const ModifiedPeptideGenerator::MapToResidueType mods = ModifiedPeptideGenerator::getModifications(modifications);
    for (const auto& [mod, modified] : mods.val)
    {
      const char origin = mod->getOrigin();
      const Residue* variant = origin == '\0' ? nullptr : rdb->getResidue(static_cast<unsigned char>(origin));
      if (variant == nullptr) { continue; } // (terminal modifications without a residue)
      const EmpiricalFormula formula = variant->getFormula(Residue::Internal) + mod->getDiffFormula();
      for (const char r : residues)
      {
        if (r == origin) { continue; }
        const Residue* reference = rdb->getResidue(static_cast<unsigned char>(r));
        if (reference != nullptr && reference->getFormula(Residue::Internal) == formula)
        {
          out.push_back({origin, r, mod->getFullId()});
        }
      }
    }
    std::sort(out.begin(), out.end(), [](const EquivalentSubstitution& a, const EquivalentSubstitution& b)
              { return std::tie(a.variant, a.reference, a.modification) < std::tie(b.variant, b.reference, b.modification); });
    return out;
  }

  const std::string& ProSEPEFF::windowSeparator()
  {
    static const std::string separator = "_VAR_";
    return separator;
  }

  namespace
  {
    // Whether @p c is the letter of a residue of known mass (as the fragment index can index it)
    bool isIndexableResidue(const char c)
    {
      if (c < 'A' || c > 'Z') { return false; }
      const Residue* residue = ResidueDB::getInstance()->getResidue(static_cast<unsigned char>(c));
      return residue != nullptr && !residue->getFormula().isEmpty();
    }

    // The starts of the cleavage segments of @p sequence (ascending, starting with 0): the positions an enzyme cuts before
    void cleavageSites(const ProteaseDigestion& segments, const std::string_view sequence, std::vector<Size>& sites)
    {
      sites.clear();
      std::vector<std::pair<Size, Size>> spans;
      segments.digestUnmodified(sequence, spans, 1, 0);
      sites.reserve(spans.size() + 1);
      sites.push_back(0);
      for (const std::pair<Size, Size>& span : spans)
      {
        if (span.first > 0) { sites.push_back(span.first); }
      }
      std::sort(sites.begin(), sites.end());
      sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
    }

    // The largest site <= @p position (0 if none)
    Size siteAtOrBefore(const std::vector<Size>& sites, const Size position)
    {
      const auto it = std::upper_bound(sites.begin(), sites.end(), position);
      return it == sites.begin() ? 0 : *std::prev(it);
    }

    // The smallest site >= @p position (@p length if none)
    Size siteAtOrAfter(const std::vector<Size>& sites, const Size position, const Size length)
    {
      const auto it = std::lower_bound(sites.begin(), sites.end(), position);
      return it == sites.end() ? length : *it;
    }

    // The smallest site > @p position (@p length if none)
    Size siteAfter(const std::vector<Size>& sites, const Size position, const Size length)
    {
      const auto it = std::upper_bound(sites.begin(), sites.end(), position);
      return it == sites.end() ? length : *it;
    }
  }

  Size ProSEPEFF::addVariantWindows(std::vector<FASTAFile::FASTAEntry>& db, const Settings& settings)
  {
    const bool whole_proteins = settings.peptide_max_length == 0 || settings.enzyme == EnzymaticDigestion::NoCleavage;
    // decoys per cleavage segment, in place (DecoyGenerator::reversePeptides(), keep_termini); 'unspecific cleavage' has
    // no segments and 'no cleavage' a single one: the decoy of an entry is then the entry reversed
    const bool segment_decoys = settings.segment_decoys && settings.enzyme != EnzymaticDigestion::UnspecificCleavage
                                && settings.enzyme != EnzymaticDigestion::NoCleavage;
    // residues beyond the changed ones: a peptide (at most peptide_max_length residues) that overlaps them ends inside
    const Size margin = settings.peptide_max_length + 1;

    ProteaseDigestion segments; // cleavage segments (no missed cleavage)
    if (segment_decoys)
    {
      segments.setEnzyme(settings.enzyme);
      segments.setMissedCleavages(0);
      segments.setSpecificity(EnzymaticDigestion::SPEC_FULL);
    }

    std::unordered_set<std::string> identifiers;
    identifiers.reserve(db.size());
    for (const FASTAFile::FASTAEntry& entry : db) { identifiers.insert(entry.identifier); }

    static const std::string simple_key = "\\VariantSimple=";
    static const std::string complex_key = "\\VariantComplex=";
    std::vector<FASTAFile::FASTAEntry> windows;
    Size skipped = 0, complex_skipped = 0;
    std::vector<Size> sites, local_sites;
    const auto record = [this](const std::string& protein, Size position, char reference, const std::string& variant,
                               VariantKind kind, const std::string& reason, const std::string& window)
    {
      VariantStatus status;
      status.protein = protein;
      status.position = position;
      status.reference = reference;
      status.variant = variant;
      status.kind = kind;
      status.searched = reason.empty();
      status.reason = reason;
      status.window = window;
      status_.push_back(std::move(status));
    };
    const Size n_entries = db.size();
    for (Size i = 0; i < n_entries; ++i)
    {
      // (fast path for entries without variants)
      if (db[i].description.find(simple_key) == std::string::npos && db[i].description.find(complex_key) == std::string::npos)
      {
        continue;
      }
      const PEFFEntry peff = PEFFFile::parseFASTAEntry(db[i]);
      const std::string& protein = db[i].sequence;
      const std::string& identifier = db[i].identifier;
      const Size L = protein.size();
      for (const PEFFVariantComplex& complex : peff.complex_variants)
      {
        // e.g. "120-122|AG", "120-122|" (a deletion)
        record(identifier, complex.start_position, '\0', StringUtils::toStr(complex.start_position) + "-"
               + StringUtils::toStr(complex.end_position) + "|" + complex.replacement, VariantKind::COMPLEX, "complex_not_supported", "");
        ++complex_skipped;
      }
      if (peff.simple_variants.empty()) { continue; }
      if (segment_decoys) { cleavageSites(segments, protein, sites); }

      std::set<std::pair<Size, char>> seen;
      for (const PEFFVariantSimple& simple : peff.simple_variants)
      {
        const char alt = static_cast<char>(std::toupper(static_cast<unsigned char>(simple.variant_aa)));
        const bool stop = alt == '*';
        const bool in_range = simple.position > 0 && simple.position <= L;
        const char ref = in_range ? static_cast<char>(std::toupper(static_cast<unsigned char>(protein[simple.position - 1]))) : '?';
        // residue (as in the database), 1-based protein position, variant residue ('*': stop codon), e.g. "R175H", "Q13*"
        const std::string variant = (in_range ? std::string(1, protein[simple.position - 1]) : std::string())
                                    + StringUtils::toStr(simple.position) + std::string(1, alt);
        if (!seen.emplace(simple.position, alt).second) { continue; } // repeated
        std::string reason;
        if (!in_range) { reason = "out_of_range"; }
        else if (!stop && !isIndexableResidue(alt)) { reason = "unknown_residue"; }
        else if (alt == ref) { reason = "same_as_reference"; }
        else if (stop && simple.position == 1) { reason = "stop_at_first_residue"; }
        if (!reason.empty())
        {
          record(identifier, simple.position, in_range ? ref : '\0', variant, stop ? VariantKind::STOP_GAIN : VariantKind::SIMPLE,
                 reason, "");
          ++skipped;
          continue;
        }
        const Size v = simple.position - 1; // 0-based
        // the variant protein ends before a stop codon: its new C-terminus is a protein terminus
        const Size variant_length = stop ? v : L;

        // The residues the variant can change in the target. A substitution: the variant, the residue before and the
        // two after it. A cleavage site depends on the residues around it (ProteaseDB rules look at most two residues
        // back and one ahead), so the variant may add or remove the sites before and after it and the one after the next
        // residue: the peptides that end or start there overlap these residues. A stop codon: the last residue, i.e. the
        // peptides that end at the new C-terminus (the protein is the same up to there).
        const Size r0 = v > 0 ? v - 1 : 0;
        const Size r1 = stop ? v : std::min(L, v + 3);

        // The protein positions the decoy of the window may change. Decoys per segment keep the residues of each segment
        // in place, reordered within it: the segments that hold the changed residues, with and without the variant (a stop
        // codon: the last segment, which reversePeptides() reverses as a whole). Otherwise the decoy is the window
        // reversed, and the changed residues mirrored (below).
        std::pair<Size, Size> decoy_changed{r0, r1};
        if (segment_decoys && stop)
        {
          decoy_changed = {siteAtOrBefore(sites, r0), v}; // (the sites before the stop are those of the protein)
        }
        else if (segment_decoys)
        {
          // The cleavage sites with the variant: those of the protein, except next to the variant, where they are taken
          // from a digest of the stretch around it (sites v..v+2 can change; with a margin of two on each side)
          const Size lo = v >= 6 ? v - 6 : 0;
          const Size hi = std::min(L, v + 8);
          std::string stretch = protein.substr(lo, hi - lo);
          stretch[v - lo] = alt;
          cleavageSites(segments, stretch, local_sites);
          const Size changed_first = v >= 2 ? v - 2 : 0, changed_last = v + 4;
          // with the variant: the last site <= r0 and the first site > r1 - 1 (r0, r1 - 1 lie in [changed_first, changed_last])
          Size before = changed_first > 0 ? siteAtOrBefore(sites, changed_first - 1) : 0;
          Size after = siteAfter(sites, changed_last, L);
          for (Size site : local_sites)
          {
            const Size p = site + lo;
            if (p < changed_first || p > changed_last || (p == lo && lo > 0) || p >= hi) { continue; }
            if (p <= r0) { before = std::max(before, p); }
            if (p >= r1) { after = std::min(after, p); }
          }
          decoy_changed.first = std::min(siteAtOrBefore(sites, r0), before);
          decoy_changed.second = std::max(siteAfter(sites, r1 - 1, L), after);
        }

        Size a = 0, b = variant_length; // the window [a, b) of the protein
        if (!whole_proteins && segment_decoys)
        {
          // the margin beyond the changed residues, at cleavage sites, and one more segment at the end (a stop codon:
          // the window ends at the new C-terminus)
          const Size changed_begin = std::min(r0, decoy_changed.first);
          const Size changed_end = std::max(r1, decoy_changed.second);
          a = changed_begin > margin ? siteAtOrBefore(sites, changed_begin - margin) : 0;
          if (!stop)
          {
            b = siteAtOrAfter(sites, changed_end + margin, L);
            if (b < L) { b = siteAfter(sites, b, L); }
          }
          while (a > 0 && protein[a] == 'M') { a = siteAtOrBefore(sites, a - 1); }
        }
        else if (!whole_proteins)
        {
          a = r0 > margin ? r0 - margin : 0;
          if (!stop) { b = std::min(L, r1 + margin); }
          while (a > 0 && protein[a] == 'M') { --a; }
        }

        FASTAFile::FASTAEntry window;
        window.identifier = identifier + windowSeparator() + variant;
        for (Size n = 2; identifiers.count(window.identifier) > 0; ++n)
        {
          window.identifier = identifier + windowSeparator() + variant + "_" + StringUtils::toStr(n);
        }
        identifiers.insert(window.identifier);
        window.sequence = protein.substr(a, b - a);
        if (!stop) { window.sequence[v - a] = alt; }
        window.description = "variant " + variant + " of " + identifier + " (residues " + StringUtils::toStr(a + 1)
                             + "-" + StringUtils::toStr(b) + (stop ? ", truncated)" : ")");

        Window info;
        info.parent = identifier;
        info.shift = a;
        info.length = b - a;
        info.region = {r0 - a, r1 - a};
        info.variant = variant;
        (*regions_)[window.identifier] = info.region;
        // the decoy's region, kept until addDecoys() (by the target's identifier)
        Window decoy = info;
        decoy.decoy = true;
        if (segment_decoys)
        {
          // in place: the residues at the protein decoy's positions
          decoy.region = {decoy_changed.first - a, decoy_changed.second - a};
        }
        else
        {
          // the decoy of the window is the window reversed: a stretch of the protein's reversed decoy
          // (a decoy that keeps the protein's initial M in place has the reversed rest after it; a window that does not
          // start the protein does not start with M, so it keeps nothing)
          const bool protein_keeps_met = settings.clip_nterm_methionine && L > 1 && protein[0] == 'M';
          decoy.shift = L - b + ((protein_keeps_met && a > 0) ? 1 : 0);
          decoy.region = {info.length - info.region.second, info.length - info.region.first};
        }
        record(identifier, simple.position, ref, variant, stop ? VariantKind::STOP_GAIN : VariantKind::SIMPLE, "", window.identifier);
        pending_decoys_.emplace(window.identifier, std::move(decoy));
        windows_.emplace(window.identifier, std::move(info));
        windows.push_back(std::move(window));
      }
    }
    if (skipped > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] database:peff_variants: " << skipped << " \\VariantSimple annotations were not searched (a "
                      << "position outside the protein, the residue it replaces, a stop codon at the first residue or a "
                      << "letter without a residue of known mass; see the variant status, -out_variant_status)." << std::endl;
    }
    if (complex_skipped > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] database:peff_variants: " << complex_skipped << " \\VariantComplex annotations (insertions, "
                      << "deletions) were not searched: not supported." << std::endl;
    }
    const Size added = windows.size();
    db.insert(db.end(), std::make_move_iterator(windows.begin()), std::make_move_iterator(windows.end()));
    return added;
  }

  void ProSEPEFF::addDecoys(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& decoy_prefix)
  {
    if (pending_decoys_.empty() || decoy_prefix.empty()) { return; }
    for (const FASTAFile::FASTAEntry& entry : db)
    {
      if (!StringUtils::hasPrefix(entry.identifier, decoy_prefix)) { continue; }
      const auto it = pending_decoys_.find(entry.identifier.substr(decoy_prefix.size()));
      if (it == pending_decoys_.end()) { continue; }
      Window decoy = it->second;
      decoy.parent = decoy_prefix + decoy.parent;
      (*regions_)[entry.identifier] = decoy.region;
      windows_[entry.identifier] = std::move(decoy);
    }
  }

  const ProSEPEFF::Window* ProSEPEFF::find(const std::string& identifier) const
  {
    const auto it = windows_.find(identifier);
    return it == windows_.end() ? nullptr : &it->second;
  }

  void ProSEPEFF::mapToParents(std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) const
  {
    if (windows_.empty()) { return; }
    std::vector<PeptideEvidence> mapped, moved_evidences;
    std::set<std::string> variants;
    for (PeptideIdentification& id : peptide_ids)
    {
      for (PeptideHit& hit : id.getHits())
      {
        const std::vector<PeptideEvidence>& evidences = hit.getPeptideEvidences();
        if (std::none_of(evidences.begin(), evidences.end(),
                         [this](const PeptideEvidence& e) { return windows_.count(e.getProteinAccession()) > 0; }))
        {
          continue;
        }
        mapped.clear();
        moved_evidences.clear();
        variants.clear();
        bool other_entries = false; // an evidence of an entry that is no window
        for (const PeptideEvidence& evidence : evidences)
        {
          const Window* window = find(evidence.getProteinAccession());
          if (window == nullptr)
          {
            other_entries = true;
            mapped.push_back(evidence);
            continue;
          }
          const Int start = evidence.getStart(), end = evidence.getEnd();
          const bool positions = start >= 0 && end >= start;
          if (positions && (static_cast<Size>(start) >= window->region.second || static_cast<Size>(end) < window->region.first))
          {
            continue; // a peptide of the protein (or of its decoy) outside the residues the variant changes
          }
          PeptideEvidence moved = evidence;
          moved.setProteinAccession(window->parent);
          if (positions)
          {
            moved.setStart(static_cast<Int>(start + window->shift));
            moved.setEnd(static_cast<Int>(end + window->shift));
          }
          moved_evidences.push_back(std::move(moved));
          variants.insert(window->parent + ":" + window->variant); // (a decoy: its decoy protein, e.g. "DECOY_sp:P01:G10A")
        }
        // The evidences moved to a protein, once per protein position: one at the position of an evidence of the protein
        // itself is that one (e.g. the peptide that ends the truncated protein of a stop codon, which the protein has with
        // the residue after it); the residues before and after the peptide are those of the protein then.
        const auto same_position = [](const PeptideEvidence& a, const PeptideEvidence& b)
        {
          return a.getProteinAccession() == b.getProteinAccession() && a.getStart() == b.getStart() && a.getEnd() == b.getEnd();
        };
        for (PeptideEvidence& moved : moved_evidences)
        {
          if (std::none_of(mapped.begin(), mapped.end(), [&](const PeptideEvidence& e) { return same_position(e, moved); }))
          {
            mapped.push_back(std::move(moved));
          }
        }
        moved_evidences.clear();
        if (mapped.empty()) { continue; } // (not expected: the hit came from a peptide in a region; keep it as it is)
        std::sort(mapped.begin(), mapped.end());
        mapped.erase(std::unique(mapped.begin(), mapped.end()), mapped.end());
        std::set<std::string> accessions;
        for (const PeptideEvidence& evidence : mapped) { accessions.insert(evidence.getProteinAccession()); }
        hit.setPeptideEvidences(mapped);
        if (hit.metaValueExists("protein_references"))
        {
          hit.setMetaValue("protein_references", accessions.size() == 1 ? "unique" : "non-unique");
        }
        if (!other_entries && !variants.empty())
        {
          std::string joined;
          for (const std::string& variant : variants) { joined += (joined.empty() ? "" : ",") + variant; }
          hit.setMetaValue(variantMetaValue(), joined);
        }
      }
    }

    for (ProteinIdentification& protein_id : protein_ids)
    {
      std::vector<ProteinHit>& hits = protein_id.getHits();
      std::unordered_set<std::string> present;
      for (const ProteinHit& hit : hits)
      {
        if (windows_.count(hit.getAccession()) == 0) { present.insert(hit.getAccession()); }
      }
      std::vector<ProteinHit> kept;
      kept.reserve(hits.size());
      std::vector<ProteinHit> added;
      for (ProteinHit& hit : hits)
      {
        const Window* window = find(hit.getAccession());
        if (window == nullptr)
        {
          kept.push_back(std::move(hit));
          continue;
        }
        if (present.insert(window->parent).second)
        {
          // the protein of the window had no hit: one like the window's (target_decoy etc.), without its residues
          ProteinHit parent = std::move(hit);
          parent.setAccession(window->parent);
          parent.setSequence("");
          parent.setDescription("");
          added.push_back(std::move(parent));
        }
      }
      kept.insert(kept.end(), std::make_move_iterator(added.begin()), std::make_move_iterator(added.end()));
      hits = std::move(kept);
    }
  }

  std::vector<ProSEPEFF::VariantCounts> ProSEPEFF::countPSMs(const std::vector<const PeptideIdentificationList*>& runs,
                                                            const double q_threshold) const
  {
    std::vector<VariantCounts> counts(status_.size());
    // the searched variants by their name in variantMetaValue() ("protein:variant"); the variants of each protein with
    // their 0-based positions (for the reference PSMs)
    std::unordered_map<std::string, Size> by_name;
    std::unordered_map<std::string, std::vector<std::pair<Size, Size>>> by_protein;
    for (Size k = 0; k < status_.size(); ++k)
    {
      const VariantStatus& status = status_[k];
      if (status.searched) { by_name.emplace(status.protein + ":" + status.variant, k); }
      if (status.position > 0) { by_protein[status.protein].emplace_back(status.position - 1, k); }
    }
    std::vector<char> in_run(status_.size());
    std::set<Size> covered;
    std::vector<std::string> names;
    for (const PeptideIdentificationList* run : runs)
    {
      if (run == nullptr) { continue; }
      std::fill(in_run.begin(), in_run.end(), 0);
      for (const PeptideIdentification& id : *run)
      {
        if (id.getHits().empty()) { continue; }
        const PeptideHit& hit = id.getHits()[0];
        if (hit.isDecoy()) { continue; }
        const bool q_values = id.getScoreType().find("q-value") != std::string::npos;
        const double q = hit.getScore();
        const bool passing = !q_values || q <= q_threshold;
        if (hit.metaValueExists(variantMetaValue()))
        {
          StringUtils::split(hit.getMetaValue(variantMetaValue()).toString(), ',', names);
          for (const std::string& name : names)
          {
            const auto it = by_name.find(name);
            if (it == by_name.end()) { continue; }
            VariantCounts& c = counts[it->second];
            ++c.psms;
            if (q_values && (c.best_q < 0.0 || q < c.best_q)) { c.best_q = q; }
            if (passing)
            {
              ++c.passing_psms;
              in_run[it->second] = 1;
            }
          }
          continue;
        }
        if (!passing) { continue; }
        covered.clear();
        const std::string residues = hit.getSequence().toUnmodifiedString();
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          const auto it = by_protein.find(evidence.getProteinAccession());
          if (it == by_protein.end() || evidence.getStart() < 0 || evidence.getEnd() < evidence.getStart()) { continue; }
          const Size start = static_cast<Size>(evidence.getStart());
          for (const auto& [position, k] : it->second)
          {
            if (position < start || position > static_cast<Size>(evidence.getEnd())) { continue; }
            // the reference residue at the position (not a variant peptide that another entry holds as well)
            const char reference = status_[k].reference;
            if (reference != '\0' && (position - start >= residues.size() || residues[position - start] != reference)) { continue; }
            covered.insert(k);
          }
        }
        for (Size k : covered) { ++counts[k].reference_passing_psms; }
      }
      for (Size k = 0; k < status_.size(); ++k) { counts[k].passing_runs += in_run[k]; }
    }
    return counts;
  }

  void ProSEPEFF::storeVariantStatus(const std::string& filename, const std::vector<const PeptideIdentificationList*>& runs,
                                     const double q_threshold) const
  {
    std::ofstream out(filename);
    if (!out)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const std::vector<VariantCounts> counts = countPSMs(runs, q_threshold);
    out << "protein\tposition\tvariant\tkind\tstatus\treason\twindow\tpsms\tpsms_q\truns_q\tbest_q\treference_psms_q\n";
    for (Size k = 0; k < status_.size(); ++k)
    {
      const VariantStatus& s = status_[k];
      const VariantCounts& c = counts[k];
      out << s.protein << '\t' << s.position << '\t' << s.variant << '\t' << kindName(s.kind) << '\t'
          << (s.searched ? "searched" : "not_searched") << '\t' << s.reason << '\t' << s.window << '\t'
          << c.psms << '\t' << c.passing_psms << '\t' << c.passing_runs << '\t';
      if (c.best_q >= 0.0) { out << c.best_q; } // (stream formatting: 0.001, not 1.0e-03)
      out << '\t' << c.reference_passing_psms << '\n';
    }
    out.close(); // (a failed write of the last block shows here)
    if (out.fail())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }
}
