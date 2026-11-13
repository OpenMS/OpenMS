// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEPEFF.h>

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/PEFFFile.h>

#include <algorithm>
#include <cctype>
#include <set>
#include <unordered_set>

namespace OpenMS
{
  const std::string& ProSEPEFF::variantMetaValue()
  {
    static const std::string name = "PEFF_variants";
    return name;
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
    // decoys per cleavage segment ('unspecific cleavage' has none: DecoyGenerator reverses the whole entry)
    const bool segment_aligned = settings.segment_decoys && settings.enzyme != EnzymaticDigestion::UnspecificCleavage
                                 && !whole_proteins;
    // residues beyond the changed ones: a peptide (at most peptide_max_length residues) that overlaps them ends inside
    const Size margin = settings.peptide_max_length + 1;

    ProteaseDigestion segments; // cleavage segments (no missed cleavage)
    if (segment_aligned)
    {
      segments.setEnzyme(settings.enzyme);
      segments.setMissedCleavages(0);
      segments.setSpecificity(EnzymaticDigestion::SPEC_FULL);
    }

    std::unordered_set<std::string> identifiers;
    identifiers.reserve(db.size());
    for (const FASTAFile::FASTAEntry& entry : db) { identifiers.insert(entry.identifier); }

    static const std::string variant_key = "\\VariantSimple=";
    std::vector<FASTAFile::FASTAEntry> windows;
    Size skipped = 0;
    std::vector<Size> sites, local_sites, variant_sites;
    const Size n_entries = db.size();
    for (Size i = 0; i < n_entries; ++i)
    {
      if (db[i].description.find(variant_key) == std::string::npos) { continue; } // (fast path for entries without variants)
      const PEFFEntry peff = PEFFFile::parseFASTAEntry(db[i]);
      if (peff.simple_variants.empty()) { continue; }
      const std::string& protein = db[i].sequence;
      const std::string& identifier = db[i].identifier;
      const Size L = protein.size();
      if (segment_aligned) { cleavageSites(segments, protein, sites); }

      std::set<std::pair<Size, char>> seen;
      for (const PEFFVariantSimple& simple : peff.simple_variants)
      {
        const char alt = static_cast<char>(std::toupper(static_cast<unsigned char>(simple.variant_aa)));
        if (simple.position == 0 || simple.position > L || !isIndexableResidue(alt)
            || alt == static_cast<char>(std::toupper(static_cast<unsigned char>(protein[simple.position - 1]))))
        {
          ++skipped;
          continue;
        }
        if (!seen.emplace(simple.position, alt).second) { continue; } // repeated
        const Size v = simple.position - 1; // 0-based

        // The residues the variant can change in the target: the variant, the residue before and the two after it. A
        // cleavage site depends on the residues around it (ProteaseDB rules look at most two residues back and one
        // ahead), so the variant may add or remove the sites before and after it and the one after the next residue:
        // the peptides that end or start there overlap these residues.
        const Size r0 = v > 0 ? v - 1 : 0;
        const Size r1 = std::min(L, v + 3);

        Size a = 0, b = L;                         // the window [a, b) of the protein
        std::pair<Size, Size> decoy_changed{r0, r1}; // protein positions the decoy of the window may change
        if (!whole_proteins)
        {
          if (segment_aligned)
          {
            // The cleavage sites with the variant: those of the protein, except next to the variant, where they are
            // taken from a digest of the stretch around it (sites v..v+2 can change; with a margin of two on each side)
            const Size lo = v >= 6 ? v - 6 : 0;
            const Size hi = std::min(L, v + 8);
            std::string stretch = protein.substr(lo, hi - lo);
            stretch[v - lo] = alt;
            cleavageSites(segments, stretch, local_sites);
            const Size changed_first = v >= 2 ? v - 2 : 0, changed_last = v + 4;
            variant_sites.clear();
            for (Size site : sites)
            {
              if (site < changed_first || site > changed_last) { variant_sites.push_back(site); }
            }
            for (Size site : local_sites)
            {
              const Size p = site + lo;
              if (p >= changed_first && p <= changed_last && (p > lo || lo == 0) && p < hi) { variant_sites.push_back(p); }
            }
            std::sort(variant_sites.begin(), variant_sites.end());
            variant_sites.erase(std::unique(variant_sites.begin(), variant_sites.end()), variant_sites.end());
            // The decoy keeps the residues of each segment in place, reordered within it: the segments that hold the
            // changed residues (with and without the variant) may change in the decoy
            decoy_changed.first = std::min(siteAtOrBefore(sites, r0), siteAtOrBefore(variant_sites, r0));
            decoy_changed.second = std::max(siteAfter(sites, r1 - 1, L), siteAfter(variant_sites, r1 - 1, L));
            // window: the margin beyond them, at cleavage sites, and one more segment at the end
            const Size changed_begin = std::min(r0, decoy_changed.first);
            const Size changed_end = std::max(r1, decoy_changed.second);
            a = changed_begin > margin ? siteAtOrBefore(sites, changed_begin - margin) : 0;
            b = siteAtOrAfter(sites, changed_end + margin, L);
            if (b < L) { b = siteAfter(sites, b, L); }
            while (a > 0 && protein[a] == 'M') { a = siteAtOrBefore(sites, a - 1); }
          }
          else
          {
            a = r0 > margin ? r0 - margin : 0;
            b = std::min(L, r1 + margin);
            while (a > 0 && protein[a] == 'M') { --a; }
          }
        }

        FASTAFile::FASTAEntry window;
        const std::string variant = std::string(1, protein[v]) + StringUtils::toStr(simple.position) + alt;
        window.identifier = identifier + windowSeparator() + variant;
        for (Size n = 2; identifiers.count(window.identifier) > 0; ++n)
        {
          window.identifier = identifier + windowSeparator() + variant + "_" + StringUtils::toStr(n);
        }
        identifiers.insert(window.identifier);
        window.sequence = protein.substr(a, b - a);
        window.sequence[v - a] = alt;
        window.description = "variant " + variant + " of " + identifier + " (residues " + StringUtils::toStr(a + 1)
                             + "-" + StringUtils::toStr(b) + ")";

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
        if (whole_proteins)
        {
          // the decoy of the window is that of the protein with the variant: every peptide (those the variant does not
          // change are the protein decoy's, indexed once)
          decoy.region = {0, L};
        }
        else if (segment_aligned)
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
        pending_decoys_.emplace(window.identifier, std::move(decoy));
        windows_.emplace(window.identifier, std::move(info));
        windows.push_back(std::move(window));
      }
    }
    if (skipped > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] database:peff_variants: " << skipped << " \\VariantSimple annotations were skipped (a "
                      << "position outside the protein, the residue it replaces, a stop codon or a letter without a residue "
                      << "of known mass)." << std::endl;
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
    std::vector<PeptideEvidence> mapped;
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
          mapped.push_back(std::move(moved));
          variants.insert(window->parent + ":" + window->variant); // (a decoy: its decoy protein, e.g. "DECOY_sp:P01:G10A")
        }
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
}
