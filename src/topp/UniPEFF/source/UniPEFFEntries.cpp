// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "UniPEFFEntries.h"

#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <cctype>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS::UniPEFFEntries
{
  /// PEFF VariantSimple newAA must be an amino-acid letter (or '*').
  bool isResidueCode(char c)
  {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '*';
  }

  /// 0 → "?", else decimal.
  std::string positionOrUnknown(int p)
  {
    return p == 0 ? std::string("?") : std::to_string(p);
  }


namespace
{
  /// "TrEMBL" → "tr"; everything else (Swiss-Prot/null/unknown) → "sp".
  std::string prefixForDataset(const std::string& dataset)
  {
    return (dataset == "TrEMBL") ? "tr" : "sp";
  }

  /// Strip the "(Microbial infection) " prefix and any trailing "; ..." qualifier.
  /// Matches the C# CleanPtmDescription verbatim.
  std::string cleanPtmDescription(std::string description)
  {
    static const std::string kMicrobialInfection = "(Microbial infection) ";
    if (description.rfind(kMicrobialInfection, 0) == 0)
    {
      description = description.substr(kMicrobialInfection.size());
    }
    const auto semi = description.find(';');
    if (semi != std::string::npos)
    {
      description = description.substr(0, semi);
    }
    return description;
  }

  // ──────────────────────────────────────────────────────────────────
  // Isoform reconstruction (UniProt "alternative products")
  // ──────────────────────────────────────────────────────────────────

  /// Split a UniProt @c ref attribute ("VSP_1 VSP_2") on whitespace.
  std::vector<std::string> splitRefs(const std::string& ref)
  {
    std::vector<std::string> out;
    std::string cur;
    for (char c : ref)
    {
      if (std::isspace(static_cast<unsigned char>(c)))
      {
        if (!cur.empty()) { out.push_back(cur); cur.clear(); }
      }
      else cur.push_back(c);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
  }

  /// One splice variant applied by buildIsoformSequence(): the canonical residues
  /// @c begin..end (1-based) are replaced by @c replacement_length residues (0 = deleted).
  struct SpliceEdit
  {
    int begin;
    int end;
    int replacement_length;
  };

  /// Apply the splice-variant features referenced by @p iso to @p canonical.
  /// A referenced range without a <variation> is deleted; one with a <variation>
  /// is replaced (a longer replacement extends the sequence). Features are
  /// applied C- to N-terminal so an application cannot shift the coordinates of
  /// the ones still pending. The applied splice variants are returned in @p edits,
  /// sorted N- to C-terminal. On any inconsistency (unknown ref, unknown or
  /// out-of-bounds coordinates, overlapping features, <original> text not
  /// matching the canonical sequence) @p error is set and "" returned.
  std::string buildIsoformSequence(const UniProtIsoform& iso, const std::string& canonical,
                                   const std::unordered_map<std::string, const UniProtFeature*>& vsp_by_id,
                                   std::vector<SpliceEdit>& edits, std::string& error)
  {
    edits.clear();
    struct Span { int begin; int end; const UniProtFeature* feature; };
    std::vector<Span> spans;
    for (const std::string& r : splitRefs(iso.sequence_ref))
    {
      auto it = vsp_by_id.find(r);
      if (it == vsp_by_id.end())
      {
        error = "references unknown splice variant feature '" + r + "'";
        return {};
      }
      const UniProtFeature* f = it->second;
      const int begin = f->has_range ? f->begin : f->position;
      const int end = f->has_range ? f->end : f->position;
      if (begin < 1 || end < begin || end > static_cast<int>(canonical.size()))
      {
        error = "splice variant '" + r + "' has invalid range " + std::to_string(begin) + "-" + std::to_string(end);
        return {};
      }
      if (!f->original.empty()
          && canonical.compare(begin - 1, static_cast<size_t>(end - begin + 1), f->original) != 0)
      {
        error = "splice variant '" + r + "' <original> does not match the canonical sequence at "
                + std::to_string(begin) + "-" + std::to_string(end);
        return {};
      }
      spans.push_back(Span{begin, end, f});
    }
    if (spans.empty())
    {
      error = "references no splice variant features";
      return {};
    }
    std::sort(spans.begin(), spans.end(), [](const Span& x, const Span& y) { return x.begin > y.begin; });
    for (size_t i = 1; i < spans.size(); ++i)
    {
      if (spans[i].end >= spans[i - 1].begin)
      {
        error = "splice variants '" + spans[i].feature->id + "' and '" + spans[i - 1].feature->id + "' overlap";
        return {};
      }
    }
    std::string seq = canonical;
    for (const Span& s : spans)
    {
      seq.replace(static_cast<size_t>(s.begin - 1), static_cast<size_t>(s.end - s.begin + 1), s.feature->variation);
    }
    if (seq.empty())
    {
      error = "reconstructed sequence is empty";
      return {};
    }
    for (auto it = spans.rbegin(); it != spans.rend(); ++it)
    {
      edits.push_back(SpliceEdit{it->begin, it->end, static_cast<int>(it->feature->variation.size())});
    }
    return seq;
  }

  // ──────────────────────────────────────────────────────────────────
  // Canonical → isoform coordinates (-map_canonical_features)
  // ──────────────────────────────────────────────────────────────────

  /// Isoform position of canonical residue @p p under @p edits (sorted N- to C-terminal),
  /// or 0 if @p p is unknown (0), beyond @p canonical_length or inside a splice variant
  /// (deleted or replaced, so the isoform has no such residue).
  int mapPosition(int p, const std::vector<SpliceEdit>& edits, int canonical_length)
  {
    if (p < 1 || p > canonical_length) return 0;
    int shift = 0;
    for (const SpliceEdit& e : edits)
    {
      if (p < e.begin) break;
      if (p <= e.end) return 0;
      shift += e.replacement_length - (e.end - e.begin + 1);
    }
    return p + shift;
  }

  /// UniProt feature types that become \\Processed annotations (see classifyAndAppend()).
  bool isProcessingFeature(const std::string& type)
  {
    return type == "chain" || type == "initiator methionine" || type == "propeptide"
           || type == "signal peptide" || type == "transit peptide";
  }

  /// Map the coordinates of the canonical feature @p f onto an isoform built by @p edits; returns
  /// false (and leaves @p mapped unspecified) if the isoform does not contain it:
  ///  - a single position, and each of the two linked residues of a disulfide bond or cross-link,
  ///    must lie outside every splice variant;
  ///  - a processed region (chain, propeptide, signal/transit peptide) is clipped to the canonical
  ///    residues it keeps; one that starts at the canonical N-terminus or ends at the canonical
  ///    C-terminus still does so in the isoform. It is dropped if it keeps no residue;
  ///  - any other region (e.g. a sequence variant range) must not overlap a splice variant.
  /// Features with unknown (0) coordinates are dropped, since they cannot be placed.
  bool mapFeature(const UniProtFeature& f, const std::vector<SpliceEdit>& edits, int canonical_length,
                  int isoform_length, UniProtFeature& mapped)
  {
    mapped = f;
    if (f.has_position)
    {
      mapped.position = mapPosition(f.position, edits, canonical_length);
      return mapped.position != 0;
    }
    if (!f.has_range || f.begin < 1 || f.begin > f.end || f.end > canonical_length) return false;

    if (f.type == "disulfide bond" || f.type == "cross-link")
    {
      // two linked residues, not a region: splice variants between them do not matter
      mapped.begin = mapPosition(f.begin, edits, canonical_length);
      mapped.end = mapPosition(f.end, edits, canonical_length);
      return mapped.begin != 0 && mapped.end != 0;
    }

    if (isProcessingFeature(f.type))
    {
      // first and last canonical residue of the region that no splice variant removes or replaces
      int first = f.begin;
      for (const SpliceEdit& e : edits)
      {
        if (first >= e.begin && first <= e.end) first = e.end + 1;
      }
      int last = f.end;
      for (auto it = edits.rbegin(); it != edits.rend(); ++it)
      {
        if (last >= it->begin && last <= it->end) last = it->begin - 1;
      }
      if (first > last) return false;
      mapped.begin = (f.begin == 1) ? 1 : mapPosition(first, edits, canonical_length);
      mapped.end = (f.end == canonical_length) ? isoform_length : mapPosition(last, edits, canonical_length);
      return true;
    }

    const bool overlaps = std::any_of(edits.begin(), edits.end(),
                                      [&f](const SpliceEdit& e) { return e.begin <= f.end && e.end >= f.begin; });
    if (overlaps) return false;
    mapped.begin = mapPosition(f.begin, edits, canonical_length);
    mapped.end = mapPosition(f.end, edits, canonical_length);
    return true;
  }

  /// First and last residue of @p f: a position p is p..p, so it equals a one-residue range; unknown is 0..0.
  std::pair<int, int> residueSpan(const UniProtFeature& f)
  {
    if (f.has_position) return {f.position, f.position};
    if (f.has_range) return {f.begin, f.end};
    return {0, 0};
  }

  /// True if @p a and @p b become the same annotation (see classifyAndAppend()): same type, residues
  /// (a position and a one-residue range at it are the same) and variation and, for modifications, the
  /// same cleaned description. Other descriptions (e.g. of a sequence variant or a chain) are not
  /// written, so features differing only in them are the same.
  bool sameAnnotationFeature(const UniProtFeature& a, const UniProtFeature& b)
  {
    if (a.type != b.type || residueSpan(a) != residueSpan(b) || a.variation != b.variation)
    {
      return false;
    }
    const bool description_written = a.type == "modified residue" || a.type == "cross-link"
                                     || a.type == "glycosylation site" || a.type == "lipid moiety-binding region";
    return !description_written || cleanPtmDescription(a.description) == cleanPtmDescription(b.description);
  }

  // ──────────────────────────────────────────────────────────────────
  // Feature classifier — mirrors C# FromUniProtXML feature switch
  // ──────────────────────────────────────────────────────────────────

  /// Lookup PTM by description in @p ptms, returns nullptr on miss.
  const PtmEntry* findPtm(const PtmMap& ptms, const std::string& description)
  {
    auto it = ptms.find(description);
    return (it == ptms.end()) ? nullptr : &it->second;
  }

  /// Classify a single <feature> and append into the right annotation list(s).
  /// @p record_processing / @p record_aa_mods / @p record_variants gate the three
  /// categories of -Omit* flags. Returns nothing; updates @p annotations in place.
  void classifyAndAppend(const UniProtFeature& f, const PtmMap& ptms, EntryAnnotations& a,
                         bool record_processing, bool record_aa_mods, bool record_variants,
                         const std::string& accession)
  {
    // Single-residue site: prefer <position>, else <begin>; unknown -> 0.
    const int mod_position = f.has_position ? f.position : (f.has_range ? f.begin : 0);

    if (record_processing)
    {
      if (f.type == "chain")
      {
        ProcessedItem p{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                        "PEFF:0001020", "mature protein"};
        a.processed.push_back(std::move(p));
        return;
      }
      if (f.type == "initiator methionine")
      {
        const int pos = f.has_position ? f.position : 0;
        a.processed.push_back(ProcessedItem{pos, pos, "PEFF:0001035", "initiator methionine"});
        return;
      }
      if (f.type == "propeptide")
      {
        a.processed.push_back(ProcessedItem{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                                            "PEFF:0001034", "propeptide"});
        return;
      }
      if (f.type == "signal peptide")
      {
        a.processed.push_back(ProcessedItem{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                                            "PEFF:0001021", "signal peptide"});
        return;
      }
      if (f.type == "transit peptide")
      {
        a.processed.push_back(ProcessedItem{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                                            "PEFF:0001022", "transit peptide"});
        return;
      }
    }

    if (record_aa_mods)
    {
      if (f.type == "cross-link")
      {
        std::string desc = cleanPtmDescription(f.description);
        if (desc.empty()) desc = "cross-link";
        // Generic ModRes (no CV accession) per UniPEFF policy.
        if (f.has_range)
        {
          a.regular_mods.push_back(ModResItem{f.begin, false, nullptr, desc});
          a.regular_mods.push_back(ModResItem{f.end,   false, nullptr, desc});
        }
        else
        {
          a.regular_mods.push_back(ModResItem{mod_position, false, nullptr, desc});
        }
        return;
      }
      if (f.type == "disulfide bond")
      {
        const PtmEntry* hc = findPtm(ptms, "Half cystine");
        if (f.has_range)
        {
          a.halfcys.push_back(ModResItem{f.begin, true, hc, hc ? hc->id : "half cystine"});
          const size_t first_idx = a.halfcys.size() - 1;
          a.halfcys.push_back(ModResItem{f.end,   true, hc, hc ? hc->id : "half cystine"});
          const size_t second_idx = a.halfcys.size() - 1;
          a.disulfides.push_back(DisulfidePairItem{first_idx, second_idx, true});
        }
        else
        {
          a.halfcys.push_back(ModResItem{mod_position, true, hc, hc ? hc->id : "half cystine"});
        }
        return;
      }
      if (f.type == "glycosylation site")
      {
        std::string desc = cleanPtmDescription(f.description);
        if (desc.empty()) desc = "glycosylation site";
        a.regular_mods.push_back(ModResItem{mod_position, false, nullptr, desc});
        return;
      }
      if (f.type == "lipid moiety-binding region")
      {
        std::string desc = cleanPtmDescription(f.description);
        if (desc.empty()) desc = "lipid moiety-binding region";
        const PtmEntry* p = findPtm(ptms, desc);
        a.regular_mods.push_back(ModResItem{mod_position, false, p, p ? p->id : desc});
        return;
      }
      if (f.type == "modified residue")
      {
        std::string desc = cleanPtmDescription(f.description);
        const PtmEntry* p = findPtm(ptms, desc);
        if (p != nullptr)
        {
          a.regular_mods.push_back(ModResItem{mod_position, false, p, p->id});
        }
        else
        {
          // Mirror C# stdout messages so users can sync their ptmlist.txt.
          OPENMS_LOG_INFO << "UniPEFF: failed to find this modified residue in ptmlist.txt: " << desc << std::endl;
        }
        return;
      }
    }

    if (record_variants && f.type == "sequence variant")
    {
      const std::string new_seq = f.variation;  // may be empty (a deletion)
      if (f.has_range)
      {
        if (f.begin == 0 || f.end == 0)
        {
          OPENMS_LOG_WARN << "UniPEFF: " << accession << " sequence variant with unknown range; omitted." << std::endl;
          return;
        }
        a.complex_variants.push_back(VariantComplexItem{f.begin, f.end, new_seq});
      }
      else if (f.has_position)
      {
        if (f.position == 0)
        {
          OPENMS_LOG_WARN << "UniPEFF: " << accession << " sequence variant with unknown position; omitted." << std::endl;
          return;
        }
        if (new_seq.size() == 1 && isResidueCode(new_seq[0]))
        {
          a.simple_variants.push_back(VariantSimpleItem{f.position, new_seq[0]});
        }
        else
        {
          a.complex_variants.push_back(VariantComplexItem{f.position, f.position, new_seq});
        }
      }
    }
  }

  /// Merge halfcys into regular_mods producing @c a.mods in UniPEFF sort order:
  /// (position == 0 sorts LAST, then by position ascending, then regular before halfcys on ties).
  /// Updates @c a.disulfides to reference indices into the merged list.
  void mergeHalfCystines(EntryAnnotations& a)
  {
    struct Indexed { ModResItem m; bool from_halfcys; size_t orig_idx; };
    std::vector<Indexed> all;
    all.reserve(a.regular_mods.size() + a.halfcys.size());
    for (size_t i = 0; i < a.regular_mods.size(); ++i)
    {
      all.push_back(Indexed{a.regular_mods[i], false, i});
    }
    for (size_t i = 0; i < a.halfcys.size(); ++i)
    {
      all.push_back(Indexed{a.halfcys[i], true, i});
    }
    // Stable sort with the UniPEFF tie-break.
    std::stable_sort(all.begin(), all.end(), [](const Indexed& x, const Indexed& y) {
      const bool xu = (x.m.position == 0);
      const bool yu = (y.m.position == 0);
      if (xu != yu) return !xu;  // !unknown < unknown -> known first
      if (x.m.position != y.m.position) return x.m.position < y.m.position;
      // Regular mods (is_halfcys=false) sort before halfcys (true).
      return !x.m.is_halfcys && y.m.is_halfcys;
    });

    // Build new merged index for each (was_halfcys, orig_idx) so disulfide pairs can be re-pointed.
    std::vector<size_t> halfcys_new_idx(a.halfcys.size(), 0);
    a.mods.clear();
    a.mods.reserve(all.size());
    for (size_t i = 0; i < all.size(); ++i)
    {
      if (all[i].from_halfcys) halfcys_new_idx[all[i].orig_idx] = a.mods.size();
      a.mods.push_back(all[i].m);
    }
    for (auto& d : a.disulfides)
    {
      d.idx_a = halfcys_new_idx[d.idx_a];
      d.idx_b = halfcys_new_idx[d.idx_b];
    }
    // After merge the source lists are no longer authoritative.
    a.regular_mods.clear();
    a.halfcys.clear();
  }

  /// Invalidate disulfide pairs with an endpoint beyond the sequence end (malformed
  /// input): the pair then gets no labels and no \DisulfideBond tuple in either mode,
  /// matching how out-of-bounds variants and processed regions are omitted. Unknown
  /// (position 0) endpoints stay valid — UniProt documents such bonds and they are
  /// emitted as '?'. An in-range half-cystine of such a pair is still written as a plain
  /// modification; the out-of-range one is omitted by the writer, like any modification
  /// beyond the sequence end.
  void invalidateOutOfRangeDisulfides(EntryAnnotations& a, const std::string& base_sequence,
                                      const std::string& accession)
  {
    const int seq_len = static_cast<int>(base_sequence.size());
    if (seq_len == 0) return;
    for (auto& d : a.disulfides)
    {
      if (!d.valid) continue;
      if (d.idx_a >= a.mods.size() || d.idx_b >= a.mods.size()) continue;
      const int pos_a = a.mods[d.idx_a].position;
      const int pos_b = a.mods[d.idx_b].position;
      if (pos_a > seq_len || pos_b > seq_len)
      {
        OPENMS_LOG_WARN << "UniPEFF: " << accession << " disulfide bond " << positionOrUnknown(pos_a)
                        << "-" << positionOrUnknown(pos_b) << " exceeds sequence length " << seq_len
                        << "; connectivity omitted.\n";
        d.valid = false;
      }
    }
  }

  // ──────────────────────────────────────────────────────────────────
  // Annotation-ID assignment
  // ──────────────────────────────────────────────────────────────────

  /// Option B (-annotation_identifiers): walk @p a in UniPEFF's emit order assigning
  /// monotonically increasing IDs 0, 1, 2, ... to every annotation tuple (ModRes split into
  /// PSI/Unimod/generic buckets, then VariantSimple, VariantComplex, Processed, finally
  /// DisulfideBond). The disulfide pair items receive IDs of their own; @c idx_a/@c idx_b
  /// already point to half-cystine ModRes items whose IDs have just been stamped.
  /// Tuples skipped here (unknown/out-of-bounds positions) keep @c kNoId — the writer
  /// skips exactly the same tuples, so every emitted tuple carries an id.
  void assignAnnotationIds(EntryAnnotations& a, const std::string& base_sequence,
                           bool emit_processed, bool emit_aa_mods, bool emit_variants)
  {
    uint32_t next_id = 0;
    const int seq_len = static_cast<int>(base_sequence.size());
    if (emit_aa_mods)
    {
      // Modifications beyond the sequence end are not written (see writePeffEntry).
      auto in_range = [seq_len](const ModResItem& m) { return !(seq_len > 0 && m.position > seq_len); };
      // PSI-MOD bucket: ModRes whose PTM has a PSI-MOD accession.
      for (auto& m : a.mods)
      {
        if (in_range(m) && m.ptm != nullptr && !m.ptm->psi_mod.empty()) m.annotation_id = next_id++;
      }
      // Unimod bucket: ModRes with Unimod but no PSI-MOD.
      for (auto& m : a.mods)
      {
        if (in_range(m) && m.ptm != nullptr && m.ptm->psi_mod.empty() && !m.ptm->unimod.empty()) m.annotation_id = next_id++;
      }
      // Generic bucket: ModRes with no CV accession (or no PTM at all).
      for (auto& m : a.mods)
      {
        if (in_range(m) && (m.ptm == nullptr || (m.ptm->psi_mod.empty() && m.ptm->unimod.empty()))) m.annotation_id = next_id++;
      }
    }

    if (emit_variants)
    {
      // VariantSimple: real simple list, then complex entries that get demoted to simple.
      for (auto& v : a.simple_variants)
      {
        if (v.position == 0) continue;
        if (seq_len > 0 && v.position > seq_len) continue;
        v.annotation_id = next_id++;
      }
      for (auto& v : a.complex_variants)
      {
        if (v.begin != 0 && v.begin == v.end && v.new_seq.size() == 1 && isResidueCode(v.new_seq[0])
            && !(seq_len > 0 && v.begin > seq_len))
        {
          v.annotation_id = next_id++;
        }
      }
      for (auto& v : a.complex_variants)
      {
        if (v.begin == 0 || v.end == 0) continue;
        if (v.begin > v.end || (seq_len > 0 && (v.begin > seq_len || v.end > seq_len))) continue;
        if (v.begin == v.end && v.new_seq.size() == 1 && isResidueCode(v.new_seq[0])) continue;  // demoted
        v.annotation_id = next_id++;
      }
    }

    if (emit_processed)
    {
      for (auto& p : a.processed)
      {
        if (p.begin == 0 || p.end == 0) continue;
        if (p.begin > p.end || (seq_len > 0 && (p.begin > seq_len || p.end > seq_len))) continue;
        p.annotation_id = next_id++;
      }
    }

    if (emit_aa_mods)
    {
      // Disulfide tuples continue the id sequence; only valid pairs whose endpoints
      // exist in the merged mod list get one (same guards as the writer). Valid pairs
      // have both endpoints in range, so both carry an id from the buckets above.
      for (auto& d : a.disulfides)
      {
        if (!d.valid) continue;
        if (d.idx_a >= a.mods.size() || d.idx_b >= a.mods.size()) continue;
        d.annotation_id = next_id++;
      }
    }
  }

  /// Default mode (no -annotation_identifiers): implement the selective labeling scheme
  /// from issue #9829 — of the K valid disulfide pairs (both half-cystines located in this
  /// sequence via <begin>/<end>: intrachain, or between chains of one precursor), pair k
  /// (counting from 0) labels its begin half-cystine 2k and its end half-cystine 2k+1, and
  /// the \DisulfideBond tuple itself is labeled 2K+k, referencing those two ids. The bonds
  /// follow the half-cystines so that the ids 0..3K-1 are unique within the entry (PEFF 1.0,
  /// section 3.4.2). All other annotations (including half-cystines from single-<position>
  /// features, whose partner lies in another molecule) stay unlabeled.
  void assignDisulfideLabels(EntryAnnotations& a)
  {
    std::vector<DisulfidePairItem*> reported;
    for (auto& d : a.disulfides)
    {
      if (!d.valid) continue;
      if (d.idx_a >= a.mods.size() || d.idx_b >= a.mods.size()) continue;
      reported.push_back(&d);
    }
    const auto num_bonds = static_cast<uint32_t>(reported.size());
    for (uint32_t k = 0; k < num_bonds; ++k)
    {
      DisulfidePairItem& d = *reported[k];
      a.mods[d.idx_a].annotation_id = 2 * k;
      a.mods[d.idx_b].annotation_id = 2 * k + 1;
      d.annotation_id = 2 * num_bonds + k;
    }
  }

  /// A feature with &lt;location sequence="..."&gt; is annotated on another isoform:
  /// its coordinates do not refer to the entry's canonical sequence.
  bool isIsoformScoped(const UniProtFeature& f)
  {
    return !f.location_sequence.empty();
  }

  /// Number of tuples classifyAndAppend() has added to @p a (valid before mergeHalfCystines()).
  size_t classifiedAnnotationCount(const EntryAnnotations& a)
  {
    return a.regular_mods.size() + a.halfcys.size() + a.simple_variants.size()
           + a.complex_variants.size() + a.processed.size();
  }
} // namespace

  PreparedEntry prepareEntry(UniProtEntry&& e, const PtmMap& ptms, const std::string& prefix_override,
                             bool option_b, bool record_processing, bool record_aa_mods,
                             bool record_variants, bool record_isoforms, bool map_canonical_features,
                             SkippedIsoformCounts& skipped_isoforms)
  {
    PreparedEntry pe;
    pe.source = std::move(e);
    pe.prefix = prefix_override.empty() ? prefixForDataset(pe.source.dataset) : prefix_override;

    std::vector<std::vector<SpliceEdit>> isoform_edits;  // parallel to pe.isoforms
    if (record_isoforms && !pe.source.isoforms.empty())
    {
      std::unordered_map<std::string, const UniProtFeature*> vsp_by_id;
      for (const auto& f : pe.source.features)
      {
        if (f.type == "splice variant" && !f.id.empty()) vsp_by_id.emplace(f.id, &f);
      }
      for (const auto& iso : pe.source.isoforms)
      {
        // The "displayed" isoform IS the canonical sequence that was just emitted.
        if (iso.sequence_type == "displayed") continue;
        if (iso.sequence_type != "described")
        {
          // "external" (sequence lives in another entry) / "not described": nothing to reconstruct.
          // Counted for the end-of-run summary; the per-isoform detail is debug output.
          ++skipped_isoforms[iso.sequence_type];
          OPENMS_LOG_DEBUG << "UniPEFF: " << pe.source.accession << " isoform "
                           << (iso.id.empty() ? std::string("<no id>") : iso.id)
                           << " has sequence type '" << iso.sequence_type << "'; not emitted.\n";
          continue;
        }
        if (iso.id.empty())
        {
          OPENMS_LOG_WARN << "UniPEFF: " << pe.source.accession
                          << " described isoform without <id>; omitted.\n";
          continue;
        }
        std::string error;
        std::vector<SpliceEdit> edits;
        std::string seq = buildIsoformSequence(iso, pe.source.sequence, vsp_by_id, edits, error);
        if (!error.empty())
        {
          OPENMS_LOG_WARN << "UniPEFF: " << pe.source.accession << " isoform " << iso.id << " "
                          << error << "; omitted.\n";
          continue;
        }
        pe.isoforms.push_back(IsoformEntry{iso.id, iso.name, std::move(seq), {}});
        isoform_edits.push_back(std::move(edits));
      }
    }

    // Isoform-scoped features go to the entry of the isoform they are annotated on. The
    // "displayed" isoform is the canonical sequence, so its features stay with the canonical entry.
    std::unordered_map<std::string, size_t> emitted_isoform;  // isoform accession -> index in pe.isoforms
    for (size_t i = 0; i < pe.isoforms.size(); ++i) emitted_isoform.emplace(pe.isoforms[i].accession, i);
    std::vector<std::string> displayed_ids;
    for (const auto& iso : pe.source.isoforms)
    {
      if (iso.sequence_type == "displayed" && !iso.id.empty()) displayed_ids.push_back(iso.id);
    }

    // With map_canonical_features: the canonical features that became annotations of the canonical
    // entry (mapped onto each isoform below), and the isoform-scoped features written to each
    // isoform (a mapped copy of one of them is not written again).
    std::vector<const UniProtFeature*> canonical_annotated;
    std::vector<std::vector<const UniProtFeature*>> isoform_scoped(pe.isoforms.size());

    for (const auto& f : pe.source.features)
    {
      if (!isIsoformScoped(f)
          || std::find(displayed_ids.begin(), displayed_ids.end(), f.location_sequence) != displayed_ids.end())
      {
        const size_t before = classifiedAnnotationCount(pe.annotations);
        classifyAndAppend(f, ptms, pe.annotations, record_processing, record_aa_mods,
                          record_variants, pe.source.accession);
        if (classifiedAnnotationCount(pe.annotations) > before) canonical_annotated.push_back(&f);
        continue;
      }
      auto target = emitted_isoform.find(f.location_sequence);
      if (target != emitted_isoform.end())
      {
        IsoformEntry& iso = pe.isoforms[target->second];
        const size_t before = classifiedAnnotationCount(iso.annotations);
        classifyAndAppend(f, ptms, iso.annotations, record_processing, record_aa_mods,
                          record_variants, iso.accession);
        if (classifiedAnnotationCount(iso.annotations) > before)
        {
          ++pe.isoform_scoped_applied;
          isoform_scoped[target->second].push_back(&f);
        }
        continue;
      }
      // The isoform is not written (external, not described, not reconstructable, or
      // -omit_isoforms). Classifying the feature into a scratch list tells whether it would
      // otherwise have become an annotation, so the summary counts only those (not e.g.
      // sequence conflicts, which UniPEFF never writes).
      EntryAnnotations not_applied;
      classifyAndAppend(f, ptms, not_applied, record_processing, record_aa_mods,
                        record_variants, f.location_sequence);
      if (classifiedAnnotationCount(not_applied) > 0)
      {
        ++pe.isoform_scoped_skipped;
        OPENMS_LOG_DEBUG << "UniPEFF: " << pe.source.accession << " " << f.type << " '" << f.description
                         << "' is annotated on isoform " << f.location_sequence
                         << ", which is not written; not applied.\n";
      }
    }

    // Canonical features mapped onto each isoform through the splice variants it is built from.
    // UniProt does not assert them for the isoform; this is opt-in.
    if (map_canonical_features)
    {
      const int canonical_length = static_cast<int>(pe.source.sequence.size());
      for (size_t i = 0; i < pe.isoforms.size(); ++i)
      {
        IsoformEntry& iso = pe.isoforms[i];
        for (const UniProtFeature* f : canonical_annotated)
        {
          UniProtFeature mapped;
          if (!mapFeature(*f, isoform_edits[i], canonical_length, static_cast<int>(iso.sequence.size()), mapped))
          {
            ++pe.canonical_unmapped;
            OPENMS_LOG_DEBUG << "UniPEFF: " << pe.source.accession << " " << f->type << " '" << f->description
                             << "' cannot be placed on isoform " << iso.accession
                             << " (in or across a splice variant, or unknown position); not mapped.\n";
            continue;
          }
          const auto& scoped = isoform_scoped[i];
          if (std::any_of(scoped.begin(), scoped.end(),
                          [&mapped](const UniProtFeature* s) { return sameAnnotationFeature(*s, mapped); }))
          {
            continue;  // UniProt annotates it on the isoform itself; already written
          }
          classifyAndAppend(mapped, ptms, iso.annotations, record_processing, record_aa_mods,
                            record_variants, iso.accession);
          ++pe.canonical_mapped;
        }
      }
    }

    // The canonical entry and each isoform entry are finished against their own sequence.
    auto finish = [&](EntryAnnotations& a, const std::string& sequence, const std::string& accession)
    {
      mergeHalfCystines(a);
      invalidateOutOfRangeDisulfides(a, sequence, accession);
      if (option_b)
      {
        assignAnnotationIds(a, sequence, record_processing, record_aa_mods, record_variants);
      }
      else
      {
        // No-op when the entry has no documented disulfide pairs
        // (including under -omit_amino_acid_modifications).
        assignDisulfideLabels(a);
      }
    };
    finish(pe.annotations, pe.source.sequence, pe.source.accession);
    for (IsoformEntry& iso : pe.isoforms) finish(iso.annotations, iso.sequence, iso.accession);
    return pe;
  }
} // namespace OpenMS::UniPEFFEntries
