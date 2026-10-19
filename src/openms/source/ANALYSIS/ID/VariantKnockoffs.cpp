// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/VariantKnockoffs.h>

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <regex>
#include <set>
#include <string_view>
#include <tuple>

namespace OpenMS
{
  const std::string VariantKnockoffs::KNOCKOFF_Q_VALUE = "knockoff_q-value";
  const std::string VariantKnockoffs::IS_KNOCKOFF = "knockoff";
  const std::string VariantKnockoffs::KNOCKOFF_STRATUM = "knockoff_stratum";

  namespace
  {
    /// substitution classes by mass, in report order; every other mass is "other"
    const std::vector<std::pair<std::string, double>> STRATA = {
      {"deamidation-like", 0.984016},
      {"methyl-like", 14.01565},
      {"oxidation-like", 15.994915}};
    const std::string OTHER = "other";
    constexpr double STRATUM_TOL_DA = 0.02;

    /// prior tiers, most reliable first
    const std::vector<std::string> TIERS = {"sample", "common", "curated", "catalogue"};

    Size tierRank(const std::string& tier)
    {
      auto it = std::find(TIERS.begin(), TIERS.end(), tier);
      return static_cast<Size>(it - TIERS.begin()); // TIERS.size() for none
    }

    /// 64-bit FNV-1a hash
    UInt64 fnv1a(const std::string& s)
    {
      UInt64 h = 0xcbf29ce484222325ULL;
      for (unsigned char c : s)
      {
        h ^= c;
        h *= 0x00000100000001b3ULL;
      }
      return h;
    }

    /// next value of a SplitMix64 generator
    UInt64 splitmix64(UInt64& state)
    {
      state += 0x9e3779b97f4a7c15ULL;
      UInt64 z = state;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      return z ^ (z >> 31);
    }

    /// class a knockoff site's next residue must share with the real site's
    char residueClass(char r)
    {
      switch (r)
      {
        case 'A': case 'V': case 'L': case 'I': case 'M': return 'h';
        case 'F': case 'W': case 'Y': return 'a';
        case 'S': case 'T': case 'N': case 'Q': return 'p';
        case 'K': case 'R': case 'H': return '+';
        case 'D': case 'E': return '-';
        default: return r;
      }
    }

    /// whether residue i is the N of an N-x-S/T sequon (x not P)
    bool sequon(const std::string& seq, Size i)
    {
      return seq[i] == 'N' && i + 2 < seq.size() && seq[i + 1] != 'P' && (seq[i + 2] == 'S' || seq[i + 2] == 'T');
    }

    /// monoisotopic mass of a residue; throws for an unknown or ambiguous one
    double residueMass(char c)
    {
      const Residue* r = ResidueDB::getInstance()->getResidue(static_cast<unsigned char>(c));
      if (r == nullptr || r->getFormula().isEmpty()) // unknown, or ambiguous (B, J, X, Z)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid", std::string(1, c));
      }
      return r->getMonoWeight(Residue::Internal);
    }

    std::string changeLabel(char ref, Size position, char alt)
    {
      return std::string(1, ref) + std::to_string(position + 1) + std::string(1, alt);
    }

    /// value of a "key=value" token of a FASTA description, or empty
    std::string descriptionToken(const std::string& description, const std::string& key)
    {
      const std::string k = key + "=";
      Size pos = 0;
      while ((pos = description.find(k, pos)) != std::string::npos)
      {
        if (pos == 0 || description[pos - 1] == ' ')
        {
          const Size start = pos + k.size();
          const Size end = description.find(' ', start);
          return description.substr(start, end == std::string::npos ? std::string::npos : end - start);
        }
        pos += k.size();
      }
      return {};
    }

    /// the accession a variant refers to: the identifier, or the middle field of 'db|ACC|NAME'
    std::vector<std::string> accessionKeys(const std::string& identifier)
    {
      std::vector<std::string> keys = {identifier};
      const Size a = identifier.find('|');
      if (a != std::string::npos)
      {
        const Size b = identifier.find('|', a + 1);
        if (b != std::string::npos && b > a + 1) keys.push_back(identifier.substr(a + 1, b - a - 1));
      }
      return keys;
    }
  }

  double VariantKnockoffs::StratumCount::falseShare() const
  {
    if (tested_before == 0) return std::numeric_limits<double>::quiet_NaN();
    return knockoff_weighted / static_cast<double>(tested_before);
  }

  std::string VariantKnockoffs::substitutionClass(double delta_mass)
  {
    for (const auto& [name, mass] : STRATA)
    {
      if (std::fabs(std::fabs(delta_mass) - mass) < STRATUM_TOL_DA) return name;
    }
    return OTHER;
  }

  std::string VariantKnockoffs::variantAccession(const std::string& accession, char ref, Size position, char alt)
  {
    return accession + "_VAR_" + changeLabel(ref, position, alt);
  }

  std::string VariantKnockoffs::knockoffAccession(const std::string& accession, char ref, Size knockoff_position, Size i, Size position, char alt)
  {
    return accession + "_KO" + std::to_string(i) + "_" + changeLabel(ref, knockoff_position, alt) + "_VAR_" + changeLabel(ref, position, alt);
  }

  bool VariantKnockoffs::parseEntry(const std::string& accession, const std::string& description, EntryInfo& info)
  {
    static const std::regex knockoff_re("^(.+)_KO[0-9]+_[A-Z][0-9]+[A-Z]_VAR_(([A-Z])[0-9]+([A-Z]))$");
    static const std::regex variant_re("^(.+)_VAR_(([A-Z])[0-9]+([A-Z]))$");
    std::smatch m;
    if (std::regex_match(accession, m, knockoff_re))
    {
      info.knockoff = true;
    }
    else if (std::regex_match(accession, m, variant_re))
    {
      info.knockoff = false;
    }
    else
    {
      return false;
    }
    info.family = m[1].str();
    info.variant = m[2].str();
    try
    {
      info.delta_mass = residueMass(m[4].str()[0]) - residueMass(m[3].str()[0]);
    }
    catch (Exception::InvalidValue&)
    {
      return false;
    }
    info.tier = descriptionToken(description, "tier");
    return true;
  }

  std::vector<Size> VariantKnockoffs::knockoffPositions(const std::string& seq, Size p, const std::vector<Size>& listed,
                                                        const ProteaseDigestion& digestion, Size per_variant, UInt64 seed)
  {
    if (per_variant == 0 || p + 1 >= seq.size()) return {};

    // enzymatic peptides without missed cleavages; span_of[i] = (start, length) of the one holding residue i
    ProteaseDigestion d(digestion);
    d.setMissedCleavages(0);
    std::vector<std::string_view> peptides;
    d.digestUnmodified(std::string_view(seq), peptides, 1, 0);
    std::vector<std::pair<Size, Size>> span_of(seq.size(), {0, seq.size()});
    for (const auto& pep : peptides)
    {
      const Size start = static_cast<Size>(pep.data() - seq.data());
      for (Size i = start; i < start + pep.size() && i < seq.size(); ++i) span_of[i] = {start, pep.size()};
    }

    const char o = seq[p];
    const auto [pa, plen] = span_of[p];
    const bool p_sequon = sequon(seq, p);
    auto next_ok = [&](Size q)
    {
      const char a = seq[p + 1], b = seq[q + 1];
      const bool exact = std::string("NQKR").find(o) != std::string::npos || a == 'P' || a == 'G' || b == 'P' || b == 'G';
      return exact ? a == b : residueClass(a) == residueClass(b);
    };

    std::vector<Size> eligible;
    for (Size q = 0; q + 1 < seq.size(); ++q)
    {
      if (seq[q] != o || !next_ok(q) || sequon(seq, q) != p_sequon) continue;
      if (!std::all_of(listed.begin(), listed.end(), [q](Size r) { return (r > q ? r - q : q - r) >= 10; })) continue;
      const auto [qa, qlen] = span_of[q];
      const Size len_diff = qlen > plen ? qlen - plen : plen - qlen;
      if (qa == pa || len_diff > 3 || qlen < 6 || qlen > 50) continue;
      eligible.push_back(q);
    }

    // partial Fisher-Yates shuffle
    UInt64 state = seed;
    const Size take = std::min(per_variant, eligible.size());
    for (Size i = 0; i < take; ++i)
    {
      const Size j = i + static_cast<Size>(splitmix64(state) % static_cast<UInt64>(eligible.size() - i));
      std::swap(eligible[i], eligible[j]);
    }
    eligible.resize(take);
    return eligible;
  }

  std::vector<FASTAFile::FASTAEntry> VariantKnockoffs::buildDatabase(const std::vector<FASTAFile::FASTAEntry>& db,
                                                                     const std::vector<Variant>& variants,
                                                                     const std::string& enzyme, Size per_variant,
                                                                     DatabaseStats& stats)
  {
    stats = DatabaseStats();
    std::map<std::string, Size> protein_of; // accession key -> index into db
    for (Size i = 0; i < db.size(); ++i)
    {
      for (const auto& key : accessionKeys(db[i].identifier)) protein_of.emplace(key, i);
    }

    // valid variants per protein, in input order
    std::map<Size, std::vector<const Variant*>> per_protein;
    std::set<std::tuple<Size, Size, char>> seen;
    for (const auto& v : variants)
    {
      auto it = protein_of.find(v.accession);
      if (it == protein_of.end())
      {
        ++stats.unknown_protein;
        continue;
      }
      const std::string& seq = db[it->second].sequence;
      if (v.position >= seq.size() || seq[v.position] != v.ref || v.ref == v.alt)
      {
        ++stats.reference_mismatch;
        continue;
      }
      residueMass(v.alt); // throws for an unknown residue
      if (!seen.emplace(it->second, v.position, v.alt).second)
      {
        ++stats.duplicates;
        continue;
      }
      per_protein[it->second].push_back(&v);
    }

    ProteaseDigestion digestion;
    digestion.setEnzyme(enzyme);
    std::vector<FASTAFile::FASTAEntry> real, knockoffs;
    for (const auto& [index, vars] : per_protein)
    {
      const FASTAFile::FASTAEntry& protein = db[index];
      std::vector<Size> listed;
      for (const Variant* v : vars) listed.push_back(v->position);
      for (const Variant* v : vars)
      {
        const std::string tier_token = v->tier.empty() ? "" : " tier=" + v->tier;
        const std::string accession = variantAccession(protein.identifier, v->ref, v->position, v->alt);
        std::string seq = protein.sequence;
        seq[v->position] = v->alt;
        std::string description = "variant_of=" + protein.identifier + " variant=" + changeLabel(v->ref, v->position, v->alt);
        if (!v->id.empty()) description += " id=" + v->id;
        real.emplace_back(accession, description + tier_token, seq);
        ++stats.variant_entries;

        const auto positions = knockoffPositions(protein.sequence, v->position, listed, digestion, per_variant,
                                                 fnv1a(protein.identifier + ":" + accession));
        if (per_variant > 0 && positions.empty()) ++stats.variants_without_knockoff;
        for (Size i = 0; i < positions.size(); ++i)
        {
          std::string ko_seq = protein.sequence;
          ko_seq[positions[i]] = v->alt;
          knockoffs.emplace_back(knockoffAccession(protein.identifier, v->ref, positions[i], i + 1, v->position, v->alt),
                                 "knockoff_of=" + accession + tier_token, ko_seq);
          ++stats.knockoff_entries;
        }
      }
    }
    real.insert(real.end(), knockoffs.begin(), knockoffs.end());
    return real;
  }

  VariantKnockoffs::Index VariantKnockoffs::indexDatabase(const std::vector<FASTAFile::FASTAEntry>& db)
  {
    Index index;
    for (const auto& e : db)
    {
      EntryInfo info;
      if (!parseEntry(e.identifier, e.description, info)) continue;
      if (info.knockoff) ++index.knockoffs_per_variant[{info.family, info.variant}];
      index.entries.emplace(e.identifier, std::move(info));
    }
    for (const auto& [key, k] : index.knockoffs_per_variant) index.max_knockoffs = std::max(index.max_knockoffs, k);
    return index;
  }

  std::vector<VariantKnockoffs::StratumCount> VariantKnockoffs::applyKnockoffQValues(
    PeptideIdentificationList& ids, const Index& index, const std::function<double(const PeptideHit&)>& rank_of,
    double cutoff, Size min_tested, bool use_all_hits)
  {
    if (index.knockoffs_per_variant.empty()) return {};
    const double k_max = static_cast<double>(std::max<Size>(index.max_knockoffs, 1));

    struct Row
    {
      PeptideHit* hit;
      double rank;
      double q;      // q-value before the knockoff q-values
      bool knockoff;
      double weight; // knockoff: 1/k of the variant it copies
      bool tested;   // real: a variant of it has knockoffs
      std::string cls;
      std::string tier;
    };
    std::vector<Row> rows;

    for (auto& id : ids)
    {
      id.sort(); // best hit first
      auto& hits = id.getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        PeptideHit& hit = hits[h];
        if (hit.isDecoy()) continue;
        const auto accessions = hit.extractProteinAccessionsSet();
        if (accessions.empty()) continue;
        bool variant_only = true, all_knockoff = true, tested = false, have_real = false;
        double weight = 1.0, delta = 0.0, ko_delta = 0.0;
        std::string tier, ko_tier;
        Size best_tier = TIERS.size() + 1;
        for (const auto& acc : accessions)
        {
          auto it = index.entries.find(acc);
          if (it == index.entries.end())
          {
            variant_only = false;
            break;
          }
          const EntryInfo& e = it->second;
          if (e.knockoff)
          {
            const Size k = index.knockoffs_per_variant.at({e.family, e.variant});
            weight = std::min(weight, 1.0 / static_cast<double>(k));
            ko_delta = e.delta_mass;
            if (ko_tier.empty() || tierRank(e.tier) < tierRank(ko_tier)) ko_tier = e.tier;
          }
          else
          {
            all_knockoff = false;
            if (!have_real) delta = e.delta_mass;
            have_real = true;
            tested = tested || index.knockoffs_per_variant.count({e.family, e.variant}) > 0;
            if (tierRank(e.tier) < best_tier)
            {
              best_tier = tierRank(e.tier);
              tier = e.tier;
            }
          }
        }
        if (!variant_only) continue;
        if (all_knockoff)
        {
          if (h == 0 || use_all_hits)
          {
            rows.push_back({&hit, rank_of(hit), hit.getScore(), true, weight, false, substitutionClass(ko_delta), ko_tier});
          }
          // never accepted, whatever rank
          hit.setScore(1.0);
          hit.setMetaValue(IS_KNOCKOFF, 1);
        }
        else if (h == 0 || use_all_hits)
        {
          rows.push_back({&hit, rank_of(hit), hit.getScore(), false, 0.0, tested, substitutionClass(delta), tier});
        }
        else
        {
          // a lower-ranked hit is not estimated: no knockoff null covers it, so it is not accepted
          hit.setScore(1.0);
        }
      }
    }

    // strata: by class, and by tier when there are several tiers; a tier with few tested PSMs at the cutoff joins the
    // class's other small tiers.
    std::set<std::string> tiers;
    for (const auto& r : rows) tiers.insert(r.tier);
    const bool by_tier = tiers.size() > 1;
    std::map<std::pair<std::string, std::string>, Size> tested_at_cutoff; // (class, tier) -> tested real PSMs
    for (const auto& r : rows)
    {
      if (!r.knockoff && r.tested && r.q <= cutoff) ++tested_at_cutoff[{r.cls, r.tier}];
    }
    auto stratum_of = [&](const Row& r) -> std::string
    {
      if (!by_tier) return r.cls;
      auto it = tested_at_cutoff.find({r.cls, r.tier});
      const Size n = it == tested_at_cutoff.end() ? 0 : it->second;
      return r.cls + "/" + (n >= min_tested ? (r.tier.empty() ? std::string("none") : r.tier) : std::string("pooled"));
    };
    std::map<std::string, std::vector<Row*>> by_stratum;
    for (auto& r : rows) by_stratum[stratum_of(r)].push_back(&r);

    std::vector<StratumCount> counts;
    for (auto& [name, members] : by_stratum)
    {
      StratumCount count;
      count.name = name;
      for (const Row* r : members)
      {
        if (r->q > cutoff) continue;
        if (r->knockoff)
        {
          count.knockoff_weighted += r->weight;
        }
        else
        {
          ++count.real_before;
          if (r->tested) ++count.tested_before;
        }
      }
      // best first; at equal rank a knockoff counts first
      std::stable_sort(members.begin(), members.end(), [](const Row* a, const Row* b)
      {
        if (a->rank != b->rank) return a->rank > b->rank;
        return a->knockoff && !b->knockoff;
      });
      std::vector<double> fdp(members.size());
      double w = 0.0;
      Size t = 0;
      for (Size i = 0; i < members.size(); ++i)
      {
        if (members[i]->knockoff) w += members[i]->weight;
        else if (members[i]->tested) ++t;
        fdp[i] = std::min(1.0, (1.0 / k_max + w) / static_cast<double>(std::max<Size>(t, 1)));
      }
      for (Size i = fdp.size(); i-- > 1;) fdp[i - 1] = std::min(fdp[i - 1], fdp[i]);
      for (Size i = 0; i < members.size(); ++i)
      {
        PeptideHit& hit = *members[i]->hit;
        hit.setMetaValue(KNOCKOFF_STRATUM, name);
        if (members[i]->knockoff) continue;
        if (members[i]->tested)
        {
          hit.setMetaValue(KNOCKOFF_Q_VALUE, fdp[i]);
          hit.setScore(std::max(hit.getScore(), fdp[i]));
          if (hit.getScore() <= cutoff) ++count.real_after;
        }
        else
        {
          hit.setScore(1.0);
        }
      }
      counts.push_back(count);
    }

    // a knockoff accession on any other PSM is not reported either
    for (auto& id : ids)
    {
      for (auto& hit : id.getHits())
      {
        if (hit.metaValueExists(IS_KNOCKOFF)) continue;
        const auto& evidences = hit.getPeptideEvidences();
        std::vector<PeptideEvidence> kept;
        for (const auto& ev : evidences)
        {
          auto it = index.entries.find(ev.getProteinAccession());
          if (it == index.entries.end() || !it->second.knockoff) kept.push_back(ev);
        }
        if (!kept.empty() && kept.size() != evidences.size()) hit.setPeptideEvidences(std::move(kept));
      }
      id.sort();
    }

    std::sort(counts.begin(), counts.end(), [](const StratumCount& a, const StratumCount& b)
    {
      auto order = [](const std::string& n)
      {
        const std::string cls = n.substr(0, n.find('/'));
        for (Size i = 0; i < STRATA.size(); ++i) if (STRATA[i].first == cls) return i;
        return STRATA.size();
      };
      return std::make_pair(order(a.name), a.name) < std::make_pair(order(b.name), b.name);
    });
    return counts;
  }
} // namespace OpenMS
