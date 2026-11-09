// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHGappedTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTCandidateRanker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassLadderAlgorithm.h>
#include <OpenMS/SYSTEM/File.h>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_in_parallel() 0
#endif
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
// One matched tag is enough to propose a protein: a correct protein often has a single tag in a spectrum,
// and the protein-level score and the extension decide between the candidates.
inline const int min_tag_count = 1;

namespace
{
// FDR classes of the PrSMs (meta value FDRClass): determined mass and termini by number of blind modifications (classes 0 to
// fdr_class_max_blind, the last one with that many or more), and undetermined by whether both termini are known. The "; N blind" form is
// kept for parsers of the label. The PrSMs of precursor-ladder candidates (all determined) form a class of their own by their source,
// so that their decoys are not hidden in the large, nearly decoy-free determined classes. So do the PrSMs of each residual round of the
// chimeric search (option chimeric_rounds; meta value ChimericRound r > 0): a search of what is left of a spectrum after its best PrSM
// has its own chance rate.
constexpr int fdr_class_max_blind = 2;
constexpr int fdr_class_both_termini = fdr_class_max_blind + 1;
constexpr int fdr_class_open_terminus = fdr_class_both_termini + 1;
constexpr int fdr_class_precursor_ladder = fdr_class_open_terminus + 1;
// the residual rounds of the chimeric search: round r is class fdr_class_residual + r - 1
constexpr int max_chimeric_rounds = 2;
constexpr int fdr_class_residual = fdr_class_precursor_ladder + 1;
constexpr int fdr_class_count = fdr_class_residual + max_chimeric_rounds;
const char* const fdr_class_names[] = {"determined; 0 blind", "determined; 1 blind", "determined; 2 blind", "undetermined, both termini",
                                       "undetermined, open terminus", "precursor ladder", "residual 1", "residual 2"};
static_assert(std::size(fdr_class_names) == (Size)fdr_class_count, "one label per FDR class");
// proteoform classes (a partition of the proteoforms; the FDR below): the determined PrSM classes are one, the undetermined ones and the
// residual rounds as for PrSMs
constexpr int pf_class_determined = 0, pf_class_both_termini = 1, pf_class_open_terminus = 2, pf_class_precursor_ladder = 3, pf_class_residual = 4;
constexpr int pf_class_count = pf_class_residual + max_chimeric_rounds;
const char* const pf_class_names[] = {"determined", "undetermined, both termini", "undetermined, open terminus", "precursor ladder", "residual 1", "residual 2"};
static_assert(std::size(pf_class_names) == (Size)pf_class_count, "one label per proteoform class");
// chimeric search: with a mass matched by the best PrSM, the masses within the spectrum's tolerance of it moved by up to this many isotope
// spacings are removed as well (FLASHDeconv's monoisotopic errors of the same fragment). On the 20 native spectra of the prototype, the
// best PrSM's fragments matched 162 and 211 masses at +1 and -1 isotope, 80 and 85 at +2 and -2, against about 6 per offset by chance
// (fragments moved by +-9, +-13, +-21 Da; final_all/chimeric/isotope_offsets.tsv); the prototype removed them up to +-2.
constexpr int chimeric_isotope_offsets = 2;
// protein level: the null of a class is exponential above this quantile of its decoy scores
constexpr double null_tail_quantile = 0.8;
// protein level: with fewer decoys the null is the empirical survival function only, without the tail
constexpr Size null_min_tail_decoys = 5;
// protein level: a class with fewer decoys uses the null of the decoys of all classes pooled
constexpr Size null_min_class_decoys = 20;
// protein level: floor of a null p-value and of a local FDR before its log
constexpr double fdr_p_floor = 1e-300;
// raw-spectrum scoring: a PrSM is ranked for the FDR by its extension score plus this weight times its raw score (-log10 of the chance
// of its matched fragment envelopes). Eclip_1480 at 1%: 4846 PrSMs by the extension score, 5189 by the raw score alone, 5241 at 4 x raw,
// 5247 at 8, 5240 at 16, 5012 at 1
constexpr double raw_score_weight = 8.0;
// path alternatives (with native, 2 per candidate protein and extension mode): the next-ranked paths within this margin of the mode's best
// path score, the penalty of one blind modification
constexpr int native_path_alternatives = 2;
constexpr int path_alternative_margin = 19;
// native: a candidate is characterized when its score (-log10 of the chance probability of its matches) exceeds that of the best decoy
// candidate of its source in its spectrum by this margin (see the native rule below)
constexpr double native_margin = 2.5;

// the registry index of a meta value name used below, registered once: a lookup by name takes the registry's critical section (the map is
// read-only after its construction, so threads share it without a lock)
UInt tntKey(const char* name)
{
  static const std::unordered_map<std::string, UInt> keys = [] {
    auto& registry = MetaInfoInterface::metaRegistry();
    std::unordered_map<std::string, UInt> k;
    for (const char* n : {"BlindModCount", "DeconvMassInfo", "DecoyPrecursor", "EndPosition", "ExtensionScore", "FDRClass", "FastaIndex", "GivenMass", "Index", "Mass", "MassLostToFeatureOwner", "ModificationEnds", "ModificationStarts", "Modifications", "PrSMIndices", "ProteinQvalue", "ProteoformIndex", "ProteoformMassByFragmentMass", "RT", "RawCoverage", "RawFragmentHypotheses", "RawMatchedFragments", "RawScore", "Representative", "Scan", "StartPosition", "TagIndices", "proqvalue", "qvalue"})
      k.emplace(n, registry.registerName(n));
    return k;
  }();
  return keys.at(name);
}

// The candidate proteoforms a blind shift opens on a protein: the windows (residue i to j) of @p sequence whose mass is within @p max_shift of
// @p mass, each of which a shift would make fit; every window when the mass is unknown (<= 0). Without a shift only the windows at the
// mass itself are candidates, so the count is the multiple-testing burden of a PrSM with a shift, as TopPIC counts its candidates by type.
double shiftedCandidateCount(const std::string& sequence, double mass, double max_shift)
{
  const Size n = sequence.size();
  if (n == 0) return 1.0;
  if (mass <= 0) return (double)n * (double)(n + 1) / 2.0;
  std::vector<double> prefix(n + 1, 0.0);
  for (Size i = 0; i < n; i++)
  {
    const Residue* r = ResidueDB::getInstance()->getResidue((unsigned char)sequence[i]);
    prefix[i + 1] = prefix[i] + (r != nullptr ? r->getMonoWeight(Residue::Internal) : 0.0);
  }
  static const double water = EmpiricalFormula("H2O").getMonoWeight();
  double count = 0;
  for (Size i = 0; i < n; i++)
  {
    // windows i..j-1 with prefix[j] in [prefix[i] + mass - water - max_shift, prefix[i] + mass - water + max_shift]
    const double low = prefix[i] + mass - water - max_shift, high = prefix[i] + mass - water + max_shift;
    const auto first = std::lower_bound(prefix.begin() + (long)i + 1, prefix.end(), low);
    const auto last = std::upper_bound(first, prefix.end(), high);
    count += (double)(last - first);
  }
  return std::max(1.0, count);
}

// mass and both termini known (unknown: -1)
bool isDetermined(const ProteinHit& hit)
{
  return ! ((double)hit.getMetaValue(tntKey("Mass")) < 0 || (int)hit.getMetaValue(tntKey("StartPosition")) < 0 || (int)hit.getMetaValue(tntKey("EndPosition")) < 0);
}

// The better of two PrSMs of a spectrum, as the best PrSM per spectrum is chosen after raw-spectrum scoring: the higher score, then a
// determined mass and termini, then the higher coverage (an exact tie keeps the earlier one).
bool isBetterPrSM(const ProteinHit& a, const ProteinHit& b)
{
  auto is_determined = [](const ProteinHit& h) {
    return (double)h.getMetaValue(tntKey("Mass")) > 0 && (int)h.getMetaValue(tntKey("StartPosition")) > 0 && (int)h.getMetaValue(tntKey("EndPosition")) > 0;
  };
  if (a.getScore() != b.getScore()) return a.getScore() > b.getScore();
  if (is_determined(a) != is_determined(b)) return is_determined(a);
  return a.getCoverage() > b.getCoverage();
}

// an item of a target-decoy ranking
template<typename Key>
struct RankedItem
{
  Key key;       // equal keys (exact) are tied
  bool decoy;
  double weight; // of a decoy
};

// q-values of a target-decoy ranking (best first). After each tie group, FDR = D / (D + T) with D the decoy weight and T the target
// count so far; q is the running minimum of the FDR from the bottom. No pseudo-count: with (D + 1) / (D + T + 1), a class of fewer than
// 99 targets above its best decoy passed nothing at 1% (a single ETD SEC fraction of PXD012247: 83 determined proteoforms and no decoy,
// none reported). @p pseudo_decoys starts the decoy count above 0 (for a class whose decoys miss its errors above the first decoy).
template<typename Key>
std::vector<double> targetDecoyQvalues(const std::vector<RankedItem<Key>>& items, double pseudo_decoys = 0.0)
{
  std::vector<double> qvalues(items.size(), 1.0);
  double decoy_count = pseudo_decoys, target_count = 0;
  for (Size i = 0; i < items.size();)
  {
    Size j = i;
    for (; j < items.size() && items[j].key == items[i].key; j++)
    {
      if (items[j].decoy) decoy_count += items[j].weight;
      else target_count++;
    }
    std::fill(qvalues.begin() + i, qvalues.begin() + j, std::min(1.0, decoy_count / std::max(1.0, decoy_count + target_count)));
    i = j;
  }
  double cummin = 1.0;
  for (Size k = items.size(); k-- > 0;)
  {
    cummin = std::min(cummin, qvalues[k]);
    qvalues[k] = cummin;
  }
  return qvalues;
}

// Null distribution of the PrSM score from a sample of decoy scores: the empirical survival function with +1 and, from
// null_min_tail_decoys on, an exponential tail above the null_tail_quantile quantile (mean excess as its scale), which goes on decreasing
// beyond the best decoy.
class DecoyScoreNull
{
public:
  explicit DecoyScoreNull(std::vector<double> scores): scores_(std::move(scores))
  {
    std::sort(scores_.begin(), scores_.end());
    const Size n = scores_.size();
    if (n < null_min_tail_decoys) return;
    with_tail_ = true;
    threshold_ = scores_[(Size)(null_tail_quantile * (double)(n - 1))];
    // mean excess over the threshold; Neumaier-compensated sum, bit-identical to the reference implementation of this FDR
    double excess = 0, compensation = 0;
    Size excess_count = 0;
    for (double y : scores_)
    {
      if (y <= threshold_) continue;
      const double e = y - threshold_, t = excess + e;
      compensation += std::abs(excess) >= std::abs(e) ? (excess - t) + e : (e - t) + excess;
      excess = t;
      excess_count++;
    }
    scale_ = excess_count > 0 ? (excess + compensation) / (double)excess_count : 1.0;
    survival_at_threshold_ = (countAtLeast_(threshold_) + 1.0) / ((double)n + 1.0);
  }

  // the chance of a decoy score >= x, floored at fdr_p_floor
  double pValue(double x) const
  {
    const double survival = with_tail_ && x > threshold_ ? survival_at_threshold_ * std::exp(-(x - threshold_) / scale_)
                                                         : (countAtLeast_(x) + 1.0) / ((double)scores_.size() + 1.0);
    return std::max(survival, fdr_p_floor);
  }

private:
  // number of decoy scores >= x
  double countAtLeast_(double x) const
  {
    return (double)(scores_.end() - std::lower_bound(scores_.begin(), scores_.end(), x));
  }

  std::vector<double> scores_; // ascending
  bool with_tail_ = false;
  double threshold_ = 0, scale_ = 1.0, survival_at_threshold_ = 1.0;
};

// Density of the targets' null p-values: the slope of the least concave majorant of their ECDF through (0,0) and (1,1) (Grenander
// estimate), non-increasing in p. A segment is left-open and right-closed: a p at a vertex takes the slope to its left.
class TargetPvalueDensity
{
public:
  explicit TargetPvalueDensity(std::vector<double> pvalues)
  {
    std::sort(pvalues.begin(), pvalues.end());
    const double n = (double)pvalues.size();
    std::vector<std::pair<double, double>> hull;
    auto add = [&hull](double x, double y) {
      while (hull.size() >= 2)
      {
        const auto [x1, y1] = hull[hull.size() - 2];
        const auto [x2, y2] = hull.back();
        if ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) >= 0) hull.pop_back(); // also drops collinear points
        else break;
      }
      hull.emplace_back(x, y);
    };
    add(0.0, 0.0);
    for (Size k = 0; k < pvalues.size(); k++) add(pvalues[k], ((double)k + 1.0) / n); // ties stay separate points
    add(1.0, 1.0);
    for (Size j = 0; j < hull.size(); j++)
    {
      xs_.push_back(hull[j].first);
      if (j + 1 == hull.size()) break;
      const double dx = hull[j + 1].first - hull[j].first;
      slopes_.push_back(dx > 0 ? (hull[j + 1].second - hull[j].second) / dx : std::numeric_limits<double>::infinity());
    }
  }

  double density(double p) const
  {
    // the hull never drops (0,0) and ends with (1,1), so there is at least one segment
    const SignedSize j = (SignedSize)(std::lower_bound(xs_.begin(), xs_.end(), p) - xs_.begin()) - 1;
    return slopes_[(Size)std::clamp(j, (SignedSize)0, (SignedSize)slopes_.size() - 1)];
  }

private:
  std::vector<double> xs_;     // hull vertices
  std::vector<double> slopes_; // of the segments between them
};

} // namespace

FLASHTnTAlgorithm::FLASHTnTAlgorithm(): DefaultParamHandler("FLASHTnTAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHTnTAlgorithm::FLASHTnTAlgorithm(const FLASHTnTAlgorithm& other): DefaultParamHandler(other), ProgressLogger(other)
{
}

FLASHTnTAlgorithm& FLASHTnTAlgorithm::operator=(const FLASHTnTAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHTnTAlgorithm::setDefaultParams_()
{
  defaults_.setValue("pro_fdr", 0.01,
                     "Specifies the proteoform-level FDR; proteins are reported at the same level. The reported PrSMs are those of the reported "
                     "proteoforms (PrSM q-values only rank the PrSMs within their FDR class; they are no filter). Also the PrSM q-value up to "
                     "which a determined PrSM is preferred as proteoform representative and counts as passing in the protein ranking (0.01 "
                     "when 1). 1 reports every proteoform with its q-value.");
  defaults_.setMinFloat("pro_fdr", 0.0);

  defaults_.setValue("discard_underdetermined", "false",
                     "Discards underdetermined proteoform IDs (e.g., those without exact precursor masses or start/end positions).");
  defaults_.setValidStrings("discard_underdetermined", {"true", "false"});

  defaults_.setValue("keep_decoy", "false", "Retains decoy hits in the results.");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaults_.setValue("protein_level_only", "false",
                     "Stops after the protein-level search, skipping proteoform characterization. The candidate proteins "
                     "and their scores are still reported, so this can be used to shortlist proteins cheaply.");
  defaults_.setValidStrings("protein_level_only", {"true", "false"});

  defaults_.setValue("max_protein_hits", 50,
                     "Number of candidate proteins kept per spectrum by the protein-level search. When proteoform search runs, only these are "
                     "characterized. With more candidates, the best of the union of the 500 best by score and by matched tags are kept by a small "
                     "model of each candidate's own evidence (FLASHTnTCandidateRanker: its score, the chance probability of the ion ladder at its "
                     "tags' offset, whether that offset is a residue boundary; not the precursor mass), trained to keep the protein of each "
                     "spectrum's best PrSM of a 500-candidate search: left out by dataset, its 10 best held it in 94.6% of the spectra (gapped tags "
                     "89.2%), the union of the 100 best by score and by matched tags in 93.8% (86.8%).");
  defaults_.setMinInt("max_protein_hits", 1);
  defaults_.addTag("max_protein_hits", "advanced");

  defaults_.setValue("min_protein_hits", 5,
                     "Number of best-scoring candidate proteins per spectrum that are always characterized. The other "
                     "candidates are characterized only if their protein, or its target-decoy partner, passes protein_fdr. Used only when the database "
                     "contains decoys (DECOY prefix); otherwise all candidates kept by max_protein_hits are characterized.");
  defaults_.setMinInt("min_protein_hits", 1);
  defaults_.addTag("min_protein_hits", "advanced");

  defaults_.setValue("protein_fdr", 0.1,
                     "Run-level protein FDR, from each protein's best tag-based score over all spectra, that a candidate "
                     "beyond min_protein_hits (or its target-decoy partner) must pass to be characterized. Needs decoys in the "
                     "database. 1 characterizes all candidates kept by max_protein_hits.");
  defaults_.setMinFloat("protein_fdr", 0.0);
  defaults_.setMaxFloat("protein_fdr", 1.0);
  defaults_.addTag("protein_fdr", "advanced");

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Specifies ion types to consider.");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  defaults_.setValue("var_mod", std::vector<std::string> {"Acetyl (Protein N-term)", "Gln->pyro-Glu (Protein N-term Q)"},
                     "Variable modifications, by name with the site as in Unimod/PSI-MOD (e.g. "
                     "'Oxidation (M)', 'Acetyl (N-term)', 'Phospho (S)') or by Unimod accession number with the site (e.g. '35 (M)', '1 (N-term)'; "
                     "'UNIMOD:35 (M)' is the same; a bare number only if it has one site). A modification of any N-terminus is restricted to the "
                     "protein N-terminus by the site 'Protein N-term <residue>' (e.g. 'Gln->pyro-Glu (Protein N-term Q)'). Searched with their own "
                     "count (max_var_mod_count), independently of blind modifications; both can occur in one proteoform. More in a file: "
                     "ex:var_mod_file. Fixed modifications: the FLASHTnT tool's fixed_mod. Terminal ones are also read by the mass ladders "
                     "(terminal_ladders, precursor_ladders, native sliding ladders). The defaults: 'Acetyl (Protein N-term)', the acetylation of the "
                     "protein's first residue, or of its second after the initiator Met; it costs no penalty and, with the initiator Met removal "
                     "(always searched, also without a penalty), gives TopPIC's default N-terminal forms (NONE, NME, NME_ACETYLATION, "
                     "M_ACETYLATION). 'Gln->pyro-Glu (Protein N-term Q)', the pyroglutamate of a protein N-terminal Gln (antibody heavy chains, "
                     "secreted proteins; in a PEFF database also at the start of a mature form), with the variable modification penalty. "
                     "'Acetyl (N-term)' is the acetylation of any proteoform N-terminus, with the variable modification penalty.");

  defaults_.setValue("max_var_mod_count", 3, "Maximum number of variable modifications per proteoform.");
  defaults_.setValue("max_blind_mod_count", 1,
                     "Maximum number of blind modifications per proteoform. With 2, two blind shifts fit the precursor mass of many chance "
                     "candidates: the best decoys of the class \"undetermined, both termini\" scored above all but about 50 of its targets, and "
                     "the class passed nothing at 1%. With 1 (TopPIC's default of one unexpected modification), Eclip_1480 at 10 ppm: 1609 "
                     "proteoforms (1448 determined, 161 undetermined with both termini) instead of 1505, 306 proteins instead of 278, B. theta "
                     "entrapment 10 of 1584 proteoforms instead of 12 of 1490, and 25% faster.");
  defaults_.setValue("max_mod_mass", 500.0, "Maximum absolute mass of blind modifications.");

  auto tparam = FLASHTaggerAlgorithm().getDefaults();
  tparam.remove("ion_type");
  tparam.remove("fixed_mod"); // unused by the tagger (FLASHIda sets it); fixed modifications are the FLASHTnT tool's fixed_mod
  // the tagger's own mass gaps: gapped tags (gapped_tags) search gaps of any residue count instead
  tparam.remove("allow_gap");
  tparam.remove("max_aa_in_gap");
  defaults_.insert("tag:", tparam);

  auto eparam = FLASHExtenderAlgorithm().getDefaults();
  eparam.remove("ion_type");
  eparam.remove("max_var_mod_count");
  eparam.remove("max_blind_mod_count");
  eparam.remove("max_mod_mass");
  eparam.remove("fixed_mod"); // unused by the extender (FLASHIda sets it)
  eparam.remove("var_mod");
  // the extender's defaults: the cleavage propensity weight and table, fitted on lysates and native proteins
  eparam.remove("cleavage_weight");
  eparam.remove("cleavage_table");
  eparam.remove("cleavage_native"); // set from -native (updateMembers_)

  defaults_.insert("ex:", eparam);

  defaults_.setValue("gapped_tags", "fallback",
                     "Candidate proteins from gapped tags (masses whose differences are sums of one or more residues): off; fallback (spectra in "
                     "which no sequence tag proposes a protein); both (every spectrum, proteins not proposed by sequence tags); only (gapped "
                     "tags instead of sequence tags). Where a gapped tag lands on a protein anchors the extension as a sequence tag does. In "
                     "the search cascade (the default) gapped tags search only the spectra the first stage did not identify (with only, the first "
                     "stage has no sequence tags): "
                     "Eclip_1480 at 10 ppm 6813 identified spectra instead of 6580 (1609 proteoforms instead of 1579), D60 human 1587 instead "
                     "of 1549, for about a third of the run time. Alone, gapped tags find about as much as sequence tags (1480 without the "
                     "cascade: 4547 spectra by gapped tags vs 4117 by sequence tags, 1538 and 1579 proteoforms, 305 proteins each).");
  defaults_.setValidStrings("gapped_tags", {"off", "fallback", "both", "only"});
  defaults_.setValue("terminal_ladders", "true",
                     "Candidate proteins from terminal ladders: deconvolved masses that are residue sums from an intact protein terminus (the N-terminus "
                     "also without the initiator Met; each terminus also with the terminal variable modifications of var_mod and ex:var_mod_file "
                     "on the residues the extension places them on, e.g. with the default the protein N-terminal acetylation), found without a "
                     "sequence tag. In every spectrum the min_protein_hits best proteins by the chance probability of their ladder (corrected for "
                     "the terminal forms tested on the protein) are characterized, anchored on the ladder.");
  defaults_.setValidStrings("terminal_ladders", {"true", "false"});
  defaults_.setValue("precursor_ladders", "true",
                     "Candidate proteins from precursor-anchored ladders: in each MSn spectrum with a precursor mass, every stretch of a protein (between "
                     "non-standard residues) whose mass equals the precursor mass (within the MS1 deconvolution tolerance of the input, or "
                     "the MS2 one if the input has no MS1 spectra; isotope errors -1 to +1) is scored by the chance probability of the "
                     "deconvolved masses on its prefix and suffix fragment ladders (a position counts once per side, whatever ion types match it). A "
                     "stretch is unmodified or carries the terminal variable modifications of var_mod and ex:var_mod_file where the extension places "
                     "them (e.g. with the default the acetylation of a stretch from the protein's first residue, or its second after the initiator "
                     "Met; the Met removal itself is a stretch from residue 2); every form is a hypothesis of the spectrum. The min_protein_hits best "
                     "stretches are characterized from tags on their ladders. Only their PrSMs of "
                     "determined mass and termini that are not the same proteoform as a determined PrSM of another source in the spectrum are kept; "
                     "they form their own FDR class (\"precursor ladder\"). Not used with protein_level_only. Without any precursor mass in the input, "
                     "results are unchanged; otherwise the PrSMs of this source can also change the others' results (the pooled protein-level null, "
                     "the best PrSM of a spectrum, proteoform clusters).");
  defaults_.setValidStrings("precursor_ladders", {"true", "false"});
  defaults_.setValue("search_cascade", "true",
                     "Search in two stages. true: stage 1 searches every spectrum with sequence tags and terminal and precursor ladders (those "
                     "that are on), which compete per spectrum as in a single stage; stage 2 searches with gapped tags only the spectra stage 1 did "
                     "not identify (no target PrSM whose proteoform passes pro_fdr, 1% without a filter). The PrSM q-values are computed per stage "
                     "(and FDR class), so the PrSMs of stage 2 compete only with the decoys of stage 2. Gapped tags then rescue spectra the other "
                     "sources leave unidentified instead of taking the place of their "
                     "candidates (Eclip_1480: 121 instead of 51 rescued spectra at the same number of proteoforms; D60 human Caco-2: 529 instead of "
                     "526 proteoforms). false: one stage searches every spectrum with every source. A single stage is also used with gapped_tags off, "
                     "with native, with protein_level_only and without decoys in the database; the chimeric rounds search with all sources.");
  defaults_.setValidStrings("search_cascade", {"true", "false"});
  defaults_.addTag("search_cascade", "advanced");
  defaults_.setValue("one_proteoform_per_feature", "true",
                     "Co-isolation: when PrSMs of different proteins report the same mass taken from MS1 (within twice the precursor tolerance, "
                     "one isotope apart allowed) within 2 min, the mass is one MS1 feature of one proteoform, and only the protein that explains it "
                     "best keeps it (by its best PrSM score, which already prices blind shifts, then fewer blind modifications, then more PrSMs at that mass); the "
                     "other proteins' PrSMs are reported without a mass. Masses from fragment pairs and decoy precursors are not affected.");
  defaults_.setValidStrings("one_proteoform_per_feature", {"true", "false"});
  defaults_.addTag("one_proteoform_per_feature", "advanced");
  defaults_.setValue("native", "false",
                     "Native and targeted top-down (one or a few spectra, often several proteins per spectrum, no MS1 survey): the candidate "
                     "proteins of each source (sequence tags; sliding ladders, i.e. masses read as residue sums from any start or to any end of a "
                     "protein, instead of terminal ladders; precursor ladders) are characterized when their candidate score "
                     "(-log10 of the chance probability of their matches) exceeds that of the best decoy candidate of the same source in the same "
                     "spectrum by 2.5 (at most max_protein_hits per source), and every such candidate's PrSM is reported, also without a "
                     "determined mass, one per protein (the better one); without native, each spectrum reports its best PrSM only. Gapped tags "
                     "are not used, with the raw spectra two more paths per candidate protein and extension mode are scored on them (the "
                     "next-ranked within one blind modification's penalty of the best), and at least one residual round of chimeric_rounds is searched (co-isolated "
                     "proteins). A PrSM beside a better PrSM of its spectrum is kept only with evidence of its own: the chance probability of its "
                     "matched masses that no better kept PrSM matched (UniqueMasses) must be below 10^-(2.5 + log10 of the target proteins) "
                     "(UniqueScore), so that a homolog of the spectrum's protein, which matches that protein's fragments, is not reported beside it. "
                     "Needs decoys in the database (DECOY prefix); without them, the candidates are chosen as without the option.");
  defaults_.setValidStrings("native", {"true", "false"});
  defaults_.setValue("decoy_precursor", "false",
                     "Decoy precursors: spectra whose PrSM took its mass from MS1 are searched again with the decoy precursor mass of their isolation "
                     "window instead of the precursor mass (with a raw mzML as -in, FLASHTnT has FLASHDeconv choose it; a deconvolved mzML must "
                     "come from FLASHDeconv -FD:report_FDR true, which writes it and its score). When that search gives a PrSM "
                     "reporting the decoy mass, the decoy precursor wins if its precursor score (FLASHDeconv's 2D Qscore) is higher than the "
                     "precursor's; the winning decoy PrSM replaces the spectrum's MS1-based PrSMs and counts as a decoy in the FDR, at the "
                     "proteoform level only against proteoforms whose mass came from MS1. Controls wrong MS1 precursor masses; PrSMs with masses "
                     "from fragment pairs or without a mass are not affected.");
  defaults_.setValidStrings("decoy_precursor", {"true", "false"});
  defaults_.setValue("chimeric_rounds", 0,
                     "Residual search of chimeric spectra: rounds after the first search (0: off). In each round, every spectrum with a PrSM (target or "
                     "decoy alike) loses the deconvolved masses its best PrSM matched (the masses on the path, also those one or two isotope "
                     "spacings away, and its proteoform mass), the 400-mass cap is applied again to what is left, and the residual is searched "
                     "again: tags, candidate proteins (sequence and gapped tags, terminal ladders; no precursor ladders) without the proteins of "
                     "the spectrum's earlier best PrSMs and their target-decoy partners, and the extension, without the precursor mass (it belongs "
                     "to the first proteoform). Every candidate's PrSM is kept as with -in_raw (scored on the raw spectrum with it), then the best "
                     "one (with native, all). Residual PrSMs are added to the scan (meta value ChimericRound r) and form their own FDR class "
                     "per round (\"residual 1\", \"residual 2\"). The evidence in a single spectrum is reported, not filtered: ResidualNullScore, the "
                     "best extension score of the other decoy PrSMs of the same residual search, ResidualScoreRatio, the extension score over it, "
                     "and MatchedMasses, the masses on the path. Without the option, results are unchanged.");
  defaults_.setMinInt("chimeric_rounds", 0);
  defaults_.setMaxInt("chimeric_rounds", max_chimeric_rounds);
  // raw-spectrum annotation (FLASHAnnotator's defaults, terminal fragments only): only its tolerance is a parameter
  const auto aparam = FLASHAnnotatorAlgorithm().getDefaults();
  defaults_.setValue("raw:tol", aparam.getValue("tol"), aparam.getDescription("tol"));
  defaults_.setMinFloat("raw:tol", 0.0);
  defaults_.setValue("scoring_model", "none",
                     "Learned PrSM scoring (needs -in_raw): a model whose linear choice score picks each spectrum's PrSM among its candidates and "
                     "whose gradient-boosted trees give the chosen PrSM's confidence, which replaces its score in the FDR (the combined score is kept "
                     "as CombinedScore). 'default': the model shipped with OpenMS (share/OpenMS/TOPDOWN/FLASHTnT_scoring_model.txt; not used with "
                     "-native); 'none' (default): no model, the PrSM by the combined score (the determined one unless the undetermined one scores over 3x); otherwise a model file.",
                     {"advanced"});
  defaults_.setValue("protein_open_terminus_model", "false",
                     "Protein q-values: score an undetermined PrSM with an open terminus with the model of the class \"undetermined, both "
                     "termini\" (its decoy null, share of false targets and target density), as precursor-ladder PrSMs and determined PrSMs "
                     "that are not the whole protein already are. PrSM and proteoform q-values are unchanged.",
                     {"advanced"});
  defaults_.setValidStrings("protein_open_terminus_model", {"true", "false"});
  defaults_.setValue("protein_ranking", "best_score",
                     "Protein q-values (without -native): 'best_score' ranks a protein by the score of its best proteoform, targets and decoys in one list; 'model' (with the learned scoring) ranks a protein by its best PrSM's confidence and "
                     "takes the FDR at a cut as the mean chance, from that confidence, that the proteins above it are wrong (no decoy counting: a "
                     "run's few decoy proteins never reached the top of the ranking); without the learned scoring it falls back to 'proteoform'. "
                     "'proteoform' ranks a protein by its best proteoform (the lowest proteoform q-value, then its best score); 'top_two_prsms' by "
                     "the sum of the scores of its two best PrSMs (one PrSM per spectrum; a protein seen in one spectrum keeps that spectrum's "
                     "score; with scoring_model, softplus of the confidence, -log P(wrong), is what is summed), targets and decoys in one ranking.",
                     {"advanced"});
  defaults_.setValidStrings("protein_ranking", {"best_score", "model", "proteoform", "top_two_prsms"});
  defaults_.setValue("shift_penalty", 0.0,
                     "With scoring_model: a determined PrSM with a blind shift loses this weight times log10 of the candidates the shift opens "
                     "(the protein's windows within max_mod_mass of its mass) from its confidence, without a floor: the multiple-testing burden "
                     "of the free shift, as TopPIC's E-value charges it. 0: off.",
                     {"advanced"});
  defaults_.setMinFloat("shift_penalty", 0.0);
  defaults_.setValue("undetermined_decoy_weight", 1.0,
                     "With the learned scoring (scoring_model): PrSM and proteoform q-values of the class \"undetermined, both termini\": each decoy counts this many times (2E/D = 1.2 "
                     "over the class on the development runs' shuffled entrapment; at 1.2 the held-out results were those of 1 within two "
                     "identifications). 1: off.",
                     {"advanced"});
  defaults_.setMinFloat("undetermined_decoy_weight", 1.0);
  defaults_.setValue("undetermined_pseudo_decoys", 1.0,
                     "With the learned scoring (scoring_model): PrSM q-values of the class \"undetermined, both termini\": the decoy count starts at this value, so that its "
                     "best PrSMs of a run do not pass at q-value 0 above the run's first decoy (where the development runs' shuffled entrapment "
                     "held its excess errors). 0: off.",
                     {"advanced"});
  defaults_.setMinFloat("undetermined_pseudo_decoys", 0.0);
  defaults_.setValue("open_terminus_prsms", "separate",
                     "With the learned scoring (scoring_model): PrSMs of the class \"undetermined, open terminus\" (one terminus and the mass unknown): 'count' gives them PrSM q-values as "
                     "every class; 'separate' reports them with q-value 1, apart from the identifications (as their proteoforms and proteins "
                     "already are): no decoy weight or pseudo-count brought the class to 1% on the development runs (3.7-7.3% shuffled entrapment).",
                     {"advanced"});
  defaults_.setValidStrings("open_terminus_prsms", {"count", "separate"});
  defaultsToParam_();
}

// Learned PrSM scoring (scoring_model): a linear choice score over a candidate PrSM's features picks the spectrum's PrSM, and gradient-boosted
// trees give the chosen PrSM's confidence, which replaces its score in the FDR. The model file is plain text, one entry per line:
//   choice <feature> <weight>        the choice score is the weighted sum of these features
//   conf_feature <feature>           the trees' features, in their order
//   conf_baseline <value>            the trees' starting value (log-odds)
//   tree <node count>                followed by one line per node: <feature> <threshold> <left> <right> <leaf 0|1> <value> <missing left 0|1>
// The features are those of the offline pipeline that fit the model (TD_survey scripts cand_features.py, add_feats.py, conf_fit.py), on the
// deconvolved spectrum's masses: residue masses without fixed modifications, fragments of the ion types at 10 ppm (+-1 isotope above 5 kDa),
// a fragment inside a modification window counted for either placement.
struct LearnedScoringModel
{
  struct Node
  {
    int feature = 0;
    double threshold = 0;
    int left = 0, right = 0;
    bool leaf = true;
    double value = 0;
    bool missing_left = true;
  };
  std::vector<std::pair<std::string, double>> choice;
  std::vector<std::string> conf_features;
  double conf_baseline = 0;
  std::vector<std::vector<Node>> trees;

  static std::shared_ptr<const LearnedScoringModel> load(const std::string& file)
  {
    std::ifstream in(file);
    if (! in) throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file);
    auto model = std::make_shared<LearnedScoringModel>();
    std::string kind;
    while (in >> kind)
    {
      if (kind == "choice")
      {
        std::string f;
        double w;
        in >> f >> w;
        model->choice.emplace_back(f, w);
      }
      else if (kind == "conf_feature")
      {
        std::string f;
        in >> f;
        model->conf_features.emplace_back(f);
      }
      else if (kind == "conf_baseline") in >> model->conf_baseline;
      else if (kind == "tree")
      {
        Size n = 0;
        in >> n;
        std::vector<Node> tree(n);
        for (auto& node : tree)
        {
          int leaf = 0, missing_left = 0;
          in >> node.feature >> node.threshold >> node.left >> node.right >> leaf >> node.value >> missing_left;
          node.leaf = leaf != 0;
          node.missing_left = missing_left != 0;
        }
        model->trees.push_back(std::move(tree));
      }
      else throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, kind, "unknown entry in the scoring model " + file);
      if (in.fail()) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, kind, "malformed entry in the scoring model " + file);
    }
    if (model->choice.empty() || model->trees.empty())
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file, "a scoring model needs choice weights and confidence trees");
    return model;
  }

  /// the trees' log-odds for the features @p x (in conf_features order)
  double confidence(const std::vector<double>& x) const
  {
    double sum = conf_baseline;
    for (const auto& tree : trees)
    {
      int n = 0;
      while (! tree[n].leaf)
      {
        const double v = x[tree[n].feature];
        n = (std::isnan(v) ? tree[n].missing_left : v <= tree[n].threshold) ? tree[n].left : tree[n].right;
      }
      sum += tree[n].value;
    }
    return sum;
  }
};

void FLASHTnTAlgorithm::updateMembers_()
{
  tagger_param_ = param_.copy("tag:", true);
  tagger_param_.setValue("ion_type", param_.getValue("ion_type"));

  extender_param_ = param_.copy("ex:", true);
  extender_param_.setValue("ion_type", param_.getValue("ion_type"));
  extender_param_.setValue("max_var_mod_count", param_.getValue("max_var_mod_count"));
  extender_param_.setValue("max_blind_mod_count", param_.getValue("max_blind_mod_count"));
  extender_param_.setValue("max_mod_mass", param_.getValue("max_mod_mass"));
  extender_param_.setValue("var_mod", param_.getValue("var_mod"));
  // c/z cleavage propensities of native proteins (ECD) with -native, of denatured ones (ETD) otherwise
  extender_param_.setValue("cleavage_native", param_.getValue("native").toString() == "true" ? "true" : "false");

  gapped_param_ = FLASHGappedTaggerAlgorithm().getDefaults();
  annotator_param_ = FLASHAnnotatorAlgorithm().getDefaults();
  annotator_param_.setValue("tol", param_.getValue("raw:tol"));
  annotator_param_.setValue("ion_type", param_.getValue("ion_type"));
  gapped_mode_ = param_.getValue("gapped_tags").toString();
  protein_open_terminus_model_ = param_.getValue("protein_open_terminus_model").toString() == "true";
  protein_top_two_ = param_.getValue("protein_ranking").toString() == "top_two_prsms";
  protein_best_score_ = param_.getValue("protein_ranking").toString() == "best_score";
  protein_model_fdr_ = param_.getValue("protein_ranking").toString() == "model";
  // scoring_model: 'default' is the model shipped in share/OpenMS/TOPDOWN (trained on denatured HCD and ETD lysates; not with -native, which
  // keeps the previous scoring), 'none' or empty switches it off, anything else is a model file
  std::string scoring_file = param_.getValue("scoring_model").toString();
  if (scoring_file == "default")
  {
    scoring_file = param_.getValue("native").toString() == "true" ? std::string("") : std::string(File::find("TOPDOWN/FLASHTnT_scoring_model.txt"));
  }
  if (scoring_file == "none") scoring_file = "";
  scoring_model_ = scoring_file.empty() ? nullptr : LearnedScoringModel::load(scoring_file);
  shift_penalty_ = param_.getValue("shift_penalty");
  undetermined_decoy_weight_ = param_.getValue("undetermined_decoy_weight");
  undetermined_pseudo_decoys_ = param_.getValue("undetermined_pseudo_decoys");
  open_terminus_separate_ = param_.getValue("open_terminus_prsms").toString() == "separate";
  terminal_ladders_ = param_.getValue("terminal_ladders").toString() == "true";
  precursor_ladders_ = param_.getValue("precursor_ladders").toString() == "true";
  search_cascade_ = param_.getValue("search_cascade").toString() == "true";
  decoy_precursor_ = param_.getValue("decoy_precursor").toString() == "true";
  one_proteoform_per_feature_ = param_.getValue("one_proteoform_per_feature").toString() == "true";
  chimeric_rounds_ = param_.getValue("chimeric_rounds");
  pro_fdr_ = param_.getValue("pro_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("discard_underdetermined").toString() == "false";
  // One PrSM per spectrum, its best (with native, every candidate's, below). Several per spectrum put homologs of the spectrum's own
  // protein beside it, which shuffled decoys do not model: three E. coli lysates with the B. theta entrapment at 1%, 5323 proteoforms
  // with 91 B. theta (73 of them among the 750 that were not their spectrum's best PrSM; of the B. theta PrSMs beside a better one, 128
  // were homologs of its protein, 92 unrelated, against 123 unrelated decoys) instead of 5083 with 39, in the same time.
  multiple_hits_per_spec_ = false;

  // the largest total mass shift of the blind modifications of a proteoform (+ 1 Da): the bound of the tag flanking masses and of the
  // candidate blind modifications annotated by the extension
  max_total_mod_mass_ = (double)param_.getValue("max_blind_mod_count") * (double)extender_param_.getValue("max_mod_mass") + 1.0;

  protein_level_only_ = param_.getValue("protein_level_only").toString() == "true";
  max_protein_hits_ = param_.getValue("max_protein_hits");
  min_protein_hits_ = param_.getValue("min_protein_hits");
  protein_fdr_ = param_.getValue("protein_fdr");
  native_ = param_.getValue("native").toString() == "true";
  // Path alternatives (native only): per candidate protein and extension mode, the next-ranked paths within path_alternative_margin of the
  // mode's best are scored on the raw spectrum too, and the PrSM is the one of the best extension + raw score. Otherwise the extension
  // keeps one path per mode, chosen among equal scores by the smaller modification mass without looking at the raw spectrum (a native
  // OMPF spectrum: 205-362 instead of 1-362 at the same score, with the weaker raw evidence).
  path_alternatives_ = native_ ? native_path_alternatives : 0;
  max_node_cntr_ = 400;
  if (native_)
  {
    gapped_mode_ = "off";
    multiple_hits_per_spec_ = true;
    // every mass of the spectrum: the cap of 400 cut ladder masses of native spectra (Q3930 DLDH: 4 of the 8 masses of its C-terminal y
    // ladder below the cap). On the 25 native and targeted test spectra: expected proteins 33 -> 34 of 43, other proteins 9 -> 9,
    // null-database calls 2 -> 2, the same time; a cap of 600 gives the same
    max_node_cntr_ = std::numeric_limits<Size>::max();
    // and at least one residual round (chimeric_rounds): on 25 native and targeted test spectra it added 4 expected proteins (three
    // antibody heavy chains, a co-isolated OMPF) without another protein called against the null databases (29 -> 33 of 43; other
    // proteins 8 -> 8; null-database calls 1 -> 1), at 4.2 instead of 3.3 s per spectrum
    chimeric_rounds_ = std::max(chimeric_rounds_, 1);
  }
}

std::set<Size> FLASHTnTAlgorithm::getProteinsPassingFDR_(const ProteoformDatabase& db) const
{
  // protein index to the protein's best tag-based score over the run, over its entries (scored: has one); the meta value read by its
  // index, without a lookup of the name in the registry (a critical section) per hit. Each protein counts once, target or decoy.
  const auto& protein_of = db.proteinOfEntry();
  std::vector<double> best_score(db.proteins.size(), 0.0);
  std::vector<bool> scored(db.proteins.size(), false);
  const UInt mi_fasta_index = MetaInfoInterface::metaRegistry().registerName("FastaIndex");
  for (const auto& hit : protein_hits_)
  {
    Size fi = protein_of[(Size)(int)hit.getMetaValue(mi_fasta_index)];
    if (! scored[fi] || best_score[fi] < hit.getScore())
    {
      best_score[fi] = hit.getScore();
      scored[fi] = true;
    }
  }

  std::vector<std::pair<double, Size>> ranked;
  for (Size fi = 0; fi < best_score.size(); fi++)
  {
    if (scored[fi]) ranked.emplace_back(best_score[fi], fi);
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

  // target-decoy FDR at each distinct score (ties together), then the running minimum from the bottom
  std::vector<double> fdr(ranked.size(), 1.0);
  double taget_count = 0, decoy_count = 0;
  for (Size i = 0; i < ranked.size();)
  {
    Size j = i;
    for (; j < ranked.size() && ranked[j].first == ranked[i].first; j++)
    {
      if (StringUtils::hasPrefix(db.proteins[ranked[j].second].identifier, "DECOY")) { decoy_count++; }
      else { taget_count++; }
    }
    double f = taget_count > 0 ? decoy_count / decoy_factor_ / taget_count : 1.0;
    std::fill(fdr.begin() + i, fdr.begin() + j, f);
    i = j;
  }

  std::set<Size> passing; // decoys included, so that the decoy hits extended stay comparable to the targets
  double qvalue = 1.0;
  for (Size k = ranked.size(); k-- > 0;)
  {
    qvalue = std::min(qvalue, fdr[k]);
    if (qvalue <= protein_fdr_) passing.insert(ranked[k].second);
  }

  // Add the target-decoy partner of every passing protein. A protein present in the sample passes through its true matches in
  // other spectra, which a decoy never has; extended in a spectrum that is not its own, it can give a false match that no decoy
  // mirrors (Eclip_1480: beyond the top candidates 237 target vs 0 decoy matches). With its decoy partner extended under the same
  // condition, that kind of false match is counted by the target-decoy competition too. The partner is that of the database (by name
  // for a FASTA file: DECOY_ or DECOY + the name, the first entry of a name).
  std::set<Size> partners;
  for (Size i : passing)
  {
    if (db.proteins[i].partner >= 0) partners.insert((Size)db.proteins[i].partner);
  }
  passing.insert(partners.begin(), partners.end());
  return passing;
}

Size FLASHTnTAlgorithm::blindModCount_(const ProteinHit& hit)
{
  if (hit.metaValueExists("BlindModCount")) return (Size)(int)hit.getMetaValue("BlindModCount");
  // without the count: the modifications that are neither variable (1) nor fixed by the database entry (2)
  if (hit.metaValueExists("ModificationVariable"))
  {
    const auto variable = hit.getMetaValue("ModificationVariable").toIntList();
    return (Size)std::count(variable.begin(), variable.end(), 0);
  }
  return hit.getMetaValue("Modifications").toDoubleList().size();
}

bool FLASHTnTAlgorithm::areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol)
{
  // Two PrSMs of one protein are counted as one proteoform when
  // - both masses are determined: the masses agree (within twice the tolerance), or differ by one isotope (1.002371 Da, within twice
  //   the tolerance), and the ranges overlap. The monoisotopic mass of a precursor is off by one isotope often enough that the same
  //   proteoform was counted twice: on Eclip_1480, 1450 determined proteoforms became 1191 with the one-isotope rule (proteins 299
  //   either way). A deamidation (+0.984 Da) cannot be told from a one-isotope error at this tolerance and is merged with it. Two
  //   isotopes apart is not merged: a reduced disulfide (+2.016 Da) is another proteoform, and the two cannot be told apart here.
  //   Where the shifts sit does not matter: a blind shift has no exact position, and the same mass placed elsewhere (+1.0 Da at
  //   residue 26, 29-30, 31 or 33 of RL32 on Eclip_1480) is no evidence of another proteoform. Requiring the same positions and
  //   names counted one proteoform up to four times (one in five proteoform rows of the lysate runs were such repeats);
  // - otherwise (an undetermined mass): every terminus known in both is the same. The modifications of an undetermined PrSM are
  //   not confirmed by a precursor mass, and an open terminus may be any position, so such a PrSM is a separate proteoform only if
  //   no other PrSM of the protein fits its known termini. Counting it apart whenever its blind shifts or its open end differed
  //   counted a determined proteoform again (232 of 511 undetermined proteoform rows of Eclip_1480 shared their known termini with a
  //   determined proteoform).
  if (a == b) return true;
  const double mass1 = a.getMetaValue(tntKey("Mass")), mass2 = b.getMetaValue(tntKey("Mass"));
  const int sp1 = a.getMetaValue(tntKey("StartPosition")), ep1 = a.getMetaValue(tntKey("EndPosition"));
  const int sp2 = b.getMetaValue(tntKey("StartPosition")), ep2 = b.getMetaValue(tntKey("EndPosition"));
  if (mass1 > 0 && mass2 > 0)
  {
    const bool overlap = (sp1 < 0 || ep2 < 0 || sp1 <= ep2) && (sp2 < 0 || ep1 < 0 || sp2 <= ep1);
    const double diff = std::abs(mass1 - mass2), max_diff = std::max(mass1, mass2) * tol * 2e-6;
    return overlap && (diff <= max_diff || std::abs(diff - Constants::ISOTOPE_MASSDIFF_55K_U) <= max_diff);
  }
  return (sp1 < 0 || sp2 < 0 || sp1 == sp2) && (ep1 < 0 || ep2 < 0 || ep1 == ep2);
}

std::vector<std::vector<Size>> FLASHTnTAlgorithm::markRepresentativeProteoformHits_(double tol, const std::vector<int>& hit_class,
                                                                                   const std::vector<double>& qvalues, double rep_q, bool by_qvalue)
{
  std::map<std::string, std::vector<int>> proteoform_map;
  const int hit_count = (int)proteoform_hits_.size();
  std::vector<bool> passing_determined(hit_count, false), with_mass(hit_count, false);

  for (int i = 0; i < hit_count; i++)
  {
    auto& hit = proteoform_hits_[i];
    hit.setMetaValue(tntKey("Index"), i);
    with_mass[i] = (double)hit.getMetaValue(tntKey("Mass")) > 0;
    if (by_qvalue)
    {
      // A PrSM in no FDR class (undetermined, with discard_underdetermined) has no q-value and is in no proteoform. Sorted as q = 0, it
      // became the representative of determined PrSMs, and their proteoform was lost when the output filter dropped it.
      if (hit_class[i] < 0) continue;
      passing_determined[i] = qvalues[i] <= rep_q && isDetermined(hit);
    }
    proteoform_map[hit.getAccession()].push_back(i);
  }

  std::vector<std::vector<Size>> clusters;
  for (auto& [acc, hit_idx] : proteoform_map)
  {
    // The representative of a proteoform is its row and is tested for consistency with the other PrSMs. It is a determined PrSM passing
    // rep_q if there is one (rather than an undetermined PrSM of another class with a lower q-value), else the PrSM with the lowest
    // q-value; then one with a mass, then the best scoring; then the order of proteoform_hits_. So a representative passes the PrSM
    // filter if one of its PrSMs does, and a proteoform is not lost while one of its PrSMs is reported.
    std::stable_sort(hit_idx.begin(), hit_idx.end(), [&](int l, int r) -> bool {
      if (passing_determined[l] != passing_determined[r]) return passing_determined[l];
      if (by_qvalue && qvalues[l] != qvalues[r]) return qvalues[l] < qvalues[r];
      if (with_mass[l] != with_mass[r]) return with_mass[l];
      return proteoform_hits_[l].getScore() > proteoform_hits_[r].getScore();
    });
    for (Size i = 0; i < hit_idx.size(); i++)
    {
      auto& hit = proteoform_hits_[hit_idx[i]];
      if (hit.metaValueExists(tntKey("Representative"))) continue;
      std::vector<int> prsm_indices {(int)hit.getMetaValue(tntKey("Index"))};
      for (Size j = i + 1; j < hit_idx.size(); j++)
      {
        auto& hit2 = proteoform_hits_[hit_idx[j]];
        if (hit2.metaValueExists(tntKey("Representative"))) continue;
        // a decoy precursor PrSM carries the target accession: merged into a target proteoform, the decoy would no longer be counted
        if (hit.metaValueExists(tntKey("DecoyPrecursor")) != hit2.metaValueExists(tntKey("DecoyPrecursor"))) continue;
        if (areConsistent_(hit, hit2, tol))
        {
          prsm_indices.push_back(hit2.getMetaValue(tntKey("Index")));
          hit2.setMetaValue(tntKey("Representative"), 0);
        }
      }
      clusters.emplace_back(prsm_indices.begin(), prsm_indices.end());
      // The row that describes the proteoform: the PrSM with the fewest blind shifts among its determined PrSMs passing rep_q (then in the
      // order above). A blind shift equal to the residues a moved terminus left out explains the same mass worse (Eclip_1480 SODM: 2-206
      // unmodified and 6-206 +464.215, the mass of residues 2-5); which one represented the proteoform followed q-values differing in the
      // fourth digit. The cluster and its q-value stay with the representative chosen above (clusters keeps it first).
      Size row = 0;
      if (passing_determined[hit_idx[i]])
      {
        for (Size k = 1; k < prsm_indices.size(); k++)
        {
          const int m = prsm_indices[k];
          if (passing_determined[m] && blindModCount_(proteoform_hits_[m]) < blindModCount_(proteoform_hits_[prsm_indices[row]])) row = k;
        }
      }
      std::swap(prsm_indices[0], prsm_indices[row]);
      for (Size k = 0; k < prsm_indices.size(); k++) proteoform_hits_[prsm_indices[k]].setMetaValue(tntKey("Representative"), k == 0 ? 1 : 0);
      proteoform_hits_[prsm_indices[0]].setMetaValue(tntKey("PrSMIndices"), prsm_indices);
    }
  }
  // in the order of the representatives in proteoform_hits_
  std::sort(clusters.begin(), clusters.end(), [](const std::vector<Size>& l, const std::vector<Size>& r) { return l[0] < r[0]; });
  return clusters;
}

void FLASHTnTAlgorithm::setProteinQvalues_(const std::vector<int>& hit_class, const std::vector<double>& qvalues, const std::vector<bool>& decoys,
                                           const std::vector<double>& decoy_weights, const std::vector<std::vector<Size>>& clusters, double rep_q,
                                           const ProteoformDatabase& db)
{
  // A protein is ranked by the evidence of all its proteoforms: the sum over them of the log of the lowest local FDR among the proteoform's
  // PrSMs, then by its best score among its passing PrSMs (rep_q; among all of them if none passes); proteins are then target-decoy counted.
  // The sum is a ranking statistic, not a probability, and is not reported. The local FDR of a PrSM comes from a two-group model of its
  // class: the null p-value of its score under the decoy scores of the class (DecoyScoreNull; the decoys of all classes pooled for a class
  // with fewer than null_min_class_decoys, as the determined classes with 0 to 3 decoys each in the E. coli lysate runs), the density of
  // the targets' p-values (TargetPvalueDensity) and the share of false targets pi0 = (decoy weight + 1) / targets. A decoy precursor
  // carries the target accession: it is a decoy of the precursor, not of the protein, and ranks no protein.
  //
  // A determined PrSM that is not the whole protein (an internal terminus: start > 2 or end before the protein's last residue), and every
  // precursor-ladder PrSM, is scored here with the model of the class "undetermined, both termini" (its null, pi0 and target density; that
  // class itself is unchanged), not with that of its determined class. Its mass agreement is closed by where the free terminus is put
  // (many windows of a long protein fit any mass), so it is no more evidence than termini without a mass; the determined classes have
  // 0 to 2 decoys, so their pi0 is about 1 / targets and their null the pooled one, and such a PrSM got a local FDR of about 0.002 at any
  // score. Over five E. coli lysate runs, 63 of 65 chance determined PrSMs (decoys and B. theta entrapment) were truncated; of the
  // determined PrSMs scoring 85 or less, whole-protein ones were 533 E. coli, 0 B. theta and 0 decoys. On Eclip_1480 against E. coli +
  // B. theta, the false proteins at 1% went from 5 to 4 (all four homologs of E. coli proteins, which shuffled decoys cannot model) and
  // the E. coli-only search from 299 to 292 proteins; on the held-out search against E. coli + human, from 2 false proteins to 0 (287 ->
  // 289 E. coli proteins). The PrSM and proteoform q-values are unchanged.
  const Size hit_count = proteoform_hits_.size();
  const int protein_class_count = fdr_class_count + 1, c3_model = fdr_class_count;
  std::vector<int> protein_class(hit_count, -1);
  for (Size i = 0; i < hit_count; i++)
  {
    const int c = hit_class[i];
    if (c < 0) continue;
    const auto& hit = proteoform_hits_[i];
    const int start = hit.getMetaValue("StartPosition"), end = hit.getMetaValue("EndPosition");
    // the whole protein: from its first residue (or the second, initiator Met removed) to its last; of a PEFF entry (EntryStart), a
    // processed form of its protein (whichever entry the PrSM came from), the second residue only of a form that starts at residue 1
    bool whole_protein = isDetermined(hit) && start >= 1 && start <= 2 && end == (int)hit.getSequence().size();
    if (hit.metaValueExists("EntryStart"))
    {
      whole_protein = false;
      for (const auto& [form_start, form_end] : db.proteins[db.proteinOfEntry()[(Size)(int)hit.getMetaValue("FastaIndex")]].forms)
      {
        if (end == form_end && (start == form_start || (form_start == 1 && start == 2))) whole_protein = isDetermined(hit);
      }
    }
    // protein_open_terminus_model: an open-terminus PrSM too (its own class has almost no decoys at the top, so its own model gave a single
    // such PrSM a local FDR near 0 and carried one-PrSM proteins: 10 of 29 shuffled-entrapment proteins at 1% on the development runs)
    protein_class[i] = (c == fdr_class_precursor_ladder || (c <= fdr_class_max_blind && ! whole_protein) || (protein_open_terminus_model_ && c == fdr_class_open_terminus))
                         ? c3_model : c;
  }

  // null of the PrSM score per class
  std::vector<std::vector<double>> class_decoy_scores(protein_class_count);
  std::vector<double> class_decoy_weight(protein_class_count, 0.0), pooled_decoy_scores;
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] < 0 || ! decoys[i]) continue;
    class_decoy_scores[protein_class[i]].push_back(proteoform_hits_[i].getScore());
    class_decoy_weight[protein_class[i]] += decoy_weights[i]; // a decoy precursor of weight 0 still enters the null's sample
    pooled_decoy_scores.push_back(proteoform_hits_[i].getScore());
  }
  std::vector<DecoyScoreNull> nulls; // nulls[0]: the pooled one
  std::vector<Size> class_null(protein_class_count, 0);
  nulls.reserve(protein_class_count + 1);
  nulls.emplace_back(std::move(pooled_decoy_scores));
  for (int c = 0; c < protein_class_count; c++)
  {
    if (class_decoy_scores[c].size() < null_min_class_decoys) continue;
    class_null[c] = nulls.size();
    nulls.emplace_back(std::move(class_decoy_scores[c]));
  }
  std::vector<double> pvalues(hit_count, 1.0);
  std::vector<std::vector<double>> class_target_pvalues(protein_class_count);
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] < 0) continue;
    pvalues[i] = nulls[class_null[protein_class[i]]].pValue(proteoform_hits_[i].getScore());
    if (! decoys[i]) class_target_pvalues[protein_class[i]].push_back(pvalues[i]);
  }

  // two-group local FDR of every PrSM, targets and decoys, in its class
  std::vector<double> pi0;
  std::vector<TargetPvalueDensity> densities;
  densities.reserve(protein_class_count);
  for (int c = 0; c < protein_class_count; c++)
  {
    pi0.push_back(std::min(1.0, (class_decoy_weight[c] + 1.0) / std::max(1.0, (double)class_target_pvalues[c].size())));
    densities.emplace_back(std::move(class_target_pvalues[c]));
  }
  // the PrSMs scored with the model of "undetermined, both termini"
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] == c3_model) pvalues[i] = nulls[class_null[fdr_class_both_termini]].pValue(proteoform_hits_[i].getScore());
  }
  std::vector<double> local_fdr(hit_count, 1.0);
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] < 0) continue;
    const int model = protein_class[i] == c3_model ? fdr_class_both_termini : protein_class[i];
    const double density = densities[model].density(pvalues[i]); // +inf gives 0
    local_fdr[i] = density <= 0 ? 1.0 : std::min(1.0, pi0[model] / density);
  }

  // protein key: (evidence, -best score), smaller is better
  std::map<std::string, std::vector<double>> evidence_terms; // accession -> log local FDR of each of its proteoforms
  for (const auto& members : clusters)
  {
    const auto& rep = proteoform_hits_[members[0]];
    if (rep.metaValueExists(tntKey("DecoyPrecursor"))) continue;
    double min_local_fdr = local_fdr[members[0]];
    for (Size i : members) min_local_fdr = std::min(min_local_fdr, local_fdr[i]);
    evidence_terms[rep.getAccession()].push_back(std::log(std::max(fdr_p_floor, min_local_fdr)));
  }
  std::map<std::string, double> best_score, best_passing_score;
  for (Size i = 0; i < hit_count; i++)
  {
    const auto& hit = proteoform_hits_[i];
    if (hit_class[i] < 0 || hit.metaValueExists(tntKey("DecoyPrecursor"))) continue;
    const auto it = best_score.find(hit.getAccession());
    if (it == best_score.end() || it->second < hit.getScore()) best_score[hit.getAccession()] = hit.getScore();
    if (qvalues[i] > rep_q) continue;
    const auto itp = best_passing_score.find(hit.getAccession());
    if (itp == best_passing_score.end() || itp->second < hit.getScore()) best_passing_score[hit.getAccession()] = hit.getScore();
  }
  std::vector<std::pair<std::pair<double, double>, std::string>> ranked;
  ranked.reserve(best_score.size());
  for (const auto& [acc, best] : best_score)
  {
    double evidence = 0;
    if (const auto ite = evidence_terms.find(acc); ite != evidence_terms.end())
    {
      auto& terms = ite->second;
      std::sort(terms.begin(), terms.end()); // summed from the most negative on, so that the sum does not depend on the proteoform order
      for (double t : terms) evidence += t;
    }
    const auto it = best_passing_score.find(acc);
    ranked.emplace_back(std::make_pair(evidence, -(it == best_passing_score.end() ? best : it->second)), acc);
  }
  std::sort(ranked.begin(), ranked.end());

  // protein q-values; an accession with decoy precursor PrSMs only is not ranked and gets 1
  std::vector<RankedItem<std::pair<double, double>>> items;
  items.reserve(ranked.size());
  const double protein_decoy_weight = 1.0 / std::max(decoy_factor_, 1e-9);
  for (const auto& [key, acc] : ranked) items.push_back({key, StringUtils::hasPrefix(acc, "DECOY"), protein_decoy_weight});
  const auto protein_qvalues = targetDecoyQvalues(items);
  std::map<std::string, double> protein_q;
  for (Size k = 0; k < ranked.size(); k++) protein_q[ranked[k].second] = protein_qvalues[k];
  for (auto& hit : proteoform_hits_)
  {
    const auto it = protein_q.find(hit.getAccession());
    hit.setMetaValue(tntKey("ProteinQvalue"), it == protein_q.end() ? 1.0 : it->second);
  }
}

void FLASHTnTAlgorithm::vectorizeProteinSequence_(const std::vector<std::string>& cleaned_protein_seqs,
                                                     std::vector<std::unordered_set<int>>& vec_pro,
                                                     std::vector<std::unordered_set<int>>& rev_vec_pro,
                                                     const std::vector<std::vector<std::pair<int, double>>>* fixed)
{
  // proteins in parallel; each set is a copy of one built by the same insertions, as push_back copied it, so that it iterates in the same
  // order (FLASHExtenderAlgorithm::getBestMassShift_ breaks ties by that order)
  const Size first = vec_pro.size(), rev_first = rev_vec_pro.size();
  vec_pro.resize(first + cleaned_protein_seqs.size());
  rev_vec_pro.resize(rev_first + cleaned_protein_seqs.size());
#pragma omp parallel for schedule(dynamic, 64)
  for (SignedSize k = 0; k < (SignedSize)cleaned_protein_seqs.size(); k++)
  { //ignore X
    const auto& seq = cleaned_protein_seqs[k];
    std::unordered_set<int> vec, rev_vec;
    double mass = 0;
    vec.insert(0);
    rev_vec.insert(0);
    // the fixed modification masses by residue (empty: none)
    std::vector<double> extra;
    if (fixed != nullptr && ! (*fixed)[k].empty())
    {
      extra.assign(seq.size(), 0.0);
      for (const auto& [p, m] : (*fixed)[k])
      {
        if (p >= 0 && p < (int)seq.size()) extra[p] += m;
      }
    }

    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = j;
      mass += ResidueDB::getInstance()->getResidue(seq[index])->getMonoWeight(Residue::Internal);
      if (! extra.empty()) mass += extra[index];
      int vindex = int(round(mass));
      vec.insert(vindex);
    }
    mass = 0;
    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = seq.size() - j - 1;
      mass += ResidueDB::getInstance()->getResidue(seq[index])->getMonoWeight(Residue::Internal);
      if (! extra.empty()) mass += extra[index];
      int vindex = int(round(mass));
      rev_vec.insert(vindex);
    }
    vec_pro[first + k] = std::unordered_set<int>(vec);
    rev_vec_pro[rev_first + k] = std::unordered_set<int>(rev_vec);
  }
}

namespace
{
  /// free the elements of @p v with all threads (each moved out and destroyed by the thread that took it), then @p v itself
  template<typename V>
  void releaseInParallel(V& v)
  {
    using T = typename V::value_type;
#pragma omp parallel for schedule(dynamic, 256)
    for (SignedSize i = 0; i < (SignedSize)v.size(); i++)
    {
      T released(std::move(v[i]));
    }
    V().swap(v);
  }
} // namespace

void FLASHTnTAlgorithm::run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& original_fasta_entry, const ProteoformDatabase* db)
{
  setLogType(CMD);

  // The database (ProteoformDatabase): without one, the FASTA entries are the proteins (the identity database), and every step below
  // that reads it is the search as before. With a PEFF database, an entry is a processed form of a protein with fixed modifications: the
  // tag hits, the residue sums and the extension add the fixed masses; the candidates are one entry per protein; the protein-level gate
  // and the chimeric exclusion work by protein; and a PrSM is decoded to its protein right after the extension.
  ProteoformDatabase identity_db;
  if (db == nullptr)
  {
    identity_db = ProteoformDatabase::fromFASTA(original_fasta_entry);
    db = &identity_db;
  }
  else if (db->entries.size() != original_fasta_entry.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "the database has another number of entries than the FASTA entries",
                                  StringUtils::toStr(db->entries.size()));
  }
  const bool by_entry = ! db->isIdentity(); // several entries per protein, or fixed modifications
  const auto& protein_of = db->proteinOfEntry();
  const bool with_fixed = db->hasFixedModifications();

  // per distinct tag sequence (tag_id_of: its position in tag_db): the proteins (index in fasta), positions and flanking masses of its hits
  struct TagDbHits
  {
    std::vector<Size> proteins, positions;
    std::vector<double> n_flanking_masses, c_flanking_masses;
  };
  std::unordered_map<std::string, Size> tag_id_of;
  std::vector<TagDbHits> tag_db;

  std::map<int, std::vector<Size>> scan_to_tag_indices;         // scan number to tag indices
  std::vector<std::string> tag_seqs, protein_seqs, cleaned_protein_seqs;
  std::vector<DeconvolvedSpectrum> dspecs;
  // per spectrum in dspecs: the tolerance of its deconvolution and, for the chimeric search only, all its masses (ascending) if the cap
  // max_node_cntr_ left some out (else empty: dspecs has them all)
  std::vector<double> spectrum_tols;
  std::vector<std::vector<PeakGroup>> uncapped_masses;

  std::vector<FASTAFile::FASTAEntry> fasta_entry;
  fasta_entry.reserve(original_fasta_entry.size());
  for (const auto& fe : original_fasta_entry)
  {
    auto seq = fe.sequence;
    std::replace(seq.begin(), seq.end(), 'I', 'L');
    fasta_entry.emplace_back(fe.identifier, fe.description, seq);
  }

  {
    protein_seqs.reserve(fasta_entry.size());
    cleaned_protein_seqs.reserve(fasta_entry.size());
    for (const auto& fe : fasta_entry)
    {
      protein_seqs.push_back(fe.sequence);
      std::string cleaned_seq;
      std::remove_copy(fe.sequence.begin(), fe.sequence.end(),
                       std::back_inserter(cleaned_seq), 'X');
      cleaned_protein_seqs.push_back(cleaned_seq);
    }
    // decoys per target by protein (a FASTA entry is a protein), so that the entries of a PEFF protein count once
    double taget_count = 0;
    double decoy_count = 0;
    for (const auto& protein : db->proteins)
    {
      if (StringUtils::hasPrefix(protein.identifier, "DECOY")) { decoy_count++; }
      else { taget_count++; }
    }

    decoy_factor_ = decoy_count / taget_count;
  }
  if (by_entry)
  {
    Size with_mods = 0, decoys = 0;
    for (const auto& entry : db->entries) with_mods += ! entry.fixed.empty();
    for (const auto& protein : db->proteins) decoys += StringUtils::hasPrefix(protein.identifier, "DECOY");
    OPENMS_LOG_INFO << "Database: " << db->proteins.size() - decoys << " target and " << decoys << " decoy proteins, " << db->entries.size()
                    << " entries (" << with_mods << " with fixed modifications)" << std::endl;
  }
  // the full-length protein sequences (I as L), to which a PrSM of an entry is decoded
  std::vector<std::string> protein_full_seqs;
  if (by_entry)
  {
    protein_full_seqs.reserve(db->proteins.size());
    for (const auto& protein : db->proteins)
    {
      std::string seq = protein.sequence;
      std::replace(seq.begin(), seq.end(), 'I', 'L');
      protein_full_seqs.push_back(seq);
    }
  }
  // the fixed modifications of each entry in its sequence (fixedMasses) and in its sequence without X (cleaned_protein_seqs)
  const std::vector<std::vector<std::pair<int, double>>>* entry_fixed = with_fixed ? &db->fixedMasses() : nullptr;
  std::vector<std::vector<std::pair<int, double>>> cleaned_fixed;
  if (with_fixed)
  {
    cleaned_fixed.resize(fasta_entry.size());
    for (Size e = 0; e < fasta_entry.size(); e++)
    {
      const auto& seq = fasta_entry[e].sequence;
      for (const auto& [p, mass] : (*entry_fixed)[e])
      {
        if (p < 0 || p >= (int)seq.size() || seq[p] == 'X') continue;
        cleaned_fixed[e].emplace_back(p - (int)std::count(seq.begin(), seq.begin() + p, 'X'), mass);
      }
    }
  }

  // the candidate blind modifications: those of ModificationsDB within the largest total blind shift
  std::map<double, std::vector<ResidueModification>> blind_mod_map;
  {
    const auto inst = ModificationsDB::getInstance();
    std::vector<std::string> mod_strs;
    inst->getAllSearchModifications(mod_strs);
    for (const auto& mod_str : mod_strs)
    {
      const auto mod = *inst->getModification(mod_str);
      if (std::abs(mod.getDiffMonoMass()) > max_total_mod_mass_) continue;
      blind_mod_map[mod.getDiffMonoMass()].push_back(mod);
    }
  }

  double precursor_tol = -1, tol, ms2_tol = -1;

  // Only MSn spectra are searched: the masses of an MS1 spectrum are whole proteoforms, not fragments, and its PrSMs (2208 of 10495 on
  // Eclip_1480, 4 of them passing 1%) only fed the decoy counts of the FDR. A file without any MSn spectrum is searched as it is.
  bool has_msn = false;
  for (const auto& spec : map)
  {
    if (spec.getMSLevel() > 1)
    {
      has_msn = true;
      break;
    }
  }

  startProgress(0, (SignedSize)map.size(), "Finding sequence tags ...");

  // Spectra in parallel, each into its own slot, collected in input order below, so that dspecs, the tolerances and tags_ are those of one
  // pass in input order. An exception (e.g. a malformed DeconvMassInfo) is rethrown for the first spectrum in input order that raised one.
  struct TaggedSpectrum
  {
    bool with_masses = false; // at least 5 masses: ms_level and tol are read
    UInt ms_level = 0;
    double tol = 0;
    bool searched = false; // dspec, all_masses and tags are set
    DeconvolvedSpectrum dspec;
    std::vector<PeakGroup> all_masses;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::exception_ptr error;
  };
  std::vector<TaggedSpectrum> tagged(map.size());
  auto tagSpectrum = [&](int index, TaggedSpectrum& out) {
    const auto& spec = map[index];
    nextProgress();

    if (spec.size() < 5) return;

    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);

    DeconvolvedSpectrum dspec(scan);
    dspec.setOriginalSpectrum(spec);
    std::string deconv_meta_str = spec.getMetaValue(tntKey("DeconvMassInfo")).toString();

    int tol_loc_s = deconv_meta_str.find("tol=") + 4;
    int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);

    const double spec_tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));
    out.with_masses = true;
    out.ms_level = spec.getMSLevel();
    out.tol = spec_tol;
    if (spec.getMSLevel() == 1 && has_msn) return;

    int q_loc_s = deconv_meta_str.find("qscore=") + 7;
    int q_loc_e = deconv_meta_str.find(";", q_loc_s);
    auto q_str = deconv_meta_str.substr(q_loc_s, q_loc_e - q_loc_s);
    Size pos = 0;
    std::vector<double> qscores;
    while (true)
    {
      Size pos_t = q_str.find(",", pos);
      if (pos_t == std::string::npos) break;
      auto token = q_str.substr(pos, pos_t - pos);
      qscores.push_back(stod(token));
      pos = pos_t + 1;
    }

    int s_loc_s = deconv_meta_str.find("snr=") + 4;
    int s_loc_e = deconv_meta_str.find(";", s_loc_s);
    auto s_str = deconv_meta_str.substr(s_loc_s, s_loc_e - s_loc_s);
    pos = 0;
    std::vector<float> snrs;
    while (true)
    {
      Size pos_t = s_str.find(",", pos);
      if (pos_t == std::string::npos) break;
      auto token = s_str.substr(pos, pos_t - pos);
      snrs.push_back(stof(token));
      pos = pos_t + 1;
    }

    int s_loc_pre_s = deconv_meta_str.find("precursorscan=") + 14;
    int s_loc_pre_e = deconv_meta_str.find(";", s_loc_pre_s);
    int precursor_scan = stoi(deconv_meta_str.substr(s_loc_pre_s, s_loc_pre_e - s_loc_pre_s));

    double precursor_mass = 0;
    if (StringUtils::hasSubstring(deconv_meta_str, "precursormass="))
    {
      int s_loc_prem_s = deconv_meta_str.find("precursormass=") + 14;
      int s_loc_prem_e = deconv_meta_str.find(";", s_loc_prem_s);
      precursor_mass = stod(deconv_meta_str.substr(s_loc_prem_s, s_loc_prem_e - s_loc_prem_s));
    }

    if (StringUtils::hasSubstring(deconv_meta_str, "decoyprecursormass="))
    {
      int s_loc_dm_s = deconv_meta_str.find("decoyprecursormass=") + 19;
      int s_loc_dm_e = deconv_meta_str.find(";", s_loc_dm_s);
      dspec.setDecoyPrecursorMass(stod(deconv_meta_str.substr(s_loc_dm_s, s_loc_dm_e == (int)std::string::npos ? std::string::npos : s_loc_dm_e - s_loc_dm_s)));
      if (StringUtils::hasSubstring(deconv_meta_str, "decoyprecursorweight="))
      {
        int s_loc_dw_s = deconv_meta_str.find("decoyprecursorweight=") + 21;
        int s_loc_dw_e = deconv_meta_str.find(";", s_loc_dw_s);
        dspec.setDecoyPrecursorWeight(stod(deconv_meta_str.substr(s_loc_dw_s, s_loc_dw_e == (int)std::string::npos ? std::string::npos : s_loc_dw_e - s_loc_dw_s)));
      }
      if (StringUtils::hasSubstring(deconv_meta_str, "decoyprecursorQscore1D="))
      {
        int s_loc_dq1_s = deconv_meta_str.find("decoyprecursorQscore1D=") + 23;
        int s_loc_dq1_e = deconv_meta_str.find(";", s_loc_dq1_s);
        dspec.setDecoyPrecursorQscore1D(stod(deconv_meta_str.substr(s_loc_dq1_s, s_loc_dq1_e == (int)std::string::npos ? std::string::npos : s_loc_dq1_e - s_loc_dq1_s)));
      }
      if (StringUtils::hasSubstring(deconv_meta_str, "decoyprecursorQscore="))
      {
        int s_loc_dq_s = deconv_meta_str.find("decoyprecursorQscore=") + 21;
        int s_loc_dq_e = deconv_meta_str.find(";", s_loc_dq_s);
        dspec.setDecoyPrecursorQscore(stod(deconv_meta_str.substr(s_loc_dq_s, s_loc_dq_e == (int)std::string::npos ? std::string::npos : s_loc_dq_e - s_loc_dq_s)));
      }
    }

    // a precursor mass may be given (e.g. by precursor_charge/precursor_mz options) even if no precursor scan is present
    if (precursor_scan > 0 || precursor_mass > 0)
    {
      PeakGroup pg;
      pg.setMonoisotopicMass(precursor_mass);
      if (StringUtils::hasSubstring(deconv_meta_str, "precursorQscore="))
      {
        int s_loc_preq_s = deconv_meta_str.find("precursorQscore=") + 16;
        int s_loc_preq_e = deconv_meta_str.find(";", s_loc_preq_s);
        double precursor_qscore = stod(deconv_meta_str.substr(s_loc_preq_s, s_loc_preq_e - s_loc_preq_s));
        pg.setQscore2D(precursor_qscore);
      }
      if (StringUtils::hasSubstring(deconv_meta_str, "precursorQscore1D="))
      {
        int s_loc_pq1_s = deconv_meta_str.find("precursorQscore1D=") + 18; // the precursor's comes before the decoy precursor's
        int s_loc_pq1_e = deconv_meta_str.find(";", s_loc_pq1_s);
        pg.setQscore(stod(deconv_meta_str.substr(s_loc_pq1_s, s_loc_pq1_e - s_loc_pq1_s)));
      }
      if (StringUtils::hasSubstring(deconv_meta_str, "precursorSNR="))
      {
        int s_loc_pres_s = deconv_meta_str.find("precursorSNR=") + 13;
        int s_loc_pres_e = deconv_meta_str.find(";", s_loc_pres_s);
        double precursor_snr = stod(deconv_meta_str.substr(s_loc_pres_s, s_loc_pres_e - s_loc_pres_s));
        pg.setSNR(precursor_snr);
      }

      dspec.setPrecursorPeakGroup(pg);
    }
    for (Size i = 0; i < spec.size(); i++)
    {
      PeakGroup peak;
      peak.setQscore(qscores[i]);
      peak.setSNR(snrs[i]);
      peak.setMonoisotopicMass(spec[i].getMZ());
      peak.setScanNumber(scan);
      dspec.push_back(peak);
    }

    std::vector<PeakGroup> all_masses; // the chimeric search applies the cap again to what is left of all masses
    if (chimeric_rounds_ > 0 && dspec.size() > max_node_cntr_) all_masses.assign(dspec.begin(), dspec.end());
    if (dspec.size() > max_node_cntr_)
    {
      dspec.sortByQscore2D();
      std::vector<PeakGroup> filtered_peaks;
      filtered_peaks.reserve(max_node_cntr_);
      for (const auto& pg : dspec)
      {
        if (filtered_peaks.size() >= max_node_cntr_) break;
        filtered_peaks.push_back(pg);
      }
      dspec.setPeakGroups(filtered_peaks);
    }

    dspec.sort();
    if (chimeric_rounds_ > 0)
    {
      std::sort(all_masses.begin(), all_masses.end()); // as dspec.sort()
      out.all_masses = std::move(all_masses);
    }
    FLASHTaggerAlgorithm tagger;
    // Run tagger
    tagger.setParameters(tagger_param_);
    tagger.run(dspec, spec_tol);
    tagger.fillTags(out.tags);
    out.dspec = std::move(dspec);
    out.searched = true;
  };
#pragma omp parallel for schedule(dynamic)
  for (int index = 0; index < (int)map.size(); index++)
  {
    try
    {
      tagSpectrum(index, tagged[index]);
    }
    catch (...)
    {
      tagged[index].error = std::current_exception();
    }
  }
  {
    Size tag_total = tags_.size();
    for (const auto& out : tagged) tag_total += out.tags.size();
    tags_.reserve(tag_total);
  }
  for (auto& out : tagged)
  {
    if (out.error) std::rethrow_exception(out.error);
    if (! out.with_masses) continue;
    tol = out.tol;
    if (out.ms_level > 1 && ms2_tol < 0) ms2_tol = tol;
    if (out.ms_level == 1 && precursor_tol < 0) precursor_tol = tol;
    if (! out.searched) continue;
    dspecs.push_back(std::move(out.dspec));
    spectrum_tols.push_back(tol);
    if (chimeric_rounds_ > 0) uncapped_masses.push_back(std::move(out.all_masses));
    for (auto& tag : out.tags) tags_.push_back(std::move(tag));
    out = TaggedSpectrum();
  }

  endProgress();

  for (Size i = 0; i < tags_.size(); i++)
  {
    const auto& seq = tags_[i].getUppercaseSequence();
    scan_to_tag_indices[tags_[i].getScan()].push_back(i);
    tags_[i].setIndex((int)i);
    if (! tag_id_of.emplace(seq, tag_db.size()).second) continue;
    tag_db.emplace_back();
    tag_seqs.push_back(seq);
  } // now we have all the tags

  // tag_seqs protein_seqs

  // Internal residue masses by one-letter code, for the flanking masses of the tag hits. AASequence::fromString() of a string of these
  // letters is one residue per letter (ResidueDB::getResidue(char)), and AASequence::getMonoWeight(Residue::Internal) adds their
  // masses to 0, from the first residue to the last; the same additions in the same order give the same doubles. A protein with any
  // other character (e.g. X, which cleaned_protein_seqs drops, so that its positions are no longer those of the hits) keeps the
  // AASequence computation.
  std::array<double, 256> aa_internal_mass {};
  std::array<bool, 256> aa_plain {};
  {
    const ResidueDB* rdb = ResidueDB::getInstance();
    const Residue* unknown = rdb->getResidue(std::string("X"));
    for (int c = 'A'; c <= 'Z'; c++)
    {
      if (c == 'X') continue;
      const Residue* r = rdb->getResidue((unsigned char)c);
      if (r == nullptr || r == unknown) continue;
      aa_plain[c] = true;
      aa_internal_mass[c] = r->getMonoWeight(Residue::Internal);
    }
  }

  // The database positions and flanking masses of the tag sequences @p seqs (their entries in tag_db exist, empty): an Aho-Corasick
  // trie of them over the proteins. Also used for the new tag sequences of the chimeric search. Proteins run in parallel in blocks (each
  // thread with its own trie state); the hits of a block are appended to tag_db in protein order, so every vector is the same as that
  // of one pass over the proteins. On an entry with fixed modifications (ProteoformDatabase), a tag (de novo: unmodified residues) does not sit
  // across a fixed residue, and the fixed masses before and after it are part of its flanking masses.
  auto searchTagSequences = [&](const std::vector<std::string>& seqs) {
    ACTrie ac_trie(0, 0); // no ambiguities, no mismatches
    ac_trie.addNeedlesAndCompress(seqs);

    // the hits of every tag sequence (tag_db does not grow during the search)
    std::vector<TagDbHits*> hits_of(seqs.size());
    for (Size j = 0; j < seqs.size(); j++) hits_of[j] = &tag_db[tag_id_of.at(seqs[j])];

    struct TagHit
    {
      OpenMS::Hit::T needle_index;
      OpenMS::Hit::T query_pos;
      bool summed; // false: the flanking masses are computed with AASequence when the hit is appended
      double nmass;
      double cmass;
      Size protein;
    };

    startProgress(0, (SignedSize)protein_seqs.size(), "Searching tags against database ...");

    const Size block_size = 4096;
    std::vector<std::vector<TagHit>> block_hits(block_size);
    Size owners = 1;
#ifdef _OPENMP
    owners = (Size)omp_get_max_threads();
#endif
    std::vector<Size> owner_offset(block_size * owners), owner_start(owners + 1);
    std::vector<TagHit> by_owner;
    for (Size block_start = 0; block_start < protein_seqs.size(); block_start += block_size)
    {
      const Size block_end = std::min(protein_seqs.size(), block_start + block_size);
#pragma omp parallel
      {
        OpenMS::ACTrieState ac_state;
        std::vector<double> prefix; // residue sums from the N-terminus
#pragma omp for schedule(dynamic, 16)
        for (SignedSize i = (SignedSize)block_start; i < (SignedSize)block_end; i++) // // Natural19WithoutI
        {
          auto& hits = block_hits[i - block_start];
          hits.clear();
          const std::string& pseq = protein_seqs[i];
          ac_state.setQuery(pseq);
          ac_trie.getAllHits(ac_state);
          if (ac_state.hits.empty()) continue;
          hits.reserve(ac_state.hits.size());

          bool plain = true;
          for (const char c : pseq)
          {
            if (! aa_plain[(unsigned char)c])
            {
              plain = false;
              break;
            }
          }
          const Size n = pseq.size();
          if (plain)
          {
            prefix.resize(n + 1);
            prefix[0] = .0;
            for (Size k = 0; k < n; k++) prefix[k + 1] = prefix[k] + aa_internal_mass[(unsigned char)pseq[k]];
          }
          const std::vector<std::pair<int, double>>* fixed = entry_fixed != nullptr && ! (*entry_fixed)[i].empty() ? &(*entry_fixed)[i] : nullptr;
          for (const auto& h : ac_state.hits)
          {
            const Size end = (Size)h.query_pos + seqs[h.needle_index].length();
            // the fixed masses before and after the tag; none inside it
            double fixed_n = 0, fixed_c = 0;
            if (fixed != nullptr)
            {
              bool inside = false;
              for (const auto& [p, mass] : *fixed)
              {
                if ((Size)p < (Size)h.query_pos) fixed_n += mass;
                else if ((Size)p >= end) fixed_c += mass;
                else inside = true;
              }
              if (inside) continue;
            }
            if (! plain || end > n)
            {
              hits.push_back({h.needle_index, h.query_pos, false, fixed_n, fixed_c, (Size)i});
              continue;
            }
            double nmass = h.query_pos <= 0 ? .0 : prefix[h.query_pos];
            double cmass = .0;
            for (Size k = end; k < n; k++) cmass += aa_internal_mass[(unsigned char)pseq[k]];
            if (fixed != nullptr)
            {
              nmass += fixed_n;
              cmass += fixed_c;
            }
            hits.push_back({h.needle_index, h.query_pos, true, nmass, cmass, (Size)i});
          }
        }
      }

      // The hits of the block go to the vectors of their sequences in protein order. Each thread appends those of a range of sequences
      // (owner): the hits are first sorted by owner (counting sort, stable in protein order) so that each thread reads only its own.
      const Size block_n = block_end - block_start;
      for (Size i = 0; i < block_n; i++) nextProgress();
      std::fill(owner_offset.begin(), owner_offset.end(), 0);
      auto ownerOf = [&](OpenMS::Hit::T needle) { return (Size)((double)needle * owners / seqs.size()); };
#pragma omp parallel for schedule(dynamic, 64)
      for (SignedSize i = 0; i < (SignedSize)block_n; i++)
      {
        for (const auto& h : block_hits[i]) owner_offset[i * owners + ownerOf(h.needle_index)]++;
      }
      Size total = 0;
      for (Size o = 0; o < owners; o++)
      {
        owner_start[o] = total;
        for (Size i = 0; i < block_n; i++)
        {
          const Size count = owner_offset[i * owners + o];
          owner_offset[i * owners + o] = total;
          total += count;
        }
      }
      owner_start[owners] = total;
      if (by_owner.size() < total) by_owner.resize(total);
#pragma omp parallel for schedule(dynamic, 64)
      for (SignedSize i = 0; i < (SignedSize)block_n; i++)
      {
        for (const auto& h : block_hits[i]) by_owner[owner_offset[i * owners + ownerOf(h.needle_index)]++] = h;
      }
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize o = 0; o < (SignedSize)owners; o++)
      {
        for (Size k = owner_start[o]; k < owner_start[o + 1]; k++)
        {
          const auto& h = by_owner[k];
          double nmass = h.nmass, cmass = h.cmass;
          if (! h.summed)
          {
            const std::string& cleanseq = cleaned_protein_seqs[h.protein];
            const std::string& seq = seqs[h.needle_index];
            nmass = h.query_pos <= 0 ? .0 : AASequence::fromString(cleanseq.substr(0, h.query_pos)).getMonoWeight(Residue::Internal);
            cmass = cleanseq.length() - seq.length() <= h.query_pos ? .0 :
              AASequence::fromString(cleanseq.substr(h.query_pos + seq.length(), cleanseq.length() - h.query_pos - seq.length())).getMonoWeight(Residue::Internal);
            if (entry_fixed != nullptr && ! (*entry_fixed)[h.protein].empty())
            {
              nmass += h.nmass; // the fixed masses before and after the tag
              cmass += h.cmass;
            }
          }
          auto& th = *hits_of[h.needle_index];
          th.proteins.push_back(h.protein);
          th.positions.push_back(h.query_pos);
          th.n_flanking_masses.push_back(nmass);
          th.c_flanking_masses.push_back(cmass);
        }
      }
    }

    endProgress();
  };
  searchTagSequences(tag_seqs);
  startProgress(0, (SignedSize)dspecs.size(), "Running candidate protein filtration ...");

  // wall time of the database indexes below (logged after the mass ladders)
  auto index_clock = std::chrono::steady_clock::now();
  auto lap = [&index_clock]() {
    const auto now = std::chrono::steady_clock::now();
    const double s = std::chrono::duration<double>(now - index_clock).count();
    index_clock = now;
    return s;
  };
  std::vector<std::unordered_set<int>> vec_pro, rev_vec_pro;
  vectorizeProteinSequence_(cleaned_protein_seqs, vec_pro, rev_vec_pro, with_fixed ? &cleaned_fixed : nullptr);
  const double t_mass_sets = lap();

  // exact residue sums from the N- and the C-terminus per protein, for the protein-level score (FLASHTaggerAlgorithm::runMatchingExact);
  // with the fixed masses of an entry, whose residues without them (pro_plain_mass; 0: none fixed) set the chance of a match
  std::vector<std::vector<double>> pro_prefix(cleaned_protein_seqs.size()), pro_suffix(cleaned_protein_seqs.size());
  std::vector<double> pro_plain_mass(with_fixed ? cleaned_protein_seqs.size() : 0, 0.0);
#pragma omp parallel for schedule(dynamic, 64)
  for (SignedSize k = 0; k < (SignedSize)cleaned_protein_seqs.size(); k++)
  {
    const auto& seq = cleaned_protein_seqs[k];
    auto& pre = pro_prefix[k];
    auto& suf = pro_suffix[k];
    pre.reserve(seq.size() + 1);
    suf.reserve(seq.size() + 1);
    pre.push_back(0);
    suf.push_back(0);
    if (with_fixed && ! cleaned_fixed[k].empty())
    {
      std::vector<double> extra(seq.size(), 0.0);
      for (const auto& [p, mass] : cleaned_fixed[k]) extra[p] += mass;
      double plain = 0; // as the sum of an entry without fixed modifications
      for (Size j = 0; j < seq.size(); j++)
      {
        const double m = ResidueDB::getInstance()->getResidue(seq[j])->getMonoWeight(Residue::Internal);
        plain += m;
        pre.push_back(pre.back() + m + extra[j]);
      }
      pro_plain_mass[k] = plain;
      for (Size j = 0; j < seq.size(); j++)
        suf.push_back(suf.back() + ResidueDB::getInstance()->getResidue(seq[seq.size() - 1 - j])->getMonoWeight(Residue::Internal) + extra[seq.size() - 1 - j]);
      continue;
    }
    for (Size j = 0; j < seq.size(); j++) pre.push_back(pre.back() + ResidueDB::getInstance()->getResidue(seq[j])->getMonoWeight(Residue::Internal));
    for (Size j = 0; j < seq.size(); j++) suf.push_back(suf.back() + ResidueDB::getInstance()->getResidue(seq[seq.size() - 1 - j])->getMonoWeight(Residue::Internal));
  }
  // their bucket indices, with which FLASHTaggerAlgorithm::scoreShifts looks the residue sum of a mass up instead of walking to it
  std::vector<std::vector<std::uint32_t>> pro_prefix_index(cleaned_protein_seqs.size()), pro_suffix_index(cleaned_protein_seqs.size());
#pragma omp parallel for schedule(dynamic, 64)
  for (int k = 0; k < (int)cleaned_protein_seqs.size(); k++)
  {
    FLASHTaggerAlgorithm::indexResidueSums(pro_prefix[k], pro_prefix_index[k]);
    FLASHTaggerAlgorithm::indexResidueSums(pro_suffix[k], pro_suffix_index[k]);
  }
  const double t_residue_sums = lap();
  const double protein_tol = ms2_tol > 0 ? ms2_tol : tol;
  // the learned candidate cut reads the ion ladders from the residue sums: c ions NH3 above them; z ions NH3 below, z+1 ions NH3 - H
  const auto ion_types = param_.getValue("ion_type").toStringVector();
  auto hasIonType = [&ion_types](const std::string& t) { return std::find(ion_types.begin(), ion_types.end(), t) != ion_types.end(); };
  const double learned_prefix_offset = hasIonType("c") ? 17.026549 : 0.0;
  const double learned_suffix_offset = hasIonType("zp1") ? -16.018724 : (hasIonType("z") ? -17.026549 : 0.0);

  // The protein-level search runs on every spectrum before any extension, because the run-level protein FDR
  // that decides which candidates are extended needs all of it. Each spectrum's candidates stay in
  // protein_hits_, in score order, at [first, second).
  std::vector<std::pair<Size, Size>> candidate_ranges(dspecs.size(), {0, 0});
  // Spectra are independent here, so they run in parallel (the shared maps are only read); each spectrum's candidates are collected
  // and appended to protein_hits_ in spectrum order afterwards, so the result does not depend on the thread schedule.
  std::vector<std::vector<ProteinHit>> spec_hits(dspecs.size());
  // candidates from terminal ladders, kept apart from the others (see the extension); their ranges in ladder_candidates
  std::vector<std::vector<ProteinHit>> spec_ladder_hits(dspecs.size());
  std::vector<ProteinHit> ladder_candidates;
  std::vector<std::pair<Size, Size>> ladder_ranges(dspecs.size(), {0, 0});
  // Meta values set by name look the name up in the global registry under a critical section; register the names once and set
  // them by index inside the parallel loop.
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_scan = registry.registerName("Scan"), mi_matched_aa = registry.registerName("MatchedAA"),
             mi_tag_indices = registry.registerName("TagIndices"), mi_tag_positions = registry.registerName("TagPositions"),
             mi_fasta_index = registry.registerName("FastaIndex"), mi_protein_rank = registry.registerName("ProteinRank"),
             mi_gapped = registry.registerName("GappedTag"), mi_terminal = registry.registerName("TerminalLadder"),
             mi_sliding = registry.registerName("SlidingLadder"),
             mi_precursor_ladder = registry.registerName("PrecursorLadder"),
             mi_pl_start = registry.registerName("PrecursorLadderStart"), mi_pl_end = registry.registerName("PrecursorLadderEnd"),
             mi_pl_isotope = registry.registerName("PrecursorLadderIsotope"), mi_entry_of_kept = registry.registerName("EntryOfKeptProtein"),
             mi_description = registry.registerName("Description");
  // gapped tags (only for spectra in which no sequence tag proposes a protein); their tags are kept per spectrum until appended to tags_
  FLASHGappedTaggerAlgorithm gapped_tagger;
  lap();
  if (gapped_mode_ != "off")
  {
    gapped_tagger.setParameters(gapped_param_);
    gapped_tagger.buildIndex(protein_seqs, entry_fixed);
    OPENMS_LOG_INFO << "Gapped tag index: " << gapped_tagger.getIndexBytes() / (1024 * 1024) << " MB" << std::endl;
  }
  std::vector<std::vector<FLASHHelperClasses::Tag>> spec_gapped_tags(dspecs.size());
  const double t_gapped = lap();

  // Mass ladders (MassLadderAlgorithm): terminal ladders from the residue sums of every protein from its termini, and precursor-anchored
  // ladders from the runs of standard residues of the database, for the window search of every spectrum with a precursor mass. The
  // precursor ladders are built only if a spectrum has a precursor mass.
  // Native mode (native): the candidates of a source in a spectrum are kept when their score exceeds that of the best decoy candidate of
  // the same source in the spectrum by native_margin (decoys thereby never), best first, at most max_protein_hits_. A spectrum is its own
  // null here: native runs have one or a few spectra, too few for an FDR over spectra. Native mode reads ladders from any start or end of a
  // protein (sliding ladders, MassLadderAlgorithm::searchSlidingLadders) instead of the intact termini only: the outer membrane proteins
  // of the targeted test spectra lack their signal peptides (OMPA from residue 22). Gapped-tag candidates were never an expected protein
  // of the native test spectra, so they are not used. Margin 2 (a hundred times the chance probability of the best decoy) was set
  // beforehand. On 25 native and targeted test spectra (OMPA, OMPF, TSX, DLDH and GapA mixtures, DLDH, AQPZ, NISTmAb 1-6, SIGMAmAb,
  // Infliximab, carbonic anhydrase, VAMP2, semiSWEET; 43 expected proteins) against E. coli + the expected proteins, native mode found as
  // many expected proteins as without it (29) and 8 other proteins instead of 27 (6 instead of 8 of them chains of another antibody or
  // porins next to OMPF); against databases without the expected proteins, it called 1 protein instead of 33. The margin is 2.5 since the
  // decoys keep their first two residues (DecoyDatabase): on these spectra with decoy seeds 1-3, null-database calls 1/1/1 -> 0/0/0,
  // expected proteins 35/33/36 -> 34/33/35 of 43, other proteins 6/7/7 -> 5/5/6, the mAb chains 20/20 unchanged. 2.3, just above the best
  // null candidate (2.27), also called no null protein (34/33/36 expected), with less room to the null.
  const bool native_rule = native_ && decoy_factor_ > 0;
  if (native_ && ! native_rule)
  {
    OPENMS_LOG_WARN << "native: the database has no decoys (DECOY prefix); the candidates are chosen as without the option." << std::endl;
  }
  bool any_precursor_mass = false;
  for (const auto& dspec : dspecs) any_precursor_mass |= dspec.getPrecursorPeakGroup().getMonoMass() > 0;
  const bool use_precursor_ladders = precursor_ladders_ && ! protein_level_only_ && any_precursor_mass;
  MassLadderAlgorithm mass_ladders;
  // the sliding ladders of native mode read the runs of standard residues of the precursor ladders
  lap();
  // The ladders read the termini with the terminal forms the extension searches: the initiator Met removal (always) and the terminal
  // variable modifications (var_mod, ex:var_mod_file at an N-terminus, the protein N-terminus or a C-terminus; one residue each), on the
  // residues the extension places them on. Sequence and gapped tags need none: their shift from the protein terminus is read from the data.
  std::vector<MassLadderAlgorithm::TerminalModification> terminal_mods;
  const int max_var_mod_count = (int)param_.getValue("max_var_mod_count");
  if (max_var_mod_count > 0)
  {
    FLASHExtenderAlgorithm extender;
    extender.setParameters(extender_param_);
    const auto names = extender.getVariableModificationNames();
    const auto residues = extender.getVariableModificationResidues();
    const auto positions = extender.getVariableModificationPositions();
    const auto mod_masses = extender.getVariableModificationMasses();
    const auto sites = extender.getVariableModificationSites();
    std::string listed;
    for (Size i = 0; i < names.size(); i++)
    {
      if (positions[i] == 0 || sites[i] != 1) continue;
      terminal_mods.push_back({mod_masses[i], positions[i], residues[i]});
      listed += std::string(listed.empty() ? "" : ", ") + names[i] + " " + StringUtils::number(mod_masses[i], 4) + " Da ("
                + (positions[i] == 3 ? "protein N-term" : positions[i] == 1 ? "N-term" : "C-term") + ", "
                + (residues[i].empty() || residues[i] == "*" ? std::string("any residue") : residues[i]) + ")";
    }
    if (! terminal_mods.empty()) OPENMS_LOG_INFO << "Mass ladders read the termini also with: " << listed << std::endl;
  }
  mass_ladders.buildIndex(protein_seqs, param_.getValue("ion_type").toStringVector(), terminal_ladders_, use_precursor_ladders || native_rule, db,
                          terminal_mods, max_var_mod_count);
  OPENMS_LOG_INFO << "Database indexes (wall s): protein mass sets " << t_mass_sets << ", residue sums " << t_residue_sums << ", gapped tags "
                  << t_gapped << ", mass ladders " << lap() << std::endl;
  // the precursor tolerance of the windows: the MS1 tolerance of the input, else (no MS1 spectra) the fragment tolerance
  const double pl_precursor_tol = precursor_tol > 0 ? precursor_tol : protein_tol;
  // per spectrum: the windows that fit its precursor mass (hypotheses), and the candidates; the same at the decoy precursor mass
  // (decoy_precursor), so that the decoy copy of a spectrum has this source as its target does
  std::vector<std::vector<ProteinHit>> spec_pl_hits(dspecs.size()), spec_decoy_pl_hits(dspecs.size());
  std::vector<ProteinHit> pl_candidates, decoy_pl_candidates;
  std::vector<std::pair<Size, Size>> pl_ranges(dspecs.size(), {0, 0}), decoy_pl_ranges(dspecs.size(), {0, 0});
  std::vector<Size> spec_pl_windows(dspecs.size(), 0);
  // Parallel over spectra only when there are at least as many spectra with tags as threads; otherwise (e.g. a targeted native run with
  // one or a few spectra) the threads go to the scoring of each spectrum's candidates (runMatchingExact).
  int spectra_with_tags = 0;
  for (const auto& dspec : dspecs)
  {
    spectra_with_tags += gapped_mode_ != "off" || terminal_ladders_ || (use_precursor_ladders && dspec.getPrecursorPeakGroup().getMonoMass() > 0)
                         || scan_to_tag_indices.count(dspec.getScanNumber()) > 0;
  }
#ifdef _OPENMP
  const bool parallel_over_spectra = spectra_with_tags >= omp_get_max_threads();
#else
  const bool parallel_over_spectra = false;
#endif
  const Size native_ladder_top = 1000; // ladder candidates ranked per spectrum and source for the decoy margin
  const std::vector<std::string> native_sources {"sequence tags", "sliding ladders", "precursor ladders"};
  std::vector<std::atomic<Size>> native_kept(native_sources.size());
  for (auto& k : native_kept) k = 0;
  // @p truncated: the ranking stopped at native_ladder_top, so a decoy missing from it scored at most the last candidate's score
  auto nativeCut = [&](auto& cands, auto score_of, auto protein_of, bool truncated, int source, std::string& note) {
    auto is_decoy = [&](const auto& c) { return StringUtils::hasPrefix(fasta_entry[protein_of(c)].identifier, "DECOY"); };
    std::stable_sort(cands.begin(), cands.end(), [&](const auto& a, const auto& b) { return score_of(a) > score_of(b); });
    double best_decoy = -1;
    for (const auto& c : cands)
    {
      if (is_decoy(c)) { best_decoy = std::max(best_decoy, (double)score_of(c)); }
    }
    if (best_decoy < 0) best_decoy = truncated && ! cands.empty() ? (double)score_of(cands.back()) : 0.0;
    const Size n = cands.size();
    std::remove_reference_t<decltype(cands)> kept;
    std::string names;
    for (auto& c : cands)
    {
      if ((int)kept.size() >= max_protein_hits_ || score_of(c) < best_decoy + native_margin) break;
      if (is_decoy(c)) continue;
      names += std::string(kept.empty() ? " " : ", ") + fasta_entry[protein_of(c)].identifier + " " + StringUtils::number(score_of(c), 1);
      kept.push_back(std::move(c));
    }
    native_kept[source] += kept.size();
    note += "; " + native_sources[source] + ": " + StringUtils::toStr(n) + " candidates, best decoy " + StringUtils::number(best_decoy, 1) + ", kept " +
            StringUtils::toStr(kept.size()) + names;
    cands.swap(kept);
  };
  // A mass-ladder candidate as a protein hit; its ladders become tags (@p gtags) spanning them, with the lowest mass as flanking mass. The
  // tag indices of the hit are -(local index + 1), made global with the gapped tags when appended to tags_.
  auto ladderCandidateHit = [&](const DeconvolvedSpectrum& dspec, const MassLadderAlgorithm::Candidate& cand,
                                std::vector<FLASHHelperClasses::Tag>& gtags) {
    const auto& pseq = protein_seqs[cand.protein];
    std::vector<int> indices, positions;
    int matched = 0;
    for (const auto& ladder : cand.ladders)
    {
      std::vector<double> tag_masses;
      std::vector<int> tag_scores;
      for (Size mi : ladder.mass_indices)
      {
        tag_masses.push_back(dspec[mi].getMonoMass());
        tag_scores.push_back(FLASHTaggerAlgorithm::getNodeScore(dspec[mi]));
      }
      matched += (int)tag_masses.size();
      gtags.emplace_back(pseq.substr(ladder.start, ladder.end - ladder.start), ladder.prefix ? tag_masses[0] : -1.0,
                         ladder.prefix ? -1.0 : tag_masses[0], tag_masses, tag_scores, dspec.getScanNumber(), (Size)(ladder.end - ladder.start));
      gtags.back().setRetentionTime((float)dspec.getOriginalSpectrum().getRT());
      indices.push_back(-(int)gtags.size());
      positions.push_back(ladder.start);
    }
    const auto& fe = fasta_entry[cand.protein];
    ProteinHit hit(cand.score, 0, fe.identifier, fe.sequence);
    hit.setMetaValue(mi_description, fe.description); // as setDescription, without its registry lookup
    hit.setMetaValue(mi_scan, dspec.getScanNumber());
    hit.setMetaValue(mi_matched_aa, matched);
    hit.setCoverage((double)matched / (double)fe.sequence.length());
    hit.setMetaValue(mi_tag_indices, indices);
    std::sort(positions.begin(), positions.end());
    hit.setMetaValue(mi_tag_positions, positions);
    hit.setMetaValue(mi_fasta_index, cand.protein);
    hit.setMetaValue(mi_protein_rank, cand.rank);
    return hit;
  };
  // The candidate proteins of one spectrum: from its sequence tags (@p tag_indices_of_scan, indices into tags_; none: nullptr), gapped tags,
  // terminal ladders and, with @p with_precursor_ladders, precursor ladders. The proteins in @p excluded (fasta indices) are not proposed
  // by any source (the chimeric search: the proteins of the spectrum's earlier best PrSMs). Writes only its outputs: the sequence- and
  // gapped-tag candidates (@p hits_out, ranked), the terminal- and precursor-ladder candidates, the gapped and ladder tags (@p gtags, local
  // indices until appended to tags_), the number of precursor-ladder windows and, with native, what the decoy margin kept (@p note).
  // @p stage: 0 every source as configured; a stage of the search cascade only its sources: 1 sequence tags and the ladders, 2 gapped
  // tags (in every spectrum searched, as with gapped_tags only).
  auto searchCandidates = [&](const DeconvolvedSpectrum& dspec, const std::vector<Size>* tag_indices_of_scan, const std::set<Size>& excluded,
                              bool with_precursor_ladders, int stage, std::vector<ProteinHit>& hits_out, std::vector<ProteinHit>& ladder_out,
                              std::vector<ProteinHit>& pl_out, std::vector<ProteinHit>& decoy_pl_out,
                              std::vector<FLASHHelperClasses::Tag>& gtags, Size& pl_windows, std::string& note,
                              const std::vector<MassLadderAlgorithm::Candidate>* sliding) {
    const bool sequence_stage = stage == 0 || stage == 1, ladder_stage = stage == 0 || stage == 1;
    with_precursor_ladders &= ladder_stage;
    // The sequence-tag candidates (proteins hit by a tag of the spectrum; thousands per spectrum) stay records until the cut below, and only
    // the kept ones become ProteinHits. Per candidate, ascending by protein: its covered positions and its tags (distinct, ascending), and
    // its shifts (protein residue sum up to the tag minus the tag flanking mass) with the tag flanking masses, in the order of the tags and
    // of their hits (as maps per protein collected them).
    struct SeqCandidate
    {
      Size protein;
      Size covered_pos;                  // distinct covered positions, in cand_positions from positions_first
      Size positions_first, tags_first;  // tags: in cand_tags
      Size tag_count;
      Size n_first, n_count, c_first, c_count; // shifts: in shift_values and flank_values
    };
    std::vector<SeqCandidate> seq_candidates;
    std::vector<Size> cand_positions, cand_tags;
    std::vector<double> shift_values, flank_values;
    std::vector<ProteinHit> hits; // the gapped-tag candidates, which come after the sequence-tag ones
    std::vector<FLASHTaggerAlgorithm::ProteinShifts> hit_shifts; // in the order of hits
    if (sequence_stage && gapped_mode_ != "only" && tag_indices_of_scan != nullptr)
    {
      const auto& tag_indices = *tag_indices_of_scan;
      // the hits of the tags in tag order, each with the slot of its protein (protein order of first hit)
      struct ShiftHit
      {
        Size slot;
        double shift, flank;
      };
      struct PositionHit
      {
        Size slot, tag, position;
      };
      std::vector<ShiftHit> n_hits, c_hits;
      std::vector<PositionHit> position_hits;
      std::vector<int> slot_of(cleaned_protein_seqs.size(), -1);
      std::vector<Size> slot_protein;
      for (const auto& tag_index : tag_indices)
      {
        const auto& tag = tags_[tag_index];
        const auto& seq = tag.getUppercaseSequence();

        const auto& db_hits = tag_db[tag_id_of.at(seq)];
        const auto& pis = db_hits.proteins;
        const auto& pps = db_hits.positions;
        const auto& n_flanking_masses = db_hits.n_flanking_masses;
        const auto& c_flanking_masses = db_hits.c_flanking_masses;

        for (Size i = 0; i < pis.size(); i++)
        {
          if (excluded.count(pis[i]) > 0) continue;
          const auto& pseq = cleaned_protein_seqs[pis[i]];
          int& slot = slot_of[pis[i]];
          if (slot < 0)
          {
            slot = (int)slot_protein.size();
            slot_protein.push_back(pis[i]);
          }

          if (tag.getNtermMass() >= 0 && pps[i] > 0 && pps[i] <= pseq.length() + seq.length())
          {
            double nmass = n_flanking_masses[i];
            if (tag.getNtermMass() > nmass + max_total_mod_mass_ + 200) continue;
            n_hits.push_back({(Size)slot, nmass - tag.getNtermMass(), tag.getNtermMass()});
          }
          if (tag.getCtermMass() >= 0 && pps[i] + seq.length() >= 0 && pseq.length() > pps[i] + seq.length())
          {
            double cmass = c_flanking_masses[i];
            if (tag.getCtermMass() > cmass + max_total_mod_mass_ + 200) continue;
            c_hits.push_back({(Size)slot, cmass - tag.getCtermMass(), tag.getCtermMass()});
          }
          position_hits.push_back({(Size)slot, tag_index, pps[i]});
        }
      }

      // the hits of each slot, in their order: positions [first[s], first[s + 1]) of the returned order
      const Size slot_count = slot_protein.size();
      auto bySlot = [slot_count](const auto& slot_hits, std::vector<Size>& first) {
        first.assign(slot_count + 1, 0);
        for (const auto& h : slot_hits) first[h.slot + 1]++;
        for (Size s = 0; s < slot_count; s++) first[s + 1] += first[s];
        std::vector<Size> next(first.begin(), first.end() - 1), order(slot_hits.size());
        for (Size k = 0; k < slot_hits.size(); k++) order[next[slot_hits[k].slot]++] = k;
        return order;
      };
      std::vector<Size> n_first, c_first, p_first;
      const auto n_order = bySlot(n_hits, n_first), c_order = bySlot(c_hits, c_first), p_order = bySlot(position_hits, p_first);

      std::vector<Size> covered_slots; // with a covered position, ascending by protein
      for (Size s = 0; s < slot_count; s++)
      {
        if (p_first[s + 1] > p_first[s]) covered_slots.push_back(s);
      }
      std::sort(covered_slots.begin(), covered_slots.end(), [&](Size a, Size b) { return slot_protein[a] < slot_protein[b]; });
      seq_candidates.reserve(covered_slots.size());
      for (Size s : covered_slots)
      {
        SeqCandidate cand;
        cand.protein = slot_protein[s];
        cand.positions_first = cand_positions.size();
        cand.tags_first = cand_tags.size();
        for (Size k = p_first[s]; k < p_first[s + 1]; k++)
        {
          cand_positions.push_back(position_hits[p_order[k]].position);
          cand_tags.push_back(position_hits[p_order[k]].tag);
        }
        std::sort(cand_positions.begin() + cand.positions_first, cand_positions.end());
        cand_positions.erase(std::unique(cand_positions.begin() + cand.positions_first, cand_positions.end()), cand_positions.end());
        std::sort(cand_tags.begin() + cand.tags_first, cand_tags.end());
        cand_tags.erase(std::unique(cand_tags.begin() + cand.tags_first, cand_tags.end()), cand_tags.end());
        cand.covered_pos = cand_positions.size() - cand.positions_first;
        cand.tag_count = cand_tags.size() - cand.tags_first;
        if (cand.covered_pos < min_tag_count) continue;
        cand.n_first = shift_values.size();
        for (Size k = n_first[s]; k < n_first[s + 1]; k++)
        {
          shift_values.push_back(n_hits[n_order[k]].shift);
          flank_values.push_back(n_hits[n_order[k]].flank);
        }
        cand.n_count = shift_values.size() - cand.n_first;
        cand.c_first = shift_values.size();
        for (Size k = c_first[s]; k < c_first[s + 1]; k++)
        {
          shift_values.push_back(c_hits[c_order[k]].shift);
          flank_values.push_back(c_hits[c_order[k]].flank);
        }
        cand.c_count = shift_values.size() - cand.c_first;
        seq_candidates.push_back(cand);
      }
    }

    // gapped tags: proteins proposed where gapped tags land. Candidates are scored and cut first; only the survivors get their
    // sequence and their landings as tags (a landing becomes a tag spanning the landed residues, with the masses at the gap
    // boundaries), so the protein score (shifts) and the extension (anchors) use it as they use a sequence tag.
    std::vector<FLASHGappedTaggerAlgorithm::Landing> landings;
    std::vector<std::vector<Size>> gapped_groups; // per gapped candidate, its landings
    std::vector<Size> gapped_tag_counts;          // per gapped candidate, its number of tags (for the cut below)
    if (stage == 2 || (stage == 0 && (gapped_mode_ == "only" || gapped_mode_ == "both" || (gapped_mode_ == "fallback" && seq_candidates.empty()))))
    {
      std::set<Size> proposed; // proteins already proposed by sequence tags (by protein: any entry of it)
      for (const auto& c : seq_candidates) proposed.insert(protein_of[c.protein]);
      gapped_tagger.search(dspec, ms2_tol > 0 ? ms2_tol : tol, landings);
      for (Size l0 = 0; l0 < landings.size();)
      {
        Size l1 = l0;
        while (l1 < landings.size() && landings[l1].protein == landings[l0].protein) ++l1; // landings are sorted by protein
        const Size pi = landings[l0].protein;
        FLASHTaggerAlgorithm::ProteinShifts sh;
        sh.fasta_index = pi;
        std::vector<Size> group;
        std::vector<int> positions;
        for (Size l = l0; l < l1 && proposed.count(protein_of[pi]) == 0 && excluded.count(pi) == 0; ++l)
        {
          const auto& la = landings[l];
          // the lowest mass borders the residues before (prefix) or after (suffix) the landed span
          const double flank = la.masses[0];
          if (flank > la.flank_mass + max_total_mod_mass_ + 200) continue;
          if (la.prefix)
          {
            sh.n_shifts.push_back(la.flank_mass - flank);
            sh.n_flanks.push_back(flank);
          }
          else
          {
            sh.c_shifts.push_back(la.flank_mass - flank);
            sh.c_flanks.push_back(flank);
          }
          group.push_back(l);
          positions.push_back(la.start);
        }
        l0 = l1;
        if (group.empty()) continue;
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
        const auto& fe = fasta_entry[pi];
        const double coverage = (double)positions.size() / (double)fe.sequence.length();
        ProteinHit hit(coverage, 0, fe.identifier, ""); // the sequence, description and tags are set if the candidate survives the cut
        hit.setMetaValue(mi_scan, dspec.getScanNumber());
        hit.setMetaValue(mi_matched_aa, (int)positions.size());
        hit.setCoverage(coverage);
        gapped_tag_counts.push_back(group.size());
        hit.setMetaValue(mi_tag_positions, positions);
        hit.setMetaValue(mi_fasta_index, pi);
        hit.setMetaValue(mi_gapped, (int)gapped_groups.size());
        gapped_groups.push_back(std::move(group));
        hits.push_back(std::move(hit));
        hit_shifts.push_back(std::move(sh));
      }
    }

    // exact masses at the deconvolution tolerance, scored by the chance probability of the match count: 1 Da integer bins and summed
    // node scores let chance matches (more of them on longer proteins) rank noise proteins above the correct one. The candidates in the
    // order they had as hits: the sequence-tag ones, then the gapped-tag ones.
    const Size seq_count = seq_candidates.size(), candidate_count = seq_count + hits.size();
    std::vector<double> score_masses {0.0}; // mass 0 stands for the terminus (FLASHTaggerAlgorithm::runMatchingExact)
    score_masses.reserve(dspec.size() + 1);
    for (const auto& pg : dspec) score_masses.push_back(pg.getMonoMass());
    std::vector<double> scores(candidate_count, 0.0);
    // by_entry: per candidate, the scores of both sides summed, and the match counts of both sides' best shifts
    std::vector<double> side_score_sums(by_entry ? candidate_count : 0, 0.0);
    std::vector<std::pair<int, int>> side_counts(by_entry ? candidate_count : 0, {0, 0});
    // a loop that is already parallel over spectra runs this serially
#pragma omp parallel for schedule(dynamic) if (! omp_in_parallel())
    for (int i = 0; i < (int)candidate_count; i++)
    {
      const Size e = (Size)i < seq_count ? seq_candidates[i].protein : hit_shifts[i - seq_count].fasta_index;
      const double plain_mass = with_fixed ? pro_plain_mass[e] : 0.0;
      int n_count = 0, c_count = 0;
      int* const n_best = by_entry ? &n_count : nullptr;
      int* const c_best = by_entry ? &c_count : nullptr;
      double n_score, c_score;
      if ((Size)i < seq_count)
      {
        const auto& c = seq_candidates[i];
        n_score = FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_prefix[c.protein], pro_prefix_index[c.protein], shift_values.data() + c.n_first,
                                                    flank_values.data() + c.n_first, c.n_count, protein_tol, plain_mass, n_best);
        c_score = FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_suffix[c.protein], pro_suffix_index[c.protein], shift_values.data() + c.c_first,
                                                    flank_values.data() + c.c_first, c.c_count, protein_tol, plain_mass, c_best);
      }
      else
      {
        const auto& sh = hit_shifts[i - seq_count];
        n_score = FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_prefix[sh.fasta_index], pro_prefix_index[sh.fasta_index], sh.n_shifts.data(),
                                                    sh.n_flanks.data(), sh.n_shifts.size(), protein_tol, plain_mass, n_best);
        c_score = FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_suffix[sh.fasta_index], pro_suffix_index[sh.fasta_index], sh.c_shifts.data(),
                                                    sh.c_flanks.data(), sh.c_shifts.size(), protein_tol, plain_mass, c_best);
      }
      scores[i] = std::max(n_score, c_score);
      if (by_entry)
      {
        side_score_sums[i] = n_score + c_score;
        side_counts[i] = {n_count, c_count};
      }
    }

    // One entry per protein (a PEFF database: several entries per protein), so that the cuts below and the extension count proteins: the
    // entry with the best scores of both sides summed (then fewer fixed modifications, then the lower index), or, of its processed form,
    // the entry with as many matches on each side and the fewest fixed modifications (then the lowest index). Both sides, because a fixed
    // mass shows in the matches of either side (between matched masses, or between a terminus, mass 0, and them), and a candidate's score
    // reads its better side only. Equal matches, because for the same matched masses an entry with fixed modifications still scores a
    // little apart (masses in range, margins), although the chance of a match is taken from its residues without fixed masses: ranked by
    // the score alone, the heaviest entry won without evidence (Eclip_1480: an annotated mass cancelled by a blind shift in 94 PrSMs at
    // 1%, now 7). The kept candidate carries the protein's best score. The choice reads neither the precursor mass nor a label: the decoy
    // precursor copy of a spectrum extends the same entries.
    std::vector<bool> active(candidate_count, true);
    // by_entry: the other candidates (entries) of each protein; the cuts below count proteins (one candidate each), and after them every
    // candidate entry of a kept protein is extended: an annotated residue is a residue of its entry, and the extension decides
    std::unordered_map<Size, std::vector<Size>> others_of;
    if (by_entry)
    {
      auto entryOf = [&](Size i) { return i < seq_count ? seq_candidates[i].protein : hit_shifts[i - seq_count].fasta_index; };
      auto fixedCount = [&](Size i) { return db->entries[entryOf(i)].fixed.size(); };
      auto simpler = [&](Size i, Size j) { return fixedCount(i) < fixedCount(j) || (fixedCount(i) == fixedCount(j) && entryOf(i) < entryOf(j)); };
      std::unordered_map<Size, std::pair<Size, double>> best_of; // protein -> candidate, best score
      for (Size i = 0; i < candidate_count; i++)
      {
        const auto [it, inserted] = best_of.emplace(protein_of[entryOf(i)], std::make_pair(i, scores[i]));
        if (inserted) continue;
        const Size j = it->second.first;
        if (side_score_sums[i] > side_score_sums[j] || (side_score_sums[i] == side_score_sums[j] && simpler(i, j))) it->second.first = i;
        it->second.second = std::max(it->second.second, scores[i]);
      }
      std::unordered_map<Size, Size> kept_of; // protein -> candidate
      for (const auto& [protein, best] : best_of) kept_of[protein] = best.first;
      for (Size i = 0; i < candidate_count; i++)
      {
        const Size e = entryOf(i), b = best_of[protein_of[e]].first;
        Size& k = kept_of[protein_of[e]];
        const auto& entry = db->entries[e];
        const auto& best_entry = db->entries[entryOf(b)];
        if (entry.offset == best_entry.offset && entry.length == best_entry.length && side_counts[i] == side_counts[b] && simpler(i, k)) k = i;
      }
      for (Size i = 0; i < candidate_count; i++)
      {
        active[i] = kept_of[protein_of[entryOf(i)]] == i;
        if (! active[i]) others_of[protein_of[entryOf(i)]].push_back(i);
      }
      for (const auto& [protein, k] : kept_of) scores[k] = best_of[protein].second;
    }

    // Ranked by score, coverage and description with std::sort. The keys are those the hits had, so every comparison and with it the
    // order, ties included, is that of sorting the hits (std::sort's moves depend on the comparisons only).
    struct RankedCandidate
    {
      double score, coverage;
      const std::string* description;
      Size index;
    };
    std::vector<RankedCandidate> ranked;
    ranked.reserve(candidate_count);
    for (Size i = 0; i < candidate_count; i++)
    {
      if (! active[i]) continue;
      if (i < seq_count)
      {
        const auto& fe = fasta_entry[seq_candidates[i].protein];
        ranked.push_back({scores[i], (double)seq_candidates[i].covered_pos / (double)fe.sequence.length(), &fe.description, i});
      }
      else { ranked.push_back({scores[i], hits[i - seq_count].getCoverage(), &fasta_entry[hit_shifts[i - seq_count].fasta_index].description, i}); }
    }
    const Size ranked_count = ranked.size();
    std::sort(ranked.begin(), ranked.end(), [](const RankedCandidate& left, const RankedCandidate& right) {
      return left.score == right.score ? (left.coverage == right.coverage ? (*left.description > *right.description) : (left.coverage > right.coverage))
                                       : (left.score > right.score);
    });

    // the candidate as a hit with its score: a sequence-tag one is built as the maps per protein built it, a gapped-tag one is moved
    auto candidateHit = [&](Size i) {
      if (i >= seq_count)
      {
        ProteinHit hit = std::move(hits[i - seq_count]);
        hit.setScore(scores[i]);
        return hit;
      }
      const auto& c = seq_candidates[i];
      const auto& fe = fasta_entry[c.protein];
      Size covered_pos = c.covered_pos;
      double coverage = (double)covered_pos / (double)fe.sequence.length();
      std::vector<int> positions(cand_positions.begin() + c.positions_first, cand_positions.begin() + (c.positions_first + covered_pos));
      ProteinHit hit(coverage, 0, fe.identifier, fe.sequence);
      hit.setMetaValue(mi_description, fe.description); // as setDescription, without its registry lookup
      hit.setMetaValue(mi_scan, dspec.getScanNumber());
      hit.setMetaValue(mi_matched_aa, covered_pos);
      hit.setCoverage(coverage);
      hit.setMetaValue(mi_tag_indices, std::vector<int>(cand_tags.begin() + c.tags_first, cand_tags.begin() + (c.tags_first + c.tag_count)));
      hit.setMetaValue(mi_tag_positions, positions);
      hit.setMetaValue(mi_fasta_index, c.protein);
      hit.setScore(scores[i]);
      return hit;
    };
    std::vector<ProteinHit> ranked_hits;
    if (native_rule)
    {
      ranked_hits.reserve(ranked_count);
      for (const auto& r : ranked) ranked_hits.push_back(candidateHit(r.index));
      nativeCut(ranked_hits, [](const ProteinHit& h) { return h.getScore(); }, [&](const ProteinHit& h) { return (Size)(int)h.getMetaValue(mi_fasta_index); },
                false, 0, note);
    }
    else
    {
      // More candidates than max_protein_hits: the max_protein_hits best by a small model of each candidate's own evidence
      // (FLASHTnTCandidateRanker: its score, the chance probability of the ion ladder at its tags' offset, whether that offset is a residue
      // boundary), among the union of the 500 best by score and by matched tags (the candidates it was trained on). The score uses the
      // masses of the whole spectrum against the protein, so a correct protein with few of them can rank low by it; the number of distinct
      // matched tags is sequence-specific evidence. Left out by dataset, the model's 10 best held the protein of the spectrum's best PrSM of
      // a 500-candidate search in 94.6% of the spectra (gapped tags 89.2%), the union of the 100 best by score and by matched tags in 93.8%
      // (86.8%).
      std::vector<bool> keep(ranked_count, (int)ranked_count <= max_protein_hits_);
      if ((int)ranked_count > max_protein_hits_)
      {
        auto entryOf = [&](Size c) { return c < seq_count ? seq_candidates[c].protein : hit_shifts[c - seq_count].fasta_index; };
        auto tagCount = [&](Size c) { return c < seq_count ? (double)seq_candidates[c].tag_count : (double)gapped_tag_counts[c - seq_count]; };
        // a side's shifts and tag flanking masses
        auto shiftsOf = [&](Size c, bool n_side, const double*& s, const double*& f, Size& n) {
          if (c < seq_count)
          {
            const auto& sc = seq_candidates[c];
            const Size first = n_side ? sc.n_first : sc.c_first;
            s = shift_values.data() + first;
            f = flank_values.data() + first;
            n = n_side ? sc.n_count : sc.c_count;
          }
          else
          {
            const auto& sh = hit_shifts[c - seq_count];
            const auto& sv = n_side ? sh.n_shifts : sh.c_shifts;
            s = sv.data();
            f = (n_side ? sh.n_flanks : sh.c_flanks).data();
            n = sv.size();
          }
        };
        constexpr Size learned_pool = 500;
        std::vector<char> in_pool(ranked_count, 0);
        std::fill(in_pool.begin(), in_pool.begin() + std::min(ranked_count, learned_pool), 1); // ranked is in score order
        if (ranked_count > learned_pool)
        {
          std::vector<std::pair<double, Size>> by_tags;
          by_tags.reserve(ranked_count);
          for (Size i = 0; i < ranked_count; i++) by_tags.emplace_back(-tagCount(ranked[i].index), i);
          std::partial_sort(by_tags.begin(), by_tags.begin() + learned_pool, by_tags.end());
          for (Size i = 0; i < learned_pool; i++) in_pool[by_tags[i].second] = 1;
        }
        FLASHTnTCandidateRanker::Spectrum spectrum;
        spectrum.masses.reserve(dspec.size());
        for (const auto& pg : dspec) spectrum.masses.push_back(pg.getMonoMass());
        std::sort(spectrum.masses.begin(), spectrum.masses.end());
        spectrum.prefix_offset = learned_prefix_offset;
        spectrum.suffix_offset = learned_suffix_offset;
        std::vector<double> key(ranked_count, -std::numeric_limits<double>::infinity());
        for (Size i = 0; i < ranked_count; i++)
        {
          if (! in_pool[i]) continue;
          const Size c = ranked[i].index, e = entryOf(c);
          FLASHTnTCandidateRanker::Candidate cand;
          cand.prefix = &pro_prefix[e];
          cand.suffix = &pro_suffix[e];
          const double* flanks;
          shiftsOf(c, true, cand.n_shifts, flanks, cand.n_count);
          shiftsOf(c, false, cand.c_shifts, flanks, cand.c_count);
          cand.score = ranked[i].score;
          key[i] = FLASHTnTCandidateRanker::score(FLASHTnTCandidateRanker::features(spectrum, cand));
        }
        // the max_protein_hits best by the key; ties in score order
        std::vector<Size> order(ranked_count);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) { return key[a] > key[b]; });
        for (int j = 0; j < max_protein_hits_; j++) keep[order[j]] = true;
      }
      for (Size i = 0; i < ranked_count; i++)
      {
        if (keep[i]) ranked_hits.push_back(candidateHit(ranked[i].index)); // in score order
      }
    }
    if (by_entry && ! others_of.empty())
    {
      std::vector<ProteinHit> with_entries;
      with_entries.reserve(ranked_hits.size());
      for (auto& kept : ranked_hits)
      {
        const Size p = protein_of[(Size)(int)kept.getMetaValue(mi_fasta_index)];
        const double kept_score = kept.getScore();
        with_entries.push_back(std::move(kept));
        auto it = others_of.find(p);
        if (it == others_of.end()) continue;
        for (Size i : it->second)
        {
          ProteinHit other = candidateHit(i);
          other.setScore(std::min(scores[i], kept_score));
          other.setMetaValue(mi_entry_of_kept, 1);
          with_entries.push_back(std::move(other));
        }
      }
      ranked_hits.swap(with_entries);
    }
    hits = std::move(ranked_hits);

    // the surviving gapped candidates: sequence and tags (indices local to this spectrum until appended to tags_ below)
    for (auto& hit : hits)
    {
      if (! hit.metaValueExists(mi_gapped)) continue;
      const Size pi = (Size)(int)hit.getMetaValue(mi_fasta_index);
      hit.setSequence(fasta_entry[pi].sequence);
      hit.setMetaValue(mi_description, fasta_entry[pi].description); // as setDescription, without its registry lookup
      std::vector<int> indices;
      for (Size l : gapped_groups[(int)hit.getMetaValue(mi_gapped)])
      {
        auto& la = landings[l];
        const double flank = la.masses[0];
        gtags.emplace_back(protein_seqs[pi].substr(la.start, la.end - la.start), la.prefix ? flank : -1.0, la.prefix ? -1.0 : flank,
                           la.masses, la.scores, dspec.getScanNumber(), (Size)(la.end - la.start));
        gtags.back().setRetentionTime((float)dspec.getOriginalSpectrum().getRT());
        indices.push_back((int)gtags.size() - 1);
      }
      hit.setMetaValue(mi_tag_indices, indices);
      hit.setMetaValue(mi_gapped, 1);
    }

    // the deconvolved masses (ascending) for the mass ladders
    std::vector<double> masses;
    if (ladder_stage && (terminal_ladders_ || with_precursor_ladders || native_rule) && ! dspec.empty())
    {
      masses.reserve(dspec.size());
      for (const auto& pg : dspec) masses.push_back(pg.getMonoMass());
    }

    // Terminal ladders (MassLadderAlgorithm::searchTerminalLadders): proteins by the chance probability of the masses landing on a residue
    // sum from a terminus. The min_protein_hits best proteins are characterized; the ladder of the best shift of each side becomes a tag
    // spanning it (as a gapped tag does), so the extension is anchored there. They are kept apart from the other candidates, also when a
    // sequence tag proposed the same protein.
    // Native mode: sliding ladders (from any start or to any end of a protein) instead, kept by the decoy margin.
    if (ladder_stage && (terminal_ladders_ || native_rule) && ! dspec.empty())
    {
      std::vector<MassLadderAlgorithm::Candidate> found;
      if (native_rule)
      {
        if (sliding != nullptr) found = *sliding; // searched for all spectra at once (searchSlidingSpectra)
        else mass_ladders.searchSlidingLadders(masses, protein_tol, native_ladder_top, excluded, found);
        nativeCut(found, [](const MassLadderAlgorithm::Candidate& c) { return c.score; }, [](const MassLadderAlgorithm::Candidate& c) { return c.protein; },
                  found.size() >= native_ladder_top, 1, note);
      }
      else { mass_ladders.searchTerminalLadders(masses, protein_tol, (Size)min_protein_hits_, excluded, found); }
      for (const auto& cand : found)
      {
        ProteinHit hit = ladderCandidateHit(dspec, cand, gtags);
        hit.setMetaValue(native_rule ? mi_sliding : mi_terminal, 1);
        ladder_out.push_back(std::move(hit));
      }
    }

    // Precursor-anchored ladders: the windows of the database (targets and decoys alike) whose unmodified mass is the precursor mass, by
    // the chance probability of the masses on their prefix and suffix ladders (MassLadderAlgorithm::searchPrecursorLadders). The
    // min_protein_hits best windows are characterized; the matched masses of the best prefix and the best suffix ion type become tags
    // spanning them (as terminal-ladder tags), from which the extension can start and end the proteoform at the window. FLASHExtender places
    // the tags by their sequence and reaches an internal start or end only within its truncation range around the protein's best mass shift
    // (getBestMassShift_), not at the window it is given: the window is stored on the candidate (PrecursorLadderStart and End, 1-based as
    // StartPosition and EndPosition, and PrecursorLadderIsotope), and the log line after the extension counts the PrSMs at its ends. Kept
    // apart from the other candidates; searched before the check for tag candidates, so spectra without any are searched too. The same at
    // the decoy precursor mass (decoy_precursor).
    if (with_precursor_ladders && ! dspec.empty() && dspec.getPrecursorPeakGroup().getMonoMass() > 0)
    {
      std::vector<double> precursor_masses {dspec.getPrecursorPeakGroup().getMonoMass()};
      if (decoy_precursor_ && dspec.getDecoyPrecursorMass() > 0) precursor_masses.push_back(dspec.getDecoyPrecursorMass());
      std::vector<MassLadderAlgorithm::PrecursorLadderResult> found;
      mass_ladders.searchPrecursorLadders(masses, precursor_masses, protein_tol, pl_precursor_tol, native_rule ? native_ladder_top : (Size)min_protein_hits_, found);
      pl_windows = found[0].windows;
      for (auto& result : found)
      {
        if (native_rule)
        {
          nativeCut(result.candidates, [](const MassLadderAlgorithm::Candidate& c) { return c.score; },
                    [](const MassLadderAlgorithm::Candidate& c) { return c.protein; }, result.candidates.size() >= native_ladder_top, 2, note);
        }
      }
      for (Size n = 0; n < found.size(); n++)
      {
        for (const auto& cand : found[n].candidates)
        {
          ProteinHit hit = ladderCandidateHit(dspec, cand, gtags);
          hit.setMetaValue(mi_precursor_ladder, 1);
          // residues i..j-1 (0-based) are i+1..j (1-based, inclusive); on the protein (the entry's offset), as the decoded PrSMs
          hit.setMetaValue(mi_pl_start, cand.window_start + 1 + db->entries[cand.protein].offset);
          hit.setMetaValue(mi_pl_end, cand.window_end + db->entries[cand.protein].offset);
          hit.setMetaValue(mi_pl_isotope, cand.isotope);
          (n == 0 ? pl_out : decoy_pl_out).push_back(std::move(hit));
        }
      }
    }

    if (hits.empty()) return;
    for (Size i = 0; i < hits.size(); i++)
    {
      hits[i].setMetaValue(mi_protein_rank, (int)i + 1);
    }
    hits_out = std::move(hits);
  };
  const std::set<Size> no_exclusion;
  // gapped and ladder tags of a spectrum join tags_; returns the index of the first, from which the tag indices of its candidates become global
  auto appendSpectrumTags = [&](std::vector<FLASHHelperClasses::Tag>& gtags) {
    const int tag_base = (int)tags_.size();
    for (auto& t : gtags)
    {
      t.setIndex((int)tags_.size());
      tags_.push_back(std::move(t));
    }
    return tag_base;
  };
  // sequence- and gapped-tag candidates: the local tag indices of gapped candidates become global (the residual rounds below)
  auto appendTagCandidates = [&](std::vector<ProteinHit>& from, std::vector<ProteinHit>& to, std::pair<Size, Size>& range, int tag_base) {
    range.first = to.size();
    for (auto& hit : from)
    {
      if (hit.metaValueExists(mi_gapped))
      {
        auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
        for (auto& i : indices) i += tag_base;
        hit.setMetaValue(mi_tag_indices, indices);
      }
      to.push_back(std::move(hit));
    }
    range.second = to.size();
  };
  // ladder candidates (terminal and precursor ladders): their tag indices -(local index + 1) become global (the residual rounds below)
  auto appendLadderCandidates = [&](std::vector<ProteinHit>& from, std::vector<ProteinHit>& to, std::pair<Size, Size>& range, int tag_base) {
    range.first = to.size();
    for (auto& hit : from)
    {
      auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
      for (auto& i : indices) i = tag_base - i - 1;
      hit.setMetaValue(mi_tag_indices, indices);
      to.push_back(std::move(hit));
    }
    range.second = to.size();
  };

  // native mode: the sliding ladders of the spectra @p which of @p ds (with the proteins @p excluded_of of each, nullptr: none) in one batched
  // pass over the database (MassLadderAlgorithm::searchSlidingLadders for many spectra), indexed by spectrum in @p out
  auto searchSlidingSpectra = [&](const std::vector<DeconvolvedSpectrum>& ds, const std::vector<int>& which, const std::vector<std::set<Size>>* excluded_of,
                                  std::vector<std::vector<MassLadderAlgorithm::Candidate>>& out) {
    out.assign(ds.size(), {});
    std::vector<std::vector<double>> masses(which.size());
    std::vector<const std::vector<double>*> masses_of;
    std::vector<const std::set<Size>*> excluded;
    for (Size k = 0; k < which.size(); k++)
    {
      for (const auto& pg : ds[which[k]]) masses[k].push_back(pg.getMonoMass()); // as searchCandidates
      masses_of.push_back(&masses[k]);
      excluded.push_back(excluded_of == nullptr ? &no_exclusion : &(*excluded_of)[which[k]]);
    }
    std::vector<std::vector<MassLadderAlgorithm::Candidate>> found;
    mass_ladders.searchSlidingLadders(masses_of, protein_tol, native_ladder_top, excluded, found);
    for (Size k = 0; k < which.size(); k++) out[which[k]] = std::move(found[k]);
  };

  // Search cascade (search_cascade): stage 1 searches every spectrum with sequence tags and the ladders, stage 2 with gapped tags only the
  // spectra stage 1 did not identify (a target PrSM whose proteoform passes rep_q, see identifiedSoFar below). Without the cascade, one
  // stage (0) searches every spectrum with every source. The stages are decided by the target-decoy FDR, so the cascade needs decoys.
  const bool cascade = search_cascade_ && gapped_mode_ != "off" && ! native_rule && ! protein_level_only_ && decoy_factor_ > 0;
  std::vector<int> stages {0};
  if (cascade) stages = {1, 2};
  // The candidates of @p spectra (indices into dspecs) from the sources of @p stage, appended to protein_hits_ and the ladder candidates
  // in spectrum order; the ranges of the other spectra are left as they are.
  auto searchStage = [&](int stage, const std::vector<int>& spectra) {
    const bool with_precursor_ladders = use_precursor_ladders && (stage == 0 || stage == 1);
    std::vector<std::string> native_notes(dspecs.size());
    // the precursor-ladder windows of all precursor masses searched below, in one pass over the database (the same windows as a sweep per mass)
    if (with_precursor_ladders)
    {
      std::vector<double> all_precursor_masses;
      for (int d : spectra)
      {
        const auto& dspec = dspecs[d];
        if (dspec.empty() || dspec.getPrecursorPeakGroup().getMonoMass() <= 0) continue;
        all_precursor_masses.push_back(dspec.getPrecursorPeakGroup().getMonoMass());
        if (decoy_precursor_ && dspec.getDecoyPrecursorMass() > 0) all_precursor_masses.push_back(dspec.getDecoyPrecursorMass());
      }
      mass_ladders.preparePrecursorWindows(all_precursor_masses, pl_precursor_tol);
    }
    std::vector<std::vector<MassLadderAlgorithm::Candidate>> sliding_of;
    if (native_rule) searchSlidingSpectra(dspecs, spectra, nullptr, sliding_of);
    // the spectra with the most masses first (each writes its own slots, so the order does not change the result): a few heavy spectra
    // at the end of the dynamic schedule were the tail of this loop
    std::vector<int> search_order(spectra);
    std::stable_sort(search_order.begin(), search_order.end(), [&](int l, int r) { return dspecs[l].size() > dspecs[r].size(); });
    const bool over_spectra = parallel_over_spectra;
#pragma omp parallel for schedule(dynamic) if (over_spectra)
    for (int k = 0; k < (int)search_order.size(); k++)
    {
      nextProgress();
      const int d = search_order[k];
      const auto tag_it = scan_to_tag_indices.find(dspecs[d].getScanNumber());
      searchCandidates(dspecs[d], tag_it == scan_to_tag_indices.end() ? nullptr : &tag_it->second, no_exclusion, with_precursor_ladders, stage,
                       spec_hits[d], spec_ladder_hits[d], spec_pl_hits[d], spec_decoy_pl_hits[d], spec_gapped_tags[d], spec_pl_windows[d],
                       native_notes[d], native_rule ? &sliding_of[d] : nullptr);
    }
    if (with_precursor_ladders) mass_ladders.clearPrecursorWindows();
    if (native_rule)
    {
      OPENMS_LOG_INFO << "Native candidates (score above the best decoy of the source in the spectrum + " << native_margin << "):";
      for (Size k = 0; k < native_sources.size(); k++)
      {
        OPENMS_LOG_INFO << (k ? ", " : " ") << native_sources[k] << " " << (Size)native_kept[k];
      }
      OPENMS_LOG_INFO << " (" << dspecs.size() << " spectra)" << std::endl;
      // the evidence per spectrum, for the few spectra of a native run
      for (Size d = 0; d < dspecs.size() && dspecs.size() <= 20; d++)
      {
        if (! native_notes[d].empty()) OPENMS_LOG_INFO << "  scan " << dspecs[d].getScanNumber() << native_notes[d] << std::endl;
      }
    }
    const Size pl_before = pl_candidates.size(), decoy_pl_before = decoy_pl_candidates.size();
    // The candidates of each spectrum are appended in spectrum order: their positions by prefix sums, then moved (and their tag indices made
    // global) by spectrum in parallel. Sequence- and gapped-tag candidates: the local tag indices of gapped candidates become global; ladder
    // candidates (terminal and precursor ladders): their tag indices -(local index + 1) become global.
    {
      Size n_tags = tags_.size();
      for (int d : spectra) n_tags += spec_gapped_tags[d].size();
      tags_.reserve(n_tags); // room for all, so that tags_ does not grow (and copy the tags) while they are appended
    }
    std::vector<int> tag_bases(dspecs.size(), 0);
    for (int d : spectra)
    {
      tag_bases[d] = appendSpectrumTags(spec_gapped_tags[d]);
      std::vector<FLASHHelperClasses::Tag>().swap(spec_gapped_tags[d]);
    }
    auto placeCandidates = [&](const std::vector<std::vector<ProteinHit>>& from, std::vector<ProteinHit>& to, std::vector<std::pair<Size, Size>>& ranges) {
      Size position = to.size();
      for (int d : spectra)
      {
        ranges[d] = {position, position + from[d].size()};
        position += from[d].size();
      }
      to.resize(position);
    };
    placeCandidates(spec_hits, protein_hits_, candidate_ranges);
    placeCandidates(spec_ladder_hits, ladder_candidates, ladder_ranges);
    placeCandidates(spec_pl_hits, pl_candidates, pl_ranges);
    placeCandidates(spec_decoy_pl_hits, decoy_pl_candidates, decoy_pl_ranges);
    auto moveLadderCandidates = [&](std::vector<ProteinHit>& from, std::vector<ProteinHit>& to, Size first, int tag_base) {
      for (Size k = 0; k < from.size(); k++)
      {
        auto& hit = from[k];
        auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
        for (auto& i : indices) i = tag_base - i - 1;
        hit.setMetaValue(mi_tag_indices, indices);
        to[first + k] = std::move(hit);
      }
      std::vector<ProteinHit>().swap(from);
    };
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < (int)spectra.size(); k++)
    {
      const int d = spectra[k];
      const int tag_base = tag_bases[d];
      for (Size j = 0; j < spec_hits[d].size(); j++)
      {
        auto& hit = spec_hits[d][j];
        if (hit.metaValueExists(mi_gapped))
        {
          auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
          for (auto& i : indices) i += tag_base;
          hit.setMetaValue(mi_tag_indices, indices);
        }
        protein_hits_[candidate_ranges[d].first + j] = std::move(hit);
      }
      std::vector<ProteinHit>().swap(spec_hits[d]);
      moveLadderCandidates(spec_ladder_hits[d], ladder_candidates, ladder_ranges[d].first, tag_base);
      moveLadderCandidates(spec_pl_hits[d], pl_candidates, pl_ranges[d].first, tag_base);
      moveLadderCandidates(spec_decoy_pl_hits[d], decoy_pl_candidates, decoy_pl_ranges[d].first, tag_base);
    }

    if (with_precursor_ladders)
    {
      std::vector<Size> windows; // per searched spectrum
      for (int d : spectra)
      {
        if (dspecs[d].getPrecursorPeakGroup().getMonoMass() > 0 && ! dspecs[d].empty()) windows.push_back(spec_pl_windows[d]);
      }
      Size total = 0;
      for (Size w : windows) total += w;
      std::nth_element(windows.begin(), windows.begin() + windows.size() / 2, windows.end());
      std::string at_decoy_masses;
      if (decoy_precursor_) at_decoy_masses = ", " + StringUtils::toStr(decoy_pl_candidates.size() - decoy_pl_before) + " at the decoy precursor masses";
      OPENMS_LOG_INFO << "Precursor ladders: " << mass_ladders.getSegmentCount() << " segments; " << windows.size() << " spectra with a precursor mass, " << total
                      << " windows (median " << (windows.empty() ? 0 : windows[windows.size() / 2]) << " per spectrum), " << pl_candidates.size() - pl_before
                      << " candidates (top " << min_protein_hits_ << ", at least " << MassLadderAlgorithm::precursor_ladder_min_k << " matched positions)" << at_decoy_masses
                      << std::endl;
    }
  };
  std::vector<int> all_spectra(dspecs.size());
  std::iota(all_spectra.begin(), all_spectra.end(), 0);
  searchStage(stages[0], all_spectra);
  endProgress();

  // the FDR below: decoys, their weights, the PrSM level of a passing PrSM and the FDR class of a PrSM; defined here for the search cascade,
  // which decides from each stage's PrSM q-values which spectra the next stage searches
  const bool use_fdr = decoy_factor_ > 0 || ! keep_underdetermined_;
  auto is_decoy = [](const ProteinHit& hit) { return StringUtils::hasPrefix(hit.getAccession(), "DECOY") || hit.metaValueExists("DecoyPrecursor"); };
  // a decoy precursor counts with its weight (1 for a signal decoy, FLASHDeconv's noise decoy weight for a noise decoy), a protein decoy
  // for 1 / decoy_factor_
  auto decoy_weight = [this](const ProteinHit& hit) {
    return hit.metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(hit.getAccession(), "DECOY") ? (double)hit.getMetaValue("DecoyPrecursor") : 1.0 / decoy_factor_;
  };
  // PrSM q-value of a passing PrSM for the representatives and the protein ranking: the reported level, or 1% without a filter
  // (pro_fdr 1, as in keep_decoy runs for offline evaluation), so that such a run has the proteoforms of a 1% run
  const double rep_q = pro_fdr_ < 1 ? pro_fdr_ : 0.01;
  const UInt mi_chimeric_round = registry.registerName("ChimericRound");
  // the FDR class of a PrSM (index into fdr_class_names); -1: in no class
  auto fdrClassOf = [&](const ProteinHit& hit) {
    const int chimeric_round = hit.metaValueExists(mi_chimeric_round) ? (int)hit.getMetaValue(mi_chimeric_round) : 0;
    if (chimeric_round > 0)
    {
      // by source: a residual round of the chimeric search; undetermined ones only with keep_underdetermined_, as below
      if (! isDetermined(hit) && ! keep_underdetermined_) return -1;
      return fdr_class_residual + std::min(chimeric_round, max_chimeric_rounds) - 1;
    }
    if (hit.metaValueExists(mi_precursor_ladder)) return (int)fdr_class_precursor_ladder; // by source; all of them are determined
    if (isDetermined(hit)) return (int)std::min(blindModCount_(hit), (Size)fdr_class_max_blind);
    if (! keep_underdetermined_) return -1;
    const bool both_termini = (int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0;
    return (int)(both_termini ? fdr_class_both_termini : fdr_class_open_terminus);
  };

  // the order of the PrSMs for the FDR: best score first, then the higher coverage, then the higher scan number
  auto byScore = [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  };
  // The PrSM and proteoform q-values of proteoform_hits_ (sorted by byScore), as described at the FDR below: per PrSM its FDR class
  // (hit_class), PrSM q-value, decoy flag and decoy weight; the proteoforms (markRepresentativeProteoformHits_) are returned and every
  // PrSM gets its proteoform's q-value (proqvalue). With @p log, the classes are logged. Used by the FDR below and, on a copy, by the
  // search cascade between its stages.
  auto proteoformQvalues = [&](std::vector<int>& hit_class, std::vector<double>& qvalues, std::vector<bool>& decoys,
                               std::vector<double>& decoy_weights, bool log) {
    const Size hit_count = proteoform_hits_.size();
    hit_class.assign(hit_count, -1);
    qvalues.assign(hit_count, 1.0);
    decoy_weights.assign(hit_count, 0.0);
    decoys.assign(hit_count, false);
    if (use_fdr)
    {
      std::vector<std::vector<Size>> class_hits(fdr_class_count); // in the order above
      for (Size i = 0; i < hit_count; i++)
      {
        const auto& hit = proteoform_hits_[i];
        decoys[i] = is_decoy(hit);
        if (decoys[i]) decoy_weights[i] = decoy_weight(hit);
        const int c = fdrClassOf(hit);
        if (c < 0) continue;
        hit_class[i] = c;
        class_hits[c].push_back(i);
      }
      // per class: targets and decoys (count and weight), and the protein-level null it uses (its own from null_min_class_decoys decoys on,
      // else the pooled one); for the precursor-ladder class also by blind modification count, for a residual class by determined mass and
      // termini. The residual classes are logged only with the chimeric search.
      for (int c = 0; log && c < fdr_class_residual + chimeric_rounds_; c++)
      {
        Size targets = 0, decoy_count = 0, blind_targets[3] = {0, 0, 0}, blind_decoys[3] = {0, 0, 0}, determined_targets = 0, determined_decoys = 0;
        double weight = 0;
        for (Size i : class_hits[c])
        {
          const Size b = std::min(blindModCount_(proteoform_hits_[i]), (Size)2);
          const bool determined = isDetermined(proteoform_hits_[i]);
          if (decoys[i])
          {
            decoy_count++;
            weight += decoy_weights[i];
            blind_decoys[b]++;
            determined_decoys += determined;
          }
          else
          {
            targets++;
            blind_targets[b]++;
            determined_targets += determined;
          }
        }
        std::string line = std::string("FDR class \"") + fdr_class_names[c] + "\": " + StringUtils::toStr(targets) + " target and " + StringUtils::toStr(decoy_count)
                      + " decoy PrSMs (decoy weight " + StringUtils::toStr(weight) + "); protein-level null: "
                      + (c == fdr_class_precursor_ladder ? std::string("that of \"undetermined, both termini\" (setProteinQvalues_)")
                         : std::string(decoy_count >= null_min_class_decoys ? "own" : "pooled")
                             + (c <= fdr_class_max_blind ? "; those not the whole protein: that of \"undetermined, both termini\"" : ""));
        if (c == fdr_class_precursor_ladder && ! class_hits[c].empty())
        {
          line += "; target/decoy PrSMs by blind modifications: 0: " + StringUtils::toStr(blind_targets[0]) + "/" + StringUtils::toStr(blind_decoys[0]) + ", 1: "
                  + StringUtils::toStr(blind_targets[1]) + "/" + StringUtils::toStr(blind_decoys[1]) + ", 2 or more: " + StringUtils::toStr(blind_targets[2]) + "/" + StringUtils::toStr(blind_decoys[2]);
        }
        if (c >= fdr_class_residual && ! class_hits[c].empty())
        {
          line += "; target/decoy PrSMs determined: " + StringUtils::toStr(determined_targets) + "/" + StringUtils::toStr(determined_decoys) + ", undetermined: "
                  + StringUtils::toStr(targets - determined_targets) + "/" + StringUtils::toStr(decoy_count - determined_decoys);
        }
        OPENMS_LOG_INFO << line << std::endl;
      }
      const UInt mi_search_stage = registry.registerName("SearchStage");
      std::map<int, std::pair<Size, Size>> stage_counts; // search cascade: target and decoy PrSMs in a class per stage
      for (int c = 0; c < fdr_class_count; c++)
      {
        // the class's PrSMs by stage (0: no search cascade, or a residual round), each in the order above
        std::map<int, std::vector<Size>> by_stage;
        for (Size i : class_hits[c])
        {
          const int stage = proteoform_hits_[i].metaValueExists(mi_search_stage) ? (int)proteoform_hits_[i].getMetaValue(mi_search_stage) : 0;
          by_stage[stage].push_back(i);
          if (stage > 0) (decoys[i] ? stage_counts[stage].second : stage_counts[stage].first)++;
        }
        // undetermined_decoy_weight / undetermined_pseudo_decoys: the both-termini class counted with heavier decoys and a pseudo-count;
        // open_terminus_prsms separate: the open-terminus class passes nothing (its rows stay, at q-value 1)
        // only with the learned scoring, which the calibration was measured with: without it (-native, no -in_raw, scoring_model none)
        // every class is counted as before
        const bool learned_scoring = scoring_model_ && raw_loader_;
        const bool both_termini_class = learned_scoring && c == (int)fdr_class_both_termini;
        const bool separate_class = learned_scoring && c == (int)fdr_class_open_terminus && open_terminus_separate_;
        for (const auto& [stage, hits] : by_stage)
        {
          std::vector<RankedItem<double>> ranked;
          ranked.reserve(hits.size());
          for (Size i : hits)
            ranked.push_back({proteoform_hits_[i].getScore(), decoys[i], decoy_weights[i] * (both_termini_class ? undetermined_decoy_weight_ : 1.0)});
          auto class_qvalues = targetDecoyQvalues(ranked, both_termini_class ? undetermined_pseudo_decoys_ : 0.0);
          if (separate_class) std::fill(class_qvalues.begin(), class_qvalues.end(), 1.0);
          for (Size k = 0; k < hits.size(); k++)
          {
            const Size i = hits[k];
            qvalues[i] = class_qvalues[k];
            proteoform_hits_[i].setMetaValue("qvalue", qvalues[i]);
            proteoform_hits_[i].setMetaValue("FDRClass", fdr_class_names[c]);
          }
        }
      }
      for (const auto& [stage, counts] : stage_counts)
      {
        if (! log) break;
        OPENMS_LOG_INFO << "Search cascade stage " << stage << ": " << counts.first << " target and " << counts.second
                        << " decoy PrSMs in FDR classes (q-values per class within the stage)" << std::endl;
      }
    }
    // the proteoforms: positions of their PrSMs in proteoform_hits_ (Index), the representative first; representatives in the order above
    auto clusters = markRepresentativeProteoformHits_(precursor_tol, hit_class, qvalues, rep_q, use_fdr);
    if (use_fdr && ! native_rule)
    {
      // Proteoform q-values per proteoform class: the classes partition the proteoforms (each in exactly one, by a fixed precedence on its
      // PrSMs, the same for targets and decoys), and each class (or group of small classes, below) is target-decoy counted on its own, so that no class
      // passes on the decoy-free margin of another. Ranked against each other in one list by their PrSM q-values, the proteoforms of
      // every class passed at 1% together (Eclip_1480: 1418), but the undetermined ones among them at 6-11% (both termini 100 targets
      // and 6 decoys, open terminus 45 and 5; in the B. theta entrapment run 7 of the 113 with both termini were B. theta); one list by
      // score passed 140 fewer determined ones instead.
      // The determined classes of the PrSMs (0, 1, 2 blind) are one proteoform class: their decoys are nearly absent (2 for 1273 targets):
      // apart, the smaller ones would be decided by no decoy. Precursor-ladder proteoforms are a class of their own (with few targets it is
      // counted together with its neighbours, below): their wrong windows are near misses of a similar sequence that shuffled decoys do not
      // show (B. theta entrapment of the three lysates: in one class with the determined proteoforms, 18 B. theta precursor-ladder
      // proteoforms at 1%, 9% estimated error; in a class of their own 4, with 9 more E. coli proteins over the runs). Precedence: a determined PrSM of the
      // main search other than a precursor-ladder one (any blind count) makes the proteoform "determined"; else a precursor-ladder PrSM
      // "precursor ladder"; else one with both termini "undetermined, both termini"; else "undetermined, open terminus"; a proteoform of
      // residual-round PrSMs only is "residual r" (its lowest round). A proteoform is ranked by the best score of the PrSMs that decide its
      // class, and its row (Representative) is one of them: of a determined or precursor-ladder proteoform the PrSM passing rep_q with the
      // fewest blind modifications, else the best scoring one. Every PrSM carries the proteoform's q-value and class (ProteoformFDRClass).
      std::vector<int> cluster_class(clusters.size(), -1);
      std::vector<double> cluster_score(clusters.size(), 0.0);
      for (Size k = 0; k < clusters.size(); k++)
      {
        auto& members = clusters[k];
        std::vector<Size> main_members, defining;
        int lowest_round = max_chimeric_rounds;
        for (Size i : members)
        {
          if (hit_class[i] < 0) continue;
          if (hit_class[i] < fdr_class_residual) main_members.push_back(i);
          else lowest_round = std::min(lowest_round, hit_class[i] - fdr_class_residual + 1);
        }
        if (! main_members.empty())
        {
          for (Size i : main_members) { if (isDetermined(proteoform_hits_[i]) && hit_class[i] != fdr_class_precursor_ladder) defining.push_back(i); }
          cluster_class[k] = pf_class_determined;
          if (defining.empty())
          {
            for (Size i : main_members) { if (hit_class[i] == fdr_class_precursor_ladder) defining.push_back(i); }
            cluster_class[k] = pf_class_precursor_ladder;
          }
          if (defining.empty())
          {
            for (Size i : main_members)
            {
              if ((int)proteoform_hits_[i].getMetaValue("StartPosition") > 0 && (int)proteoform_hits_[i].getMetaValue("EndPosition") > 0) defining.push_back(i);
            }
            cluster_class[k] = pf_class_both_termini;
          }
          if (defining.empty())
          {
            defining = main_members;
            cluster_class[k] = pf_class_open_terminus;
          }
        }
        else
        {
          for (Size i : members) { if (hit_class[i] >= 0) defining.push_back(i); }
          if (defining.empty()) continue; // in no class
          cluster_class[k] = pf_class_residual + lowest_round - 1;
        }
        Size row = defining[0];
        for (Size i : defining)
        {
          cluster_score[k] = std::max(cluster_score[k], proteoform_hits_[i].getScore());
          if (cluster_class[k] == pf_class_determined || cluster_class[k] == pf_class_precursor_ladder)
          {
            // passing first, then fewer blind modifications, then the higher score
            const bool pass_i = qvalues[i] <= rep_q, pass_row = qvalues[row] <= rep_q;
            if (pass_i != pass_row) { if (pass_i) row = i; continue; }
            const Size blind_i = blindModCount_(proteoform_hits_[i]), blind_row = blindModCount_(proteoform_hits_[row]);
            if (blind_i != blind_row) { if (blind_i < blind_row) row = i; continue; }
          }
          if (proteoform_hits_[i].getScore() > proteoform_hits_[row].getScore()) row = i;
        }
        // the row first, the other PrSMs in their order
        std::stable_partition(members.begin(), members.end(), [row](Size i) { return i == row; });
        std::vector<int> prsm_indices;
        for (Size i : members)
        {
          proteoform_hits_[i].setMetaValue("Representative", i == row ? 1 : 0);
          prsm_indices.push_back((int)i);
        }
        proteoform_hits_[row].setMetaValue("PrSMIndices", prsm_indices);
        proteoform_hits_[row].setMetaValue("ProteoformEvidenceClass", cluster_class[k]); // its own class, also when counted with others
      }
      // Classes with few targets are counted together: in the order of the classes, a group takes classes until it has min_class_targets
      // targets, and a last group with fewer joins the one before. The classes of a lysate run have hundreds of targets each and stay
      // apart; in a small run (a single SEC fraction of PXD012247: 83 determined, 20 and 90 undetermined targets) a class alone is
      // decided by its few decoys, or by none.
      // The open-terminus class (undetermined mass, one terminus unknown) is not counted: its proteoforms keep their rows, with q-value 1,
      // and no protein (below). The class holds about as many decoys as targets (B. theta entrapment of the three lysates, three seeds:
      // 722-835 target and 735-858 decoy proteoforms per run), so the 5-15 targets above its first decoy passed at q = 0: 83 E. coli and
      // 5 B. theta proteoforms over the nine runs (5.7%, against 0.5% for the other classes), and 5 of the 27 B. theta proteins at 1%.
      // Target-decoy counting with +1 passes none of them either, but in every class it would also drop the precursor-ladder class
      // (about 83 targets per run, 0.5% B. theta). Without the class: proteoforms 14989 -> 14906, proteins 2953 -> 2923 E. coli.
      constexpr Size min_class_targets = 100;
      std::vector<Size> class_targets(pf_class_count, 0), class_size(pf_class_count, 0);
      for (Size k = 0; k < clusters.size(); k++)
      {
        if (cluster_class[k] < 0) continue;
        class_size[cluster_class[k]]++;
        class_targets[cluster_class[k]] += decoys[clusters[k][0]] ? 0 : 1;
      }
      std::vector<std::vector<int>> groups;
      Size group_targets = 0;
      for (int c = 0; c < pf_class_count; c++)
      {
        if (class_size[c] == 0 || c == pf_class_open_terminus) continue;
        if (groups.empty() || group_targets >= min_class_targets)
        {
          groups.push_back({});
          group_targets = 0;
        }
        groups.back().push_back(c);
        group_targets += class_targets[c];
      }
      if (groups.size() > 1 && group_targets < min_class_targets)
      {
        groups[groups.size() - 2].insert(groups[groups.size() - 2].end(), groups.back().begin(), groups.back().end());
        groups.pop_back();
      }
      std::vector<int> group_of(pf_class_count, -1);
      std::vector<std::string> group_names;
      for (Size g = 0; g < groups.size(); g++)
      {
        std::string name;
        for (int c : groups[g])
        {
          group_of[c] = (int)g;
          name += (name.empty() ? "" : " + ") + std::string(pf_class_names[c]);
        }
        group_names.push_back(name);
      }
      std::map<int, std::array<Size, 3>> class_counts; // per group: targets, decoys, targets passing rep_q
      for (int c = 0; c < (int)groups.size(); c++)
      {
        std::vector<Size> order;
        for (Size k = 0; k < clusters.size(); k++) { if (cluster_class[k] >= 0 && group_of[cluster_class[k]] == c) order.push_back(k); }
        if (order.empty()) continue;
        std::stable_sort(order.begin(), order.end(), [&](Size l, Size r) { return cluster_score[l] > cluster_score[r]; });
        std::vector<RankedItem<double>> ranked;
        ranked.reserve(order.size());
        // the both-termini proteoforms with the class's decoy weight (as their PrSMs above). Not its pseudo-count: a group passes nothing at 1%
        // below 99 targets above its first decoy, and the class shares its group with the precursor-ladder class in most runs (held-out seed
        // 2: 42 -> 0 proteoforms of "undetermined, both termini + precursor ladder" in a CRC lysate run, 864 fewer proteoforms in all)
        for (Size k : order)
        {
          const bool bt = scoring_model_ && raw_loader_ && cluster_class[k] == pf_class_both_termini; // with the learned scoring only, as above
          ranked.push_back({-cluster_score[k], decoys[clusters[k][0]], decoy_weights[clusters[k][0]] * (bt ? undetermined_decoy_weight_ : 1.0)});
        }
        const double group_pseudo = 0.0;
        const auto proteoform_qvalues = targetDecoyQvalues(ranked, group_pseudo);
        // Decoy precursors model wrong MS1 masses, which only a proteoform whose mass came from MS1 can have: they count against those, and a
        // proteoform without a mass or with a mass from fragment pairs takes its q-value from the same ranking without them (identical when
        // there are no decoy precursors).
        std::vector<RankedItem<double>> ranked_no_dp;
        std::vector<Size> no_dp_pos(order.size(), (Size)-1);
        for (Size k = 0; k < order.size(); k++)
        {
          const Size rep = clusters[order[k]][0];
          if (decoys[rep] && proteoform_hits_[rep].metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(proteoform_hits_[rep].getAccession(), "DECOY")) continue;
          no_dp_pos[k] = ranked_no_dp.size();
          ranked_no_dp.push_back(ranked[k]);
        }
        const bool with_dp = ranked_no_dp.size() < ranked.size();
        const auto proteoform_qvalues_no_dp = with_dp ? targetDecoyQvalues(ranked_no_dp, group_pseudo) : proteoform_qvalues;
        // and by blind group, as in the competition: a proteoform with an MS1 mass and no blind modification counts only the decoy
        // precursors without one, one with a blind modification only those with one (the targets and the database decoys as they are)
        const auto blindGroup = [this, &clusters, &order](Size k) { return blindModCount_(proteoform_hits_[clusters[order[k]][0]]) > 0 ? 1 : 0; };
        std::array<std::vector<RankedItem<double>>, 2> ranked_dp;
        std::array<std::vector<Size>, 2> dp_pos {std::vector<Size>(order.size(), (Size)-1), std::vector<Size>(order.size(), (Size)-1)};
        for (Size k = 0; with_dp && k < order.size(); k++)
        {
          const bool dp = no_dp_pos[k] == (Size)-1;
          for (int g = 0; g < 2; g++)
          {
            if (dp && blindGroup(k) != g) continue;
            dp_pos[g][k] = ranked_dp[g].size();
            ranked_dp[g].push_back(ranked[k]);
          }
        }
        std::array<std::vector<double>, 2> proteoform_qvalues_dp;
        for (int g = 0; with_dp && g < 2; g++) proteoform_qvalues_dp[g] = targetDecoyQvalues(ranked_dp[g], group_pseudo);
        auto& counts = class_counts[c];
        counts = {0, 0, 0};
        for (Size k = 0; k < order.size(); k++)
        {
          const Size rep = clusters[order[k]][0];
          const bool ms1_mass = (double)proteoform_hits_[rep].getMetaValue("Mass", -1.0) > 0
                                && (int)proteoform_hits_[rep].getMetaValue("ProteoformMassByFragmentMass", 0) == 0;
          double pq = proteoform_qvalues[k]; // without decoy precursors the only ranking
          if (with_dp && ms1_mass) pq = proteoform_qvalues_dp[blindGroup(k)][dp_pos[blindGroup(k)][k]];
          else if (with_dp && no_dp_pos[k] != (Size)-1) pq = proteoform_qvalues_no_dp[no_dp_pos[k]];
          for (Size i : clusters[order[k]])
          {
            proteoform_hits_[i].setMetaValue("proqvalue", pq);
            proteoform_hits_[i].setMetaValue("ProteoformFDRClass", group_names[c]);
          }
          if (decoys[rep]) counts[1]++;
          else { counts[0]++; counts[2] += pq <= rep_q; }
        }
      }
      for (const auto& [c, counts] : class_counts)
      {
        if (! log) break;
        OPENMS_LOG_INFO << "Proteoform class \"" << group_names[c] << "\": " << counts[0] << " target and " << counts[1] << " decoy proteoforms, "
                        << counts[2] << " targets at q <= " << rep_q << std::endl;
      }
      // the open-terminus class, not counted (above)
      for (Size k = 0; k < clusters.size(); k++)
      {
        if (cluster_class[k] != pf_class_open_terminus) continue;
        for (Size i : clusters[k])
        {
          proteoform_hits_[i].setMetaValue("proqvalue", 1.0);
          proteoform_hits_[i].setMetaValue("ProteoformFDRClass", std::string(pf_class_names[pf_class_open_terminus]));
        }
      }
      if (log && class_size[pf_class_open_terminus] > 0)
      {
        OPENMS_LOG_INFO << "Proteoform class \"" << pf_class_names[pf_class_open_terminus] << "\": " << class_targets[pf_class_open_terminus]
                        << " target and " << class_size[pf_class_open_terminus] - class_targets[pf_class_open_terminus]
                        << " decoy proteoforms, not counted (q-value 1)" << std::endl;
      }
    }
    else if (use_fdr)
      // -native: all representatives (targets and decoys of every class) in one ranking by their PrSM q-value, as before the proteoform
      // classes (a native run has a few proteoforms per class at most). Every PrSM of a proteoform
      // carries the proteoform's q-value, which can be lower than its own PrSM q-value.
      {
        std::vector<Size> order(clusters.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](Size l, Size r) {
          const Size rep_l = clusters[l][0], rep_r = clusters[r][0];
          if (qvalues[rep_l] != qvalues[rep_r]) return qvalues[rep_l] < qvalues[rep_r];
          return proteoform_hits_[rep_l].getScore() > proteoform_hits_[rep_r].getScore();
        });
        std::vector<RankedItem<std::pair<double, double>>> ranked;
        ranked.reserve(order.size());
        for (Size k : order)
        {
          const Size rep = clusters[k][0];
          ranked.push_back({{qvalues[rep], -proteoform_hits_[rep].getScore()}, decoys[rep], decoy_weights[rep]});
        }
        const auto proteoform_qvalues = targetDecoyQvalues(ranked);
        // Decoy precursors model wrong MS1 masses, which only a proteoform whose mass came from MS1 can have: they count against those, and a
        // proteoform without a mass or with a mass from fragment pairs takes its q-value from the same ranking without them (identical when
        // there are no decoy precursors).
        std::vector<RankedItem<std::pair<double, double>>> ranked_no_dp;
        std::vector<Size> no_dp_pos(order.size(), (Size)-1);
        for (Size k = 0; k < order.size(); k++)
        {
          const Size rep = clusters[order[k]][0];
          if (decoys[rep] && proteoform_hits_[rep].metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(proteoform_hits_[rep].getAccession(), "DECOY")) continue;
          no_dp_pos[k] = ranked_no_dp.size();
          ranked_no_dp.push_back(ranked[k]);
        }
        const auto proteoform_qvalues_no_dp = ranked_no_dp.size() < ranked.size() ? targetDecoyQvalues(ranked_no_dp) : proteoform_qvalues;
        for (Size k = 0; k < order.size(); k++)
        {
          const Size rep = clusters[order[k]][0];
          const bool ms1_mass = (double)proteoform_hits_[rep].getMetaValue("Mass", -1.0) > 0
                                && (int)proteoform_hits_[rep].getMetaValue("ProteoformMassByFragmentMass", 0) == 0;
          const double pq = (ms1_mass || no_dp_pos[k] == (Size)-1) ? proteoform_qvalues[k] : proteoform_qvalues_no_dp[ranked_no_dp.size() < ranked.size() ? no_dp_pos[k] : k];
          for (Size i : clusters[order[k]]) proteoform_hits_[i].setMetaValue("proqvalue", pq);
        }
      }
    return clusters;
  };

  if (! protein_level_only_)
  {
    // Beyond the min_protein_hits best candidates of a spectrum (a fair target-decoy competition: they are chosen by this spectrum
    // alone), extend only those whose protein, or its target-decoy partner, passes the run-level protein FDR. That needs decoys;
    // without them every candidate is extended, as before.
    // In the search cascade, it is recomputed before each stage's extension from the tag candidates of the stages so far.
    const bool use_protein_fdr = decoy_factor_ > 0 && protein_fdr_ < 1.0;
    std::set<Size> passing_proteins; // protein indices
    if (use_protein_fdr) { passing_proteins = getProteinsPassingFDR_(*db); }

    // spectra in parallel, each with its own extender; proteoforms are appended in spectrum order afterwards
    std::vector<std::vector<ProteinHit>> spec_proteoforms(dspecs.size());
    // as above: over spectra only with enough spectra to extend, otherwise over each spectrum's candidates inside the extender
    auto hasCandidates = [&](int d) {
      return candidate_ranges[d].first != candidate_ranges[d].second || ladder_ranges[d].first != ladder_ranges[d].second
             || pl_ranges[d].first != pl_ranges[d].second;
    };
    auto overSpectra = [&](const std::vector<int>& spectra) {
      [[maybe_unused]] const auto with_candidates = std::count_if(spectra.begin(), spectra.end(), hasCandidates);
#ifdef _OPENMP
      return with_candidates >= omp_get_max_threads();
#else
      return false;
#endif
    };
    const bool keep_all_candidates = (bool)raw_loader_;
    // precursor-ladder PrSMs of the extension that is kept: kept (determined; of them at both window ends, at the window start only, at the
    // window end only), dropped as the same proteoform as a determined PrSM of another source, dropped as undetermined
    std::atomic<Size> pl_determined {0}, pl_at_both_ends {0}, pl_at_start {0}, pl_at_end {0}, pl_duplicates {0}, pl_undetermined {0};
    // A PrSM of a database entry (by_entry) is decoded to its protein right after the extension: its positions (and TagPositions) moved to
    // the full-length protein, the entry's fixed modifications inside the proteoform merged into its modifications by position
    // (ModificationVariable 2; BlindModCount and VariableModifications stay those of the search), the sequence the full-length protein, and
    // the entry's processed form as EntryStart and EntryEnd (1-based). FastaIndex stays the entry. At an open terminus (-1) the proteoform
    // runs to the entry's end, but a modification of that terminus (kind 2 or 3) is not appended: no terminus is determined there for it.
    // A proteoform with one half cystine of a bond of its entry (cut between them; the extension can truncate an entry anywhere) is
    // dropped (decode returns false): its residue masses had half the bond. The extension never sees a decoded PrSM: every pass extends
    // the candidates. The meta values are set by index (registered here, outside the parallel loops).
    const UInt mi_start_position = registry.registerName("StartPosition"), mi_end_position = registry.registerName("EndPosition"),
               mi_modifications = registry.registerName("Modifications"), mi_mod_starts = registry.registerName("ModificationStarts"),
               mi_mod_ends = registry.registerName("ModificationEnds"), mi_mod_ids = registry.registerName("ModificationIDs"),
               mi_mod_accs = registry.registerName("ModificationACCs"), mi_mod_variable = registry.registerName("ModificationVariable"),
               mi_entry_start = registry.registerName("EntryStart"), mi_entry_end = registry.registerName("EntryEnd");
    std::atomic<Size> cut_bonds {0}; // PrSMs dropped by decode
    auto decode = [&](ProteinHit& h) -> bool {
      const auto& entry = db->entries[(Size)(int)h.getMetaValue(mi_fasta_index)];
      const int offset = entry.offset;
      const int start = h.getMetaValue(mi_start_position), end = h.getMetaValue(mi_end_position);
      const int first = start > 0 ? start - 1 : 0, last = end > 0 ? end : entry.length; // the proteoform's residues [first, last) of the entry
      for (const auto& fm : entry.fixed)
      {
        if (fm.bond >= 0 && fm.position >= first && fm.position < last && (fm.bond < first || fm.bond >= last)) return false;
      }
      std::vector<double> masses = h.getMetaValue(mi_modifications);
      std::vector<int> starts = h.getMetaValue(mi_mod_starts).toIntList(), ends = h.getMetaValue(mi_mod_ends).toIntList();
      std::vector<std::string> ids = h.getMetaValue(mi_mod_ids).toStringList(), accs = h.getMetaValue(mi_mod_accs).toStringList();
      std::vector<int> variable = h.getMetaValue(mi_mod_variable).toIntList();
      for (auto& v : starts) v += offset;
      for (auto& v : ends) v += offset;
      const Size searched = masses.size();
      for (const auto& fm : entry.fixed)
      {
        if (fm.position < first || fm.position >= last) continue;
        if ((fm.kind == 2 && start <= 0) || (fm.kind == 3 && end <= 0)) continue;
        masses.push_back(fm.mass);
        starts.push_back(offset + fm.position);
        ends.push_back(offset + fm.position);
        std::string name = fm.name;
        StringUtils::substitute(name, ',', '_'); // a comma separates the alternatives of a modification ID
        ids.push_back(name);
        accs.push_back(fm.accession);
        variable.push_back(2);
      }
      if (masses.size() > searched && masses.size() == starts.size() && masses.size() == ends.size() && masses.size() == ids.size()
          && masses.size() == accs.size() && masses.size() == variable.size())
      {
        // all modifications by position (the searched ones first at a residue)
        std::vector<Size> order(masses.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&starts](Size l, Size r) { return starts[l] < starts[r]; });
        auto permuted = [&order](const auto& v) {
          std::remove_const_t<std::remove_reference_t<decltype(v)>> out;
          out.reserve(v.size());
          for (Size k : order) out.push_back(v[k]);
          return out;
        };
        masses = permuted(masses);
        starts = permuted(starts);
        ends = permuted(ends);
        ids = permuted(ids);
        accs = permuted(accs);
        variable = permuted(variable);
      }
      h.setMetaValue(mi_modifications, masses);
      h.setMetaValue(mi_mod_starts, starts);
      h.setMetaValue(mi_mod_ends, ends);
      h.setMetaValue(mi_mod_ids, ids);
      h.setMetaValue(mi_mod_accs, accs);
      h.setMetaValue(mi_mod_variable, variable);
      if (start > 0) h.setMetaValue(mi_start_position, start + offset);
      if (end > 0) h.setMetaValue(mi_end_position, end + offset);
      if (offset > 0 && h.metaValueExists(mi_tag_positions))
      {
        std::vector<int> positions = h.getMetaValue(mi_tag_positions).toIntList();
        for (auto& v : positions) v += offset;
        h.setMetaValue(mi_tag_positions, positions);
      }
      h.setSequence(protein_full_seqs[entry.protein]);
      h.setMetaValue(mi_entry_start, offset + 1);
      h.setMetaValue(mi_entry_end, offset + entry.length);
      return true;
    };
    // The extension of one spectrum's candidates with the given extension parameters: the sequence- and gapped-tag candidates
    // tag_candidates[tag_range], the terminal-ladder candidates ladder_candidates_now[ladder_range] and the precursor-ladder candidates
    // pl_candidates_now[pl_range]. With keep_all, every candidate's PrSM is kept (multiple_hits_per_spec_, raw-spectrum scoring, the
    // residual rounds of the chimeric search); otherwise the extension keeps one.
    // the raw-evidence score of path alternatives (path_alternatives; set below, with the raw spectra loaded)
    std::function<double(const ProteinHit&)> alternative_scorer;
    auto extendCandidatesInto = [&](const Param& eparam, const DeconvolvedSpectrum& dspec, const std::vector<ProteinHit>& tag_candidates,
                                    std::pair<Size, Size> tag_range, const std::vector<ProteinHit>& ladder_candidates_now, std::pair<Size, Size> ladder_range,
                                    const std::vector<ProteinHit>& pl_candidates_now, std::pair<Size, Size> pl_range, std::vector<ProteinHit>& out,
                                    bool keep_all) {
      out.clear();
      const auto [first, second] = tag_range;
      const auto [pl_first, pl_second] = pl_range;
      if (first == second && ladder_range.first == ladder_range.second && pl_first == pl_second) return;

      std::vector<ProteinHit> hits;
      // the min_protein_hits best candidates of each origin (sequence tags, gapped tags) are always extended, so gapped candidates
      // add chances without taking the places of sequence-tag candidates
      int rank_tag = 0, rank_gapped = 0;
      bool top = false; // of the protein of the last counted candidate: the other entries of a kept protein follow it (by_entry)
      for (Size i = first; i < second; i++)
      {
        const bool gapped = tag_candidates[i].metaValueExists(mi_gapped);
        // native: the decoy margin chose the candidates
        if (! tag_candidates[i].metaValueExists(mi_entry_of_kept)) top = (gapped ? rank_gapped++ : rank_tag++) < min_protein_hits_ || native_rule;
        if (use_protein_fdr && ! top && passing_proteins.count(protein_of[(Size)(int)tag_candidates[i].getMetaValue(mi_fasta_index)]) == 0) continue;
        hits.push_back(tag_candidates[i]);
      }

      std::vector<int> spec_vec;
      spec_vec.reserve(dspec.size() + 1);
      spec_vec.push_back(0);
      for (const auto& pg : dspec)
      {
        spec_vec.push_back(int(round(pg.getMonoMass())));
      }

      // Terminal ladder candidates are extended on their own, with their ladder tags only. A ladder of a few masses reads as an
      // undetermined proteoform on any protein: mixed with the others, such PrSMs (2226 decoys vs 8 without them on Eclip_1480) raised the
      // score threshold of the undetermined proteoforms from 2 to 33, and ladder tags added to the sequence tag candidates still raised it
      // to 22. So a ladder PrSM of determined mass and termini competes with the others (the better PrSM is kept as one extension would
      // keep it: the best score, the earlier candidate among equal ones), and an undetermined one is reported only for a spectrum without
      // any other PrSM, in its own FDR class.
      std::vector<ProteinHit> ladder_hits(ladder_candidates_now.begin() + ladder_range.first, ladder_candidates_now.begin() + ladder_range.second);
      // several entries of a protein extended (by_entry): its PrSM is the one of determined mass, then fewer blind modifications, then the
      // higher score, then more fixed (annotated) modifications; the others are dropped. The same rule for targets and decoys (decoy
      // entries carry the annotations moved with their residues) and for the decoy precursor copy.
      auto blindCount = [](const ProteinHit& h) { return h.metaValueExists("BlindModCount") ? (int)h.getMetaValue("BlindModCount") : 0; };
      auto fixedCount = [&](const ProteinHit& h) {
        if (! h.metaValueExists(mi_mod_variable)) return 0;
        const auto v = h.getMetaValue(mi_mod_variable).toIntList();
        return (int)std::count(v.begin(), v.end(), 2);
      };
      auto determined = [](const ProteinHit& h) { return h.metaValueExists("Mass") && (double)h.getMetaValue("Mass") > 0; };
      auto better = [&](const ProteinHit& a, const ProteinHit& b) {
        if (determined(a) != determined(b)) return determined(a);
        if (blindCount(a) != blindCount(b)) return blindCount(a) < blindCount(b);
        if (a.getScore() != b.getScore()) return a.getScore() > b.getScore();
        return fixedCount(a) > fixedCount(b);
      };
      auto onePerProtein = [&](std::vector<ProteinHit>& proteoforms) {
        std::vector<ProteinHit> one_per_protein;
        one_per_protein.reserve(proteoforms.size());
        for (auto& h : proteoforms)
        {
          auto it = std::find_if(one_per_protein.begin(), one_per_protein.end(),
                                 [&h](const ProteinHit& o) { return o.getAccession() == h.getAccession(); });
          if (it == one_per_protein.end()) one_per_protein.push_back(std::move(h));
          else if (better(h, *it)) *it = std::move(h);
        }
        proteoforms.swap(one_per_protein);
      };
      auto extendAll = [&](std::vector<ProteinHit>& candidates, std::vector<ProteinHit>& proteoforms) {
        proteoforms.clear();
        if (candidates.empty()) return;
        FLASHExtenderAlgorithm extender;
        extender.setParameters(eparam);
        extender.setCandidateBlindModificationMap(blind_mod_map);
        extender.setFixedModifications(db);
        // the chimeric search removes the masses of a spectrum's best PrSM before the next round
        extender.setReportMatchedMasses(chimeric_rounds_ > 0);
        extender.setPathAlternatives(alternative_scorer, path_alternatives_, path_alternative_margin);
        // with raw-spectrum scoring, every candidate's PrSM is kept here and the best one per spectrum is chosen by the combined score in
        // annotateOnRawSpectra_: the extension score alone (deconvolved masses) picked another protein in 1 of 15 spectra of Eclip_1480
        // where TopPIC's protein had the better raw evidence (offline: 5252 -> 5451 PrSMs, 1435 -> 1500 proteoforms at 1%)
        extender.run(candidates, dspec, spec_vec, vec_pro, rev_vec_pro, tags_, tol, keep_all);
        extender.takeProteoforms(proteoforms);
        if (by_entry)
        {
          Size kept = 0;
          for (Size k = 0; k < proteoforms.size(); k++)
          {
            if (! decode(proteoforms[k]))
            {
              cut_bonds++;
              continue;
            }
            if (kept != k) proteoforms[kept] = std::move(proteoforms[k]);
            kept++;
          }
          proteoforms.resize(kept);
          onePerProtein(proteoforms);
        }
      };
      // With several entries per protein, the protein's ranked entry is extended first, and its other entries (EntryOfKeptProtein) only
      // if that gave a PrSM another entry may better: an undetermined mass or a blind modification (an annotation can explain a shift).
      // A settled PrSM (determined, no blind modification) is not bettered by the rule above except by the score among equally settled
      // ones; a ranked entry without any PrSM rarely has an entry that has one (Eclip_1480: 1529 -> 1526 proteoforms, 54 -> 44 s when
      // these are skipped too). With one PrSM per spectrum the extension compares all entries at once.
      auto extend = [&](std::vector<ProteinHit>& candidates, std::vector<ProteinHit>& proteoforms) {
        if (! by_entry || ! keep_all)
        {
          extendAll(candidates, proteoforms);
          return;
        }
        std::vector<ProteinHit> first, later;
        for (auto& c : candidates) (c.metaValueExists(mi_entry_of_kept) ? later : first).push_back(std::move(c));
        extendAll(first, proteoforms);
        if (later.empty()) return;
        // the proteins whose ranked entry gave a PrSM that another entry may better: an undetermined mass or a blind modification
        std::unordered_set<std::string> open;
        for (const auto& h : proteoforms)
          if (! determined(h) || blindCount(h) > 0) open.insert(h.getAccession());
        std::vector<ProteinHit> rest;
        for (auto& c : later)
          if (open.count(c.getAccession()) > 0) rest.push_back(std::move(c));
        if (rest.empty()) return;
        std::vector<ProteinHit> more;
        extendAll(rest, more);
        for (auto& h : more) proteoforms.push_back(std::move(h));
        onePerProtein(proteoforms);
      };
      extend(hits, out);
      // native: every candidate's PrSM is reported, one per protein (the better one)
      auto mergeNative = [&out](std::vector<ProteinHit>& from) {
        for (auto& h : from)
        {
          auto it = std::find_if(out.begin(), out.end(), [&h](const ProteinHit& o) { return o.getAccession() == h.getAccession(); });
          if (it == out.end()) { out.push_back(std::move(h)); }
          else if (h.getScore() > it->getScore()) { *it = std::move(h); }
        }
      };
      auto is_determined = [](const ProteinHit& h) {
        return (double)h.getMetaValue(tntKey("Mass")) > 0 && (int)h.getMetaValue(tntKey("StartPosition")) > 0 && (int)h.getMetaValue(tntKey("EndPosition")) > 0;
      };
      if (! ladder_hits.empty())
      {
        std::vector<ProteinHit> ladder_out;
        extend(ladder_hits, ladder_out);
        std::vector<ProteinHit> determined, undetermined;
        for (auto& h : ladder_out)
        {
          h.setMetaValue(native_rule ? mi_sliding : mi_terminal, 1);
          (is_determined(h) ? determined : undetermined).push_back(h);
        }
        const bool fill = out.empty();
        if (native_rule)
        {
          mergeNative(determined);
          mergeNative(undetermined);
        }
        else if (keep_all) { out.insert(out.end(), determined.begin(), determined.end()); }
        else if (! determined.empty() && (out.empty() || ! is_determined(out[0]) || determined[0].getScore() > out[0].getScore()))
        {
          // a ladder PrSM that explains the precursor mass also replaces an undetermined one of a higher score
          out.swap(determined);
        }
        if (native_rule) {}
        else if (fill && out.empty()) out.swap(undetermined);
        else if (fill && keep_all) { out.insert(out.end(), undetermined.begin(), undetermined.end()); }
      }
      // Precursor-ladder candidates: extended on their own, as the terminal-ladder ones. Only their PrSMs of determined mass and termini
      // are kept (the same rule for targets and decoys): the source proposes a window that explains the precursor mass, and an
      // undetermined PrSM from it means the extension could not; it would only add chance PrSMs to a class meant to be determined.
      // A kept PrSM that is the same proteoform (same accession, areConsistent_) as a determined PrSM another source found in this spectrum
      // is dropped too: the proteoform stays with that source and its class, so the precursor-ladder class holds what this source adds,
      // and correct duplicates neither pad the class meant to show this source's own chance rate nor count twice where several PrSMs per
      // spectrum are kept. The rule looks at no label, so it treats targets and decoys alike.
      // With one PrSM per spectrum, a precursor-ladder PrSM replaces the spectrum's PrSM only with a higher score. A terminal-ladder PrSM
      // also replaces an undetermined one of a higher score, because an intact-terminus ladder is strong evidence; a window with a few
      // matched positions among the windows that fit a precursor mass (about 4e4 per spectrum on a human-size database) is not, and such
      // a rule would replace most undetermined
      // PrSMs of the spectra with a precursor mass by chance windows (harness: rank-1 windows in 92 of 190 E. coli spectra, 26 of them
      // decoys). Being last, they lose the ties (here, and in the best PrSM per spectrum by the combined score: determined, then coverage,
      // then the earlier PrSM), so a PrSM is in their class only if their source gave the better one.
      if (pl_first != pl_second)
      {
        std::vector<ProteinHit> pl_hits(pl_candidates_now.begin() + pl_first, pl_candidates_now.begin() + pl_second), pl_out, pl_kept;
        extend(pl_hits, pl_out);
        Size undetermined = 0, duplicates = 0, at_both_ends = 0, at_start = 0, at_end = 0;
        for (auto& h : pl_out)
        {
          h.setMetaValue(mi_precursor_ladder, 1);
          if (! is_determined(h))
          {
            undetermined++;
            continue;
          }
          const bool duplicate = std::any_of(out.begin(), out.end(), [&](const ProteinHit& o) {
            return is_determined(o) && o.getAccession() == h.getAccession() && areConsistent_(o, h, pl_precursor_tol);
          });
          if (duplicate)
          {
            duplicates++;
            continue;
          }
          const bool start_ok = (int)h.getMetaValue(tntKey("StartPosition")) == (int)h.getMetaValue(mi_pl_start);
          const bool end_ok = (int)h.getMetaValue(tntKey("EndPosition")) == (int)h.getMetaValue(mi_pl_end);
          at_both_ends += start_ok && end_ok;
          at_start += start_ok && ! end_ok;
          at_end += end_ok && ! start_ok;
          pl_kept.push_back(std::move(h));
        }
        pl_determined += pl_kept.size();
        pl_at_both_ends += at_both_ends;
        pl_at_start += at_start;
        pl_at_end += at_end;
        pl_duplicates += duplicates;
        pl_undetermined += undetermined;
        if (native_rule) { mergeNative(pl_kept); }
        else if (keep_all) { out.insert(out.end(), pl_kept.begin(), pl_kept.end()); }
        else if (! pl_kept.empty() && (out.empty() || pl_kept[0].getScore() > out[0].getScore())) { out.swap(pl_kept); }
      }
    };
    // The extension of spectrum d's candidates. With decoy_mass, the spectrum is the decoy copy at the decoy precursor mass
    // (decoy_precursor), and its precursor-ladder candidates are those of the windows that fit that mass.
    auto extendSpectrumInto = [&](int d, const Param& eparam, const DeconvolvedSpectrum& dspec, std::vector<ProteinHit>& out, bool decoy_mass) {
      extendCandidatesInto(eparam, dspec, protein_hits_, candidate_ranges[d], ladder_candidates, ladder_ranges[d],
                           decoy_mass ? decoy_pl_candidates : pl_candidates, decoy_mass ? decoy_pl_ranges[d] : pl_ranges[d], out,
                           multiple_hits_per_spec_ || keep_all_candidates);
    };
    // Extends the given spectra with the given parameters. Each spectrum is extended on its own and writes only its own slot, so the result
    // depends neither on the thread count and schedule nor on which other spectra are extended in the same call.
    auto extendSpectra = [&](const std::vector<int>& spectra, const Param& eparam) {
      const bool over_spectra = overSpectra(spectra);
#pragma omp parallel for schedule(dynamic) if (over_spectra)
      for (int k = 0; k < (int)spectra.size(); k++)
      {
        nextProgress();
        const int d = spectra[k];
        extendSpectrumInto(d, eparam, dspecs[d], spec_proteoforms[d], false);
      }
    };

    // path alternatives are chosen on the raw spectra: all of them are loaded before the first extension (annotateOnRawSpectra_ uses them)
    if (raw_loader_ && path_alternatives_ > 0)
    {
      std::set<int> scans;
      for (const auto& dspec : dspecs) scans.insert(dspec.getScanNumber());
      loadRawSpectra_(scans);
      alternative_scorer = alternativeRawScorer_();
    }

    // An undetermined PrSM (no mass) with an open terminus has no precursor mass to constrain a blind shift: the shift takes up whatever
    // the sequence of a similar protein differs by. On the B. theta entrapment of the three lysates, 7 of the 8 B. theta proteoforms of the
    // class at 1% had one, against 11 of the 58 E. coli ones; without such PrSMs the class held 1 B. theta of 58, and the three lysates gave
    // 85 more proteoforms (the spectra's next PrSMs). Targets and decoys alike.
    auto dropOpenTerminusShifts = [&](std::vector<ProteinHit>& hits) {
      const auto shifted = [&](const ProteinHit& h) {
        if ((double)h.getMetaValue(tntKey("Mass"), -1.0) > 0) return false;
        if ((int)h.getMetaValue(tntKey("StartPosition"), -1) > 0 && (int)h.getMetaValue(tntKey("EndPosition"), -1) > 0) return false;
        return blindModCount_(h) > 0;
      };
      const Size before = hits.size();
      hits.erase(std::remove_if(hits.begin(), hits.end(), shifted), hits.end());
      if (hits.size() < before) OPENMS_LOG_INFO << "Undetermined PrSMs with an open terminus and a blind shift: " << before - hits.size() << " dropped" << std::endl;
    };

    // The PrSMs of the spectra of one stage (@p spectra; the first stage: every spectrum), appended to proteoform_hits_ in spectrum order.
    const UInt mi_stage = registry.registerName("SearchStage");
    auto extendStage = [&](int stage, const std::vector<int>& spectra) {
      const bool first_stage = stage == stages[0];
      startProgress(0, (SignedSize)spectra.size(),
                    first_stage ? std::string("Running extension algorithm ...") : "Running extension algorithm (search cascade stage " + StringUtils::toStr(stage) + ") ...");
      extendSpectra(spectra, extender_param_);
      const bool extend_over_spectra = overSpectra(spectra);
      if (use_precursor_ladders && (stage == 0 || stage == 1))
      {
        OPENMS_LOG_INFO << "Precursor-ladder PrSMs: " << pl_determined.load() << " of determined mass and termini kept (" << pl_at_both_ends.load()
                        << " at both window ends, " << pl_at_start.load() << " at the window start only, " << pl_at_end.load()
                        << " at the window end only), " << pl_duplicates.load() << " dropped as the same proteoform as a determined PrSM of another source, "
                        << pl_undetermined.load() << " undetermined dropped" << std::endl;
      }
      // Decoy precursors: a spectrum whose PrSM took its mass from MS1 is searched again as a decoy spectrum whose precursor mass is the
      // decoy precursor mass (noise or signal decoy) of its isolation window, as FLASHDeconv chose it with report_FDR. A PrSM of that copy
      // that reports the decoy mass (not a mass from fragment pairs) is a decoy of the precursor: it competes with the spectrum's MS1-based
      // PrSMs (the higher score stays) and counts as a decoy in the FDR. PrSMs with masses from fragment pairs or without a mass (-1) are
      // not affected.
      if (decoy_precursor_)
      {
        const double iso = Constants::C13C12_MASSDIFF_U;
        auto ms1Based = [](const ProteinHit& h) {
          return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("ProteoformMassByFragmentMass") == 0;
        };
        const Param& eparam_now = extender_param_;
        std::vector<std::vector<ProteinHit>> decoy_proteoforms(dspecs.size());
        int copies = 0, without_decoy = 0;
        std::vector<int> todo;
        for (int d : spectra)
        {
          if (dspecs[d].getPrecursorPeakGroup().getMonoMass() <= 0) continue;
          if (std::none_of(spec_proteoforms[d].begin(), spec_proteoforms[d].end(), ms1Based)) continue;
          if (dspecs[d].getDecoyPrecursorMass() <= 0)
          {
            without_decoy++;
            continue;
          }
          todo.push_back(d);
        }
        copies = (int)todo.size();
        if (copies == 0 && first_stage)
        {
          OPENMS_LOG_WARN << "Decoy precursors: the input has no decoy precursor masses; run FLASHDeconv with -FD:report_FDR true" << std::endl;
        }
  #pragma omp parallel for schedule(dynamic) if (extend_over_spectra)
        for (int k = 0; k < (int)todo.size(); k++)
        {
          const int d = todo[k];
          DeconvolvedSpectrum decoy_spec = dspecs[d];
          auto pg = decoy_spec.getPrecursorPeakGroup();
          const double decoy_mass = decoy_spec.getDecoyPrecursorMass();
          pg.setMonoisotopicMass(decoy_mass);
          decoy_spec.setPrecursorPeakGroup(pg);
          std::vector<ProteinHit> out;
          extendSpectrumInto(d, eparam_now, decoy_spec, out, true);
          for (auto& h : out)
          {
            if (! ms1Based(h)) continue;
            if (std::abs((double)h.getMetaValue("Mass") - decoy_mass) > 3 * iso + 0.1) continue; // the moved mass, up to isotope errors
            h.setMetaValue("DecoyPrecursor", decoy_spec.getDecoyPrecursorWeight()); // its weight as a decoy (FLASHDeconv's noise decoy weight)
            decoy_proteoforms[d].push_back(std::move(h));
          }
        }
        // Per spectrum, the decoy precursor competes with the spectrum's PrSMs that took their mass from MS1, on the precursor: the decoy
        // spectrum has to give a PrSM that reports the decoy mass (the try), and that decoy wins when its precursor score (FLASHDeconv's 2D
        // Qscore, the score its own FDR compares target and decoy masses with) is higher than the precursor's. The fragments cannot decide
        // the mass: a related or wrong mass explains nearly the same ladders, only the MS1 signal tells them apart. Inputs from FLASHDeconv
        // versions that do not write the decoy precursor's score fall back to the fragment score (the higher PrSM score wins).
        // The competition is held per blind group (no blind modification, or one or more): a decoy PrSM replaces only the spectrum's
        // MS1-based PrSMs of its own group, and counts against the proteoforms of that group (the proteoform FDR below). A decoy mass that
        // fits a sequence only with a blind modification says nothing against a PrSM that matches the unmodified sequence mass, and the other
        // way round (the three lysates with decoy precursors and the MS1 noise decoy weight: of the proteoforms lost against the MS2 weight,
        // those without a blind modification matched their unmodified sequence mass within 10 ppm in 88-92% and were found in another lysate
        // in 58-72%, the ones with a blind modification in 22-31%).
        int decoy_wins = 0, decoy_tries = 0, decoy_by_fragments = 0;
        std::array<int, 2> group_wins {0, 0}; // by blind group
        const auto blindGroup = [](const ProteinHit& h) { return blindModCount_(h) > 0 ? 1 : 0; };
        for (int d : todo)
        {
          if (decoy_proteoforms[d].empty()) continue;
          decoy_tries++;
          // on the Qscores in the survey scan alone (1D) when the input has them: a decoy is rarely traced over several scans, so the 2D
          // Qscore favours the precursor (Eclip_1480: of 137 decoy spectra giving a PrSM with the decoy mass, 12 decoys beat the precursor on
          // the 1D Qscore, 6 on the 2D one). The precursor decides for both groups; without its scores, the fragments within each group.
          const double decoy_qscore1d = dspecs[d].getDecoyPrecursorQscore1D(), decoy_qscore = dspecs[d].getDecoyPrecursorQscore();
          int by_precursor = -1; // 1: the decoy precursor beats the precursor, 0: it does not, -1: no decoy precursor score
          if (decoy_qscore1d >= 0) by_precursor = decoy_qscore1d > dspecs[d].getPrecursorPeakGroup().getQscore() ? 1 : 0;
          else if (decoy_qscore >= 0) by_precursor = decoy_qscore > dspecs[d].getPrecursorPeakGroup().getQscore2D() ? 1 : 0;
          else decoy_by_fragments++;
          if (by_precursor == 0) continue;
          std::array<double, 2> best_decoy {-1e10, -1e10};
          for (const auto& h : decoy_proteoforms[d])
          {
            auto& best = best_decoy[blindGroup(h)];
            best = std::max(best, h.getScore());
          }
          std::array<bool, 2> wins {false, false};
          for (int g = 0; g < 2; g++)
          {
            if (best_decoy[g] == -1e10) continue; // no decoy PrSM of the group
            if (by_precursor < 0)
            {
              double best_own = -1e10;
              for (const auto& h : spec_proteoforms[d])
              {
                if (ms1Based(h) && blindGroup(h) == g) best_own = std::max(best_own, h.getScore());
              }
              if (best_decoy[g] <= best_own) continue;
            }
            wins[g] = true;
            group_wins[g]++;
          }
          if (! wins[0] && ! wins[1]) continue;
          decoy_wins++;
          std::vector<ProteinHit> kept;
          for (auto& h : spec_proteoforms[d])
          {
            if (! ms1Based(h) || ! wins[blindGroup(h)]) kept.push_back(std::move(h));
          }
          for (auto& h : decoy_proteoforms[d])
          {
            const int g = blindGroup(h);
            if (wins[g] && h.getScore() == best_decoy[g]) kept.push_back(std::move(h));
          }
          spec_proteoforms[d].swap(kept);
        }
        OPENMS_LOG_INFO << "Decoy precursors: " << copies << " spectra with an MS1-based PrSM and an MS1 decoy mass (" << without_decoy
                        << " without one); the decoy spectrum gave a PrSM with the decoy mass in " << decoy_tries << ", the decoy precursor wins in "
                        << decoy_wins << " (against PrSMs without a blind modification in " << group_wins[0] << ", with one in " << group_wins[1] << ")"
                        << (decoy_by_fragments > 0 ? " (" + StringUtils::toStr(decoy_by_fragments) + " decided on the fragments: no decoy precursor score in the input)" : std::string(""))
                        << std::endl;
      }

      std::vector<ProteinHit> stage_hits;
      for (int d : spectra)
      {
        for (auto& h : spec_proteoforms[d]) stage_hits.push_back(std::move(h));
        std::vector<ProteinHit>().swap(spec_proteoforms[d]);
      }
      endProgress();

      // raw-spectrum evidence of every PrSM, targets and decoys alike, before the FDR
      if (raw_loader_ && ! stage_hits.empty())
      {
        std::set<int> scans;
        for (const auto& hit : stage_hits) scans.insert((int)hit.getMetaValue("Scan"));
        loadRawSpectra_(scans);
        annotateOnRawSpectra_(stage_hits);
        dropOpenTerminusShifts(stage_hits);
        if (const char* trace_file = std::getenv("FLASH_TRACE_FEATURES"))
        {
          try { traceCandidates_(stage_hits, trace_file); }
          catch (const Exception::BaseException& e) { OPENMS_LOG_WARN << "FLASH_TRACE_FEATURES failed: " << e.what() << " (" << e.getName() << ")" << std::endl; }
          catch (const std::exception& e) { OPENMS_LOG_WARN << "FLASH_TRACE_FEATURES failed: " << e.what() << std::endl; }
        }
        // one PrSM per spectrum: the best by the combined score, or with a scoring model by its choice score (its score then its confidence);
        // the extension kept every candidate's PrSM for this choice
        if (! multiple_hits_per_spec_)
        {
          if (scoring_model_)
          {
            std::map<int, const MSSpectrum*> spectra;
            for (Size index = 0; index < map.size(); index++)
              if (map[index].getMSLevel() > 1) spectra[FLASHDeconvAlgorithm::getScanNumber(map, index)] = &map[index];
            applyLearnedScoring_(stage_hits, spectra);
          }
          else keepBestPerSpectrum_(stage_hits, stage == 0 ? std::string("") : " (search cascade stage " + StringUtils::toStr(stage) + ")", "combined score");
        }
      }
      else { dropOpenTerminusShifts(stage_hits); }
      for (auto& h : stage_hits)
      {
        if (stage > 0) h.setMetaValue(mi_stage, stage);
        proteoform_hits_.push_back(std::move(h));
      }
    };

    // The spectra identified after a stage of the search cascade: whose best PrSM (by score, of any stage so far) is a target whose
    // proteoform q-value is at most rep_q, the level at which proteoforms are reported. A passing target below a better decoy or a better
    // failing target does not identify its spectrum (survey benchmark: 187 gapped-tag identifications had a decoy best PrSM in stage 1 and
    // were never searched with gapped tags). The q-values are those the FDR below gives the PrSMs so far (per class and stage, the
    // proteoforms over all of them), computed on a copy; the one proteoform per MS1 feature rule is applied at the end only. An identified
    // spectrum is searched by no later stage.
    std::map<int, int> spectrum_of_scan_number;
    for (int d = 0; d < (int)dspecs.size(); d++) spectrum_of_scan_number[dspecs[d].getScanNumber()] = d;
    auto identifiedSoFar = [&]() {
      std::vector<ProteinHit> kept = proteoform_hits_;
      std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), byScore);
      std::vector<int> hit_class;
      std::vector<double> qvalues, decoy_weights;
      std::vector<bool> decoys;
      proteoformQvalues(hit_class, qvalues, decoys, decoy_weights, false);
      std::set<int> identified, seen;
      for (Size i = 0; i < proteoform_hits_.size(); i++)
      {
        const auto& hit = proteoform_hits_[i];
        const auto it = spectrum_of_scan_number.find((int)hit.getMetaValue(mi_scan));
        // the spectrum's best PrSM only (proteoform_hits_ is in byScore order here)
        if (it == spectrum_of_scan_number.end() || ! seen.insert(it->second).second) continue;
        if (decoys[i] || hit_class[i] < 0 || (double)hit.getMetaValue("proqvalue", 1.0) > rep_q) continue;
        identified.insert(it->second);
      }
      proteoform_hits_.swap(kept);
      return identified;
    };

    std::vector<int> stage_spectra = all_spectra;
    for (Size s = 0; s < stages.size(); s++)
    {
      const int stage = stages[s];
      if (s > 0)
      {
        startProgress(0, (SignedSize)stage_spectra.size(), "Running candidate protein filtration (search cascade stage " + StringUtils::toStr(stage) + ") ...");
        searchStage(stage, stage_spectra);
        endProgress();
        if (use_protein_fdr) passing_proteins = getProteinsPassingFDR_(*db);
      }
      const Size from = proteoform_hits_.size();
      extendStage(stage, stage_spectra);
      if (stage == 0) continue;
      if (s + 1 == stages.size())
      {
        OPENMS_LOG_INFO << "Search cascade stage " << stage << " (" << (stage == 1 ? "sequence tags and ladders" : "gapped tags")
                        << "): " << stage_spectra.size() << " spectra searched, " << proteoform_hits_.size() - from << " PrSMs" << std::endl;
        continue;
      }
      const auto identified = identifiedSoFar();
      std::vector<int> rest;
      for (int d : stage_spectra)
      {
        if (identified.count(d) == 0) rest.push_back(d);
      }
      OPENMS_LOG_INFO << "Search cascade stage " << stage << " (" << (stage == 1 ? "sequence tags and ladders" : "gapped tags") << "): " << stage_spectra.size()
                      << " spectra searched, " << proteoform_hits_.size() - from << " PrSMs, " << stage_spectra.size() - rest.size()
                      << " of the spectra identified (a target PrSM at proteoform q-value <= " << rep_q << "); the next stage searches the other "
                      << rest.size() << std::endl;
      stage_spectra.swap(rest);
    }
    if (cut_bonds.load() > 0)
    {
      OPENMS_LOG_INFO << "Database entries: " << cut_bonds.load() << " PrSMs dropped with one half cystine of a disulfide bond of their entry"
                      << std::endl;
    }
    // the chimeric rounds score their PrSMs on the raw spectra of every stage
    if (cascade && raw_loader_ && chimeric_rounds_ > 0 && ! proteoform_hits_.empty())
    {
      std::set<int> scans;
      for (const auto& hit : proteoform_hits_) scans.insert((int)hit.getMetaValue(tntKey("Scan")));
      loadRawSpectra_(scans);
    }

    // Chimeric spectra (chimeric_rounds): a spectrum of two or more co-isolated proteoforms is searched again without the masses its best
    // PrSM explains. In the full spectrum, a second proteoform competes with chance PrSMs that take many of their masses from the first one
    // (native prototype, final_all/chimeric: 55% of the masses of chance PrSMs were the first proteoform's; 6% in the residual), so its
    // per-spectrum null is lower in the residual. Every spectrum with a PrSM is searched again, whether its best PrSM is a target or a
    // decoy, so that the residual PrSMs of targets and decoys come from the same procedure; they form an FDR class per round. Per round:
    // - removed: the masses on the path of the spectrum's best PrSM of the last round (the extension's MatchedMassIndices; with
    //   one PrSM per spectrum its only PrSM, else (native) the best by score, then a determined mass and termini, then coverage), each also up to
    //   chimeric_isotope_offsets isotope spacings away, and its proteoform mass, all within the spectrum's tolerance. They are removed from
    //   all masses of the spectrum, and the cap max_node_cntr_ is applied again to the rest, as to the input.
    // - searched: tags of the residual (a tag found before in the scan keeps its index; new tag sequences are looked up in the database),
    //   the candidates of searchCandidates without precursor ladders and without the proteins of the spectrum's earlier best PrSMs and
    //   their target-decoy partners, and the extension without a precursor mass (the MS1 mass is the first proteoform's; a mass can still
    //   come from complementary fragment pairs). Every candidate's PrSM is kept, scored on the raw spectrum with -in_raw, and the best
    //   one per spectrum is reported (by the combined score with -in_raw, as in the first search; native: every candidate's).
    // - reported per residual PrSM: ChimericRound, RemovedMasses (removed from the spectrum before this round), ResidualNullScore (the
    //   best extension score of the other decoy PrSMs of the same residual search; -1: none) and ResidualScoreRatio (extension score /
    //   ResidualNullScore; -1: no null). They are not a filter: the FDR is the target-decoy competition of the round's class, and the
    //   native prototype's rule (ratio >= 2.5, >= 8 matched masses) was calibrated with a stronger null (also three shuffled-only
    //   databases), so it is applied by the user, e.g. to a single native spectrum, where no class has enough PrSMs for a q-value.
    if (chimeric_rounds_ > 0 && ! proteoform_hits_.empty())
    {
      const Param& eparam_now = extender_param_;
      const UInt mi_round = registry.registerName("ChimericRound"), mi_matched_masses = registry.registerName("MatchedMassIndices"),
                 mi_removed = registry.registerName("RemovedMasses"), mi_null = registry.registerName("ResidualNullScore"),
                 mi_ratio = registry.registerName("ResidualScoreRatio"), mi_mass = registry.registerName("Mass");
      for (auto& hit : proteoform_hits_) hit.setMetaValue(mi_round, 0);
      std::map<int, int> spectrum_of_scan;
      for (int d = 0; d < (int)dspecs.size(); d++) spectrum_of_scan[dspecs[d].getScanNumber()] = d;
      // the entries of a protein (by an entry of it) and of its target-decoy partner (the database's; by name for a FASTA file)
      auto addProteinEntries = [&](Size fasta_index, std::set<Size>& to) {
        const auto& protein = db->proteins[protein_of[fasta_index]];
        to.insert(protein.entries.begin(), protein.entries.end());
        if (protein.partner >= 0) to.insert(db->proteins[(Size)protein.partner].entries.begin(), db->proteins[(Size)protein.partner].entries.end());
      };
      const double isotope_spacing = Constants::ISOTOPE_MASSDIFF_55K_U; // FLASHDeconv's isotope spacing
      // per spectrum: all masses (ascending) and the searched masses (after the cap) of its last residual, the proteins excluded and the
      // number of masses removed so far
      std::vector<std::vector<PeakGroup>> last_all(dspecs.size());
      std::vector<DeconvolvedSpectrum> last_searched(dspecs.size());
      std::vector<std::set<Size>> excluded(dspecs.size());
      std::vector<Size> removed_count(dspecs.size(), 0);
      // the sequence tags found so far in each scan; a residual tag equal to one of them keeps its index
      std::map<int, std::vector<Size>> tags_of_scan = scan_to_tag_indices;
      const std::vector<ProteinHit> no_candidates;
      Size last_first = 0, last_end = proteoform_hits_.size(); // the PrSMs of the last round in proteoform_hits_
      for (int round = 1; round <= chimeric_rounds_; round++)
      {
        // the best PrSM of each spectrum in the last round, targets and decoys alike
        std::vector<long> best(dspecs.size(), -1);
        for (Size i = last_first; i < last_end; i++)
        {
          const auto it = spectrum_of_scan.find((int)proteoform_hits_[i].getMetaValue(mi_scan));
          if (it == spectrum_of_scan.end()) continue;
          long& b = best[it->second];
          if (b < 0 || isBetterPrSM(proteoform_hits_[i], proteoform_hits_[(Size)b])) b = (long)i;
        }
        // the residual spectra
        std::vector<int> searched; // ascending
        std::vector<DeconvolvedSpectrum> residual(dspecs.size());
        std::vector<std::vector<PeakGroup>> residual_all(dspecs.size());
        Size removed_total = 0;
        for (int d = 0; d < (int)dspecs.size(); d++)
        {
          if (best[d] < 0) continue;
          const auto& best_hit = proteoform_hits_[(Size)best[d]];
          const DeconvolvedSpectrum& last = round == 1 ? dspecs[d] : last_searched[d];
          std::vector<PeakGroup> all; // the last round's spectrum with the masses its cap left out
          if (round > 1) all = std::move(last_all[d]);
          else if (! uncapped_masses[d].empty()) all = std::move(uncapped_masses[d]);
          else all.assign(last.begin(), last.end());
          std::vector<double> centers; // the masses on the path of the best PrSM, and its proteoform mass
          if (best_hit.metaValueExists(mi_matched_masses))
          {
            for (int i : best_hit.getMetaValue(mi_matched_masses).toIntList())
            {
              if (i >= 0 && i < (int)last.size()) centers.push_back(last[i].getMonoMass());
            }
          }
          const double proteoform_mass = best_hit.getMetaValue(mi_mass);
          if (proteoform_mass > 0) centers.push_back(proteoform_mass);
          std::vector<double> removed_masses;
          removed_masses.reserve(centers.size() * (2 * chimeric_isotope_offsets + 1));
          for (double c : centers)
          {
            for (int k = -chimeric_isotope_offsets; k <= chimeric_isotope_offsets; k++) removed_masses.push_back(c + k * isotope_spacing);
          }
          std::sort(removed_masses.begin(), removed_masses.end());
          std::vector<PeakGroup> kept;
          kept.reserve(all.size());
          for (const auto& pg : all)
          {
            const double m = pg.getMonoMass(), margin = spectrum_tols[d] * 1e-6 * m;
            const auto it = std::lower_bound(removed_masses.begin(), removed_masses.end(), m - margin);
            if (it != removed_masses.end() && *it <= m + margin) continue;
            kept.push_back(pg);
          }
          removed_count[d] += all.size() - kept.size();
          removed_total += all.size() - kept.size();
          // the protein of the best PrSM and its target-decoy partner are not proposed again in this spectrum; in native mode, where every
          // kept candidate's PrSM is reported, none of the proteins reported in the spectrum (one PrSM per protein)
          std::vector<Size> reported {(Size)(int)best_hit.getMetaValue(mi_fasta_index)};
          if (native_)
          {
            for (Size i = last_first; i < last_end; i++)
            {
              const auto& h = proteoform_hits_[i];
              const auto it = spectrum_of_scan.find((int)h.getMetaValue(mi_scan));
              if (it != spectrum_of_scan.end() && it->second == d && h.metaValueExists(mi_fasta_index)) reported.push_back((Size)(int)h.getMetaValue(mi_fasta_index));
            }
          }
          for (Size fasta_index : reported) addProteinEntries(fasta_index, excluded[d]);
          // the residual: no precursor mass, and of the original spectrum only what the search reads (retention time, MS level)
          DeconvolvedSpectrum res(dspecs[d].getScanNumber());
          MSSpectrum original;
          original.setRT(dspecs[d].getOriginalSpectrum().getRT());
          original.setMSLevel(dspecs[d].getOriginalSpectrum().getMSLevel());
          res.setOriginalSpectrum(original);
          res.setPeakGroups(kept);
          if (res.size() > max_node_cntr_) // the cap of the input
          {
            res.sortByQscore2D();
            std::vector<PeakGroup> filtered_peaks;
            filtered_peaks.reserve(max_node_cntr_);
            for (const auto& pg : res)
            {
              if (filtered_peaks.size() >= max_node_cntr_) break;
              filtered_peaks.push_back(pg);
            }
            res.setPeakGroups(filtered_peaks);
          }
          res.sort();
          residual[d] = std::move(res);
          residual_all[d] = std::move(kept);
          searched.push_back(d);
        }

        // tags of the residual spectra, in spectrum order; new tag sequences are looked up in the database
        std::vector<std::vector<Size>> residual_tags(dspecs.size());
        std::vector<std::string> new_seqs;
        Size tag_count = 0, new_tag_count = 0;
        // the tagging of the residual spectra in parallel (each spectrum its own tagger and slot), the tags then taken in spectrum order
        std::vector<std::vector<FLASHHelperClasses::Tag>> found_of(dspecs.size());
#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < (int)searched.size(); k++)
        {
          const int d = searched[k];
          FLASHTaggerAlgorithm tagger;
          tagger.setParameters(tagger_param_);
          tagger.run(residual[d], spectrum_tols[d]);
          tagger.fillTags(found_of[d]);
        }
        for (int d : searched)
        {
          auto& found = found_of[d];
          tag_count += found.size();
          auto& known = tags_of_scan[residual[d].getScanNumber()];
          for (auto& tag : found)
          {
            Size index = tags_.size();
            for (Size ti : known)
            {
              if (tags_[ti] == tag && tags_[ti].getMzs() == tag.getMzs())
              {
                index = ti;
                break;
              }
            }
            if (index == tags_.size())
            {
              const std::string seq = tag.getUppercaseSequence();
              if (tag_id_of.emplace(seq, tag_db.size()).second)
              {
                tag_db.emplace_back();
                new_seqs.push_back(seq);
              }
              tag.setIndex((int)index);
              tags_.push_back(std::move(tag));
              known.push_back(index);
              new_tag_count++;
            }
            residual_tags[d].push_back(index);
          }
          std::vector<FLASHHelperClasses::Tag>().swap(found);
        }
        if (! new_seqs.empty()) searchTagSequences(new_seqs);

        // candidates of the residual spectra (spectra in parallel, each writing its own slots), appended in spectrum order
        std::vector<std::vector<ProteinHit>> r_hits(dspecs.size()), r_ladder_hits(dspecs.size()), r_pl_hits(dspecs.size()), r_decoy_pl_hits(dspecs.size());
        std::vector<std::vector<FLASHHelperClasses::Tag>> r_gapped_tags(dspecs.size());
        std::vector<Size> r_windows(dspecs.size(), 0);
        std::vector<std::vector<MassLadderAlgorithm::Candidate>> r_sliding;
        if (native_rule) searchSlidingSpectra(residual, searched, &excluded, r_sliding);
#ifdef _OPENMP
        const bool residual_over_spectra = (int)searched.size() >= omp_get_max_threads();
#else
        const bool residual_over_spectra = false;
#endif
#pragma omp parallel for schedule(dynamic) if (residual_over_spectra)
        for (int k = 0; k < (int)searched.size(); k++)
        {
          const int d = searched[k];
          std::string note; // the native note of the first search is the one logged
          searchCandidates(residual[d], &residual_tags[d], excluded[d], false, 0, r_hits[d], r_ladder_hits[d], r_pl_hits[d], r_decoy_pl_hits[d],
                           r_gapped_tags[d], r_windows[d], note, native_rule ? &r_sliding[d] : nullptr);
        }
        std::vector<ProteinHit> r_candidates, r_ladder_candidates;
        std::vector<std::pair<Size, Size>> r_ranges(dspecs.size(), {0, 0}), r_ladder_ranges(dspecs.size(), {0, 0});
        int with_candidates = 0;
        for (int d : searched)
        {
          const int tag_base = appendSpectrumTags(r_gapped_tags[d]);
          appendTagCandidates(r_hits[d], r_candidates, r_ranges[d], tag_base);
          appendLadderCandidates(r_ladder_hits[d], r_ladder_candidates, r_ladder_ranges[d], tag_base);
          with_candidates += r_ranges[d].first != r_ranges[d].second || r_ladder_ranges[d].first != r_ladder_ranges[d].second;
        }

        // the extension of the residual spectra, every candidate's PrSM kept
        std::vector<std::vector<ProteinHit>> r_prsms(dspecs.size());
#ifdef _OPENMP
        const bool extend_residual_over_spectra = with_candidates >= omp_get_max_threads();
#else
        const bool extend_residual_over_spectra = false;
#endif
#pragma omp parallel for schedule(dynamic) if (extend_residual_over_spectra)
        for (int k = 0; k < (int)searched.size(); k++)
        {
          const int d = searched[k];
          extendCandidatesInto(eparam_now, residual[d], r_candidates, r_ranges[d], r_ladder_candidates, r_ladder_ranges[d], no_candidates,
                               std::pair<Size, Size>(0, 0), r_prsms[d], true);
        }

        // the null of each PrSM: the best extension score of the other decoy PrSMs of the same residual search
        std::vector<ProteinHit> round_hits;
        for (int d : searched)
        {
          auto& prsms = r_prsms[d];
          double best_decoy = -1, second_decoy = -1;
          long best_decoy_at = -1;
          for (Size i = 0; i < prsms.size(); i++)
          {
            if (! StringUtils::hasPrefix(prsms[i].getAccession(), "DECOY")) continue;
            const double score = prsms[i].getScore();
            if (score > best_decoy)
            {
              second_decoy = best_decoy;
              best_decoy = score;
              best_decoy_at = (long)i;
            }
            else if (score > second_decoy) { second_decoy = score; }
          }
          for (Size i = 0; i < prsms.size(); i++)
          {
            auto& h = prsms[i];
            h.setMetaValue(mi_round, round);
            h.setMetaValue(mi_removed, (int)removed_count[d]);
            h.setMetaValue(mi_null, (long)i == best_decoy_at ? second_decoy : best_decoy);
            round_hits.push_back(std::move(h));
          }
        }
        const Size all_prsms = round_hits.size();
        if (raw_loader_ && ! round_hits.empty()) annotateOnRawSpectra_(round_hits);
        if (! multiple_hits_per_spec_) keepBestPerSpectrum_(round_hits, " (chimeric round " + StringUtils::toStr(round) + ")", raw_loader_ ? "combined score" : "score");
        Size targets = 0, decoys = 0;
        for (auto& h : round_hits)
        {
          const double extension_score = h.metaValueExists(tntKey("ExtensionScore")) ? (double)h.getMetaValue(tntKey("ExtensionScore")) : h.getScore();
          const double null_score = h.getMetaValue(mi_null);
          h.setMetaValue(mi_ratio, null_score > 0 ? extension_score / null_score : -1.0);
          (StringUtils::hasPrefix(h.getAccession(), "DECOY") ? decoys : targets)++;
        }
        OPENMS_LOG_INFO << "Chimeric round " << round << ": " << searched.size() << " spectra searched without the masses of their best PrSM ("
                        << removed_total << " masses removed), " << tag_count << " tags (" << new_tag_count << " new, " << new_seqs.size()
                        << " new sequences), " << r_candidates.size() << " tag and " << r_ladder_candidates.size() << " ladder candidates, "
                        << all_prsms << " PrSMs of all candidates; reported: " << targets << " target and " << decoys << " decoy PrSMs" << std::endl;

        last_first = proteoform_hits_.size();
        for (auto& h : round_hits) proteoform_hits_.push_back(std::move(h));
        last_end = proteoform_hits_.size();
        for (int d : searched)
        {
          last_searched[d] = std::move(residual[d]);
          last_all[d] = std::move(residual_all[d]);
        }
      }
    }
    // Native mode: a PrSM beside a better PrSM of its spectrum needs evidence of its own. Native mode reports every candidate above its
    // decoy margin, and a homolog of the spectrum's protein passes on that protein's fragments, which shuffled decoys do not model (the
    // same leak as several PrSMs per spectrum in a lysate, above). The PrSMs of a spectrum are taken in the order native mode keeps them
    // (round, then isBetterPrSM); the first is kept, and each later one only if its matched masses that no kept PrSM of the spectrum
    // matched (UniqueMasses; the spectrum's tolerance, +-1 isotope spacing) are unlikely by chance: k ~ Poisson(lambda), lambda = 2 tol x
    // its fragment positions (stretch length - 1, per ion type) / the residue mass of its stretch x the summed masses of the spectrum below
    // that mass that no kept PrSM explains (as the terminal ladders, MassLadderAlgorithm::searchTerminalLadders), UniqueScore = -log10 P,
    // at least native_margin + log10 of the target proteins (the second protein is chosen from all of them). Measured before the rule
    // (25 native test spectra, PXD047341 carbonic anhydrase and bovine hemoglobin): the real second proteins of a spectrum (DLDH, mAb
    // chains, the other hemoglobin chain) had 197 such PrSMs, score median 19, lowest 9.7 apart from a second database entry of the same
    // protein (score 0, same masses); human CA2 beside bovine CA 70 (median 1.7, 1 above the threshold), human hemoglobin chains beside
    // the bovine ones 117 (2 above); no expected protein of a file was lost.
    if (native_rule)
    {
      static const std::map<char, double> residue_mass {{'G', 57.02146}, {'A', 71.03711}, {'S', 87.03203}, {'P', 97.05276}, {'V', 99.06841},
        {'T', 101.04768}, {'C', 103.00919}, {'L', 113.08406}, {'I', 113.08406}, {'N', 114.04293}, {'D', 115.02694}, {'Q', 128.05858},
        {'K', 128.09496}, {'E', 129.04259}, {'M', 131.04049}, {'H', 137.05891}, {'F', 147.06841}, {'R', 156.10111}, {'Y', 163.06333},
        {'W', 186.07931}, {'U', 150.95364}, {'O', 237.14773}};
      const double ion_series = (double)param_.getValue("ion_type").toStringVector().size();
      Size target_proteins = 0;
      if (db != nullptr) { for (const auto& p : db->proteins) target_proteins += p.decoy ? 0 : 1; }
      else { for (const auto& e : fasta_entry) target_proteins += StringUtils::hasPrefix(e.identifier, "DECOY") ? 0 : 1; }
      const double unique_min = native_margin + std::log10((double)std::max(target_proteins, (Size)1));
      const double isotope_spacing = Constants::ISOTOPE_MASSDIFF_55K_U;
      const UInt mi_values = MetaInfoInterface::metaRegistry().registerName("MatchedMassValues"),
                 mi_unique = MetaInfoInterface::metaRegistry().registerName("UniqueMasses"),
                 mi_unique_score = MetaInfoInterface::metaRegistry().registerName("UniqueScore");
      std::map<int, int> spectrum_of;
      for (int d = 0; d < (int)dspecs.size(); d++) spectrum_of[dspecs[d].getScanNumber()] = d;
      std::map<int, std::vector<Size>> prsms_of;
      for (Size i = 0; i < proteoform_hits_.size(); i++)
      {
        const auto it = spectrum_of.find((int)proteoform_hits_[i].getMetaValue(tntKey("Scan")));
        if (it != spectrum_of.end()) prsms_of[it->second].push_back(i);
      }
      std::vector<bool> drop(proteoform_hits_.size(), false);
      Size tested = 0, dropped = 0;
      for (auto& [d, list] : prsms_of)
      {
        std::sort(list.begin(), list.end(), [&](Size a, Size b) {
          const int ra = proteoform_hits_[a].getMetaValue("ChimericRound", 0), rb = proteoform_hits_[b].getMetaValue("ChimericRound", 0);
          if (ra != rb) return ra < rb;
          return isBetterPrSM(proteoform_hits_[a], proteoform_hits_[b]);
        });
        const double tol = spectrum_tols[d];
        std::vector<double> explained; // the masses of the kept PrSMs, +-1 isotope spacing, ascending
        auto isExplained = [&](double m) {
          const double margin = tol * 1e-6 * m;
          const auto it = std::lower_bound(explained.begin(), explained.end(), m - margin);
          return it != explained.end() && *it <= m + margin;
        };
        for (Size r = 0; r < list.size(); r++)
        {
          auto& hit = proteoform_hits_[list[r]];
          std::vector<double> values;
          if (hit.metaValueExists(mi_values)) values = hit.getMetaValue(mi_values).toDoubleList();
          if (r > 0)
          {
            Size k = 0;
            for (double m : values) k += isExplained(m) ? 0 : 1;
            const std::string& seq = hit.getSequence();
            int s = hit.getMetaValue("StartPosition", -1), e = hit.getMetaValue("EndPosition", -1);
            if (s < 1) s = 1;
            if (e < 1 || e > (int)seq.size()) e = (int)seq.size();
            double spread = 0;
            for (int p = s - 1; p < e; p++)
            {
              const auto rm = residue_mass.find(seq[p]);
              spread += rm == residue_mass.end() ? 111.1 : rm->second;
            }
            double summed = 0;
            for (const auto& pg : dspecs[d])
            {
              const double m = pg.getMonoMass();
              if (m < spread && ! isExplained(m)) summed += m;
            }
            const double lambda = spread > 0 ? 2.0 * tol * 1e-6 * (double)std::max(e - s, 0) * ion_series / spread * summed : 0.0;
            const double score = k > 0 ? -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300)) : 0.0;
            hit.setMetaValue(mi_unique, (int)k);
            hit.setMetaValue(mi_unique_score, score);
            tested++;
            if (score < unique_min)
            {
              drop[list[r]] = true;
              dropped++;
              continue;
            }
          }
          for (double m : values)
          {
            for (int iso = -1; iso <= 1; iso++) explained.push_back(m + iso * isotope_spacing);
          }
          std::sort(explained.begin(), explained.end());
        }
      }
      Size w = 0;
      for (Size i = 0; i < proteoform_hits_.size(); i++)
      {
        if (drop[i]) continue;
        if (w != i) proteoform_hits_[w] = std::move(proteoform_hits_[i]);
        w++;
      }
      proteoform_hits_.resize(w);
      OPENMS_LOG_INFO << "Native: " << tested << " PrSMs beside a better PrSM of their spectrum, " << dropped
                      << " removed without evidence of their own (unique-mass score below " << StringUtils::number(unique_min, 2) << ")" << std::endl;
    }
    // the raw spectra stay loaded for the combined coverage after the FDR (setCombinedCoverage_ releases them)
  }

  // the order below and every ranking of the FDR need finite scores: a NaN breaks the sorts and never closes a tie group, and an infinite
  // decoy score leaves the null of the protein level undefined
  for (const auto& hit : proteoform_hits_)
  {
    if (std::isfinite(hit.getScore())) continue;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "PrSM score must be finite (scan " + hit.getMetaValue(tntKey("Scan")).toString() + ")", std::to_string(hit.getScore()));
  }
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), byScore);

  // One MS1 feature, one proteoform: a co-isolated precursor is one real mass of one proteoform. When PrSMs of different proteins report
  // the same MS1-based mass in the same elution (within 2 min), FLASHExtender fitted at least one of them to it by truncation or blind
  // shifts (Eclip_1480: 55 PrSMs at 1% reported a mass of another protein's 1% proteoform nearby and never their own, against 8 with the
  // retention times moved by 20 min; e.g. TATA 41-89 +16.02 and MALE 328-373 -90.1 +473.1 both at 5549.62 in scans 2325 and 2326). Per
  // group of such PrSMs, each protein is ranked by its best score, then by its fewest blind modifications, then by its number of PrSMs at
  // the mass; the first keeps the mass, the others' PrSMs lose it (Mass -1: undetermined, their own FDR class). The score already prices
  // a blind shift (raw score minus log10 of the windows it opens); ranked by blind count first, low-scoring PrSMs that fit the mass by
  // truncation alone (mostly decoys and PrSMs far below 1%) took the masses of real one-shift proteoforms in the same feature (Eclip_1480:
  // 544 PrSMs lost the mass, 625 with the B. theta entrapment database; 323 of the 342 E. coli proteoforms that only the entrapment search
  // left undetermined had lost the mass to a PrSM outside the 1% list).
  if (one_proteoform_per_feature_ && precursor_tol > 0)
  {
    const double iso = Constants::ISOTOPE_MASSDIFF_55K_U, rt_window = 120.0;
    std::vector<Size> idx;
    for (Size i = 0; i < proteoform_hits_.size(); i++)
    {
      const auto& h = proteoform_hits_[i];
      if ((double)h.getMetaValue(tntKey("Mass")) <= 0 || (int)h.getMetaValue(tntKey("ProteoformMassByFragmentMass"), 0) != 0 || h.metaValueExists(tntKey("DecoyPrecursor"))) continue;
      idx.push_back(i);
    }
    std::sort(idx.begin(), idx.end(), [&](Size l, Size r) { return (double)proteoform_hits_[l].getMetaValue(tntKey("Mass")) < (double)proteoform_hits_[r].getMetaValue(tntKey("Mass")); });
    // conflicts: pairs of PrSMs of different proteins at the same mass (or one isotope apart) and elution; union-find over them
    std::vector<Size> parent(idx.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<Size(Size)> root = [&](Size a) { return parent[a] == a ? a : parent[a] = root(parent[a]); };
    std::vector<bool> in_conflict(idx.size(), false);
    for (Size a = 0; a < idx.size(); a++)
    {
      const auto& ha = proteoform_hits_[idx[a]];
      const double ma = ha.getMetaValue(tntKey("Mass")), ra = ha.getMetaValue(tntKey("RT")), max_diff = ma * precursor_tol * 2e-6;
      for (Size b = a + 1; b < idx.size(); b++)
      {
        const auto& hb = proteoform_hits_[idx[b]];
        const double mb = hb.getMetaValue(tntKey("Mass"));
        if (mb - ma > iso + max_diff) break;
        if (std::abs((double)hb.getMetaValue(tntKey("RT")) - ra) > rt_window) continue;
        const double diff = mb - ma;
        if (std::abs(diff) > max_diff && std::abs(diff - iso) > max_diff) continue;
        // the same protein is the same proteoform here; a target and its decoy are different proteins
        if (ha.getAccession() == hb.getAccession()) { parent[root(a)] = root(b); continue; }
        in_conflict[a] = in_conflict[b] = true;
        parent[root(a)] = root(b);
      }
    }
    std::map<Size, std::vector<Size>> groups;
    for (Size a = 0; a < idx.size(); a++) groups[root(a)].push_back(a);
    Size demoted = 0, groups_with_conflict = 0;
    for (const auto& [g, members] : groups)
    {
      if (std::none_of(members.begin(), members.end(), [&](Size a) { return (bool)in_conflict[a]; })) continue;
      groups_with_conflict++;
      // per protein: best score, fewest blind modifications, PrSM count
      std::map<std::string, std::tuple<Size, double, Size>> by_protein;
      for (Size a : members)
      {
        const auto& h = proteoform_hits_[idx[a]];
        auto it = by_protein.find(h.getAccession());
        const Size blind = blindModCount_(h);
        if (it == by_protein.end()) { by_protein[h.getAccession()] = {blind, h.getScore(), 1}; continue; }
        auto& [b0, s0, n0] = it->second;
        b0 = std::min(b0, blind);
        s0 = std::max(s0, h.getScore());
        n0++;
      }
      std::string owner;
      std::tuple<Size, double, Size> best {};
      bool first = true;
      for (const auto& [acc, key] : by_protein)
      {
        const auto& [b, s, n] = key;
        const auto& [bb, bs, bn] = best;
        if (first || s > bs || (s == bs && (b < bb || (b == bb && n > bn))))
        {
          owner = acc;
          best = key;
          first = false;
        }
      }
      for (Size a : members)
      {
        auto& h = proteoform_hits_[idx[a]];
        if (h.getAccession() == owner) continue;
        h.setMetaValue(tntKey("Mass"), -1.0);
        h.setMetaValue(tntKey("MassLostToFeatureOwner"), owner);
        demoted++;
      }
    }
    OPENMS_LOG_INFO << "One proteoform per MS1 feature: " << groups_with_conflict << " masses claimed by PrSMs of several proteins; " << demoted
                    << " PrSMs of the other proteins lost the mass" << std::endl;
  }

  // FDR by target-decoy counting (targetDecoyQvalues, no pseudo-count) at three levels:
  // - PrSMs, per class (fdr_class_names): determined mass and termini with 0, 1, or 2 or more blind modifications (variable modifications
  //   do not count), and undetermined with both termini known (the mass is then open only inside) or with an open terminus. Their decoy
  //   rates differ widely (Eclip_1480 with decoys kept: 2 decoys for 3770 determined targets, 97 for 1070 undetermined with both
  //   termini, 1030 for 2228 with an open terminus). The PrSMs of precursor-ladder candidates (meta value PrecursorLadder; only
  //   determined ones that no other source found in the spectrum are kept) are a class of their own, assigned by that source before the
  //   rules above: a source that searches every window fitting the precursor mass has its own chance rate, which the nearly decoy-free
  //   determined classes would hide. The class does not split by blind modifications: a window fits the unmodified mass, so such a PrSM
  //   with blind shifts left its window (the log line of the class counts them, targets and decoys, to decide whether it should be split).
  //   At the protein level the class uses the pooled null while it has fewer than null_min_class_decoys decoys, and its decoys join the
  //   pooled null (so the null of the sparse determined classes changes with the option; not yet measured). The PrSMs of each residual
  //   round of the chimeric search (meta value ChimericRound r > 0) are class "residual r", also by source and before the rules above;
  //   determined and undetermined ones share it (the log line of the class counts both), and its decoys join the pooled null like those
  //   of the precursor-ladder class. Each PrSM is in exactly one class, so no PrSM passes twice. With discard_underdetermined the undetermined PrSMs are in no class and take part in nothing below,
  //   the proteoforms included.
  // - proteoforms: per proteoform class (a partition: determined, undetermined with both termini, with an open terminus, residual r),
  //   ranked by the best score of the PrSMs that decide the class (proteoformQvalues; -native: one ranking by PrSM q-value);
  // - proteins: ranked by their best proteoform (the lowest proteoform q-value, then the best score of its PrSMs) and target-decoy counted,
  //   so that a protein is identified by a proteoform; -native: setProteinQvalues_ (local-FDR evidence of all PrSMs).
  // Without decoys in the database and without discard_underdetermined no q-value is set.
  //   With search_cascade, the PrSM q-values are computed per class and stage (meta value SearchStage), so that the gapped-tag PrSMs of
  //   stage 2, from the spectra stage 1 did not identify, compete only with the decoys of stage 2.
  std::vector<int> hit_class; // index into fdr_class_names; -1: in no class
  std::vector<double> qvalues, decoy_weights;
  std::vector<bool> decoys;
  // the proteoforms: positions of their PrSMs in proteoform_hits_ (Index), the representative first; representatives in the order above
  const auto clusters = proteoformQvalues(hit_class, qvalues, decoys, decoy_weights, true);
  const Size hit_count = proteoform_hits_.size();
  setCombinedCoverage_(clusters, qvalues, rep_q, use_fdr);

  if (use_fdr && native_rule) { setProteinQvalues_(hit_class, qvalues, decoys, decoy_weights, clusters, rep_q, *db); }
  else if (use_fdr && protein_model_fdr_ && scoring_model_ && raw_loader_)
  {
    // Protein q-values from the learned confidence (protein_ranking model): each target protein is ranked by the confidence c of its best PrSM
    // (one PrSM per spectrum, as below), a log-odds whose P(the PrSM is wrong) = 1 / (1 + e^c); the FDR at a cut is the mean of that chance
    // over the proteins above it, q its running minimum from the bottom (ties take the value at the end of their group). No decoy is
    // counted: per run the decoy proteins never reached the top of a ranking by proteoform q-value, so the proteins above a run's first
    // decoy passed at q = 0 however many were false. Held-out benchmark (seed 2, 31 runs, proteins per run, upper bound 2E/N of shuffled
    // entrapment): 2.69 % with the proteoform ranking, 1.34 % with this one (TopPIC 1.38 %); foreign entrapment 1.53 % (TopPIC 2.31 %).
    // Decoy proteins keep q-value 1. It relies on the confidence being calibrated, which the learned model was trained for.
    std::map<int, Size> best_of_scan;
    for (Size i = 0; i < hit_count; i++)
    {
      if (hit_class[i] < 0) continue;
      const auto& hit = proteoform_hits_[i];
      if (hit.metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(hit.getAccession(), "DECOY")) continue;
      auto [it, fresh] = best_of_scan.try_emplace((int)hit.getMetaValue(tntKey("Scan")), i);
      if (! fresh && hit.getScore() > proteoform_hits_[it->second].getScore()) it->second = i;
    }
    std::map<std::string, double> best_confidence; // target proteins
    for (const auto& [scan, i] : best_of_scan)
    {
      const auto& hit = proteoform_hits_[i];
      if (StringUtils::hasPrefix(hit.getAccession(), "DECOY")) continue;
      auto [it, fresh] = best_confidence.try_emplace(hit.getAccession(), hit.getScore());
      if (! fresh) it->second = std::max(it->second, hit.getScore());
    }
    std::vector<std::pair<double, std::string>> ranked_proteins; // -confidence, accession
    for (const auto& [acc, c] : best_confidence) ranked_proteins.emplace_back(-c, acc);
    std::sort(ranked_proteins.begin(), ranked_proteins.end());
    const Size n = ranked_proteins.size();
    std::vector<double> fdr(n, 1.0);
    double sum_wrong = 0;
    for (Size k = 0; k < n;)
    {
      Size j = k;
      for (; j < n && ranked_proteins[j].first == ranked_proteins[k].first; j++) sum_wrong += 1.0 / (1.0 + std::exp(-ranked_proteins[j].first));
      std::fill(fdr.begin() + k, fdr.begin() + j, std::min(1.0, sum_wrong / (double)j));
      k = j;
    }
    for (Size k = n; k-- > 1;) fdr[k - 1] = std::min(fdr[k - 1], fdr[k]);
    std::map<std::string, double> protein_q;
    Size passing = 0;
    for (Size k = 0; k < n; k++)
    {
      protein_q[ranked_proteins[k].second] = fdr[k];
      passing += fdr[k] <= rep_q;
    }
    OPENMS_LOG_INFO << "Proteins by the learned confidence of their best PrSM: " << n << " target proteins, " << passing << " at q <= " << rep_q
                    << std::endl;
    for (auto& hit : proteoform_hits_)
    {
      const auto it = protein_q.find(hit.getAccession());
      hit.setMetaValue("ProteinQvalue", it == protein_q.end() ? 1.0 : it->second);
    }
  }
  else if (use_fdr && protein_best_score_)
  {
    // Protein q-values by the rule of the PrSMs and proteoforms (protein_ranking best_score): a protein's score is the score of its best
    // proteoform, and the proteins, targets and decoys alike, are target-decoy counted in one ranking. A decoy precursor of a target and a
    // proteoform of the open-terminus class (not counted, proteoformQvalues) rank no protein.
    std::map<std::string, double> best_score;
    for (const auto& members : clusters)
    {
      const auto& rep = proteoform_hits_[members[0]];
      if (rep.metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(rep.getAccession(), "DECOY")) continue;
      if (rep.metaValueExists("ProteoformEvidenceClass") && (int)rep.getMetaValue("ProteoformEvidenceClass") == pf_class_open_terminus) continue;
      double score = std::numeric_limits<double>::lowest();
      for (Size i : members) score = std::max(score, proteoform_hits_[i].getScore());
      auto [it, fresh] = best_score.try_emplace(rep.getAccession(), score);
      if (! fresh) it->second = std::max(it->second, score);
    }
    std::vector<std::pair<double, std::string>> ranked_proteins; // -score, accession
    for (const auto& [acc, score] : best_score) ranked_proteins.emplace_back(-score, acc);
    std::sort(ranked_proteins.begin(), ranked_proteins.end());
    const double protein_decoy_weight = 1.0 / std::max(decoy_factor_, 1e-9);
    std::vector<RankedItem<double>> items;
    items.reserve(ranked_proteins.size());
    for (const auto& [key, acc] : ranked_proteins) items.push_back({key, StringUtils::hasPrefix(acc, "DECOY"), protein_decoy_weight});
    const auto protein_qvalues = targetDecoyQvalues(items);
    std::map<std::string, double> protein_q;
    Size targets = 0, decoy_count = 0, passing = 0;
    for (Size k = 0; k < ranked_proteins.size(); k++)
    {
      protein_q[ranked_proteins[k].second] = protein_qvalues[k];
      if (items[k].decoy) decoy_count++;
      else { targets++; passing += protein_qvalues[k] <= rep_q; }
    }
    OPENMS_LOG_INFO << "Proteins by the score of their best proteoform: " << targets << " target and " << decoy_count << " decoy proteins, " << passing
                    << " targets at q <= " << rep_q << std::endl;
    for (auto& hit : proteoform_hits_)
    {
      const auto it = protein_q.find(hit.getAccession());
      hit.setMetaValue("ProteinQvalue", it == protein_q.end() ? 1.0 : it->second);
    }
  }
  else if (use_fdr && protein_top_two_)
  {
    // Protein q-values from the two best PrSMs of each protein (protein_ranking top_two_prsms): a protein's score is the sum of the scores
    // of its two best PrSMs (a protein seen in one spectrum keeps that spectrum's score), and the proteins, targets and decoys alike, are
    // target-decoy counted in one ranking. A learned score (scoring_model) is a log-odds: softplus(score) = -log P(the PrSM is wrong) is
    // what adds over independent PrSMs. A single strong PrSM, often an undetermined one with an open C-terminus, no longer carries a
    // protein alone: on the development runs (shuffled entrapment, seed 1, first run of 11 data sets) 812 target and 2 entrapment
    // proteins pass at 1%, against 905 and 27 ranked by the best proteoform. A decoy precursor of a target ranks no protein.
    // one PrSM per spectrum: a spectrum can keep two PrSMs of one protein (e.g. one with both termini and one with an open terminus), and
    // its best alone is evidence; counted twice, one spectrum would pass a protein as two
    const double ninf = -std::numeric_limits<double>::infinity();
    std::map<int, Size> best_of_scan;
    for (Size i = 0; i < hit_count; i++)
    {
      if (hit_class[i] < 0) continue;
      const auto& hit = proteoform_hits_[i];
      if (hit.metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(hit.getAccession(), "DECOY")) continue;
      auto [it, fresh] = best_of_scan.try_emplace((int)hit.getMetaValue(tntKey("Scan")), i);
      if (! fresh && hit.getScore() > proteoform_hits_[it->second].getScore()) it->second = i;
    }
    std::map<std::string, std::array<double, 2>> top_two;
    for (const auto& [scan, i] : best_of_scan)
    {
      const auto& hit = proteoform_hits_[i];
      const double s = hit.getScore();
      const double evidence = scoring_model_ ? (s > 30 ? s : std::log1p(std::exp(s))) : s;
      auto& t = top_two.try_emplace(hit.getAccession(), std::array<double, 2> {ninf, ninf}).first->second;
      if (evidence > t[0]) { t[1] = t[0]; t[0] = evidence; }
      else if (evidence > t[1]) t[1] = evidence;
    }
    std::vector<std::pair<double, std::string>> ranked_proteins; // -score, accession
    for (const auto& [acc, t] : top_two) ranked_proteins.emplace_back(-(t[0] + (t[1] > ninf ? t[1] : 0.0)), acc);
    std::sort(ranked_proteins.begin(), ranked_proteins.end());
    const double protein_decoy_weight = 1.0 / std::max(decoy_factor_, 1e-9);
    std::vector<RankedItem<double>> items;
    items.reserve(ranked_proteins.size());
    for (const auto& [key, acc] : ranked_proteins) items.push_back({key, StringUtils::hasPrefix(acc, "DECOY"), protein_decoy_weight});
    const auto protein_qvalues = targetDecoyQvalues(items);
    std::map<std::string, double> protein_q;
    Size targets = 0, decoy_count = 0, passing = 0;
    for (Size k = 0; k < ranked_proteins.size(); k++)
    {
      protein_q[ranked_proteins[k].second] = protein_qvalues[k];
      if (items[k].decoy) decoy_count++;
      else { targets++; passing += protein_qvalues[k] <= rep_q; }
    }
    OPENMS_LOG_INFO << "Proteins by their two best PrSMs: " << targets << " target and " << decoy_count << " decoy proteins, " << passing
                    << " targets at q <= " << rep_q << std::endl;
    for (auto& hit : proteoform_hits_)
    {
      const auto it = protein_q.find(hit.getAccession());
      hit.setMetaValue("ProteinQvalue", it == protein_q.end() ? 1.0 : it->second);
    }
  }
  else if (use_fdr)
  {
    // Protein q-values from the proteoforms: a protein is ranked by its best proteoform (the lowest proteoform q-value, then the best
    // score of the proteoform's PrSMs), targets and decoys alike, and the proteins are target-decoy counted. The local-FDR evidence of
    // all PrSMs (setProteinQvalues_) also counted PrSMs of no 1% proteoform and was not calibrated against the proteoform classes: on
    // Eclip_1480, 16 proteins passed without a 1% proteoform and 18 with one did not. A decoy precursor carries the target accession: a
    // decoy of the mass, not of the protein, it ranks no protein.
    // The proteins are counted in two classes by their best proteoform: one that reaches a terminus of the protein (its first residue,
    // or the second after an initiator Met, or its last; of a PEFF entry, those of a processed form) and one inside it (an internal
    // stretch). An internal stretch of a long protein fits a mass at many windows, so it is weaker evidence of the protein: counted
    // with the others, it passed B. theta entrapment proteins on Eclip_1480 (3 of 6, single stretches of about 2.9 kDa); apart (73 target
    // and 18 decoy proteins), 5 of 325 proteins were B. theta instead of 6 of 321, with the same 322 E. coli proteins of the E. coli search.
    struct ProteinBest
    {
      std::pair<double, double> key {2.0, 0.0}; // lowest proteoform q-value, -best score
      bool internal = false;
    };
    // A proteoform of the open-terminus class (not counted, proteoformQvalues) ranks no protein: with its q-value 1 it would still rank
    // a protein below all others, and pass with the decoy-free margin above it (the proteins that only such a proteoform identified were
    // 5 B. theta of 35 at 1% when the class was counted).
    std::map<std::string, ProteinBest> best_of;
    for (const auto& members : clusters)
    {
      const auto& rep = proteoform_hits_[members[0]];
      if (rep.metaValueExists("DecoyPrecursor") && ! StringUtils::hasPrefix(rep.getAccession(), "DECOY")) continue;
      if (rep.metaValueExists("ProteoformEvidenceClass") && (int)rep.getMetaValue("ProteoformEvidenceClass") == pf_class_open_terminus) continue;
      const double q = rep.getMetaValue("proqvalue", 1.0);
      double score = 0;
      for (Size i : members) score = std::max(score, proteoform_hits_[i].getScore());
      const std::pair<double, double> key {q, -score};
      auto& b = best_of[rep.getAccession()];
      if (! (key < b.key)) continue;
      b.key = key;
      const int start = rep.getMetaValue("StartPosition"), end = rep.getMetaValue("EndPosition");
      bool terminal = (start >= 1 && start <= 2) || (end > 0 && end == (int)rep.getSequence().size());
      if (rep.metaValueExists("EntryStart"))
      {
        terminal = false;
        for (const auto& [form_start, form_end] : db->proteins[db->proteinOfEntry()[(Size)(int)rep.getMetaValue("FastaIndex")]].forms)
        {
          terminal = terminal || (start > 0 && (start == form_start || (form_start == 1 && start == 2))) || (end > 0 && end == form_end);
        }
      }
      b.internal = ! terminal;
    }
    std::map<std::string, double> protein_q;
    const double protein_decoy_weight = 1.0 / std::max(decoy_factor_, 1e-9);
    for (const bool internal : {false, true})
    {
      std::vector<std::pair<std::pair<double, double>, std::string>> ranked_proteins;
      for (const auto& [acc, b] : best_of)
      {
        if (b.internal == internal) ranked_proteins.emplace_back(b.key, acc);
      }
      std::sort(ranked_proteins.begin(), ranked_proteins.end());
      std::vector<RankedItem<std::pair<double, double>>> items;
      items.reserve(ranked_proteins.size());
      for (const auto& [key, acc] : ranked_proteins) items.push_back({key, StringUtils::hasPrefix(acc, "DECOY"), protein_decoy_weight});
      const auto protein_qvalues = targetDecoyQvalues(items);
      Size targets = 0, decoy_count = 0, passing = 0;
      for (Size k = 0; k < ranked_proteins.size(); k++)
      {
        protein_q[ranked_proteins[k].second] = protein_qvalues[k];
        if (items[k].decoy) decoy_count++;
        else { targets++; passing += protein_qvalues[k] <= rep_q; }
      }
      OPENMS_LOG_INFO << "Protein class \"" << (internal ? "internal stretch" : "reaches a terminus") << "\": " << targets << " target and " << decoy_count
                      << " decoy proteins, " << passing << " targets at q <= " << rep_q << std::endl;
    }
    for (auto& hit : proteoform_hits_)
    {
      const auto it = protein_q.find(hit.getAccession());
      hit.setMetaValue("ProteinQvalue", it == protein_q.end() ? 1.0 : it->second);
    }
  }
  if (use_fdr)
  {
    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits_.size());
    for (Size i = 0; i < hit_count; i++)
    {
      if (hit_class[i] < 0) continue; // in no class, without a q-value (undetermined, with keep_underdetermined_ false)
      if (! keep_decoy_ && decoys[i]) continue;
      // no PrSM-level filter: a PrSM is reported when its proteoform passes (every member carries its cluster's proteoform q-value)
      if (pro_fdr_ < 1 && (double)proteoform_hits_[i].getMetaValue(tntKey("proqvalue"), 1.0) > pro_fdr_) continue;
      filtered_proteoform_hits.push_back(std::move(proteoform_hits_[i]));
    }
    proteoform_hits_.swap(filtered_proteoform_hits);
  }
  {
    // by RT: std::sort of the positions with the RT values taken out once (the same comparisons give std::sort the same permutation as
    // sorting the hits, without a registry lookup per comparison), then the hits moved into that order
    const UInt mi_rt = MetaInfoInterface::metaRegistry().registerName("RT");
    std::vector<DataValue> rts(proteoform_hits_.size());
    std::vector<Size> order(proteoform_hits_.size());
    for (Size i = 0; i < proteoform_hits_.size(); i++) rts[i] = proteoform_hits_[i].getMetaValue(mi_rt);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&rts](Size left, Size right) { return rts[left] < rts[right]; });
    std::vector<ProteinHit> sorted;
    sorted.reserve(order.size());
    for (Size i : order) sorted.push_back(std::move(proteoform_hits_[i]));
    proteoform_hits_.swap(sorted);
  }


  // the candidates of a database entry: the length of its protein (their sequence is the entry's processed form)
  if (by_entry)
  {
    for (auto& hit : protein_hits_)
    {
      if (hit.metaValueExists("FastaIndex"))
        hit.setMetaValue("ProteinLength", (int)db->proteins[protein_of[(Size)(int)hit.getMetaValue("FastaIndex")]].sequence.size());
    }
  }

  //  define prsm index and tag to prsm indices.
  int prsm_index = 0;
  std::map<int, int> prsm_index_map;

  for (auto& hit : proteoform_hits_)
  {
    int fasta_index = hit.getMetaValue("FastaIndex");
    hit.setSequence(db->proteins[protein_of[fasta_index]].sequence); // the full-length protein as given (with I)

    prsm_index_map[hit.getMetaValue(tntKey("Index"))] = prsm_index; // from previous to new index
    hit.setMetaValue(tntKey("Index"), prsm_index);
    for (int tag_index : (std::vector<int>)hit.getMetaValue(tntKey("TagIndices")).toIntList())
    {
      matching_hits_indices_[tag_index].push_back(prsm_index);
    }
    prsm_index++;
  }

  // redefine prsm indices for representative proteoforms
  for (auto& hit : proteoform_hits_)
  {
    if ((int)hit.getMetaValue(tntKey("Representative")) == 0) continue;
    auto prsm_indices = hit.getMetaValue(tntKey("PrSMIndices")).toIntList();
    std::vector<int> new_prsm_indices;
    for (int i : prsm_indices)
    {
      // PrSMs of the proteoform that did not pass the PrSM-level FDR are not reported (operator[] made them PrSM 0 of the output)
      const auto it = prsm_index_map.find(i);
      if (it != prsm_index_map.end()) new_prsm_indices.push_back(it->second);
    }
    std::sort(new_prsm_indices.begin(), new_prsm_indices.end());
    hit.setMetaValue(tntKey("PrSMIndices"), new_prsm_indices);
  }
  // the index of the reported proteoform (as FLASHTnTFile::writeProteoforms numbers them) on every PrSM of it; -1: not reported
  {
    const double pro_fdr = param_.getValue("pro_fdr");
    std::map<int, int> prsm_to_proteoform;
    int proteoform_index = 0;
    for (const auto& hit : proteoform_hits_)
    {
      if ((int)hit.getMetaValue(tntKey("Representative")) == 0) continue;
      if (hit.metaValueExists(tntKey("proqvalue")) && (double)hit.getMetaValue(tntKey("proqvalue")) > pro_fdr) continue;
      for (int i : hit.getMetaValue(tntKey("PrSMIndices")).toIntList()) prsm_to_proteoform[i] = proteoform_index;
      prsm_to_proteoform[(int)hit.getMetaValue(tntKey("Index"))] = proteoform_index;
      proteoform_index++;
    }
    for (auto& hit : proteoform_hits_)
    {
      const auto it = prsm_to_proteoform.find((int)hit.getMetaValue(tntKey("Index")));
      hit.setMetaValue(tntKey("ProteoformIndex"), it == prsm_to_proteoform.end() ? -1 : it->second);
    }
  }
  if (native_rule) setNativeTermini_();
  setProteoformLevels_(*db);
  // the large locals freed with all threads instead of by their destructors on return
  releaseInParallel(dspecs);
  releaseInParallel(uncapped_masses);
  releaseInParallel(tag_db);
  releaseInParallel(vec_pro);
  releaseInParallel(rev_vec_pro);
  releaseInParallel(pro_prefix);
  releaseInParallel(pro_suffix);
  releaseInParallel(ladder_candidates);
  releaseInParallel(pl_candidates);
  releaseInParallel(decoy_pl_candidates);
  // TODO
  // per scan, get the scan number and the precursor mass - maybe given or maybe deconvolved.
  // get all protoeforms - mass and score.
  // get all combinations with a certain tolerance.
  // record the best ranking combination(s) along with summed score (can score be non overlapping?) and mass difference
}

// The proteoform level of every hit, by the five-level classification of proteoform identifications (Smith et al., Nat. Methods 16,
// 939-940, 2019) with the labels of its reference code (mzLib ProteoformLevelClassifier): four aspects, each known or ambiguous; level 1
// all known; 2A, 2B, 2C or 2D one ambiguous (A PTM localization, B PTM identity, C amino acid sequence, D gene of origin); 3, 4 and 5
// two, three and four ambiguous. The reference code compares the tied candidates of a spectrum; a PrSM here is one proteoform, so each
// aspect is judged by what the PrSM itself leaves open:
// - PTM localization: known when every modification sits on one residue; a range of residues (ModStart < ModEnd) leaves it open.
// - PTM identity: known when every modification is named: a variable modification with one name (or a sum of named ones, "A+B"), or a
//   fixed one of the database entry (PEFF). A blind shift is a mass, not an identity, even where Unimod entries of its mass are listed
//   (ModID), as the reference code counts a mass shift ("+15.99") against the identity.
// - Amino acid sequence: known when both termini are known and the residues have no ambiguous letter (B, J, X, Z). A blind shift does not
//   open the sequence; as in the reference code, a mass shift counts against the PTM identity only.
// - Gene: known when no other gene's protein of the same kind (target or decoy) contains the proteoform's residues, I and L counted as
//   one (the same mass). At an open terminus the residues are those of the ProteinSequence column (up to the end of the entry's form).
//   The gene of a protein is GN= of its UniProt description or \GName= of its PEFF description, else its accession without an isoform
//   suffix (-n), so that the isoforms of a gene are one gene.
// A PrSM without a proteoform mass (undetermined) has a mass difference that no modification accounts for: it may be a PTM of any identity
// at any position, so both PTM aspects are open; with an open terminus the sequence is open as well.
void FLASHTnTAlgorithm::setProteoformLevels_(const ProteoformDatabase& db)
{
  if (proteoform_hits_.empty()) return;
  const auto& protein_of = db.proteinOfEntry();

  // the gene of every protein, as an index
  std::vector<Size> gene_of(db.proteins.size());
  {
    std::unordered_map<std::string, Size> gene_index;
    auto value_after = [](const std::string& text, const std::string& key) -> std::string {
      for (Size k = text.find(key); k != std::string::npos; k = text.find(key, k + 1))
      {
        if (k > 0 && text[k - 1] != ' ' && text[k - 1] != '\\') continue; // a key starts a word (GN=) or follows its backslash (\GName=)
        const Size from = k + key.size();
        const Size to = text.find_first_of(" \\", from);
        std::string value = text.substr(from, to == std::string::npos ? std::string::npos : to - from);
        if (! value.empty()) return value;
      }
      return "";
    };
    for (Size p = 0; p < db.proteins.size(); p++)
    {
      const auto& protein = db.proteins[p];
      std::string gene = value_after(protein.description, "GN=");
      if (gene.empty()) gene = value_after(protein.description, "\\GName=");
      if (! gene.empty()) gene = "GN " + gene;
      else
      {
        // the accession: the middle field of db|accession|name, without an isoform suffix (-n)
        std::vector<std::string> fields;
        StringUtils::split(protein.identifier, '|', fields);
        gene = fields.size() >= 3 ? fields[1] : protein.identifier;
        const Size dash = gene.rfind('-');
        if (dash != std::string::npos && dash + 1 < gene.size() &&
            std::all_of(gene.begin() + (long)dash + 1, gene.end(), [](char c) { return std::isdigit((unsigned char)c) != 0; }))
        {
          gene.erase(dash);
        }
        gene = "AC " + gene;
      }
      gene_of[p] = gene_index.emplace(gene, gene_index.size()).first->second;
    }
  }

  // the residues of every hit (I as L), the distinct ones as needles
  std::vector<std::string> needles;
  std::unordered_map<std::string, Size> needle_index;
  std::vector<Size> needle_of_hit(proteoform_hits_.size());
  std::vector<bool> ambiguous_letters(proteoform_hits_.size(), false);
  for (Size i = 0; i < proteoform_hits_.size(); i++)
  {
    const auto& hit = proteoform_hits_[i];
    const int start = hit.getMetaValue(tntKey("StartPosition")), end = hit.getMetaValue(tntKey("EndPosition"));
    const std::string& sequence = hit.getSequence();
    const int size = (int)sequence.size();
    // as FLASHTnTFile's ProteinSequence column: an open terminus bounded by the entry's processed form
    const int entry_first = hit.metaValueExists("EntryStart") ? (int)hit.getMetaValue("EntryStart") - 1 : 0;
    const int entry_end = hit.metaValueExists("EntryEnd") ? (int)hit.getMetaValue("EntryEnd") : size;
    const int first = std::max(0, start < 0 ? entry_first : start - 1);
    const int last = std::min(size, end < 0 ? entry_end : end);
    std::string residues = last > first ? sequence.substr((Size)first, (Size)(last - first)) : std::string();
    std::replace(residues.begin(), residues.end(), 'I', 'L');
    ambiguous_letters[i] = residues.find_first_of("BJXZ") != std::string::npos;
    needle_of_hit[i] = needle_index.emplace(residues, needles.size()).first->second;
    if (needle_of_hit[i] == needles.size()) needles.push_back(residues);
  }

  // the proteins that contain each needle: an Aho-Corasick trie of the needles of plain residues over all proteins (I as L), in parallel;
  // the few needles with other letters (U, O, ambiguous ones) by a plain search
  static const std::string plain_letters = "ACDEFGHKLMNPQRSTVWY";
  std::vector<bool> plain(needles.size());
  std::vector<std::string> trie_needles;
  std::vector<Size> trie_to_needle;
  for (Size n = 0; n < needles.size(); n++)
  {
    plain[n] = ! needles[n].empty() && needles[n].find_first_not_of(plain_letters) == std::string::npos;
    if (! plain[n]) continue;
    trie_needles.push_back(needles[n]);
    trie_to_needle.push_back(n);
  }
  std::vector<std::string> haystacks(db.proteins.size());
  for (Size p = 0; p < db.proteins.size(); p++)
  {
    haystacks[p] = db.proteins[p].sequence;
    std::replace(haystacks[p].begin(), haystacks[p].end(), 'I', 'L');
  }
  std::vector<std::vector<Size>> proteins_of_needle(needles.size());
  if (! trie_needles.empty())
  {
    ACTrie ac_trie(0, 0); // no ambiguities, no mismatches
    ac_trie.addNeedlesAndCompress(trie_needles);
    std::vector<std::pair<Size, Size>> found; // needle, protein
#pragma omp parallel
    {
      OpenMS::ACTrieState ac_state;
      std::vector<std::pair<Size, Size>> local;
#pragma omp for schedule(dynamic, 64) nowait
      for (SignedSize p = 0; p < (SignedSize)haystacks.size(); p++)
      {
        ac_state.setQuery(haystacks[p]);
        ac_trie.getAllHits(ac_state);
        for (const auto& h : ac_state.hits) local.emplace_back(trie_to_needle[h.needle_index], (Size)p);
      }
#pragma omp critical(proteoform_levels)
      found.insert(found.end(), local.begin(), local.end());
    }
    for (const auto& [n, p] : found) proteins_of_needle[n].push_back(p);
  }
  for (Size n = 0; n < needles.size(); n++)
  {
    if (plain[n] || needles[n].empty()) continue;
    for (Size p = 0; p < haystacks.size(); p++)
    {
      if (haystacks[p].find(needles[n]) != std::string::npos) proteins_of_needle[n].push_back(p);
    }
  }

  // the four aspects and the level of every hit
  const char* const aspect_names[] = {"PTM localization", "PTM identity", "sequence", "gene"};
  std::map<std::string, Size> level_counts;
  std::array<Size, 4> open_counts {};
  Size counted = 0;
  const double pro_fdr = param_.getValue("pro_fdr");
  for (Size i = 0; i < proteoform_hits_.size(); i++)
  {
    auto& hit = proteoform_hits_[i];
    const bool determined = isDetermined(hit);
    const int start = hit.getMetaValue(tntKey("StartPosition")), end = hit.getMetaValue(tntKey("EndPosition"));
    bool localized = determined, identified = determined;
    const bool sequence_known = start > 0 && end > 0 && ! ambiguous_letters[i];
    if (determined && hit.metaValueExists(tntKey("Modifications")))
    {
      const std::vector<double> masses = hit.getMetaValue(tntKey("Modifications"));
      const std::vector<int> starts = hit.getMetaValue(tntKey("ModificationStarts"));
      const std::vector<int> ends = hit.getMetaValue(tntKey("ModificationEnds"));
      const std::vector<std::string> ids = hit.metaValueExists("ModificationIDs") ? hit.getMetaValue("ModificationIDs").toStringList() : std::vector<std::string>();
      const std::vector<int> variable = hit.metaValueExists("ModificationVariable") ? (std::vector<int>)hit.getMetaValue("ModificationVariable").toIntList() : std::vector<int>();
      for (Size m = 0; m < masses.size(); m++)
      {
        if (m >= starts.size() || m >= ends.size() || starts[m] != ends[m]) localized = false;
        const int kind = m < variable.size() ? variable[m] : 0; // 1 variable, 2 fixed of the entry, 0 blind
        if (kind == 2) continue;
        std::string names = m < ids.size() ? ids[m] : std::string();
        while (StringUtils::hasSuffix(names, ",")) names.pop_back();
        // one name (or a sum of named ones); alternatives "A,B" leave the identity open, as the ProForma column writes them as a mass
        if (kind != 1 || names.empty() || StringUtils::has(names, ',')) identified = false;
      }
    }
    const Size own_protein = protein_of[(Size)(int)hit.getMetaValue(tntKey("FastaIndex"))];
    bool gene_known = true;
    for (Size p : proteins_of_needle[needle_of_hit[i]])
    {
      if (db.proteins[p].decoy == db.proteins[own_protein].decoy && gene_of[p] != gene_of[own_protein])
      {
        gene_known = false;
        break;
      }
    }
    const bool known[] = {localized, identified, sequence_known, gene_known};
    const int known_count = (int)localized + (int)identified + (int)sequence_known + (int)gene_known;
    std::string level;
    if (known_count == 3) level = ! localized ? "2A" : ! identified ? "2B" : ! sequence_known ? "2C" : "2D";
    else level = std::to_string(5 - known_count);
    hit.setMetaValue("ProteoformLevel", level);

    // the log: the reported target proteoforms (as FLASHTnTFile::writeProteoforms)
    if (! hit.metaValueExists(tntKey("Representative")) || (int)hit.getMetaValue(tntKey("Representative")) == 0) continue;
    if (hit.metaValueExists(tntKey("proqvalue")) && (double)hit.getMetaValue(tntKey("proqvalue")) > pro_fdr) continue;
    if (db.proteins[own_protein].decoy) continue;
    counted++;
    level_counts[level]++;
    for (Size a = 0; a < 4; a++) open_counts[a] += ! known[a];
  }
  std::string line;
  for (const auto& [level, count] : level_counts) line += (line.empty() ? "" : ", ") + level + " " + std::to_string(count);
  std::string open_line;
  for (Size a = 0; a < 4; a++) open_line += (a == 0 ? "" : ", ") + std::string(aspect_names[a]) + " " + std::to_string(open_counts[a]);
  OPENMS_LOG_INFO << "Proteoform levels (Smith et al. 2019) of the " << counted << " reported target proteoforms: " << line
                  << "; ambiguous: " << open_line << std::endl;
}

void FLASHTnTAlgorithm::setNativeTermini_()
{
  // Native mode reads ladders and tags from any start or end of a protein, and its spectra often have no precursor mass. A PrSM without a
  // proteoform mass has no mass that places its termini, so each terminus rests on its own side's fragments (N: a b c; C: x y z), yet
  // FLASHExtender reports one without them too: the protein's own terminus with its free forms (initiator Met removal, protein N-terminal
  // acetyl). Such a terminus is reported undetermined (-1) when its side's terminal fragment envelopes on the raw spectrum (the envelope
  // rule of the raw score) are as likely by chance: P(k or more) > 0.01, k ~ Poisson(lambda), lambda = the side's envelopes per decoy
  // fragment set (decoy_shifts). Measured with FLASHAnnotator and the same rule: PXD047341 bacteriorhodopsin 3.4 N-side envelopes per PrSM
  // against 2.3 per decoy set (17 of its 21 PrSMs started at residue 1 or 2 with acetyl), bovine carbonic anhydrase the N side of 101 of
  // 103 PrSMs, the heavy-chain C-termini of the native mAb spectra; ADH none (its acetylated Ser2 has its fragments). In the lysates the
  // rule removed no error: 1% proteoforms of "undetermined, both termini" with such a terminus were 580 E. coli and 5 B. theta (0.85%),
  // the others 705 and 4 (0.56%) (B. theta entrapment, nine runs), so it applies to native mode, where termini are reported, not counted.
  // It runs after the FDR, the proteins and the combined coverage, which take the termini as the search found them.
  if (! raw_loader_) return;
  std::vector<Size> hits;
  std::set<int> scans;
  for (Size i = 0; i < proteoform_hits_.size(); i++)
  {
    const auto& h = proteoform_hits_[i];
    if ((double)h.getMetaValue(tntKey("Mass"), -1.0) > 0) continue;
    if ((int)h.getMetaValue(tntKey("StartPosition")) <= 0 && (int)h.getMetaValue(tntKey("EndPosition")) <= 0) continue;
    hits.push_back(i);
    scans.insert((int)h.getMetaValue(tntKey("Scan")));
  }
  if (hits.empty()) return;
  loadRawSpectra_(scans);
  FLASHAnnotatorAlgorithm annotator;
  Param aparam = annotator_param_;
  aparam.setValue("internal_ions", "false");
  aparam.setValue("report_decoy_envelopes", "true");
  annotator.setParameters(aparam);
  const std::vector<double> shifts = aparam.getValue("decoy_shifts").toDoubleVector();
  if (shifts.empty())
  {
    releaseRawSpectra_();
    return;
  }
  std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms;
  std::vector<double> precursor_masses;
  std::vector<Size> hit_of; // proteoform_hits_ index of prsms[j]
  double max_mass = 0;
  for (Size i : hits)
  {
    FLASHAnnotatorAlgorithm::PrSM prsm;
    double precursor_mass = -1.0;
    if (! toRawPrSM_(proteoform_hits_[i], prsm, precursor_mass)) continue;
    max_mass = std::max({max_mass, precursor_mass, FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)});
    prsms.push_back(std::move(prsm));
    precursor_masses.push_back(precursor_mass);
    hit_of.push_back(i);
  }
  std::vector<double> iso_mzs, iso_intensities;
  annotator.getIsotopeEnvelope(max_mass + 200.0, 1, iso_mzs, iso_intensities);
  for (Size j = 0; j < prsms.size(); j++)
  {
    int reason = 4;
    prsms[j].precursor_charge = annotator.precursorChargeCap(raw_spectra_[prsms[j].spectrum_index], precursor_masses[j], reason);
  }
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  annotator.annotate(raw_peaks_, prsms, results);
  const auto& ion_types = annotator.getIonTypes();
  Size undetermined[2] = {0, 0};
  for (Size j = 0; j < prsms.size(); j++)
  {
    if (! results[j].annotated) continue;
    double k[2] = {0, 0}, d[2] = {0, 0}; // per side: target envelopes, decoy envelopes
    for (const auto& e : results[j].envelopes)
    {
      if (e.first >= 0 || e.ion_type < 0 || e.ion_type >= (int)ion_types.size() || ion_types[e.ion_type].empty()) continue; // terminal only
      const int side = (ion_types[e.ion_type][0] == 'a' || ion_types[e.ion_type][0] == 'b' || ion_types[e.ion_type][0] == 'c') ? 0 : 1;
      (e.decoy_shift == 0 ? k : d)[side] += 1;
    }
    auto& hit = proteoform_hits_[hit_of[j]];
    for (int side = 0; side < 2; side++)
    {
      const auto key = tntKey(side == 0 ? "StartPosition" : "EndPosition");
      if ((int)hit.getMetaValue(key) <= 0) continue;
      const double lambda = std::max(d[side] / (double)shifts.size(), 1e-12);
      const double p = k[side] > 0 ? boost::math::gamma_p(k[side], lambda) : 1.0;
      if (p <= 0.01) continue;
      hit.setMetaValue(key, -1);
      undetermined[side]++;
    }
  }
  releaseRawSpectra_();
  OPENMS_LOG_INFO << "Native termini: " << prsms.size() << " PrSMs without a proteoform mass; reported undetermined for fragments at chance: "
                  << undetermined[0] << " N-termini, " << undetermined[1] << " C-termini" << std::endl;
}

void FLASHTnTAlgorithm::loadRawSpectra_(const std::set<int>& scans)
{
  // Only the spectra not loaded yet are read, and they join those loaded: the stages, the chimeric rounds and the combined coverage of the
  // proteoforms (setCombinedCoverage_) reuse what an earlier step read. A raw file given with -in_raw is read in full for every call that
  // needs a spectrum not loaded yet (FLASHAnnotatorAlgorithm::loadSpectra).
  std::set<int> missing;
  for (int scan : scans)
  {
    if (raw_slot_.count(scan) == 0) missing.insert(scan);
  }
  if (missing.empty()) return;
  std::map<int, MSSpectrum> raw;
  raw_loader_(missing, raw);
  const Size first = raw_spectra_.size();
  raw_spectra_.reserve(first + raw.size());
  for (auto& [scan, spec] : raw)
  {
    raw_slot_[scan] = (int)raw_spectra_.size();
    raw_spectra_.push_back(std::move(spec));
  }
  raw_peaks_.resize(raw_spectra_.size());
  int profile_count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : profile_count)
  for (int i = (int)first; i < (int)raw_spectra_.size(); i++)
  {
    profile_count += FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(raw_spectra_[i], raw_peaks_[i]) ? 1 : 0;
    raw_spectra_[i].clear(false); // the precursor information stays
    raw_spectra_[i].shrink_to_fit();
  }
  raw_profile_count_ += profile_count;
}

void FLASHTnTAlgorithm::releaseRawSpectra_()
{
  std::vector<MSSpectrum>().swap(raw_spectra_);
  std::vector<FLASHAnnotatorAlgorithm::PeakArrays>().swap(raw_peaks_);
  raw_slot_.clear();
  raw_profile_count_ = 0;
}

bool FLASHTnTAlgorithm::toRawPrSM_(const ProteinHit& hit, FLASHAnnotatorAlgorithm::PrSM& prsm, double& precursor_mass) const
{
  const UInt mi_entry_start = MetaInfoInterface::metaRegistry().registerName("EntryStart"),
             mi_entry_end = MetaInfoInterface::metaRegistry().registerName("EntryEnd");
  const std::string& seq = hit.getSequence();
  const int start = hit.getMetaValue("StartPosition"), end = hit.getMetaValue("EndPosition");
  // an open terminus is bounded by the database entry's processed form (EntryStart, EntryEnd; a FASTA protein: its ends)
  const int entry_first = hit.metaValueExists(mi_entry_start) ? (int)hit.getMetaValue(mi_entry_start) - 1 : 0;
  const int entry_end = hit.metaValueExists(mi_entry_end) ? (int)hit.getMetaValue(mi_entry_end) : (int)seq.size();
  const int start_in_seq = start <= 0 ? std::min((int)seq.size(), entry_first) : std::min((int)seq.size(), start - 1);
  const int end_in_seq = end <= 0 ? std::max(start_in_seq, std::min((int)seq.size(), entry_end)) : std::max(start_in_seq, std::min((int)seq.size(), end));
  prsm.sequence = seq.substr(start_in_seq, end_in_seq - start_in_seq);
  prsm.known_n_term = start > 0;
  prsm.known_c_term = end > 0;
  if (hit.metaValueExists("Modifications"))
  {
    const std::vector<double> masses = hit.getMetaValue("Modifications");
    const std::vector<int> starts = hit.getMetaValue("ModificationStarts"), ends = hit.getMetaValue("ModificationEnds");
    for (Size m = 0; m < masses.size() && m < starts.size() && m < ends.size(); m++)
    {
      FLASHAnnotatorAlgorithm::Modification mod;
      mod.mass = masses[m];
      mod.start = starts[m] - start_in_seq;
      mod.end = ends[m] - start_in_seq;
      prsm.modifications.push_back(mod);
    }
  }
  precursor_mass = hit.metaValueExists("GivenMass") ? (double)hit.getMetaValue("GivenMass") : -1.0;
  const auto it = raw_slot_.find((int)hit.getMetaValue("Scan"));
  if (it == raw_slot_.end()) return false;
  prsm.spectrum_index = it->second;
  return true;
}

std::function<double(const ProteinHit&)> FLASHTnTAlgorithm::alternativeRawScorer_() const
{
  // as annotateOnRawSpectra_, one PrSM at a time; serial, as the scores of FLASHAnnotatorAlgorithm::annotate are (std::lgamma)
  auto annotator = std::make_shared<FLASHAnnotatorAlgorithm>();
  Param aparam = annotator_param_;
  aparam.setValue("internal_ions", "false");
  annotator->setParameters(aparam);
  return [this, annotator](const ProteinHit& hit) {
    FLASHAnnotatorAlgorithm::PrSM prsm;
    double precursor_mass = -1.0;
    if (! toRawPrSM_(hit, prsm, precursor_mass)) return 0.0;
    std::vector<FLASHAnnotatorAlgorithm::Result> results;
#pragma omp critical(flashtnt_alternative_raw_score)
    {
      std::vector<double> iso_mzs, iso_intensities;
      annotator->getIsotopeEnvelope(std::max(precursor_mass, FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)) + 200.0, 1,
                                    iso_mzs, iso_intensities);
      int reason = 4;
      prsm.precursor_charge = annotator->precursorChargeCap(raw_spectra_[prsm.spectrum_index], precursor_mass, reason);
      annotator->annotate(raw_peaks_, {prsm}, results);
    }
    return raw_score_weight * std::max(0.0, results[0].score);
  };
}

void FLASHTnTAlgorithm::annotateOnRawSpectra_(std::vector<ProteinHit>& hits) const
{
  // The extension score counts deconvolved masses on the path; the raw score (FLASHAnnotatorAlgorithm) is the chance of the fragment isotope
  // envelopes found in the raw spectrum, at every charge, with the chance rate measured on mass-shifted decoy fragments of the same spectrum.
  // Decoy PrSMs are annotated the same way, so the FDR of the combined score stays target-decoy controlled. Eclip_1480: at 1%, 4846 PrSMs
  // by the extension score, 5189 by the raw score alone, 5241 by extension + 4 x raw, 5247 at 8 x raw (5240 at 16, 5012 at 1); the
  // three lysates: 4278 proteoforms / 991 proteins at 4, 4293 / 1002 at 8. The raw spectra are those of loadRawSpectra_.
  FLASHAnnotatorAlgorithm annotator;
  // only the terminal results are used below (RawScore, its fragment counts, RawCoverage): internal fragments are not annotated
  Param aparam = annotator_param_;
  aparam.setValue("internal_ions", "false");
  annotator.setParameters(aparam);
  std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms(hits.size());
  std::vector<double> precursor_masses(hits.size(), -1.0);
  double max_mass = 0;
  Size without_spectrum = 0;
  for (Size i = 0; i < hits.size(); i++)
  {
    auto& prsm = prsms[i];
    if (! toRawPrSM_(hits[i], prsm, precursor_masses[i])) without_spectrum++;
    max_mass = std::max({max_mass, precursor_masses[i], FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)});
  }
  // the isotope table once, beyond the heaviest PrSM; then the precursor charge caps
  std::vector<double> iso_mzs, iso_intensities;
  annotator.getIsotopeEnvelope(max_mass + 200.0, 1, iso_mzs, iso_intensities);
  for (Size i = 0; i < prsms.size(); i++)
  {
    if (prsms[i].spectrum_index < 0) continue;
    int reason = 4;
    prsms[i].precursor_charge = annotator.precursorChargeCap(raw_spectra_[prsms[i].spectrum_index], precursor_masses[i], reason);
  }
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  annotator.annotate(raw_peaks_, prsms, results);

  Size annotated = 0;
  for (Size i = 0; i < hits.size(); i++)
  {
    auto& hit = hits[i];
    const auto& res = results[i];
    annotated += res.annotated ? 1 : 0;
    hit.setMetaValue(tntKey("ExtensionScore"), hit.getScore());
    hit.setMetaValue(tntKey("RawScore"), res.score);
    hit.setMetaValue(tntKey("RawMatchedFragments"), res.matched_fragments);
    hit.setMetaValue(tntKey("RawFragmentHypotheses"), res.fragment_hypotheses);
    hit.setMetaValue(tntKey("RawCoverage"), res.coverage);
    // with a blind shift, the raw score's chance probability is multiplied by the candidates the shift opens (shiftedCandidateCount): a
    // shift fits the precursor mass for any window within max_total_mod_mass_ of it, and absorbed the sequence differences of a similar
    // protein. Before, a fixed extension penalty stood for it: on the B. theta entrapment the determined proteoforms with one shift held
    // 15 B. theta of 2188, without one 6 of 1786 (and preferring PrSMs with a mass over those without took in 25 more, all with a shift).
    double raw = std::max(0.0, res.score);
    if (blindModCount_(hit) > 0)
    {
      const double mass = (double)hit.getMetaValue(tntKey("Mass"), -1.0);
      raw = std::max(0.0, raw - std::log10(shiftedCandidateCount(hit.getSequence(), mass, max_total_mod_mass_)));
    }
    hit.setScore(hit.getScore() + raw_score_weight * raw);
  }
  OPENMS_LOG_INFO << "Raw spectra: " << raw_peaks_.size() << " loaded (" << raw_profile_count_ << " centroided); PrSMs annotated: " << annotated
                  << " of " << hits.size() << " (" << without_spectrum << " without their scan in the raw file)" << std::endl;
}

void FLASHTnTAlgorithm::setCombinedCoverage_(const std::vector<std::vector<Size>>& clusters, const std::vector<double>& qvalues, double rep_q, bool use_fdr)
{
  // A proteoform identified by several spectra has the raw evidence of all of them. Its combined coverage counts the terminal bonds covered
  // by its best PrSM (the most covered bonds; the representative among equal ones) and those covered in at least two of its PrSMs, by the
  // envelope rule of the raw score (FLASHAnnotatorAlgorithm, the same for every spectrum). Its chance is the same rule on the decoy fragment
  // sets of the raw score (decoy_shifts), each set combined over the PrSMs as the targets are, averaged over the sets. The PrSMs combined
  // are those of the proteoform at a PrSM q-value of at most rep_q, of one spectrum each, not of a residual round, with the termini and the
  // modification masses (to 0.05 Da) of its representative; the representative always takes part. Prototype on the f5c72eda25 outputs
  // (proteoforms with two or more such PrSMs; mean coverage %, chance in brackets): best PrSM, this rule, the union of all PrSMs, one
  // spectrum averaged from the calibrated raw spectra:
  //   Eclip_1480  35.9 (1.6), 38.4 (2.1), 43.0 (4.1), 40.4 (3.8)
  //   Caco-2      29.3 (1.9), 32.3 (3.0), 37.8 (6.7), 34.5 (5.3)
  //   ETD PXD012247 (14 fractions) 33.4 (0.9), 34.4 (1.0), 37.8 (1.6), 35.8 (1.4)
  // The union and the averaged spectrum gain coverage at the price of their chance: decoy bonds recur in the spectra of one proteoform,
  // and an averaged spectrum holds the peaks of all of them. Targets and decoys alike, so a later use in a score stays target-decoy fair.
  if (! raw_loader_) return;
  if (clusters.empty())
  {
    releaseRawSpectra_();
    return;
  }
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_chimeric = registry.registerName("ChimericRound"), mi_prsms = registry.registerName("CombinedPrSMs"),
             mi_coverage = registry.registerName("CombinedCoverage"), mi_chance = registry.registerName("CombinedCoverageChance");
  auto modKey = [](const ProteinHit& h) {
    std::vector<long> k;
    if (h.metaValueExists("Modifications"))
    {
      for (double m : h.getMetaValue("Modifications").toDoubleList()) k.push_back(std::lround(m / 0.05));
    }
    std::sort(k.begin(), k.end());
    return k;
  };
  std::vector<std::vector<Size>> members(clusters.size());
  std::set<int> scans;
  for (Size c = 0; c < clusters.size(); c++)
  {
    if (clusters[c].empty()) continue;
    const auto& rep = proteoform_hits_[clusters[c][0]];
    const int start = rep.getMetaValue(tntKey("StartPosition")), end = rep.getMetaValue(tntKey("EndPosition"));
    const auto key = modKey(rep);
    std::set<int> seen;
    for (Size i : clusters[c])
    {
      const auto& h = proteoform_hits_[i];
      if (i != clusters[c][0])
      {
        if (use_fdr && qvalues[i] > rep_q) continue;
        if (h.metaValueExists(mi_chimeric) && (int)h.getMetaValue(mi_chimeric) > 0) continue;
        if ((int)h.getMetaValue(tntKey("StartPosition")) != start || (int)h.getMetaValue(tntKey("EndPosition")) != end || modKey(h) != key) continue;
      }
      const int scan = h.getMetaValue(tntKey("Scan"));
      if (! seen.insert(scan).second) continue;
      members[c].push_back(i);
      scans.insert(scan);
    }
  }
  loadRawSpectra_(scans);

  FLASHAnnotatorAlgorithm annotator;
  Param aparam = annotator_param_;
  aparam.setValue("internal_ions", "false");
  aparam.setValue("report_decoy_envelopes", "true");
  annotator.setParameters(aparam);
  const std::vector<double> shifts = aparam.getValue("decoy_shifts").toDoubleVector();
  std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms;
  std::vector<double> precursor_masses;
  std::vector<Size> first_of(clusters.size() + 1, 0); // the PrSMs of cluster c: prsms[first_of[c], first_of[c + 1])
  double max_mass = 0;
  for (Size c = 0; c < clusters.size(); c++)
  {
    first_of[c] = prsms.size();
    for (Size i : members[c])
    {
      FLASHAnnotatorAlgorithm::PrSM prsm;
      double precursor_mass = -1.0;
      if (! toRawPrSM_(proteoform_hits_[i], prsm, precursor_mass)) continue;
      max_mass = std::max({max_mass, precursor_mass, FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)});
      prsms.push_back(std::move(prsm));
      precursor_masses.push_back(precursor_mass);
    }
  }
  first_of[clusters.size()] = prsms.size();
  std::vector<double> iso_mzs, iso_intensities;
  annotator.getIsotopeEnvelope(max_mass + 200.0, 1, iso_mzs, iso_intensities);
  for (Size j = 0; j < prsms.size(); j++)
  {
    int reason = 4;
    prsms[j].precursor_charge = annotator.precursorChargeCap(raw_spectra_[prsms[j].spectrum_index], precursor_masses[j], reason);
  }
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  annotator.annotate(raw_peaks_, prsms, results);

  Size combined = 0;
  double sum_best = 0, sum_combined = 0, sum_chance = 0;
  for (Size c = 0; c < clusters.size(); c++)
  {
    const Size first = first_of[c], n = first_of[c + 1] - first;
    if (n == 0 || members[c].empty() || ! results[first].annotated) continue;
    const int bonds = (int)prsms[first].sequence.size() - 1;
    if (bonds <= 0) continue;
    // covered[v][j][b]: bond b covered in PrSM j by the target fragments (v = 0) or by decoy set v - 1
    std::vector<std::vector<std::vector<char>>> covered(shifts.size() + 1, std::vector<std::vector<char>>(n, std::vector<char>(bonds + 1, 0)));
    for (Size j = 0; j < n; j++)
    {
      for (const auto& e : results[first + j].envelopes)
      {
        if (e.first >= 0 || e.bond < 1 || e.bond > bonds) continue; // terminal fragments only
        Size v = 0;
        if (e.decoy_shift != 0)
        {
          v = std::find_if(shifts.begin(), shifts.end(), [&e](double s) { return std::abs(s - e.decoy_shift) < 1e-9; }) - shifts.begin() + 1;
          if (v > shifts.size()) continue;
        }
        covered[v][j][e.bond] = 1;
      }
    }
    Size best = 0;
    int best_count = -1;
    for (Size j = 0; j < n; j++)
    {
      const int count = (int)std::count(covered[0][j].begin(), covered[0][j].end(), (char)1);
      if (count > best_count)
      {
        best_count = count;
        best = j;
      }
    }
    auto combinedCount = [&](Size v) {
      int count = 0;
      for (int b = 1; b <= bonds; b++)
      {
        int in = 0;
        for (Size j = 0; j < n; j++) in += covered[v][j][b];
        count += (covered[v][best][b] || in >= 2) ? 1 : 0;
      }
      return count;
    };
    const double coverage = 100.0 * combinedCount(0) / bonds;
    double chance = 0;
    for (Size v = 1; v <= shifts.size(); v++) chance += 100.0 * combinedCount(v) / bonds;
    chance = shifts.empty() ? -1.0 : chance / (double)shifts.size();
    auto& rep = proteoform_hits_[clusters[c][0]];
    rep.setMetaValue(mi_prsms, (int)n);
    rep.setMetaValue(mi_coverage, coverage);
    rep.setMetaValue(mi_chance, chance);
    if (n >= 2 && ! StringUtils::hasPrefix(rep.getAccession(), "DECOY"))
    {
      combined++;
      sum_best += 100.0 * best_count / bonds;
      sum_combined += coverage;
      sum_chance += chance;
    }
  }
  releaseRawSpectra_();
  if (combined > 0)
  {
    OPENMS_LOG_INFO << "Combined coverage: " << combined << " target proteoforms with two or more PrSMs combined; their mean coverage " << sum_best / combined
                    << "% by the best PrSM, " << sum_combined / combined << "% combined (chance " << sum_chance / combined << "%)" << std::endl;
  }
}

namespace
{
// residue masses of the offline features (no fixed modifications: the model was fit on them)
double scoringResidueMass(char c)
{
  switch (c)
  {
    case 'G': return 57.02146;
    case 'A': return 71.03711;
    case 'S': return 87.03203;
    case 'P': return 97.05276;
    case 'V': return 99.06841;
    case 'T': return 101.04768;
    case 'C': return 103.00919;
    case 'L':
    case 'I': return 113.08406;
    case 'N': return 114.04293;
    case 'D': return 115.02694;
    case 'Q': return 128.05858;
    case 'K': return 128.09496;
    case 'E': return 129.04259;
    case 'M': return 131.04049;
    case 'H': return 137.05891;
    case 'F': return 147.06841;
    case 'R': return 156.10111;
    case 'Y': return 163.06333;
    case 'W': return 186.07931;
    case 'U': return 150.95364;
    case 'O': return 237.14773;
    default: return -1;
  }
}

// common modification masses (common_shift: every shift within 0.05 Da of one of them)
constexpr double scoring_common_shifts[] = {15.994915, 31.98983,   42.010565, 14.01565,  28.0313,    42.04695,   79.966331, 27.994915, 0.984016,  21.981943,
                                            37.955882, 57.021464, -18.010565, -17.026549, 43.005814, 100.016044, -2.01565,  -131.040485, -89.02992, 114.042927,
                                            15.010899, 47.984744, 79.956815,  72.021129,  86.000394, 162.052824, 203.079373, 304.20536};

// fragment-evidence features of a proteoform in a deconvolved spectrum (masses ascending, intensities aligned); empty if a residue is unknown
std::map<std::string, double> scoringEvidence(const std::string& seq, const std::vector<FLASHAnnotatorAlgorithm::Modification>& mods, bool known_n, bool known_c,
                                         const std::vector<double>& masses, const std::vector<double>& intensities, const std::vector<double>& n_offsets,
                                         const std::vector<double>& c_offsets)
{
  const int L = (int)seq.size();
  if (L < 2) return {};
  std::vector<double> prefix(L + 1, 0.0);
  for (int i = 0; i < L; i++)
  {
    const double m = scoringResidueMass(seq[i]);
    if (m < 0) return {};
    prefix[i + 1] = prefix[i] + m;
  }
  double total_intensity = std::accumulate(intensities.begin(), intensities.end(), 0.0);
  if (total_intensity == 0) total_intensity = 1.0;
  double m_any = 0, m_exact = 0, n_side = 0, c_side = 0, matched_intensity = 0;
  std::set<Size> used;
  std::vector<double> errors;
  std::vector<bool> cut_hit(L - 1, false);
  constexpr double iso = 1.0033548;
  auto find = [&](double x) -> long { // index of a mass within 10 ppm of x (the first at or above x - 10 ppm), -1 if none
    const auto it = std::lower_bound(masses.begin(), masses.end(), x - 10e-6 * x);
    return it != masses.end() && *it <= x + 10e-6 * x ? (long)(it - masses.begin()) : -1;
  };
  for (int cut = 1; cut < L; cut++)
  {
    double left = 0, right = 0, inside = 0;
    for (const auto& mod : mods)
    {
      if (mod.end < cut) left += mod.mass;
      if (mod.start >= cut) right += mod.mass;
      if (mod.start < cut && cut <= mod.end) inside += mod.mass;
    }
    struct Candidate { bool n; double m1, m2; };
    std::vector<Candidate> candidates;
    if (known_n)
      for (double o : n_offsets) candidates.push_back({true, prefix[cut] + left + o, prefix[cut] + left + inside + o});
    if (known_c)
      for (double o : c_offsets) candidates.push_back({false, prefix[L] - prefix[cut] + right + o, prefix[L] - prefix[cut] + right + inside + o});
    for (const auto& cand : candidates)
    {
      long index = -1;
      int iso_off = 0;
      double ppm = 0;
      for (int which = 0; which < (inside == 0 ? 1 : 2) && index < 0; which++)
      {
        const double m = which == 0 ? cand.m1 : cand.m2;
        for (int k : {0, -1, 1})
        {
          if (k != 0 && m < 5000) continue;
          const double x = m + k * iso;
          index = find(x);
          if (index >= 0)
          {
            iso_off = k;
            ppm = (masses[index] - x) / x * 1e6;
            break;
          }
        }
      }
      if (index < 0) continue;
      m_any++;
      m_exact += iso_off == 0;
      cut_hit[cut - 1] = true;
      (cand.n ? n_side : c_side)++;
      if (used.insert((Size)index).second) matched_intensity += intensities[index];
      errors.push_back(std::fabs(ppm));
    }
  }
  int longest = 0, run = 0, hit_cuts = 0;
  for (bool h : cut_hit)
  {
    run = h ? run + 1 : 0;
    longest = std::max(longest, run);
    hit_cuts += h;
  }
  double err_med = 10.0;
  if (! errors.empty())
  {
    std::sort(errors.begin(), errors.end());
    const Size n = errors.size();
    err_med = n % 2 ? errors[n / 2] : (errors[n / 2 - 1] + errors[n / 2]) / 2;
  }
  return {{"m_any", m_any}, {"m_exact", m_exact}, {"cuts_hit", (double)hit_cuts}, {"cut_frac", (double)hit_cuts / (L - 1)},
          {"longest_run", (double)longest}, {"n_side", n_side}, {"c_side", c_side}, {"int_frac", matched_intensity / total_intensity},
          {"peaks_used", (double)used.size()}, {"err_med", err_med}};
}

// sub_like (a shift equal to replacing one residue of its window by another amino acid, 0.02 Da) and wlen (the longest shift window)
std::pair<double, double> scoringShiftShape(const std::string& seq, const std::vector<FLASHAnnotatorAlgorithm::Modification>& mods)
{
  static const std::string residues = "GASPVTCLNDQKEMHFRYW";
  double sub = 0, wlen = 0;
  for (const auto& mod : mods)
  {
    const int s = std::max(mod.start, 0), e = mod.end;
    wlen = std::max(wlen, (double)(e - s + 1));
    if (sub > 0) continue;
    for (int i = s; i <= e && i < (int)seq.size() && sub == 0; i++)
    {
      const char r = seq[i] == 'I' ? 'L' : seq[i];
      if (residues.find(r) == std::string::npos) continue;
      for (char a : residues)
      {
        const double d = scoringResidueMass(a) - scoringResidueMass(r);
        if (a != r && std::fabs(d) > 0.05 && std::fabs(d - mod.mass) <= 0.02)
        {
          sub = 1;
          break;
        }
      }
    }
  }
  return {sub, wlen};
}
} // namespace

void FLASHTnTAlgorithm::applyLearnedScoring_(std::vector<ProteinHit>& hits, const std::map<int, const MSSpectrum*>& spectra) const
{
  const auto& model = *scoring_model_;
  constexpr double h2o = 18.010565; // as the offline features
  std::vector<double> n_offsets, c_offsets;
  for (const auto& t : param_.getValue("ion_type").toStringVector())
  {
    if (t == "a") n_offsets.push_back(-27.994915);
    else if (t == "b") n_offsets.push_back(0.0);
    else if (t == "c") n_offsets.push_back(17.026549);
    else if (t == "x") c_offsets.push_back(h2o + 25.979265);
    else if (t == "y") c_offsets.push_back(h2o);
    else if (t == "z") c_offsets.push_back(h2o - 17.026549);
    else if (t == "zp1") c_offsets.push_back(h2o - 16.018724);
    else if (t == "zp2") c_offsets.push_back(h2o - 15.010899);
  }
  if (n_offsets.empty()) n_offsets.push_back(0.0);
  if (c_offsets.empty()) c_offsets.push_back(h2o);

  // per hit: its features, its choice score, and a key that makes duplicate candidates count once (as the offline features did)
  const Size n = hits.size();
  std::vector<std::map<std::string, double>> features(n);
  std::vector<double> choice(n, 0.0);
  std::vector<std::string> keys(n);
#pragma omp parallel for schedule(dynamic)
  for (int64_t i = 0; i < (int64_t)n; i++)
  {
    const auto& hit = hits[i];
    FLASHAnnotatorAlgorithm::PrSM prsm;
    double precursor = -1;
    toRawPrSM_(hit, prsm, precursor);
    auto& f = features[i];
    f["score"] = hit.getScore();
    f["ext"] = (double)hit.getMetaValue(tntKey("ExtensionScore"), 0.0);
    f["raw"] = (double)hit.getMetaValue(tntKey("RawScore"), 0.0);
    f["raw_k"] = (int)hit.getMetaValue(tntKey("RawMatchedFragments"), 0);
    f["det"] = isDetermined(hit) ? 1 : 0;
    f["blind"] = (double)blindModCount_(hit);
    f["ladder_p"] = hit.metaValueExists("PrecursorLadder") ? 1 : 0;
    f["ladder_t"] = ! hit.metaValueExists("PrecursorLadder") && hit.metaValueExists("TerminalLadder") ? 1 : 0;
    double shift_abs = 0;
    bool common = ! prsm.modifications.empty();
    std::string mods;
    for (const auto& m : prsm.modifications)
    {
      shift_abs += std::fabs(m.mass);
      double best = std::numeric_limits<double>::max();
      for (double c : scoring_common_shifts) best = std::min(best, std::fabs(c - m.mass));
      common = common && best <= 0.05;
      mods += StringUtils::toStr(m.mass) + "@" + StringUtils::toStr(m.start) + "-" + StringUtils::toStr(m.end) + ";";
    }
    f["shift_abs"] = std::log1p(shift_abs);
    f["common_shift"] = common ? 1 : 0;
    f["len"] = (double)prsm.sequence.size();
    const auto [sub, wlen] = scoringShiftShape(prsm.sequence, prsm.modifications);
    f["sub_like"] = sub;
    f["wlen"] = wlen;
    const int scan = (int)hit.getMetaValue(tntKey("Scan"));
    const auto sit = spectra.find(scan);
    if (sit != spectra.end())
    {
      std::vector<std::pair<double, double>> peaks;
      for (const auto& p : *sit->second) peaks.emplace_back(p.getMZ(), p.getIntensity());
      std::sort(peaks.begin(), peaks.end());
      std::vector<double> masses, intensities;
      for (const auto& [m, in] : peaks)
      {
        masses.push_back(m);
        intensities.push_back(in);
      }
      for (const auto& [k, v] : scoringEvidence(prsm.sequence, prsm.modifications, prsm.known_n_term, prsm.known_c_term, masses, intensities, n_offsets, c_offsets))
        f[k] = v;
    }
    for (const auto& [name, w] : model.choice)
    {
      const auto it = f.find(name);
      choice[i] += w * (it == f.end() ? 0.0 : it->second);
    }
    keys[i] = hit.getAccession() + "|" + StringUtils::toStr((int)hit.getMetaValue("StartPosition", 0)) + "|" + StringUtils::toStr((int)hit.getMetaValue("EndPosition", 0)) + "|"
              + StringUtils::toStr((double)hit.getMetaValue(tntKey("Mass"), -1.0)) + "|" + mods;
  }

  if (const char* trace_file = std::getenv("FLASH_TRACE_SCORING")) // diagnostics: every candidate's features and choice score
  {
    std::ofstream out(trace_file, std::ios::app);
    out << std::setprecision(10);
    for (Size i = 0; i < n; i++)
    {
      out << (int)hits[i].getMetaValue(tntKey("Scan")) << "\t" << keys[i] << "\t" << choice[i];
      for (const auto& [k, v] : features[i]) out << "\t" << k << "=" << v;
      out << "\n";
    }
  }
  // per spectrum: the chosen PrSM (the best choice score), its spectrum context, its confidence
  std::map<int, std::vector<Size>> by_scan;
  for (Size i = 0; i < n; i++) by_scan[(int)hits[i].getMetaValue(tntKey("Scan"))].push_back(i);
  std::vector<bool> keep(n, false);
  // the candidates arrive in the order the parallel search found them: a total order by content first, so that of duplicates (the same
  // proteoform from a tag and from a ladder, which fall into different FDR classes) the same one is kept and ties of the choice score go the
  // same way whatever the number of threads
  const auto source_rank = [](const ProteinHit& h) { return h.metaValueExists("PrecursorLadder") ? 2 : h.metaValueExists("TerminalLadder") ? 1 : 0; };
  const auto content_order = [&](Size a, Size b) {
    if (keys[a] != keys[b]) return keys[a] < keys[b];
    if (isBetterPrSM(hits[a], hits[b])) return true;
    if (isBetterPrSM(hits[b], hits[a])) return false;
    if (source_rank(hits[a]) != source_rank(hits[b])) return source_rank(hits[a]) < source_rank(hits[b]);
    const double ea = hits[a].getMetaValue(tntKey("ExtensionScore"), 0.0), eb = hits[b].getMetaValue(tntKey("ExtensionScore"), 0.0);
    if (ea != eb) return ea > eb;
    const double ra = hits[a].getMetaValue(tntKey("RawScore"), 0.0), rb = hits[b].getMetaValue(tntKey("RawScore"), 0.0);
    if (ra != rb) return ra > rb;
    return hits[a].getMetaValue(tntKey("TagIndices"), std::string()).toString() < hits[b].getMetaValue(tntKey("TagIndices"), std::string()).toString();
  };
  for (auto& [scan, idx] : by_scan)
  {
    std::stable_sort(idx.begin(), idx.end(), content_order);
    std::vector<Size> unique;
    std::set<std::string> seen;
    for (Size i : idx)
      if (seen.insert(keys[i]).second) unique.push_back(i);
    Size best = unique.front();
    for (Size i : unique)
      if (choice[i] > choice[best]) best = i;
    double other = std::numeric_limits<double>::lowest(), other_protein = std::numeric_limits<double>::lowest();
    for (Size i : unique)
    {
      if (i == best) continue;
      other = std::max(other, hits[i].getScore());
      if (hits[i].getAccession() != hits[best].getAccession()) other_protein = std::max(other_protein, hits[i].getScore());
    }
    auto& f = features[best];
    f["n_cand"] = std::log1p((double)unique.size());
    f["gap1"] = hits[best].getScore() - (other == std::numeric_limits<double>::lowest() ? 0.0 : other);
    f["gap_oth"] = hits[best].getScore() - (other_protein == std::numeric_limits<double>::lowest() ? 0.0 : other_protein);
    std::vector<double> x;
    for (const auto& name : model.conf_features)
    {
      const auto it = f.find(name);
      x.push_back(it == f.end() ? 0.0 : it->second);
    }
    auto& hit = hits[best];
    hit.setMetaValue("CombinedScore", hit.getScore());
    hit.setMetaValue("ChoiceScore", choice[best]);
    double confidence = model.confidence(x);
    // the free shift's multiple-testing burden (shift_penalty): a determined PrSM with a blind shift could have fit any of the protein's windows
    // within max_mod_mass of its mass; uncapped, unlike the raw score's. Held-out seed 2 (0.5): shuffled-entrapment proteoforms 1.33 -> 1.19 %,
    // foreign 1.41 -> 1.29 %, for under 1 % of the proteoforms (emulated on the class's proteoform q-values)
    if (shift_penalty_ > 0 && isDetermined(hit) && blindModCount_(hit) > 0)
    {
      confidence -= shift_penalty_ * std::log10(shiftedCandidateCount(hit.getSequence(), (double)hit.getMetaValue(tntKey("Mass"), -1.0), max_total_mod_mass_));
    }
    hit.setScore(confidence);
    keep[best] = true;
  }
  // no self-move: with libc++ (macOS) a defaulted move assignment onto itself empties the hit's short strings (its accession)
  Size w = 0;
  for (Size i = 0; i < n; i++)
  {
    if (! keep[i]) continue;
    if (w != i) hits[w] = std::move(hits[i]);
    w++;
  }
  hits.resize(w);
}

void FLASHTnTAlgorithm::traceCandidates_(const std::vector<ProteinHit>& hits, const char* file) const
{
  // every candidate PrSM before the spectrum's choice, with what a scoring function needs to be rebuilt offline: its proteoform (residues and
  // shifts), the termini, the mass and the scores. Appended (one call per stage of the search cascade). Diagnostics only.
  std::ofstream out(file, std::ios::app);
  out << std::setprecision(10);
  for (const auto& hit : hits)
  {
    FLASHAnnotatorAlgorithm::PrSM prsm;
    double precursor = -1;
    toRawPrSM_(hit, prsm, precursor);
    std::string mods;
    for (const auto& m : prsm.modifications) mods += StringUtils::toStr(m.mass) + "@" + StringUtils::toStr(m.start) + "-" + StringUtils::toStr(m.end) + ";";
    const char* source = hit.metaValueExists("PrecursorLadder") ? "precursor_ladder" : hit.metaValueExists("TerminalLadder") ? "terminal_ladder" : "tag";
    out << (int)hit.getMetaValue(tntKey("Scan")) << "\t" << hit.getAccession() << "\t" << (int)hit.getMetaValue("StartPosition", 0) << "\t"
        << (int)hit.getMetaValue("EndPosition", 0) << "\t" << (prsm.known_n_term ? 1 : 0) << (prsm.known_c_term ? 1 : 0) << "\t" << prsm.sequence << "\t" << mods << "\t"
        << (double)hit.getMetaValue(tntKey("Mass"), -1.0) << "\t" << precursor << "\t" << hit.getScore() << "\t" << (double)hit.getMetaValue(tntKey("ExtensionScore"), 0.0) << "\t"
        << (double)hit.getMetaValue(tntKey("RawScore"), 0.0) << "\t" << (int)hit.getMetaValue(tntKey("RawMatchedFragments"), 0) << "\t"
        << (int)hit.getMetaValue(tntKey("RawFragmentHypotheses"), 0) << "\t" << (isDetermined(hit) ? 1 : 0) << "\t" << blindModCount_(hit) << "\t" << source << "\n";
  }
}

void FLASHTnTAlgorithm::keepBestPerSpectrum_(std::vector<ProteinHit>& hits, const std::string& label, const std::string& score_name)
{
  // Of a spectrum's PrSMs, the best one whose proteoform fits the precursor mass (determined) is reported, unless the best one without a
  // mass scores more than three times higher; a PrSM without a mass that wins is then replaced by the best determined PrSM of its protein,
  // if there is one (the protein by the fragments, the proteoform by the mass).
  // A higher score without a mass is weak evidence: the terminal paths are scored apart, so the b and y ions of one bond count twice there
  // and once on the precursor-mass path, which is also penalised for an internal end (Eclip_1480: RL3 180-209 at the right MS1 mass scored
  // 116, the same range without a mass 169). And a blind shift on a PrSM without a mass can take up the sequence differences of a similar
  // protein (B. theta entrapment, Eclip_1480 scan 9171: KPRS_ECOLI 243-315 at the MS1 mass 78, Q8A535_BACTN 1-799 without a mass 112).
  // A mass is no strong evidence for its protein either: one blind shift fits any precursor mass, and a weak PrSM with one (Eclip_1479
  // scan 10849: Q8A047_BACTN 44) was reported instead of the right protein without a mass (AHPC_ECOLI 449); as the spectrum then counted
  // as identified, the gapped tags of the search cascade did not see it. Where the protein changed against always keeping the determined
  // PrSM (B. theta entrapment of the Eclip lysates), the PrSM without a mass was more often a B. theta or decoy protein than the determined
  // one below three times its score, and less often above.
  constexpr double undetermined_score_factor = 3.0;
  std::map<int, Size> best_determined, best_undetermined;
  for (Size i = 0; i < hits.size(); i++)
  {
    auto& best_of_kind = isDetermined(hits[i]) ? best_determined : best_undetermined;
    const int scan = hits[i].getMetaValue(tntKey("Scan"));
    auto it = best_of_kind.find(scan);
    if (it == best_of_kind.end()) best_of_kind[scan] = i;
    else if (isBetterPrSM(hits[i], hits[it->second])) it->second = i;
  }
  std::map<int, Size> best = best_determined;
  for (const auto& [scan, u] : best_undetermined)
  {
    auto it = best.find(scan);
    if (it == best.end() || hits[u].getScore() > undetermined_score_factor * hits[it->second].getScore()) best[scan] = u;
  }
  {
    std::map<int, Size> best_with_mass;
    for (Size i = 0; i < hits.size(); i++)
    {
      const auto& h = hits[i];
      const int scan = h.getMetaValue(tntKey("Scan"));
      const auto& b = hits[best[scan]];
      if (isDetermined(b) || ! isDetermined(h) || h.getAccession() != b.getAccession()) continue;
      auto it = best_with_mass.find(scan);
      if (it == best_with_mass.end() || isBetterPrSM(h, hits[it->second])) best_with_mass[scan] = i;
    }
    for (const auto& [scan, i] : best_with_mass) best[scan] = i;
    if (! best_with_mass.empty())
      OPENMS_LOG_INFO << "PrSMs" << label << ": " << best_with_mass.size() << " spectra report the protein's PrSM with a mass instead of the one without"
                      << std::endl;
  }
  std::vector<ProteinHit> kept;
  kept.reserve(best.size());
  for (Size i = 0; i < hits.size(); i++)
  {
    if (best[(int)hits[i].getMetaValue(tntKey("Scan"))] == i) kept.push_back(std::move(hits[i]));
  }
  OPENMS_LOG_INFO << "PrSMs" << label << ": " << hits.size() << " of all candidates, " << kept.size() << " best per spectrum by the " << score_name
                  << std::endl;
  hits.swap(kept);
}

void FLASHTnTAlgorithm::getProteoformHitsMatchedBy(const FLASHHelperClasses::Tag& tag, std::vector<ProteinHit>& hits) const
{
  int index = tag.getIndex();

  if (index < 0 || matching_hits_indices_.find(index) == matching_hits_indices_.end()) return;

  for (auto i : matching_hits_indices_.at(index))
  {
    hits.push_back(proteoform_hits_[i]);
  }
}

void FLASHTnTAlgorithm::releaseResults()
{
  // hundreds of thousands of hits with sequences and meta values: their destructors in parallel instead of serially at exit
  releaseInParallel(protein_hits_);
  releaseInParallel(proteoform_hits_);
  releaseInParallel(tags_);
  matching_hits_indices_.clear();
  releaseRawSpectra_();
}

void FLASHTnTAlgorithm::getTags(std::vector<FLASHHelperClasses::Tag>& tags) const
{
  for (const auto& tag : tags_)
  {
    tags.push_back(tag);
  }
}

} // namespace OpenMS