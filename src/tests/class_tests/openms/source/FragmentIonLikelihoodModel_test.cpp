// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FragmentIonLikelihoodModel.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

using Model = FragmentIonLikelihoodModel;

// b/y ions (b1 included) with ion names and charges, as ProSE generates them
static MSSpectrum theoretical_(const AASequence& peptide, int max_charge = 1)
{
  TheoreticalSpectrumGenerator tsg;
  Param p = tsg.getParameters();
  p.setValue("add_metainfo", "true");
  p.setValue("add_first_prefix_ion", "true");
  tsg.setParameters(p);
  MSSpectrum theo;
  tsg.getSpectrum(theo, peptide, 1, max_charge);
  return theo;
}

// The run's fragmentation pattern: every y ion intense, b ions of odd ordinal weak, b ions of even
// ordinal missing, plus five noise peaks far above the fragments. Peaks deviate by @p shift_ppm.
static MSSpectrum observed_(const AASequence& peptide, double shift_ppm = 0.0)
{
  const MSSpectrum theo = theoretical_(peptide);
  const auto& names = theo.getStringDataArrays()[0];
  MSSpectrum spec;
  for (Size i = 0; i < theo.size(); ++i)
  {
    bool prefix = false;
    Size ordinal = 0;
    Model::parseIonName(names[i], prefix, ordinal);
    if (prefix && ordinal % 2 == 0) continue;
    spec.emplace_back(theo[i].getMZ() * (1.0 + shift_ppm * 1e-6), prefix ? 10.0f : 1000.0f);
  }
  for (int k = 0; k < 5; ++k) spec.emplace_back(3000.0 + 3.0 * k, 1.0f);
  spec.sortByPosition();
  return spec;
}

// All but the C-terminal residue reversed: the noise hypothesis of a confident PSM
static AASequence reversed_(const AASequence& peptide)
{
  std::string s = peptide.toUnmodifiedString();
  std::reverse(s.begin(), s.end() - 1);
  return AASequence::fromString(s);
}

static const std::vector<std::string> training_peptides_ = {
  "DFPIANGER", "TVMENFVAFVDK", "LVNELTEFAK", "YLYEIAR", "HPYFYAPELLFFAK", "AEFVEVTK",
  "DDSPDLPK", "LGEYGFQNALIVR", "QTALVELLK", "SLHTLFGDELCK", "VPQVSTPTLVEVSR", "ETYGDMADCCEK"};

// The ions of a PSM: @p peptide against the spectrum of @p spectrum_peptide, as the only peak list
static std::vector<Model::Ion> ions_(const Model& model, const AASequence& spectrum_peptide, const AASequence& peptide,
                                     int precursor_charge = 2, double tolerance = 20.0, bool ppm = true, double shift_ppm = 0.0)
{
  Model::PeakLists peaks;
  peaks.reset(1);
  peaks.assign(0, observed_(spectrum_peptide, shift_ppm));
  std::vector<Model::Ion> ions;
  model.matchIons(peaks, 0, theoretical_(peptide), peptide, precursor_charge, tolerance, ppm, ions);
  return ions;
}

// Signal: the pattern above for each training peptide; noise: the reversed peptides against the same spectra.
static Model trained_model_(Model::ContextSet context_set = Model::ContextSet::RICH, double pseudo_count = 20.0, bool reverse_order = false)
{
  Model model(context_set, pseudo_count);
  std::vector<std::string> peptides = training_peptides_;
  if (reverse_order) std::reverse(peptides.begin(), peptides.end());
  for (const std::string& s : peptides)
  {
    const AASequence peptide = AASequence::fromString(s);
    model.addObservations(ions_(model, peptide, peptide), false);
    model.addObservations(ions_(model, peptide, reversed_(peptide)), true);
  }
  model.finalize();
  return model;
}

// Flat context index of the RICH context set
static Size rich_context_(Size series, Size bucket, Size charge, Size position, Size site, Size complement)
{
  return ((((series * Model::PRECURSOR_BUCKETS + bucket) * Model::FRAGMENT_CHARGES + charge) * Model::POSITION_BINS + position)
          * Model::CLEAVAGE_SITES + site) * Model::COMPLEMENT_STATES + complement;
}

// Flat context index of the BASIC context set (#10378)
static Size basic_context_(Size series, Size bucket, Size charge, Size position)
{
  return ((series * Model::PRECURSOR_BUCKETS + bucket) * Model::FRAGMENT_CHARGES + charge) * Model::POSITION_BINS + position;
}

START_TEST(FragmentIonLikelihoodModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Model* ptr = nullptr;
Model* null_ptr = nullptr;

START_SECTION((explicit FragmentIonLikelihoodModel(ContextSet context_set = ContextSet::RICH, double pseudo_count = 20.0)))
{
  ptr = new Model();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->contextSet() == Model::ContextSet::RICH, true)
  TEST_REAL_SIMILAR(ptr->pseudoCount(), 20.0)
  TEST_EQUAL(ptr->isTrained(), false)
  TEST_EQUAL(ptr->signalPsms(), 0)
  TEST_EQUAL(ptr->noisePsms(), 0)
  // 2 series x 3 precursor buckets x 2 fragment charges x 10 positions x 3 cleavage sites x 2 complement states
  TEST_EQUAL(ptr->contexts(), 720)
  // 7 rank bins x 3 error bins, absent
  TEST_EQUAL(ptr->outcomes(), 22)
  TEST_EQUAL(ptr->absentOutcome(), 21)
  delete ptr;

  Model basic(Model::ContextSet::BASIC, 5.0);
  TEST_EQUAL(basic.contextSet() == Model::ContextSet::BASIC, true)
  TEST_REAL_SIMILAR(basic.pseudoCount(), 5.0)
  TEST_EQUAL(basic.contexts(), 120)
  TEST_EQUAL(basic.outcomes(), 8)
  TEST_EQUAL(basic.absentOutcome(), 7)

  TEST_EXCEPTION(Exception::InvalidParameter, Model(Model::ContextSet::RICH, 0.0))
  TEST_EXCEPTION(Exception::InvalidParameter, Model(Model::ContextSet::RICH, -1.0))
  TEST_EXCEPTION(Exception::InvalidParameter, Model(Model::ContextSet::BASIC, std::numeric_limits<double>::quiet_NaN()))
  TEST_EXCEPTION(Exception::InvalidParameter, Model(Model::ContextSet::BASIC, std::numeric_limits<double>::infinity()))
}
END_SECTION

START_SECTION((static Size rankBin(Size rank)))
{
  TEST_EQUAL(Model::rankBin(1), 0)
  TEST_EQUAL(Model::rankBin(2), 0)
  TEST_EQUAL(Model::rankBin(3), 1)
  TEST_EQUAL(Model::rankBin(5), 1)
  TEST_EQUAL(Model::rankBin(6), 2)
  TEST_EQUAL(Model::rankBin(10), 2)
  TEST_EQUAL(Model::rankBin(11), 3)
  TEST_EQUAL(Model::rankBin(20), 3)
  TEST_EQUAL(Model::rankBin(21), 4)
  TEST_EQUAL(Model::rankBin(40), 4)
  TEST_EQUAL(Model::rankBin(41), 5)
  TEST_EQUAL(Model::rankBin(80), 5)
  TEST_EQUAL(Model::rankBin(81), 6)
  TEST_EQUAL(Model::rankBin(100000), Model::RANK_BINS - 1)
}
END_SECTION

START_SECTION((static Size errorBin(double relative_error)))
{
  TEST_EQUAL(Model::errorBin(0.0), 0)
  TEST_EQUAL(Model::errorBin(0.2499), 0)
  TEST_EQUAL(Model::errorBin(0.25), 1)
  TEST_EQUAL(Model::errorBin(0.4999), 1)
  TEST_EQUAL(Model::errorBin(0.5), 2)
  TEST_EQUAL(Model::errorBin(1.0), Model::ERROR_BINS - 1)
}
END_SECTION

START_SECTION((static bool parseIonName(const std::string& name, bool& prefix, Size& ordinal)))
{
  bool prefix = false;
  Size ordinal = 0;
  TEST_EQUAL(Model::parseIonName("b5+", prefix, ordinal), true)
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 5)
  TEST_EQUAL(Model::parseIonName("y3++", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 3)
  TEST_EQUAL(Model::parseIonName("z.4+", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 4)
  TEST_EQUAL(Model::parseIonName("z'12+", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 12)
  TEST_EQUAL(Model::parseIonName("c7+", prefix, ordinal), true)
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 7)
  TEST_EQUAL(Model::parseIonName("a1+", prefix, ordinal), true)
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 1)
  TEST_EQUAL(Model::parseIonName("x2+", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 2)
  TEST_EQUAL(Model::parseIonName("b3-H2O1+", prefix, ordinal), true) // neutral loss of b3
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 3)
  TEST_EQUAL(Model::parseIonName("PEPTIDE$y6+", prefix, ordinal), true) // cross-link annotation
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 6)
  for (const std::string& name : {"", "b", "b+", "b0+", "[M+H]+", "iP+", "y-H2O+", "Q", "$"})
  {
    TEST_EQUAL(Model::parseIonName(name, prefix, ordinal), false)
  }
}
END_SECTION

START_SECTION((PeakLists: void reset(Size spectra), void clear(), Size size() const, void assign(Size index, const MSSpectrum& spectrum), Size peaks(Size index) const, double mz(Size index, Size peak) const, Size rankBin(Size index, Size peak) const, Size totalPeaks() const, Size memoryUsage() const))
{
  Model::PeakLists lists;
  TEST_EQUAL(lists.size(), 0)
  TEST_EQUAL(lists.totalPeaks(), 0)
  lists.reset(3);
  TEST_EQUAL(lists.size(), 3)
  TEST_EQUAL(lists.peaks(0), 0)

  // ranks: 1 = most intense, ties by ascending m/z
  MSSpectrum spec;
  spec.emplace_back(100.0, 10.0f);
  spec.emplace_back(200.0, 30.0f);
  spec.emplace_back(300.0, 20.0f);
  spec.emplace_back(400.0, 30.0f);
  lists.assign(1, spec);
  TEST_EQUAL(lists.peaks(1), 4)
  TEST_REAL_SIMILAR(lists.mz(1, 2), 300.0)
  TEST_EQUAL(lists.rankBin(1, 1), 0) // rank 1
  TEST_EQUAL(lists.rankBin(1, 3), 0) // rank 2
  TEST_EQUAL(lists.rankBin(1, 2), 1) // rank 3
  TEST_EQUAL(lists.rankBin(1, 0), 1) // rank 4
  TEST_EQUAL(lists.totalPeaks(), 4)
  TEST_TRUE(lists.memoryUsage() >= 4 * (sizeof(float) + 1))

  // The rank bins of a large spectrum equal those of a full stable sort by intensity, ties by peak order.
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> mz_dist(100.0, 2000.0);
  std::uniform_int_distribution<int> intensity_dist(1, 50); // many ties
  MSSpectrum large;
  for (int i = 0; i < 500; ++i) large.emplace_back(mz_dist(rng), static_cast<float>(intensity_dist(rng)));
  large.sortByPosition();
  lists.assign(2, large);
  std::vector<Size> order(large.size());
  std::iota(order.begin(), order.end(), Size(0));
  std::stable_sort(order.begin(), order.end(), [&large](Size a, Size b) { return large[a].getIntensity() > large[b].getIntensity(); });
  Size mismatches = 0;
  for (Size rank = 0; rank < order.size(); ++rank)
  {
    if (lists.rankBin(2, order[rank]) != Model::rankBin(rank + 1)) ++mismatches;
  }
  TEST_EQUAL(mismatches, 0)
  TEST_EQUAL(lists.totalPeaks(), 504)

  // reassigning replaces the list
  lists.assign(1, MSSpectrum());
  TEST_EQUAL(lists.peaks(1), 0)
  TEST_EQUAL(lists.totalPeaks(), 500)

  // invalid input
  TEST_EXCEPTION(Exception::IndexOverflow, lists.assign(3, spec))
  MSSpectrum unsorted;
  unsorted.emplace_back(200.0, 1.0f);
  unsorted.emplace_back(100.0, 1.0f);
  TEST_EXCEPTION(Exception::IllegalArgument, lists.assign(0, unsorted))

  lists.clear();
  TEST_EQUAL(lists.size(), 0)
  TEST_EQUAL(lists.totalPeaks(), 0)
}
END_SECTION

START_SECTION((PeakLists: Int findNearest(Size index, double mz, double tolerance) const))
{
  // Same rule as MSSpectrum::findNearest(mz, tolerance), the matching of HyperScore.
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> mz_dist(150.0, 1500.0);
  MSSpectrum spec;
  for (int i = 0; i < 300; ++i) spec.emplace_back(mz_dist(rng), 1.0f);
  spec.emplace_back(700.0, 1.0f);
  spec.emplace_back(700.02, 1.0f); // 700.01 lies exactly between these two
  spec.sortByPosition();
  // float m/z of the lists: compare against the spectrum with float-rounded m/z
  MSSpectrum rounded = spec;
  for (auto& peak : rounded) peak.setMZ(static_cast<double>(static_cast<float>(peak.getMZ())));
  Model::PeakLists lists;
  lists.reset(2);
  lists.assign(0, spec);
  Size mismatches = 0;
  for (int i = 0; i < 5000; ++i)
  {
    const double mz = 100.0 + 1500.0 * static_cast<double>(i) / 5000.0;
    for (const double tolerance : {0.001, 0.02, 0.5})
    {
      if (lists.findNearest(0, mz, tolerance) != rounded.findNearest(mz, tolerance)) ++mismatches;
    }
  }
  TEST_EQUAL(mismatches, 0)
  TEST_EQUAL(lists.findNearest(0, 700.01, 0.05), rounded.findNearest(700.01, 0.05))
  TEST_EQUAL(lists.findNearest(0, 1.0, 0.5), -1)
  TEST_EQUAL(lists.findNearest(1, 700.0, 0.5), -1) // empty list
}
END_SECTION

START_SECTION((void matchIons(const PeakLists& peaks, Size index, const MSSpectrum& theoretical, const AASequence& peptide, int precursor_charge, double tolerance, bool ppm, std::vector<Ion>& ions) const))
{
  // DFPIANGER: D0 F1 P2 I3 A4 N5 G6 E7 R8; b1-b8 and y1-y8, every y ion and b1, b3, b5, b7 in the spectrum
  const AASequence peptide = AASequence::fromString("DFPIANGER");
  const Model rich;
  const std::vector<Model::Ion> ions = ions_(rich, peptide, peptide, 2);
  TEST_EQUAL(ions.size(), 16)
  Size matched = 0;
  for (const Model::Ion& ion : ions)
  {
    TEST_TRUE(ion.context < rich.contexts())
    TEST_TRUE(ion.outcome < rich.outcomes())
    if (ion.outcome != rich.absentOutcome()) ++matched;
  }
  TEST_EQUAL(matched, 12)

  // the contexts of all ions, by name
  const MSSpectrum theo = theoretical_(peptide);
  const auto& names = theo.getStringDataArrays()[0];
  ABORT_IF(names.size() != ions.size())
  auto context_of = [&](const std::string& name) -> Size {
    for (Size i = 0; i < names.size(); ++i) if (names[i] == name) return ions[i].context;
    return std::numeric_limits<Size>::max();
  };
  auto outcome_of = [&](const std::string& name) -> Size {
    for (Size i = 0; i < names.size(); ++i) if (names[i] == name) return ions[i].outcome;
    return std::numeric_limits<Size>::max();
  };
  // b1 (D|F: after D), complementary y8 present; position 10 * 1 / 9 = 1
  TEST_EQUAL(context_of("b1+"), rich_context_(0, 0, 0, 1, 2, 1))
  // b2 (F|P: before P) is missing, its complement y7 present
  TEST_EQUAL(context_of("b2+"), rich_context_(0, 0, 0, 2, 1, 1))
  TEST_EQUAL(outcome_of("b2+"), rich.absentOutcome())
  // y7 (F|P), complement b2 missing; position 10 * 7 / 9 = 7
  TEST_EQUAL(context_of("y7+"), rich_context_(1, 0, 0, 7, 1, 0))
  // y3 (N|G), complement b6 missing
  TEST_EQUAL(context_of("y3+"), rich_context_(1, 0, 0, 3, 0, 0))
  // y2 (G|E), complement b7 present
  TEST_EQUAL(context_of("y2+"), rich_context_(1, 0, 0, 2, 0, 1))
  // b8 (E|R: after E), complement y1 present; position 10 * 8 / 9 = 8
  TEST_EQUAL(context_of("b8+"), rich_context_(0, 0, 0, 8, 2, 1))
  // outcomes: y ions are the most intense peaks (rank bins of 8 y ions: ranks 1-8), exact m/z: error bin 0
  TEST_EQUAL(outcome_of("y1+") % Model::ERROR_BINS, 0)
  TEST_TRUE(outcome_of("y1+") / Model::ERROR_BINS <= 2)
  // b ions (intensity 10) rank below the y ions: ranks 9-12, bins 2-3
  TEST_TRUE(outcome_of("b3+") / Model::ERROR_BINS >= 2)

  // precursor charge buckets
  TEST_EQUAL(ions_(rich, peptide, peptide, 3)[0].context / (Model::FRAGMENT_CHARGES * Model::POSITION_BINS * Model::CLEAVAGE_SITES * Model::COMPLEMENT_STATES) % Model::PRECURSOR_BUCKETS, 1)
  TEST_EQUAL(ions_(rich, peptide, peptide, 5)[0].context / (Model::FRAGMENT_CHARGES * Model::POSITION_BINS * Model::CLEAVAGE_SITES * Model::COMPLEMENT_STATES) % Model::PRECURSOR_BUCKETS, 2)

  // mass error bins relative to the tolerance: 8 ppm off at 20 ppm is 0.4 of the tolerance, 12 ppm 0.6
  const std::vector<Model::Ion> shifted = ions_(rich, peptide, peptide, 2, 20.0, true, 8.0);
  const std::vector<Model::Ion> shifted_more = ions_(rich, peptide, peptide, 2, 20.0, true, 12.0);
  const std::vector<Model::Ion> outside = ions_(rich, peptide, peptide, 2, 20.0, true, 25.0);
  for (Size i = 0; i < ions.size(); ++i)
  {
    if (ions[i].outcome == rich.absentOutcome()) continue;
    TEST_EQUAL(shifted[i].outcome, ions[i].outcome + 1)
    TEST_EQUAL(shifted_more[i].outcome, ions[i].outcome + 2)
    TEST_EQUAL(outside[i].outcome, rich.absentOutcome())
    TEST_EQUAL(shifted[i].context, ions[i].context)
  }

  // BASIC: the contexts and rank outcomes of #10378
  const Model basic(Model::ContextSet::BASIC);
  const std::vector<Model::Ion> basic_ions = ions_(basic, peptide, peptide, 2);
  ABORT_IF(basic_ions.size() != ions.size())
  for (Size i = 0; i < ions.size(); ++i)
  {
    bool prefix = false;
    Size ordinal = 0;
    Model::parseIonName(names[i], prefix, ordinal);
    TEST_EQUAL(basic_ions[i].context, basic_context_(prefix ? 0 : 1, 0, 0, std::min<Size>(9, 10 * ordinal / 9)))
    TEST_EQUAL(basic_ions[i].outcome, ions[i].outcome == rich.absentOutcome() ? basic.absentOutcome() : ions[i].outcome / Model::ERROR_BINS)
  }

  // doubly charged fragments get their own charge bin
  Model::PeakLists peaks;
  peaks.reset(1);
  peaks.assign(0, observed_(peptide));
  std::vector<Model::Ion> two_charges;
  rich.matchIons(peaks, 0, theoretical_(peptide, 2), peptide, 3, 20.0, true, two_charges);
  TEST_EQUAL(two_charges.size(), 32)
  Size charge2 = 0;
  for (const Model::Ion& ion : two_charges)
  {
    if (ion.context / (Model::POSITION_BINS * Model::CLEAVAGE_SITES * Model::COMPLEMENT_STATES) % Model::FRAGMENT_CHARGES == 1) ++charge2;
  }
  TEST_EQUAL(charge2, 16)

  // invalid input
  std::vector<Model::Ion> out;
  MSSpectrum unnamed = theo;
  unnamed.getStringDataArrays().clear();
  TEST_EXCEPTION(Exception::InvalidValue, rich.matchIons(peaks, 0, unnamed, peptide, 2, 20.0, true, out))
  TEST_EXCEPTION(Exception::InvalidParameter, rich.matchIons(peaks, 0, theo, peptide, 2, 0.0, true, out))
  TEST_EXCEPTION(Exception::InvalidParameter, rich.matchIons(peaks, 0, theo, peptide, 2, -1.0, false, out))
  TEST_EXCEPTION(Exception::InvalidParameter, rich.matchIons(peaks, 0, theo, peptide, 2, std::numeric_limits<double>::infinity(), true, out))

  // an empty peak list leaves every ion absent
  Model::PeakLists empty;
  empty.reset(1);
  rich.matchIons(empty, 0, theo, peptide, 2, 20.0, true, out);
  TEST_EQUAL(out.size(), 16)
  TEST_EQUAL(std::count_if(out.begin(), out.end(), [&rich](const Model::Ion& ion) { return ion.outcome == rich.absentOutcome(); }), 16)
}
END_SECTION

START_SECTION((void addObservations(const std::vector<Ion>& ions, bool noise)))
{
  const AASequence peptide = AASequence::fromString("DFPIANGER");
  Model model;
  model.addObservations(ions_(model, peptide, peptide), false);
  TEST_EQUAL(model.signalPsms(), 1)
  TEST_EQUAL(model.noisePsms(), 0)
  TEST_EQUAL(model.isTrained(), false)
  model.addObservations(ions_(model, peptide, reversed_(peptide)), true);
  TEST_EQUAL(model.signalPsms(), 1)
  TEST_EQUAL(model.noisePsms(), 1)
  model.finalize();
  TEST_EQUAL(model.isTrained(), true)
  // later observations need another finalize()
  model.addObservations(ions_(model, peptide, peptide), false);
  TEST_EQUAL(model.isTrained(), false)
  TEST_EQUAL(model.signalPsms(), 2)
}
END_SECTION

START_SECTION((void finalize()))
{
  // Without observations every outcome is equally likely under signal and noise: the flat prior.
  for (const Model::ContextSet context_set : {Model::ContextSet::BASIC, Model::ContextSet::RICH})
  {
    Model model(context_set);
    model.finalize();
    TEST_EQUAL(model.isTrained(), true)
    for (Size outcome = 0; outcome < model.outcomes(); ++outcome)
    {
      TEST_REAL_SIMILAR(model.logLikelihoodRatio(3, outcome), 0.0)
    }
    TEST_REAL_SIMILAR(model.presenceProbability(3), 1.0 - 1.0 / static_cast<double>(model.outcomes()))
  }

  // The counts are integers: the model does not depend on the order of the PSMs.
  const Model forward = trained_model_();
  const Model backward = trained_model_(Model::ContextSet::RICH, 20.0, true);
  Size differences = 0;
  for (Size context = 0; context < forward.contexts(); ++context)
  {
    for (Size outcome = 0; outcome < forward.outcomes(); ++outcome)
    {
      if (forward.logLikelihoodRatio(context, outcome) != backward.logLikelihoodRatio(context, outcome)) ++differences;
    }
  }
  TEST_EQUAL(differences, 0)
}
END_SECTION

START_SECTION((bool isTrained() const))
{
  NOT_TESTABLE // see addObservations() and finalize()
}
END_SECTION

START_SECTION((Size signalPsms() const))
{
  NOT_TESTABLE // see addObservations()
}
END_SECTION

START_SECTION((Size noisePsms() const))
{
  NOT_TESTABLE // see addObservations()
}
END_SECTION

START_SECTION((double pseudoCount() const))
{
  NOT_TESTABLE // see the constructor
}
END_SECTION

START_SECTION((double logLikelihoodRatio(Size context, Size outcome) const))
{
  // y7 and b2 of a 9-mer at precursor charge 2, singly charged, bond before P
  const Size suffix = rich_context_(1, 0, 0, 7, 1, 0);
  const Size prefix = rich_context_(0, 0, 0, 2, 1, 1);
  Model untrained;
  TEST_EXCEPTION(Exception::Precondition, untrained.logLikelihoodRatio(suffix, 0))

  const Model model = trained_model_();
  // y ions are always present and intense in the training spectra, rarely in their reversed sequences
  TEST_TRUE(model.logLikelihoodRatio(suffix, 0) > 0.0)
  TEST_TRUE(model.logLikelihoodRatio(suffix, model.absentOutcome()) < 0.0)
  // b ions are weak when present and absent half of the time: a missing b ion says less than a missing y ion
  TEST_TRUE(model.logLikelihoodRatio(prefix, model.absentOutcome()) > model.logLikelihoodRatio(suffix, model.absentOutcome()))
  // outcomes beyond the table count as absent
  TEST_REAL_SIMILAR(model.logLikelihoodRatio(suffix, 1000), model.logLikelihoodRatio(suffix, model.absentOutcome()))
  // unseen contexts (precursor charge 5, doubly charged) back off to their parents and stay finite
  const Size unseen = rich_context_(1, 2, 1, 7, 0, 1);
  for (Size outcome = 0; outcome < model.outcomes(); ++outcome)
  {
    TEST_TRUE(std::isfinite(model.logLikelihoodRatio(unseen, outcome)))
  }

  // the same for the BASIC context set
  const Model basic = trained_model_(Model::ContextSet::BASIC);
  const Size basic_suffix = basic_context_(1, 0, 0, 7);
  TEST_TRUE(basic.logLikelihoodRatio(basic_suffix, 0) > 0.0)
  TEST_TRUE(basic.logLikelihoodRatio(basic_suffix, basic.absentOutcome()) < 0.0)
  const Size basic_unseen = basic_context_(1, 2, 1, 7);
  TEST_TRUE(basic.logLikelihoodRatio(basic_unseen, 0) > 0.0)
  TEST_TRUE(basic.logLikelihoodRatio(basic_unseen, basic.absentOutcome()) < 0.0)
}
END_SECTION

START_SECTION((double presenceProbability(Size context) const))
{
  const Size suffix = rich_context_(1, 0, 0, 4, 0, 0);
  const Size prefix = rich_context_(0, 0, 0, 4, 0, 1);
  Model untrained;
  TEST_EXCEPTION(Exception::Precondition, untrained.presenceProbability(suffix))

  const Model model = trained_model_();
  TEST_TRUE(model.presenceProbability(suffix) > 0.8)
  TEST_TRUE(model.presenceProbability(suffix) <= 1.0)
  TEST_TRUE(model.presenceProbability(prefix) > 0.1)
  TEST_TRUE(model.presenceProbability(prefix) < 0.9)
  TEST_TRUE(model.presenceProbability(suffix) > model.presenceProbability(prefix))
  // a smaller pseudo-count follows the counts more closely
  const Model sharp = trained_model_(Model::ContextSet::RICH, 1.0);
  TEST_TRUE(sharp.presenceProbability(suffix) > model.presenceProbability(suffix))
}
END_SECTION

START_SECTION((Features score(const std::vector<Ion>& ions, Size top_k = 6) const))
{
  const AASequence peptide = AASequence::fromString("DFPIANGER");
  Model untrained;
  TEST_EXCEPTION(Exception::Precondition, untrained.score(ions_(untrained, peptide, peptide)))

  for (const Model::ContextSet context_set : {Model::ContextSet::RICH, Model::ContextSet::BASIC})
  {
    const Model model = trained_model_(context_set);
    const Model::Features features = model.score(ions_(model, peptide, peptide));
    TEST_EQUAL(features.theoretical_ions, 16) // b1-b8, y1-y8
    TEST_EQUAL(features.matched_ions, 12)     // every y ion and b1, b3, b5, b7
    TEST_TRUE(features.log_likelihood_ratio > 0.0)
    TEST_TRUE(features.explained_presence > 0.6)
    TEST_TRUE(features.explained_presence <= 1.0)
    TEST_REAL_SIMILAR(features.top_predicted_observed, 1.0) // the six most likely ions are y ions, all present

    // The reversed sequence explains the same spectrum badly.
    const Model::Features noise = model.score(ions_(model, peptide, reversed_(peptide)));
    TEST_EQUAL(noise.theoretical_ions, 16)
    TEST_TRUE(noise.matched_ions < features.matched_ions)
    TEST_TRUE(noise.log_likelihood_ratio < 0.0)
    TEST_TRUE(noise.explained_presence < features.explained_presence)
    TEST_TRUE(noise.top_predicted_observed < features.top_predicted_observed)

    // A peptide the model never saw, at a precursor charge it never saw, is scored by the run's pattern.
    const AASequence unseen = AASequence::fromString("GLSDGEWQQVLNVWGK");
    const Model::Features unseen_features = model.score(ions_(model, unseen, unseen, 3));
    const Model::Features unseen_noise = model.score(ions_(model, unseen, reversed_(unseen), 3));
    TEST_TRUE(unseen_features.log_likelihood_ratio > 0.0)
    TEST_TRUE(unseen_features.log_likelihood_ratio > unseen_noise.log_likelihood_ratio)

    // Da and ppm tolerances match the same peaks here.
    const Model::Features dalton = model.score(ions_(model, peptide, peptide, 2, 0.02, false));
    TEST_EQUAL(dalton.matched_ions, features.matched_ions)

    // top_k bounds the third feature; zero disables it
    TEST_REAL_SIMILAR(model.score(ions_(model, peptide, peptide), 0).top_predicted_observed, 0.0)
    TEST_REAL_SIMILAR(model.score(ions_(model, peptide, peptide), 100).top_predicted_observed, 12.0 / 16.0)

    // nothing to score
    const Model::Features none = model.score(std::vector<Model::Ion>());
    TEST_EQUAL(none.theoretical_ions, 0)
    TEST_EQUAL(none.matched_ions, 0)
    TEST_REAL_SIMILAR(none.log_likelihood_ratio, 0.0)
    TEST_REAL_SIMILAR(none.explained_presence, 0.0)
    TEST_REAL_SIMILAR(none.top_predicted_observed, 0.0)
  }

  // RICH: matches far off in mass are worth less than exact ones (the training matches are exact)
  const Model model = trained_model_();
  const Model::Features exact = model.score(ions_(model, peptide, peptide));
  const Model::Features off = model.score(ions_(model, peptide, peptide, 2, 20.0, true, 15.0));
  TEST_EQUAL(off.matched_ions, exact.matched_ions)
  TEST_TRUE(off.log_likelihood_ratio < exact.log_likelihood_ratio)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
