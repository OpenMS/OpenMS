// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
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
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

// b/y ions (1+, b1 included) with ion names and charges, as ProSE generates them
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
// ordinal missing, plus five noise peaks far above the fragments.
static MSSpectrum observed_(const AASequence& peptide)
{
  const MSSpectrum theo = theoretical_(peptide);
  const auto& names = theo.getStringDataArrays()[0];
  MSSpectrum spec;
  for (Size i = 0; i < theo.size(); ++i)
  {
    bool prefix = false;
    Size ordinal = 0;
    FragmentIonLikelihoodModel::parseIonName(names[i], prefix, ordinal);
    if (prefix && ordinal % 2 == 0) continue;
    spec.emplace_back(theo[i].getMZ(), prefix ? 10.0f : 1000.0f);
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

// Signal: the pattern above for each training peptide; noise: the reversed peptides against the same spectra.
static FragmentIonLikelihoodModel trained_model_(double pseudo_count = 20.0, int precursor_charge = 2)
{
  FragmentIonLikelihoodModel model(pseudo_count);
  for (const std::string& s : training_peptides_)
  {
    const AASequence peptide = AASequence::fromString(s);
    const MSSpectrum spec = observed_(peptide);
    const std::vector<Size> ranks = FragmentIonLikelihoodModel::intensityRanks(spec);
    model.addObservations(spec, ranks, theoretical_(peptide), peptide.size(), precursor_charge, 20.0, true, false);
    model.addObservations(spec, ranks, theoretical_(reversed_(peptide)), peptide.size(), precursor_charge, 20.0, true, true);
  }
  model.finalize();
  return model;
}

START_TEST(FragmentIonLikelihoodModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FragmentIonLikelihoodModel* ptr = nullptr;
FragmentIonLikelihoodModel* null_ptr = nullptr;

START_SECTION(FragmentIonLikelihoodModel())
{
  ptr = new FragmentIonLikelihoodModel();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_REAL_SIMILAR(ptr->pseudoCount(), 20.0)
  TEST_EQUAL(ptr->isTrained(), false)
  TEST_EQUAL(ptr->signalPsms(), 0)
  TEST_EQUAL(ptr->noisePsms(), 0)
}
END_SECTION

START_SECTION(~FragmentIonLikelihoodModel())
{
  delete ptr;
}
END_SECTION

START_SECTION(FragmentIonLikelihoodModel(double pseudo_count))
{
  FragmentIonLikelihoodModel model(5.0);
  TEST_REAL_SIMILAR(model.pseudoCount(), 5.0)
  TEST_EXCEPTION(Exception::InvalidParameter, FragmentIonLikelihoodModel(0.0))
  TEST_EXCEPTION(Exception::InvalidParameter, FragmentIonLikelihoodModel(-1.0))
  TEST_EXCEPTION(Exception::InvalidParameter, FragmentIonLikelihoodModel{std::numeric_limits<double>::quiet_NaN()})
  TEST_EXCEPTION(Exception::InvalidParameter, FragmentIonLikelihoodModel{std::numeric_limits<double>::infinity()})
}
END_SECTION

START_SECTION((static std::vector<Size> intensityRanks(const MSSpectrum& spectrum)))
{
  MSSpectrum spec;
  spec.emplace_back(100.0, 10.0f);
  spec.emplace_back(200.0, 30.0f);
  spec.emplace_back(300.0, 20.0f);
  spec.emplace_back(400.0, 30.0f);
  const std::vector<Size> ranks = FragmentIonLikelihoodModel::intensityRanks(spec);
  ABORT_IF(ranks.size() != 4)
  TEST_EQUAL(ranks[0], 4)
  TEST_EQUAL(ranks[1], 1) // ties keep peak order
  TEST_EQUAL(ranks[2], 3)
  TEST_EQUAL(ranks[3], 2)
  TEST_EQUAL(FragmentIonLikelihoodModel::intensityRanks(MSSpectrum()).empty(), true)
}
END_SECTION

START_SECTION((static Size rankOutcome(Size rank)))
{
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(1), 0)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(2), 0)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(3), 1)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(5), 1)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(6), 2)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(10), 2)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(11), 3)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(20), 3)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(21), 4)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(40), 4)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(41), 5)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(80), 5)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(81), 6)
  TEST_EQUAL(FragmentIonLikelihoodModel::rankOutcome(100000), 6)
  TEST_TRUE(FragmentIonLikelihoodModel::rankOutcome(100000) < FragmentIonLikelihoodModel::ABSENT)
  TEST_EQUAL(FragmentIonLikelihoodModel::ABSENT + 1, FragmentIonLikelihoodModel::OUTCOMES)
}
END_SECTION

START_SECTION((static Context contextOf(bool prefix, int precursor_charge, int fragment_charge, Size fragment_length, Size peptide_length)))
{
  using Model = FragmentIonLikelihoodModel;
  TEST_EQUAL(Model::contextOf(true, 2, 1, 3, 10).series, 0)
  TEST_EQUAL(Model::contextOf(false, 2, 1, 3, 10).series, 1)
  TEST_EQUAL(Model::contextOf(true, 1, 1, 3, 10).precursor_bucket, 0)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 3, 10).precursor_bucket, 0)
  TEST_EQUAL(Model::contextOf(true, 3, 1, 3, 10).precursor_bucket, 1)
  TEST_EQUAL(Model::contextOf(true, 4, 1, 3, 10).precursor_bucket, 2)
  TEST_EQUAL(Model::contextOf(true, 7, 1, 3, 10).precursor_bucket, 2)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 3, 10).fragment_charge, 0)
  TEST_EQUAL(Model::contextOf(true, 2, 2, 3, 10).fragment_charge, 1)
  TEST_EQUAL(Model::contextOf(true, 2, 3, 3, 10).fragment_charge, 1)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 1, 10).position_bin, 1)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 5, 10).position_bin, 5)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 9, 10).position_bin, 9)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 10, 10).position_bin, Model::POSITION_BINS - 1) // capped
  TEST_EQUAL(Model::contextOf(true, 2, 1, 1, 30).position_bin, 0)
  TEST_EQUAL(Model::contextOf(true, 2, 1, 1, 0).position_bin, 0)
}
END_SECTION

START_SECTION((static bool parseIonName(const std::string& name, bool& prefix, Size& ordinal)))
{
  bool prefix = false;
  Size ordinal = 0;
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("b5+", prefix, ordinal), true)
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 5)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("y3++", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 3)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("z.4+", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 4)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("z'12+", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 12)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("c7+", prefix, ordinal), true)
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 7)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("a1+", prefix, ordinal), true)
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 1)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("x2+", prefix, ordinal), true)
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 2)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("b3-H2O1+", prefix, ordinal), true) // neutral loss of b3
  TEST_EQUAL(prefix, true)
  TEST_EQUAL(ordinal, 3)
  TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName("PEPTIDE$y6+", prefix, ordinal), true) // cross-link annotation
  TEST_EQUAL(prefix, false)
  TEST_EQUAL(ordinal, 6)
  for (const std::string& name : {"", "b", "b+", "b0+", "[M+H]+", "iP+", "y-H2O+", "Q", "$"})
  {
    TEST_EQUAL(FragmentIonLikelihoodModel::parseIonName(name, prefix, ordinal), false)
  }
}
END_SECTION

START_SECTION((void addObservations(const MSSpectrum& spectrum, const std::vector<Size>& ranks, const MSSpectrum& theoretical, Size peptide_length, int precursor_charge, double tolerance, bool ppm, bool noise)))
{
  const AASequence peptide = AASequence::fromString("DFPIANGER");
  const MSSpectrum spec = observed_(peptide);
  const MSSpectrum theo = theoretical_(peptide);
  const std::vector<Size> ranks = FragmentIonLikelihoodModel::intensityRanks(spec);

  FragmentIonLikelihoodModel model;
  model.addObservations(spec, ranks, theo, peptide.size(), 2, 20.0, true, false);
  TEST_EQUAL(model.signalPsms(), 1)
  TEST_EQUAL(model.noisePsms(), 0)
  TEST_EQUAL(model.isTrained(), false)
  model.addObservations(spec, ranks, theoretical_(reversed_(peptide)), peptide.size(), 2, 20.0, true, true);
  TEST_EQUAL(model.signalPsms(), 1)
  TEST_EQUAL(model.noisePsms(), 1)
  model.finalize();
  TEST_EQUAL(model.isTrained(), true)
  // later observations need another finalize()
  model.addObservations(spec, ranks, theo, peptide.size(), 2, 20.0, true, false);
  TEST_EQUAL(model.isTrained(), false)
  TEST_EQUAL(model.signalPsms(), 2)

  // invalid input
  MSSpectrum unnamed = theo;
  unnamed.getStringDataArrays().clear();
  TEST_EXCEPTION(Exception::InvalidValue, model.addObservations(spec, ranks, unnamed, peptide.size(), 2, 20.0, true, false))
  TEST_EXCEPTION(Exception::InvalidParameter, model.addObservations(spec, ranks, theo, peptide.size(), 2, 0.0, true, false))
  TEST_EXCEPTION(Exception::InvalidParameter, model.addObservations(spec, ranks, theo, peptide.size(), 2, -1.0, false, false))
  TEST_EXCEPTION(Exception::InvalidParameter, model.addObservations(spec, ranks, theo, peptide.size(), 2, std::numeric_limits<double>::infinity(), true, false))
  TEST_EXCEPTION(Exception::InvalidValue, model.addObservations(spec, std::vector<Size>(ranks.size() + 1, 1), theo, peptide.size(), 2, 20.0, true, false))
  TEST_EQUAL(model.signalPsms(), 2)
  // an empty spectrum leaves every ion absent
  model.addObservations(MSSpectrum(), std::vector<Size>(), theo, peptide.size(), 2, 20.0, true, true);
  TEST_EQUAL(model.noisePsms(), 2)
}
END_SECTION

START_SECTION((void finalize()))
{
  // Without observations every outcome is equally likely under signal and noise: the flat prior.
  FragmentIonLikelihoodModel model;
  model.finalize();
  TEST_EQUAL(model.isTrained(), true)
  const FragmentIonLikelihoodModel::Context context = FragmentIonLikelihoodModel::contextOf(false, 2, 1, 3, 10);
  for (Size outcome = 0; outcome < FragmentIonLikelihoodModel::OUTCOMES; ++outcome)
  {
    TEST_REAL_SIMILAR(model.logLikelihoodRatio(context, outcome), 0.0)
  }
  TEST_REAL_SIMILAR(model.presenceProbability(context), 1.0 - 1.0 / static_cast<double>(FragmentIonLikelihoodModel::OUTCOMES))
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
  NOT_TESTABLE // see the constructors
}
END_SECTION

START_SECTION((double logLikelihoodRatio(const Context& context, Size outcome) const))
{
  const FragmentIonLikelihoodModel::Context suffix = FragmentIonLikelihoodModel::contextOf(false, 2, 1, 4, 9);
  const FragmentIonLikelihoodModel::Context prefix = FragmentIonLikelihoodModel::contextOf(true, 2, 1, 4, 9);
  FragmentIonLikelihoodModel untrained;
  TEST_EXCEPTION(Exception::Precondition, untrained.logLikelihoodRatio(suffix, 0))

  const FragmentIonLikelihoodModel model = trained_model_();
  // y ions are always present and intense in the training spectra, rarely in their reversed sequences
  TEST_TRUE(model.logLikelihoodRatio(suffix, 0) > 0.0)
  TEST_TRUE(model.logLikelihoodRatio(suffix, FragmentIonLikelihoodModel::ABSENT) < 0.0)
  // b ions are weak when present and absent half of the time: a missing b ion says less than a missing y ion
  TEST_TRUE(model.logLikelihoodRatio(prefix, FragmentIonLikelihoodModel::ABSENT) > model.logLikelihoodRatio(suffix, FragmentIonLikelihoodModel::ABSENT))
  // outcomes beyond the table count as absent
  TEST_REAL_SIMILAR(model.logLikelihoodRatio(suffix, 1000), model.logLikelihoodRatio(suffix, FragmentIonLikelihoodModel::ABSENT))
  // unseen contexts back off to their series and stay finite
  const FragmentIonLikelihoodModel::Context unseen = FragmentIonLikelihoodModel::contextOf(false, 5, 3, 4, 9);
  for (Size outcome = 0; outcome < FragmentIonLikelihoodModel::OUTCOMES; ++outcome)
  {
    TEST_TRUE(std::isfinite(model.logLikelihoodRatio(unseen, outcome)))
  }
  TEST_TRUE(model.logLikelihoodRatio(unseen, 0) > 0.0)
  TEST_TRUE(model.logLikelihoodRatio(unseen, FragmentIonLikelihoodModel::ABSENT) < 0.0)
}
END_SECTION

START_SECTION((double presenceProbability(const Context& context) const))
{
  const FragmentIonLikelihoodModel::Context suffix = FragmentIonLikelihoodModel::contextOf(false, 2, 1, 4, 9);
  const FragmentIonLikelihoodModel::Context prefix = FragmentIonLikelihoodModel::contextOf(true, 2, 1, 4, 9);
  FragmentIonLikelihoodModel untrained;
  TEST_EXCEPTION(Exception::Precondition, untrained.presenceProbability(suffix))

  const FragmentIonLikelihoodModel model = trained_model_();
  TEST_TRUE(model.presenceProbability(suffix) > 0.9)
  TEST_TRUE(model.presenceProbability(suffix) <= 1.0)
  TEST_TRUE(model.presenceProbability(prefix) > 0.1)
  TEST_TRUE(model.presenceProbability(prefix) < 0.9)
  TEST_TRUE(model.presenceProbability(suffix) > model.presenceProbability(prefix))
  // a smaller pseudo-count follows the counts more closely
  const FragmentIonLikelihoodModel sharp = trained_model_(1.0);
  TEST_TRUE(sharp.presenceProbability(suffix) > model.presenceProbability(suffix))
}
END_SECTION

START_SECTION((Features score(const MSSpectrum& spectrum, const std::vector<Size>& ranks, const MSSpectrum& theoretical, Size peptide_length, int precursor_charge, double tolerance, bool ppm, Size top_k = 6) const))
{
  const AASequence peptide = AASequence::fromString("DFPIANGER");
  const MSSpectrum spec = observed_(peptide);
  const MSSpectrum theo = theoretical_(peptide);
  const std::vector<Size> ranks = FragmentIonLikelihoodModel::intensityRanks(spec);

  FragmentIonLikelihoodModel untrained;
  TEST_EXCEPTION(Exception::Precondition, untrained.score(spec, ranks, theo, peptide.size(), 2, 20.0, true))

  const FragmentIonLikelihoodModel model = trained_model_();
  const FragmentIonLikelihoodModel::Features features = model.score(spec, ranks, theo, peptide.size(), 2, 20.0, true);
  TEST_EQUAL(features.theoretical_ions, 16) // b1-b8, y1-y8
  TEST_EQUAL(features.matched_ions, 12)     // every y ion and b1, b3, b5, b7
  TEST_TRUE(features.log_likelihood_ratio > 0.0)
  TEST_TRUE(features.explained_presence > 0.6)
  TEST_TRUE(features.explained_presence <= 1.0)
  TEST_REAL_SIMILAR(features.top_predicted_observed, 1.0) // the six most likely ions are y ions, all present

  // The reversed sequence explains the same spectrum badly.
  const FragmentIonLikelihoodModel::Features noise = model.score(spec, ranks, theoretical_(reversed_(peptide)), peptide.size(), 2, 20.0, true);
  TEST_EQUAL(noise.theoretical_ions, 16)
  TEST_TRUE(noise.matched_ions < features.matched_ions)
  TEST_TRUE(noise.log_likelihood_ratio < 0.0)
  TEST_TRUE(noise.explained_presence < features.explained_presence)
  TEST_TRUE(noise.top_predicted_observed < features.top_predicted_observed)

  // A peptide the model never saw, at a precursor charge it never saw, is scored by the run's pattern.
  const AASequence unseen = AASequence::fromString("GLSDGEWQQVLNVWGK");
  const MSSpectrum unseen_spec = observed_(unseen);
  const std::vector<Size> unseen_ranks = FragmentIonLikelihoodModel::intensityRanks(unseen_spec);
  const FragmentIonLikelihoodModel::Features unseen_features = model.score(unseen_spec, unseen_ranks, theoretical_(unseen), unseen.size(), 3, 20.0, true);
  const FragmentIonLikelihoodModel::Features unseen_noise = model.score(unseen_spec, unseen_ranks, theoretical_(reversed_(unseen)), unseen.size(), 3, 20.0, true);
  TEST_TRUE(unseen_features.log_likelihood_ratio > 0.0)
  TEST_TRUE(unseen_features.log_likelihood_ratio > unseen_noise.log_likelihood_ratio)

  // Da and ppm tolerances match the same peaks here.
  const FragmentIonLikelihoodModel::Features dalton = model.score(spec, ranks, theo, peptide.size(), 2, 0.02, false);
  TEST_REAL_SIMILAR(dalton.log_likelihood_ratio, features.log_likelihood_ratio)
  TEST_EQUAL(dalton.matched_ions, features.matched_ions)

  // top_k bounds the third feature; zero disables it
  TEST_REAL_SIMILAR(model.score(spec, ranks, theo, peptide.size(), 2, 20.0, true, 0).top_predicted_observed, 0.0)
  TEST_REAL_SIMILAR(model.score(spec, ranks, theo, peptide.size(), 2, 20.0, true, 100).top_predicted_observed, 12.0 / 16.0)

  // nothing to score
  MSSpectrum empty_theo;
  empty_theo.getStringDataArrays().emplace_back();
  const FragmentIonLikelihoodModel::Features none = model.score(spec, ranks, empty_theo, peptide.size(), 2, 20.0, true);
  TEST_EQUAL(none.theoretical_ions, 0)
  TEST_EQUAL(none.matched_ions, 0)
  TEST_REAL_SIMILAR(none.log_likelihood_ratio, 0.0)
  TEST_REAL_SIMILAR(none.explained_presence, 0.0)
  TEST_REAL_SIMILAR(none.top_predicted_observed, 0.0)

  // an empty spectrum leaves every ion absent
  const FragmentIonLikelihoodModel::Features absent = model.score(MSSpectrum(), std::vector<Size>(), theo, peptide.size(), 2, 20.0, true);
  TEST_EQUAL(absent.matched_ions, 0)
  TEST_TRUE(absent.log_likelihood_ratio < noise.log_likelihood_ratio)
  TEST_REAL_SIMILAR(absent.explained_presence, 0.0)
  TEST_REAL_SIMILAR(absent.top_predicted_observed, 0.0)

  // invalid input
  TEST_EXCEPTION(Exception::InvalidValue, model.score(spec, ranks, MSSpectrum(), peptide.size(), 2, 20.0, true))
  TEST_EXCEPTION(Exception::InvalidParameter, model.score(spec, ranks, theo, peptide.size(), 2, 0.0, true))
  TEST_EXCEPTION(Exception::InvalidValue, model.score(spec, std::vector<Size>(), theo, peptide.size(), 2, 20.0, true))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
