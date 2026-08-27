// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepModEncoder.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>

using namespace OpenMS;
using namespace OpenMS::ML;
using namespace std;

/// Element indices in AlphaPeptDeep's mod_elements list, used by the expectations below.
static const size_t IDX_C = 0;
static const size_t IDX_H = 1;
static const size_t IDX_N = 2;
static const size_t IDX_O = 3;
static const size_t IDX_P = 4;

/// mod_x value at (row, element) for a peptide block of `rows` rows.
static float at(const vector<float>& mod_x, size_t row, size_t element)
{
  return mod_x[row * PeptDeepModEncoder::MOD_FEATURE_SIZE + element];
}

START_TEST(PeptDeepModEncoder, "$Id$")

START_SECTION(static const std::vector<std::string>& elementOrder())
{
  const auto& els = PeptDeepModEncoder::elementOrder();
  // The order is the model's; any reordering silently changes every prediction.
  TEST_EQUAL(els.size(), PeptDeepModEncoder::MOD_FEATURE_SIZE)
  TEST_EQUAL(els[IDX_C], "C")
  TEST_EQUAL(els[IDX_H], "H")
  TEST_EQUAL(els[IDX_N], "N")
  TEST_EQUAL(els[IDX_O], "O")
  TEST_EQUAL(els[IDX_P], "P")
  TEST_EQUAL(els[5], "S")
  TEST_EQUAL(els.back(), "?")
}
END_SECTION

START_SECTION(static size_t elementIndex(const std::string& symbol))
{
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("C"), IDX_C)
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("O"), IDX_O)
  // OpenMS spells isotopes "(13)C"; PeptDeep spells them "13C". Both must land in one slot.
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("(13)C"), PeptDeepModEncoder::elementIndex("13C"))
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("(2)H"), PeptDeepModEncoder::elementIndex("2H"))
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("(15)N"), PeptDeepModEncoder::elementIndex("15N"))
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("(18)O"), PeptDeepModEncoder::elementIndex("18O"))
  // Unmodelled elements fold into "?" instead of being dropped.
  TEST_EQUAL(PeptDeepModEncoder::elementIndex("Uue"), PeptDeepModEncoder::unknownElementIndex())
  TEST_EQUAL(PeptDeepModEncoder::unknownElementIndex(), PeptDeepModEncoder::MOD_FEATURE_SIZE - 1)
}
END_SECTION

START_SECTION([EXTRA] encode() reproduces AlphaPeptDeep parse_mod_feature)
{
  // Expectations generated with peptdeep.model.featurize.parse_mod_feature; the
  // feature of a modification is the atom count vector of its diff formula.
  {
    // Carbamidomethyl (H3C2NO) on the C at residue 8 -> row 8.
    AASequence seq = AASequence::fromString("PEPTIDEC(Carbamidomethyl)K");
    const size_t rows = seq.size() + 2;
    vector<float> mod_x(rows * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
    PeptDeepModEncoder::encode(seq, mod_x.data(), rows);
    TEST_REAL_SIMILAR(at(mod_x, 8, IDX_C), 2.0)
    TEST_REAL_SIMILAR(at(mod_x, 8, IDX_H), 3.0)
    TEST_REAL_SIMILAR(at(mod_x, 8, IDX_N), 1.0)
    TEST_REAL_SIMILAR(at(mod_x, 8, IDX_O), 1.0)
    // Everything else stays zero, including the terminal rows.
    TEST_REAL_SIMILAR(at(mod_x, 0, IDX_C), 0.0)
    TEST_REAL_SIMILAR(at(mod_x, 7, IDX_C), 0.0)
    TEST_REAL_SIMILAR(at(mod_x, rows - 1, IDX_O), 0.0)
  }
  {
    // Oxidation (O) on the M at residue 4 -> row 4.
    AASequence seq = AASequence::fromString("PEPM(Oxidation)TIDEK");
    const size_t rows = seq.size() + 2;
    vector<float> mod_x(rows * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
    PeptDeepModEncoder::encode(seq, mod_x.data(), rows);
    TEST_REAL_SIMILAR(at(mod_x, 4, IDX_O), 1.0)
    TEST_REAL_SIMILAR(at(mod_x, 4, IDX_C), 0.0)
  }
  {
    // Phospho (HO3P) on the S at residue 1 -> row 1.
    AASequence seq = AASequence::fromString("S(Phospho)PEPTIDEK");
    const size_t rows = seq.size() + 2;
    vector<float> mod_x(rows * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
    PeptDeepModEncoder::encode(seq, mod_x.data(), rows);
    TEST_REAL_SIMILAR(at(mod_x, 1, IDX_H), 1.0)
    TEST_REAL_SIMILAR(at(mod_x, 1, IDX_O), 3.0)
    TEST_REAL_SIMILAR(at(mod_x, 1, IDX_P), 1.0)
  }
  {
    // Two modifications on one peptide occupy their own rows.
    AASequence seq = AASequence::fromString("PEPTIDEC(Carbamidomethyl)M(Oxidation)K");
    const size_t rows = seq.size() + 2;
    vector<float> mod_x(rows * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
    PeptDeepModEncoder::encode(seq, mod_x.data(), rows);
    TEST_REAL_SIMILAR(at(mod_x, 8, IDX_C), 2.0)
    TEST_REAL_SIMILAR(at(mod_x, 8, IDX_H), 3.0)
    TEST_REAL_SIMILAR(at(mod_x, 9, IDX_O), 1.0)
    TEST_REAL_SIMILAR(at(mod_x, 9, IDX_C), 0.0)
  }
  {
    // An N-terminal modification lands on row 0, not on residue 1.
    AASequence seq = AASequence::fromString(".(Acetyl)PEPTIDEK");
    const size_t rows = seq.size() + 2;
    vector<float> mod_x(rows * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
    PeptDeepModEncoder::encode(seq, mod_x.data(), rows);
    TEST_REAL_SIMILAR(at(mod_x, 0, IDX_C), 2.0)
    TEST_REAL_SIMILAR(at(mod_x, 0, IDX_H), 2.0)
    TEST_REAL_SIMILAR(at(mod_x, 0, IDX_O), 1.0)
    TEST_REAL_SIMILAR(at(mod_x, 1, IDX_C), 0.0)
  }
}
END_SECTION

START_SECTION([EXTRA] encode() leaves unmodified peptides and padding untouched)
{
  AASequence seq = AASequence::fromString("PEPTIDEK");
  const size_t rows = seq.size() + 6; // deliberately padded beyond length+2
  vector<float> mod_x(rows * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
  PeptDeepModEncoder::encode(seq, mod_x.data(), rows);
  double sum = 0.0;
  for (float v : mod_x) { sum += v; }
  TEST_REAL_SIMILAR(sum, 0.0)
}
END_SECTION

START_SECTION([EXTRA] encode() rejects an undersized mod_x block)
{
  AASequence seq = AASequence::fromString("PEPTIDEK");
  vector<float> mod_x(seq.size() * PeptDeepModEncoder::MOD_FEATURE_SIZE, 0.0f);
  TEST_EXCEPTION(Exception::InvalidValue, PeptDeepModEncoder::encode(seq, mod_x.data(), seq.size()))
}
END_SECTION

START_SECTION([EXTRA] buildModifiedPeptideBatch encodes modifications into the batch)
{
  vector<AASequence> peptides = {
    AASequence::fromString("PEPTIDEC(Carbamidomethyl)K"),
    AASequence::fromString("PEPTIDEKK")};
  PeptDeepInputBatch batch = PeptDeepInputBuilder::buildModifiedPeptideBatch(peptides);

  TEST_EQUAL(batch.batch_size, 2)
  TEST_EQUAL(batch.sequence_length, 11) // 9 residues + 2 terminal tokens
  TEST_EQUAL(batch.mod_x.size(), batch.batch_size * batch.sequence_length * PeptDeepModEncoder::MOD_FEATURE_SIZE)

  const size_t stride = batch.sequence_length * PeptDeepModEncoder::MOD_FEATURE_SIZE;
  // First peptide carries the modification ...
  TEST_REAL_SIMILAR(batch.mod_x[8 * PeptDeepModEncoder::MOD_FEATURE_SIZE + IDX_C], 2.0)
  // ... the second is unmodified, so its whole block stays zero.
  double sum_second = 0.0;
  for (size_t i = 0; i < stride; ++i) { sum_second += batch.mod_x[stride + i]; }
  TEST_REAL_SIMILAR(sum_second, 0.0)

  // The token stream must be the unmodified backbone in both cases.
  TEST_EQUAL(batch.aa_indices[0], 0)                 // leading terminal token
  TEST_EQUAL(batch.aa_indices[1], 'P' - 'A' + 1)
}
END_SECTION

END_TEST
