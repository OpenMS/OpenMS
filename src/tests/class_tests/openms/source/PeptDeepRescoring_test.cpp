// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeptDeepRescoring, "$Id$")

/////////////////////////////////////////////////////////////

PeptDeepRescoring* ptr = nullptr;
PeptDeepRescoring* null_ptr = nullptr;

START_SECTION((PeptDeepRescoring()))
  ptr = new PeptDeepRescoring();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION((~PeptDeepRescoring()))
  delete ptr;
END_SECTION

START_SECTION((static StringList featureNames()))
{
  StringList names = PeptDeepRescoring::featureNames();
  TEST_EQUAL(names.size(), 5)
  // The order is part of the contract: it is what gets appended to extra_features.
  TEST_STRING_EQUAL(names[0], Constants::UserParam::MS2_COSINE)
  TEST_STRING_EQUAL(names[1], Constants::UserParam::MS2_SPECTRAL_ANGLE)
  TEST_STRING_EQUAL(names[2], Constants::UserParam::MS2_PEARSON)
  TEST_STRING_EQUAL(names[3], Constants::UserParam::MS2_FRAC_PRED_FOUND)
  TEST_STRING_EQUAL(names[4], Constants::UserParam::RT_ABS_ERROR)
}
END_SECTION

START_SECTION([EXTRA] default parameters)
{
  PeptDeepRescoring r;
  Param p = r.getParameters();
  // Models have no sensible default, so they start empty and annotate() refuses to run.
  TEST_STRING_EQUAL(p.getValue("ms2_model").toString(), "")
  TEST_STRING_EQUAL(p.getValue("rt_model").toString(), "")
  // A negative NCE means "select it from the data".
  TEST_REAL_SIMILAR((double)p.getValue("nce"), -1.0)
  TEST_STRING_EQUAL(p.getValue("rt_model_type").toString(), "b_spline")
  TEST_STRING_EQUAL(p.getValue("instrument").toString(), "QE")
  // Nothing has been run yet.
  TEST_REAL_SIMILAR(r.getUsedNCE(), -1.0)
}
END_SECTION

START_SECTION([EXTRA] annotate() refuses to run without models)
{
  PeptDeepRescoring r;
  PeakMap exp;
  vector<ProteinIdentification> prot(1);
  PeptideIdentificationList peps;
  PeptideIdentification pi;
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  pi.setHits({hit});
  peps.push_back(pi);

  // Both model paths empty: the features cannot be computed at all.
  TEST_EXCEPTION(Exception::MissingInformation, r.annotate(exp, prot, peps))

  // A configured but absent model must fail loudly rather than silently emit zeros.
  Param p = r.getParameters();
  p.setValue("ms2_model", "does_not_exist_ms2.onnx");
  p.setValue("rt_model", "does_not_exist_rt.onnx");
  r.setParameters(p);
  TEST_EXCEPTION(Exception::FileNotFound, r.annotate(exp, prot, peps))
}
END_SECTION

START_SECTION([EXTRA] annotate() is a no-op on empty input)
{
  PeptDeepRescoring r;
  Param p = r.getParameters();
  p.setValue("ms2_model", OPENMS_GET_TEST_DATA_PATH("PeptDeepRescoring_test.txt"));
  p.setValue("rt_model", OPENMS_GET_TEST_DATA_PATH("PeptDeepRescoring_test.txt"));
  r.setParameters(p);

  PeakMap exp;
  vector<ProteinIdentification> prot(1);
  PeptideIdentificationList peps;
  // No PSMs: returns before touching the ONNX runtime, so the placeholder files above
  // are never actually loaded.
  r.annotate(exp, prot, peps);
  TEST_EQUAL(peps.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
