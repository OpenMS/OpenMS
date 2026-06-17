// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/CometNativeIDRemapper.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CometNativeIDRemapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((Size rewriteToIndex(MSExperiment & exp)))
{
  // synthetic Bruker DDA experiment: non-monotonic "frame=F scan=S" native IDs that trigger
  // Comet's mzParser sort-UB; rewriteToIndex must replace them with monotonic "index=N".
  MSExperiment exp;
  MSSpectrum s0;
  s0.setNativeID("frame=1 scan=100 precursor=1");
  s0.setMSLevel(2);
  MSSpectrum s1;
  s1.setNativeID("frame=1 scan=105 precursor=2");
  s1.setMSLevel(2);
  MSSpectrum s2;
  s2.setNativeID("frame=2 scan=200 precursor=1");
  s2.setMSLevel(2);
  exp.addSpectrum(s0);
  exp.addSpectrum(s1);
  exp.addSpectrum(s2);

  Size n = CometNativeIDRemapper::rewriteToIndex(exp);
  TEST_EQUAL(n, 3)

  // native IDs are now monotonic "index=N"
  TEST_EQUAL(exp[0].getNativeID(), "index=0")
  TEST_EQUAL(exp[1].getNativeID(), "index=1")
  TEST_EQUAL(exp[2].getNativeID(), "index=2")

  // the original native IDs are preserved under the ORIGINAL_NATIVE_ID MetaValue
  TEST_EQUAL(exp[0].metaValueExists(CometNativeIDRemapper::ORIGINAL_NATIVE_ID), true)
  TEST_EQUAL(exp[0].getMetaValue(CometNativeIDRemapper::ORIGINAL_NATIVE_ID).toString(), "frame=1 scan=100 precursor=1")
  TEST_EQUAL(exp[1].getMetaValue(CometNativeIDRemapper::ORIGINAL_NATIVE_ID).toString(), "frame=1 scan=105 precursor=2")
  TEST_EQUAL(exp[2].getMetaValue(CometNativeIDRemapper::ORIGINAL_NATIVE_ID).toString(), "frame=2 scan=200 precursor=1")

  // empty experiment -> nothing rewritten
  MSExperiment empty_exp;
  TEST_EQUAL(CometNativeIDRemapper::rewriteToIndex(empty_exp), 0)
}
END_SECTION

START_SECTION((void translateReferencesBack(const MSExperiment& exp, std::vector<PeptideIdentification>& pids)))
{
  // mixed "frame=" and plain "scan=" native IDs; rewrite, then translate PSM references back.
  MSExperiment exp;
  MSSpectrum s0;
  s0.setNativeID("frame=1 scan=100");
  MSSpectrum s1;
  s1.setNativeID("scan=4242"); // a plain scan= native ID
  MSSpectrum s2;
  s2.setNativeID("frame=2 scan=200");
  exp.addSpectrum(s0);
  exp.addSpectrum(s1);
  exp.addSpectrum(s2);
  CometNativeIDRemapper::rewriteToIndex(exp); // -> index=0/1/2, originals stashed

  // PSMs as Comet would emit them: referencing the rewritten "index=N" IDs
  PeptideHit h;
  h.setSequence(AASequence::fromString("PEPTIDE"));
  h.setScore(1.0);
  h.setCharge(2);
  PeptideIdentificationList pids;
  for (const std::string& ref : {std::string("index=0"), std::string("index=1"), std::string("index=2")})
  {
    PeptideIdentification p;
    p.setSpectrumReference(ref);
    p.insertHit(h);
    pids.push_back(p);
  }

  CometNativeIDRemapper::translateReferencesBack(exp, pids);

  // every PSM reference is restored to its ORIGINAL native ID
  TEST_EQUAL(pids[0].getSpectrumReference(), "frame=1 scan=100")
  TEST_EQUAL(pids[1].getSpectrumReference(), "scan=4242")
  TEST_EQUAL(pids[2].getSpectrumReference(), "frame=2 scan=200")

  // the hits themselves are untouched by the translation
  TEST_EQUAL(pids[0].getHits().size(), 1)
  TEST_EQUAL(pids[0].getHits()[0].getSequence().toString(), "PEPTIDE")
}
END_SECTION

START_SECTION([EXTRA] translateReferencesBack guard + unmatched references)
{
  // (a) empty experiment -> no-op, no crash
  {
    MSExperiment empty_exp;
    PeptideIdentification p;
    p.setSpectrumReference("index=0");
    PeptideIdentificationList pids;
    pids.push_back(p);
    CometNativeIDRemapper::translateReferencesBack(empty_exp, pids);
    TEST_EQUAL(pids[0].getSpectrumReference(), "index=0") // unchanged
  }

  // (b) experiment never rewritten (no ORIGINAL_NATIVE_ID) -> no-op
  {
    MSExperiment exp;
    MSSpectrum s;
    s.setNativeID("scan=1");
    exp.addSpectrum(s);
    PeptideIdentification p;
    p.setSpectrumReference("scan=1");
    PeptideIdentificationList pids;
    pids.push_back(p);
    CometNativeIDRemapper::translateReferencesBack(exp, pids);
    TEST_EQUAL(pids[0].getSpectrumReference(), "scan=1") // unchanged
  }

  // (c) a reference that matches no rewritten spectrum is left untouched
  {
    MSExperiment exp;
    MSSpectrum s;
    s.setNativeID("frame=1 scan=9");
    exp.addSpectrum(s);
    CometNativeIDRemapper::rewriteToIndex(exp); // -> index=0
    PeptideIdentification p;
    p.setSpectrumReference("index=999"); // no such spectrum
    PeptideIdentificationList pids;
    pids.push_back(p);
    CometNativeIDRemapper::translateReferencesBack(exp, pids);
    TEST_EQUAL(pids[0].getSpectrumReference(), "index=999") // unchanged
  }

  // (d) a later rewritten spectrum is still used even if the first spectrum has no ORIGINAL_NATIVE_ID
  {
    MSExperiment exp;
    MSSpectrum s0;
    s0.setNativeID("scan=not-rewritten");
    exp.addSpectrum(s0);
    MSSpectrum s1;
    s1.setNativeID("index=7");
    s1.setMetaValue(CometNativeIDRemapper::ORIGINAL_NATIVE_ID, "frame=7 scan=42");
    exp.addSpectrum(s1);
    PeptideIdentification p;
    p.setSpectrumReference("index=7");
    PeptideIdentificationList pids;
    pids.push_back(p);
    CometNativeIDRemapper::translateReferencesBack(exp, pids);
    TEST_EQUAL(pids[0].getSpectrumReference(), "frame=7 scan=42")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
