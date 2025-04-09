// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <boost/assign/std/vector.hpp>

#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMMapping, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMMapping* ptr = nullptr;
MRMMapping* nullPointer = nullptr;

START_SECTION(MRMMapping())
{
	ptr = new MRMMapping();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMMapping())
{
  delete ptr;
}
END_SECTION

START_SECTION(void mapExperiment(const OpenMS::PeakMap& input_chromatograms, const OpenMS::TargetedExperiment& targeted_exp, OpenMS::PeakMap& output)) 
{
  MRMMapping m;

  MSExperiment exp;
  exp.setComment("comment1");
  MSChromatogram c; 
  exp.addChromatogram(c);

  TEST_EQUAL(exp.getNrChromatograms(), 1)

  TargetedExperiment targ;
  ReactionMonitoringTransition t;
  targ.addTransition(t);
  MSExperiment out;

  m.mapExperiment(exp, targ, out);
  TEST_EQUAL(out.getNrChromatograms(), 0)

  {
    Param p = m.getDefaults();
    p.setValue("map_multiple_assays", "true");
    m.setParameters(p);

    m.mapExperiment(exp, targ, out);
    TEST_EQUAL(out.getNrChromatograms(), 1) // both transition and chromatogram have zero m/z
    TEST_EQUAL(out.getComment(), "comment1") // should preserve the meta data
  }

  exp.setComment("comment2");
  {
    Param p = m.getDefaults();
    p.setValue("map_multiple_assays", "true");
    p.setValue("precursor_tolerance", 9999.0);
    p.setValue("product_tolerance", 9999.0);
    m.setParameters(p);

    m.mapExperiment(exp, targ, out);
    TEST_EQUAL(out.getNrChromatograms(), 1)
    TEST_EQUAL(out.getComment(), "comment2") // should preserve the meta data
  }

  // Now set some precursor and fragment ion values, and check whether we can map one chromatogram to two transitions
  exp.getChromatograms()[0].getPrecursor().setMZ(500);
  exp.getChromatograms()[0].getProduct().setMZ(500);

  targ.addTransition(t);
  auto tr = targ.getTransitions();
  tr[0].setPrecursorMZ(500);
  tr[0].setProductMZ(500.1);
  tr[0].setNativeID("tr1");
  tr[1].setPrecursorMZ(500);
  tr[1].setProductMZ(500.0);
  tr[1].setNativeID("tr2");
  targ.setTransitions(tr);

  {
    Param p = m.getDefaults();
    p.setValue("map_multiple_assays", "true");
    p.setValue("precursor_tolerance", 1.0);
    p.setValue("product_tolerance", 1.0);
    m.setParameters(p);

    m.mapExperiment(exp, targ, out);
    TEST_EQUAL(exp.getNrChromatograms(), 1)
    TEST_EQUAL(out.getNrChromatograms(), 2)
    TEST_EQUAL(out.getChromatograms()[0].getNativeID(), "tr1")
    TEST_EQUAL(out.getChromatograms()[1].getNativeID(), "tr2")
  }

  // test that we cannot map when we don't allow multiple assays per chromatogram
  {
    Param p = m.getDefaults();
    p.setValue("map_multiple_assays", "false");
    p.setValue("precursor_tolerance", 1.0);
    p.setValue("product_tolerance", 1.0);
    m.setParameters(p);

    TEST_EXCEPTION(Exception::IllegalArgument, m.mapExperiment(exp, targ, out) )
  }

  // with a smaller mapping tolerance, we should only see a 1:1 mapping
  {
    Param p = m.getDefaults();
    p.setValue("map_multiple_assays", "true");
    p.setValue("precursor_tolerance", 0.05);
    p.setValue("product_tolerance", 0.05);
    m.setParameters(p);

    MSExperiment out2;
    m.mapExperiment(exp, targ, out2);
    TEST_EQUAL(exp.getNrChromatograms(), 1)
    TEST_EQUAL(out2.getNrChromatograms(), 1)
    TEST_EQUAL(out2.getChromatograms()[0].getNativeID(), "tr2")
  }

  // test error on unmapped chromatograms
  exp.getChromatograms()[0].getPrecursor().setMZ(600);
  exp.getChromatograms()[0].getProduct().setMZ(700);
  {
    Param p = m.getDefaults();
    p.setValue("map_multiple_assays", "true");
    p.setValue("precursor_tolerance", 1.0);
    p.setValue("product_tolerance", 1.0);
    m.setParameters(p);

    // that should still work
    MSExperiment out2;
    m.mapExperiment(exp, targ, out2);

    // not this
    p.setValue("error_on_unmapped", "true");
    m.setParameters(p);
    TEST_EXCEPTION(Exception::IllegalArgument, m.mapExperiment(exp, targ, out2) )
  }

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

