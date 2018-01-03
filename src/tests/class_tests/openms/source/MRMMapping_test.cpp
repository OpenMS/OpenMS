// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

