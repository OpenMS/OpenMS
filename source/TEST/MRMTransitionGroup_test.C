// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

///////////////////////////

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>

using namespace OpenMS;
using namespace std;

typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
typedef OpenMS::ReactionMonitoringTransition TransitionType;
typedef MRMTransitionGroup<RichPeakChromatogram, TransitionType> MRMTransitionGroupType;

///////////////////////////

START_TEST(MRMTransitionGroup, "$Id$")

/////////////////////////////////////////////////////////////

MRMTransitionGroupType* ptr = 0;
MRMTransitionGroupType* nullPointer = 0;

START_SECTION(MRMTransitionGroup())
{
	ptr = new MRMTransitionGroupType();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMTransitionGroup())
{
  delete ptr;
}
END_SECTION


RichPeakChromatogram chrom1;
RichPeakChromatogram chrom2;
TransitionType trans1;
TransitionType trans2;
MRMFeature feature1;
MRMFeature feature2;

START_SECTION(MRMTransitionGroup(const MRMTransitionGroup &rhs))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");

  MRMTransitionGroupType tmp(mrmtrgroup);
  TEST_EQUAL(mrmtrgroup.size(), tmp.size() )
}
END_SECTION

    
START_SECTION( MRMTransitionGroup& operator=(const MRMTransitionGroup &rhs) )
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");

  MRMTransitionGroupType tmp = mrmtrgroup;
  TEST_EQUAL(mrmtrgroup.size(), tmp.size() )
}
END_SECTION

START_SECTION (Size size() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.size(), 1)
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.size(), 2)
}
END_SECTION

START_SECTION (  const String & getTransitionGroupID() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.setTransitionGroupID("some_id");
  TEST_EQUAL(mrmtrgroup.getTransitionGroupID(), "some_id")
}
END_SECTION

START_SECTION (  void setTransitionGroupID(const String & tr_gr_id))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( std::vector<TransitionType>& getTransitionsMuteable())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  mrmtrgroup.addTransition(trans2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getTransitionsMuteable().size(), 2)
}
END_SECTION

START_SECTION ( void addTransition(const TransitionType &transition, String key))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( const TransitionType& getTransition(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  trans1.setLibraryIntensity(42);
  mrmtrgroup.addTransition(trans1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getTransition("dummy1").getLibraryIntensity(), 42)
}
END_SECTION

START_SECTION (  const std::vector<TransitionType>& getTransitions() const )
{
  MRMTransitionGroupType mrmtrgroup;
  trans1.setLibraryIntensity(42);
  mrmtrgroup.addTransition(trans1, "dummy1");
  trans2.setLibraryIntensity(-2);
  mrmtrgroup.addTransition(trans2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getTransitions()[0].getLibraryIntensity(), 42)
  TEST_EQUAL(mrmtrgroup.getTransitions()[1].getLibraryIntensity(), -2)
}
END_SECTION

START_SECTION (  bool hasTransition(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasTransition("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasTransition("dummy2"), false)
}
END_SECTION

START_SECTION (  const std::vector<SpectrumType>& getChromatograms() const ) 
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getChromatograms().size(), 2)
}
END_SECTION

START_SECTION (  std::vector<SpectrumType>& getChromatograms())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getChromatograms().size(), 2)
}
END_SECTION

START_SECTION (  void addChromatogram(SpectrumType &chromatogram, String key)) 
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION (  SpectrumType& getChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  chrom1.setMetaValue("some_value", 1);
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getChromatogram("dummy1").getMetaValue("some_value"), 1)
}
END_SECTION

START_SECTION (  bool hasChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasChromatogram("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasChromatogram("dummy2"), false)
}
END_SECTION

START_SECTION (  const std::vector<MRMFeature> & getFeatures() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addFeature(feature1);
  mrmtrgroup.addFeature(feature2);
  TEST_EQUAL(mrmtrgroup.getFeatures().size(), 2)
}
END_SECTION

START_SECTION (  std::vector<MRMFeature> & getFeaturesMuteable())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addFeature(feature1);
  mrmtrgroup.addFeature(feature2);
  TEST_EQUAL(mrmtrgroup.getFeaturesMuteable().size(), 2)
}
END_SECTION

START_SECTION (  void addFeature(MRMFeature & feature))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( void getLibraryIntensity(std::vector<double> & result) const)
{
  TransitionType new_trans1;
  TransitionType new_trans2;
  MRMTransitionGroupType mrmtrgroup;
  new_trans1.setLibraryIntensity(3);
  new_trans2.setLibraryIntensity(-2);
  mrmtrgroup.addTransition(new_trans1, "dummy1");
  mrmtrgroup.addTransition(new_trans2, "dummy2");
  std::vector< double > result;
  mrmtrgroup.getLibraryIntensity(result);
  TEST_EQUAL(result.size(), 2)
  TEST_REAL_SIMILAR(result[0], 3)
  TEST_REAL_SIMILAR(result[1], 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST

