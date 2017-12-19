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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/DigestSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DigestSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DigestSimulation* ptr = nullptr;
DigestSimulation* nullPointer = nullptr;
START_SECTION(DigestSimulation())
{
  ptr = new DigestSimulation();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DigestSimulation())
{
  delete ptr;
}
END_SECTION

START_SECTION((DigestSimulation(const DigestSimulation &source)))
{
  DigestSimulation a;
  Param p = a.getParameters();
  p.setValue("enzyme","no cleavage");
  a.setParameters(p);
  DigestSimulation b(a);

  TEST_EQUAL(b.getParameters(),a.getParameters());
}
END_SECTION

START_SECTION((DigestSimulation& operator=(const DigestSimulation &source)))
{
  DigestSimulation a,b;
  Param p = a.getParameters();
  p.setValue("enzyme","no cleavage");
  a.setParameters(p);

  TEST_NOT_EQUAL(b.getParameters(),a.getParameters());
  b = a;
  TEST_EQUAL(b.getParameters(),a.getParameters());
}
END_SECTION


START_SECTION((void digest(SimTypes::FeatureMapSim & feature_map)))
{
  SimTypes::FeatureMapSim fm;
  ProteinIdentification protIdent;
  // add new ProteinHit to ProteinIdentification
  {
  ProteinHit protHit(0.0, 1, "Hit1", "ACDKDDLDDFRLNN");
  protHit.setMetaValue("description", "desc 1");
  protHit.setMetaValue("intensity", 100);
  protIdent.insertHit(protHit);
  }
  {
  ProteinHit protHit(0.0, 1, "Hit2", "ACDKDDLASSRL");
  protHit.setMetaValue("description", "desc 1");
  protHit.setMetaValue("intensity", 50);
  protIdent.insertHit(protHit);
  }
  std::vector<ProteinIdentification> vec_protIdent;
  vec_protIdent.push_back(protIdent);
  fm.setProteinIdentifications(vec_protIdent);

  DigestSimulation a;
  Param p;
  p.setValue("model", "naive");
  a.setParameters(p);
  a.digest(fm);

  TEST_EQUAL(fm.size(), 8)
  ABORT_IF(fm.size() != 8)

  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("LNN"))
  TEST_EQUAL(fm[0].getIntensity(), 72)

  TEST_EQUAL(fm[1].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("ACDK"))
  TEST_EQUAL(fm[1].getIntensity(), 108)

  TEST_EQUAL(fm[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("DDLASSR"))
  TEST_EQUAL(fm[2].getIntensity(), 36)

  TEST_EQUAL(fm[3].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("DDLDDFR"))
  TEST_EQUAL(fm[3].getIntensity(), 72)

  TEST_EQUAL(fm[4].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("DDLASSRL"))
  TEST_EQUAL(fm[4].getIntensity(), 36)

  TEST_EQUAL(fm[5].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("DDLDDFRLNN"))
  TEST_EQUAL(fm[5].getIntensity(), 72)

  TEST_EQUAL(fm[6].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("ACDKDDLASSR"))
  TEST_EQUAL(fm[6].getIntensity(), 36)

  TEST_EQUAL(fm[7].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("ACDKDDLDDFR"))
  TEST_EQUAL(fm[7].getIntensity(), 72)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

