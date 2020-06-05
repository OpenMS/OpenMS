// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

using namespace OpenMS;
using namespace std;

START_TEST(IDConflictResolverAlgorithm, "$Id$")

START_SECTION(resolveBetweenFeatures())
{
  FeatureMap map;
  Feature f1;
  Feature f2;
  Feature f3;
  Feature f4;

  PeptideHit hit;
  hit.setScore(23);
  hit.setSequence(AASequence::fromString("MORRISSEY"));
  PeptideIdentification id;
  id.insertHit(hit);
  std::vector<PeptideIdentification> ids;
  ids.push_back(id);
  
  PeptideHit hit2;
  hit2.setScore(23);
  hit2.setSequence(AASequence::fromString("M(Oxidation)ORRISSEY"));
  PeptideIdentification id2;
  id2.insertHit(hit2);
  std::vector<PeptideIdentification> ids2;
  ids2.push_back(id2);
  
  f1.setRT(1600.5);
  f1.setMZ(400.7);
  f1.setIntensity(1000.0);
  f1.setCharge(2);
  f1.setOverallQuality(1.0);
  f1.setPeptideIdentifications(ids);
  
  f2.setRT(1600.5);
  f2.setMZ(400.7);
  f2.setIntensity(10000.0);
  f2.setCharge(2);
  f2.setOverallQuality(1.0);
  f2.setPeptideIdentifications(ids);
  
  f3.setRT(1600.5);
  f3.setMZ(400.7);
  f3.setIntensity(1000.0);
  f3.setCharge(3);
  f3.setOverallQuality(1.0);
  f3.setPeptideIdentifications(ids);
  
  f4.setRT(1600.5);
  f4.setMZ(400.7);
  f4.setIntensity(1001.0);
  f4.setCharge(2);
  f4.setOverallQuality(1.0);
  f4.setPeptideIdentifications(ids2);
  
  map.push_back(f1);
  map.push_back(f2);
  
  IDConflictResolverAlgorithm::resolveBetweenFeatures(map);
  
  for (FeatureMap::ConstIterator it = map.begin(); it != map.end(); ++it)
  {
    
    if ((it->getIntensity() == 1000.0) && (it->getCharge() == 2))
    {
      // This identification was removed by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), true)
    }
      
    if ((it->getIntensity() == 10000.0) && (it->getCharge() == 2))
    {
      // This identification remains unchanged by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), false)
    }
    
    if ((it->getIntensity() == 1000.0) && (it->getCharge() == 3))
    {
      // This identification remains unchanged by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), false)
    }
    
    if ((it->getIntensity() == 1001.0) && (it->getCharge() == 2))
    {
      // This identification remains unchanged by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), false)
    }

  }
      
}
END_SECTION

END_TEST
