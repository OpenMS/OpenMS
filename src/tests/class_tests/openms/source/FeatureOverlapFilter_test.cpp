// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FILTERING/DATAREDUCTION/FeatureOverlapFilter.h>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(FeatureOverlapFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION((Filter FeatureMap))
  //feature with convex hulls
  Feature feature1;
  feature1.getPosition()[0] = 5.25;
  feature1.getPosition()[1] = 1.5;
  feature1.setIntensity(0.5f);
  feature1.setOverallQuality(8);
  std::vector< ConvexHull2D > hulls(1);
  hulls[0].addPoint(DPosition<2>(-1.0,2.0));
  hulls[0].addPoint(DPosition<2>(4.0,1.2));
  hulls[0].addPoint(DPosition<2>(5.0,3.123));
  feature1.setConvexHulls(hulls);

  Feature feature2;
  feature2.getPosition()[0] = 5.25;
  feature2.getPosition()[1] = 1.5;
  feature2.setIntensity(0.5f);
  feature2.setOverallQuality(10);
  std::vector< ConvexHull2D > hulls2(1);
  hulls2[0].addPoint(DPosition<2>(-1.0,2.0));
  hulls2[0].addPoint(DPosition<2>(4.0,1.2));
  hulls2[0].addPoint(DPosition<2>(5.5,3.123));
  feature2.setConvexHulls(hulls2);

  Feature feature3;
  feature3.getPosition()[0] = 5.25;
  feature3.getPosition()[1] = 1.5;
  feature3.setIntensity(0.5f);
  feature3.setOverallQuality(7);
  std::vector< ConvexHull2D > hulls3(1);
  hulls3[0].addPoint(DPosition<2>(4.5,2.0));
  hulls3[0].addPoint(DPosition<2>(10,1.2));
  hulls3[0].addPoint(DPosition<2>(10,3.123));
  feature3.setConvexHulls(hulls3);

  Feature feature4;
  feature4.getPosition()[0] = 20.;
  feature4.getPosition()[1] = 10.;
  feature4.setIntensity(0.5f);
  feature4.setOverallQuality(7);
  std::vector< ConvexHull2D > hulls4(1);
  hulls4[0].addPoint(DPosition<2>(20,5));
  hulls4[0].addPoint(DPosition<2>(22,10));
  hulls4[0].addPoint(DPosition<2>(22,14));
  feature4.setConvexHulls(hulls4);

  Feature feature5;
  feature5.getPosition()[0] = 20.;
  feature5.getPosition()[1] = 11.;
  feature5.setIntensity(0.5f);
  feature5.setOverallQuality(0.);
  std::vector< ConvexHull2D > hulls5(1);
  hulls5[0].addPoint(DPosition<2>(20,12.));
  hulls5[0].addPoint(DPosition<2>(21,16.));
  hulls5[0].addPoint(DPosition<2>(21,18.));
  feature5.setConvexHulls(hulls5);

  FeatureMap fmap;
  fmap.emplace_back(feature1);
  fmap.emplace_back(feature2);
  fmap.emplace_back(feature3);
  fmap.emplace_back(feature4);
  fmap.emplace_back(feature5);

  fmap.updateRanges();
  for (auto& f : fmap)
  {
    f.ensureUniqueId();
  }
  
  FeatureOverlapFilter::filter(fmap, 
    [](const Feature& left, const Feature& right){ return left.getOverallQuality() > right.getOverallQuality(); }, 
    [](const Feature&, const Feature&) { return true; },
    false);

  TEST_EQUAL(fmap[0].getOverallQuality(), 10)
  TEST_EQUAL(fmap[1].getOverallQuality(), 7)

END_SECTION

END_TEST
