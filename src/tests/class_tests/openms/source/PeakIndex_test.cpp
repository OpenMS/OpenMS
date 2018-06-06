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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/PeakIndex.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakIndex, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakIndex* ptr = nullptr;
PeakIndex* nullPointer = nullptr;
START_SECTION((PeakIndex()))
{
	ptr = new PeakIndex();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~PeakIndex()))
{
	delete ptr;
}
END_SECTION

START_SECTION((PeakIndex(Size peak)))
{
  PeakIndex i(17);
	TEST_EQUAL(i.peak,17)
}
END_SECTION

START_SECTION((PeakIndex(Size spectrum, Size peak)))
{
  PeakIndex i(5,17);
	TEST_EQUAL(i.peak,17)
  TEST_EQUAL(i.spectrum,5)
}
END_SECTION

START_SECTION((bool isValid() const))
{
  PeakIndex i;
	TEST_EQUAL(i.isValid(),false)
  i.peak = 5;
	i.spectrum = 17;
	TEST_EQUAL(i.isValid(),true)
}
END_SECTION

START_SECTION((void clear()))
{
  PeakIndex i(5,17);
	TEST_EQUAL(i.isValid(),true)
	i.clear();
	TEST_EQUAL(i.isValid(),false)
	TEST_NOT_EQUAL(i.peak,17)
	TEST_NOT_EQUAL(i.spectrum,5)
}
END_SECTION

START_SECTION((bool operator==(const PeakIndex &rhs) const))
{
  PeakIndex i1, i2;
	TEST_EQUAL(i1==i2, true)
	i1.peak = 1;
	TEST_EQUAL(i1==i2, false)
	i2.peak = 1;
	TEST_EQUAL(i1==i2, true)
	i1.spectrum = 2;
	TEST_EQUAL(i1==i2, false)
	i2.spectrum = 2;
	TEST_EQUAL(i1==i2, true)
}
END_SECTION

START_SECTION((bool operator!=(const PeakIndex &rhs) const))
{
  PeakIndex i1, i2;
	TEST_EQUAL(i1!=i2, false)
	i1.peak = 1;
	TEST_EQUAL(i1!=i2, true)
	i2.peak = 1;
	TEST_EQUAL(i1!=i2, false)
	i1.spectrum = 2;
	TEST_EQUAL(i1!=i2, true)
	i2.spectrum = 2;
	TEST_EQUAL(i1!=i2, false)
}
END_SECTION

FeatureMap map;
map.resize(5);
map[0].setMZ(1);
map[1].setMZ(2);
map[2].setMZ(3);
map[3].setMZ(4);
map[4].setMZ(5);

ConsensusMap c_map;
c_map.resize(5);
c_map[0].setMZ(1.1);
c_map[1].setMZ(2.1);
c_map[2].setMZ(3.1);
c_map[3].setMZ(4.1);
c_map[4].setMZ(5.1);

START_SECTION((template <typename FeatureMapType> const FeatureMapType::value_type& getFeature(const FeatureMapType &map) const ))
{
  PeakIndex i;
	TEST_PRECONDITION_VIOLATED(i.getFeature(map))
	i.peak = 4;
	TEST_REAL_SIMILAR(i.getFeature(map).getMZ(),5.0)
	TEST_REAL_SIMILAR(i.getFeature(c_map).getMZ(),5.1)
	i.peak = 0;
	TEST_REAL_SIMILAR(i.getFeature(map).getMZ(),1.0)
	TEST_REAL_SIMILAR(i.getFeature(c_map).getMZ(),1.1)
	i.peak = 5;
	TEST_PRECONDITION_VIOLATED(i.getFeature(map))
}
END_SECTION

PeakMap exp;
exp.resize(3);
exp[0].setRT(1);
exp[0].resize(15);
exp[1].setRT(2);
exp[2].setRT(3);
exp[2].resize(3);
exp[2][0].setMZ(1.0);
exp[2][1].setMZ(2.0);
exp[2][2].setMZ(3.0);


START_SECTION((template <typename PeakMapType> const PeakMapType::SpectrumType& getSpectrum(const PeakMapType &map) const ))
{
  PeakIndex i;
	TEST_PRECONDITION_VIOLATED(i.getSpectrum(exp))
	i.spectrum = 0;
	TEST_REAL_SIMILAR(i.getSpectrum(exp).getRT(),1.0)
	i.spectrum = 2;
	TEST_REAL_SIMILAR(i.getSpectrum(exp).getRT(),3.0)
	i.spectrum = 3;
	TEST_PRECONDITION_VIOLATED(i.getSpectrum(exp))
}
END_SECTION

START_SECTION((template <typename PeakMapType> const PeakMapType::PeakType& getPeak(const PeakMapType &map) const ))
{
  PeakIndex i;
	TEST_PRECONDITION_VIOLATED(i.getPeak(exp))
	i.peak = 0;
	i.spectrum = 0;
	TEST_REAL_SIMILAR(i.getPeak(exp).getMZ(),0.0)
	i.peak = 0;
	i.spectrum = 2;
	TEST_REAL_SIMILAR(i.getPeak(exp).getMZ(),1.0)
	i.peak = 2;
	TEST_REAL_SIMILAR(i.getPeak(exp).getMZ(),3.0)
	i.peak = 16;
	TEST_PRECONDITION_VIOLATED(i.getPeak(exp))
	i.peak = 0;
	i.spectrum = 3;
	TEST_PRECONDITION_VIOLATED(i.getPeak(exp))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



