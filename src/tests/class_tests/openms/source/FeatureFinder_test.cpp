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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
///////////////////////////


#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinder* ptr = nullptr;
FeatureFinder* nullPointer = nullptr;
START_SECTION((FeatureFinder()))
	ptr = new FeatureFinder();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureFinder()))
	delete ptr;
END_SECTION

START_SECTION((template <class PeakType, class FeatureType> void run(const String &algorithm_name, MSExperiment< PeakType > const &input_map, FeatureMap< FeatureType > &features, const Param &param, const FeatureMap<FeatureType>& seeds)))
	FeatureFinder ff;
	FeatureMap features;

	//empty map works -> nothing to do
	PeakMap map;
	ff.run("none", map, features, Param(), FeatureMap());

	//no updateRanges -> exception
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	TEST_EXCEPTION(Exception::IllegalArgument, ff.run("none", map, features, Param(), FeatureMap()))

	//updateRanges -> it works again
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap());

	//MS2 scans -> exception
	map[0].setMSLevel(1);
	map[0].setMSLevel(2);
	map.updateRanges();
	TEST_EXCEPTION(Exception::IllegalArgument, ff.run("none", map, features, Param(), FeatureMap()))
END_SECTION

START_SECTION((const Flag& getPeakFlag(const IndexPair& index) const))
	FeatureFinder ff;
	FeatureMap features;
	PeakMap map;
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap());
	TEST_EQUAL(ff.getPeakFlag(make_pair(0,0)),FeatureFinderDefs::UNUSED)
	TEST_EQUAL(ff.getPeakFlag(make_pair(1,0)),FeatureFinderDefs::UNUSED)
END_SECTION

START_SECTION((Flag& getPeakFlag(const IndexPair& index)))
	FeatureFinder ff;
	FeatureMap features;
	PeakMap map;
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap());
	ff.getPeakFlag(make_pair(0,0)) = FeatureFinderDefs::USED;
	TEST_EQUAL(ff.getPeakFlag(make_pair(0,0)),FeatureFinderDefs::USED)
	TEST_EQUAL(ff.getPeakFlag(make_pair(1,0)),FeatureFinderDefs::UNUSED)
END_SECTION

START_SECTION((Param getParameters(const String& algorithm_name) const))
	FeatureFinder ff;
	TEST_EQUAL(ff.getParameters("none")==Param(),true)
	TEST_EQUAL(ff.getParameters("centroided")==Param(),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



