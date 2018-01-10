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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(FeaFiModule, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeaFiModule<Peak1D>* ptr = nullptr;
FeaFiModule<Peak1D>* nullPointer = nullptr;
START_SECTION((FeaFiModule(const MSExperiment<PeakType>* map, FeatureMap* features, FeatureFinder* ff)))
	ptr = new FeaFiModule<Peak1D>(nullptr,nullptr,nullptr);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeaFiModule()))
	delete ptr;
END_SECTION

//create dummy MSExperiment
PeakMap exp;
exp.resize(2);
exp[0].setMSLevel(1);
exp[0].setRT(1.1);
exp[1].setMSLevel(1);
exp[1].setRT(2.2);
//scan 1
Peak1D p;
p.setMZ(500.0);
p.setIntensity(501.0);
exp[0].push_back(p);
p.setMZ(700.0);
p.setIntensity(701.0);
exp[0].push_back(p);
p.setMZ(900.0);
p.setIntensity(901.0);
exp[0].push_back(p);
//scan 2
p.setMZ(600.0);
p.setIntensity(601.0);
exp[1].push_back(p);
p.setMZ(1000.0);
p.setIntensity(1001.0);
exp[1].push_back(p);

START_SECTION(IntensityType getPeakIntensity(const FeatureFinderDefs::IndexPair& index) const)
	FeaFiModule<Peak1D> t(&exp,nullptr,nullptr);
	TEST_REAL_SIMILAR(t.getPeakIntensity(make_pair(0,0)), 501.0)
	TEST_REAL_SIMILAR(t.getPeakIntensity(make_pair(0,1)), 701.0)
	TEST_REAL_SIMILAR(t.getPeakIntensity(make_pair(0,2)), 901.0)
	TEST_REAL_SIMILAR(t.getPeakIntensity(make_pair(1,0)), 601.0)
	TEST_REAL_SIMILAR(t.getPeakIntensity(make_pair(1,1)), 1001.0)
END_SECTION

START_SECTION(CoordinateType getPeakMz(const FeatureFinderDefs::IndexPair& index) const)
	FeaFiModule<Peak1D> t(&exp,nullptr,nullptr);
	TEST_REAL_SIMILAR(t.getPeakMz(make_pair(0,0)), 500.0)
	TEST_REAL_SIMILAR(t.getPeakMz(make_pair(0,1)), 700.0)
	TEST_REAL_SIMILAR(t.getPeakMz(make_pair(0,2)), 900.0)
	TEST_REAL_SIMILAR(t.getPeakMz(make_pair(1,0)), 600.0)
	TEST_REAL_SIMILAR(t.getPeakMz(make_pair(1,1)), 1000.0)
END_SECTION

START_SECTION(CoordinateType getPeakRt(const FeatureFinderDefs::IndexPair& index) const)
	FeaFiModule<Peak1D> t(&exp,nullptr,nullptr);
	TEST_REAL_SIMILAR(t.getPeakRt(make_pair(0,0)), 1.1)
	TEST_REAL_SIMILAR(t.getPeakRt(make_pair(0,1)), 1.1)
	TEST_REAL_SIMILAR(t.getPeakRt(make_pair(0,2)), 1.1)
	TEST_REAL_SIMILAR(t.getPeakRt(make_pair(1,0)), 2.2)
	TEST_REAL_SIMILAR(t.getPeakRt(make_pair(1,1)), 2.2)
END_SECTION

START_SECTION(void getNextMz(FeatureFinderDefs::IndexPair& index) const )
	FeaFiModule<Peak1D> t(&exp,nullptr,nullptr);
	//scan one
	FeatureFinderDefs::IndexPair i = make_pair(0,0);
	t.getNextMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,1)
	t.getNextMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,2)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextMz(i));

	//scan two
	i = make_pair(1,0);
	t.getNextMz(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,1)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextMz(i));

	//test for corrupt index
	i = make_pair(5,0);
	TEST_PRECONDITION_VIOLATED(t.getNextMz(i));
	i = make_pair(1,5);
	TEST_PRECONDITION_VIOLATED(t.getNextMz(i));
END_SECTION

START_SECTION(void getPrevMz(FeatureFinderDefs::IndexPair& index) const )
	FeaFiModule<Peak1D> t(&exp,nullptr,nullptr);
	//scan one
	FeatureFinderDefs::IndexPair i = make_pair(0,2);
	t.getPrevMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,1)
	t.getPrevMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getPrevMz(i));
	//scan two
	i = make_pair(1,1);
	t.getPrevMz(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getPrevMz(i));

	//test for corrupt index
	i = make_pair(5,0);
	TEST_PRECONDITION_VIOLATED(t.getPrevMz(i));
	i = make_pair(1,5);
	TEST_PRECONDITION_VIOLATED(t.getPrevMz(i));
END_SECTION

START_SECTION(void getNextRt(FeatureFinderDefs::IndexPair& index) )

	PeakMap exp2 = exp;
	exp2.resize(3);
	exp2[2].resize(1);
	exp2[2][0].setMZ(800.0);
	exp2[0].resize(5);
	exp2[0][2].setMZ(799.0);
	exp2[0][3].setMZ(801.0);
	exp2[0][4].setMZ(900.0);

	FeaFiModule<Peak1D> t(&exp2,nullptr,nullptr);

	FeatureFinderDefs::IndexPair i;

	std::cout << "peak one" << std::endl;
	i = make_pair(0,0);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextRt(i));

	std::cout << "peak two" << std::endl;
	i = make_pair(0,1);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextRt(i));

	std::cout << "peak three" << std::endl;
	i = make_pair(0,2);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextRt(i));

	std::cout << "peak four" << std::endl;
	i = make_pair(0,3);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,1)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextRt(i));

	std::cout << "peak five" << std::endl;
	i = make_pair(0,4);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,1)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getNextRt(i));

#ifdef OPENMS_DEBUG
	std::cout << "test for corrupt index" << std::endl;
	i = make_pair(5,0);
	TEST_PRECONDITION_VIOLATED(t.getNextRt(i));
	i = make_pair(1,5);
	TEST_PRECONDITION_VIOLATED(t.getNextRt(i));
#endif
END_SECTION

START_SECTION(void getPrevRt(FeatureFinderDefs::IndexPair& index) )
	PeakMap exp2 = exp;
	exp2[1].resize(4);
	exp2[1][0].setMZ(599.0);
	exp2[1][1].setMZ(799.0);
	exp2[1][2].setMZ(801.0);
	exp2[1][3].setMZ(1000.0);

	FeaFiModule<Peak1D> t(&exp2,nullptr,nullptr);
	FeatureFinderDefs::IndexPair i;

	//peak one
	i = make_pair(1,0);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getPrevRt(i));

	//peak two
	i = make_pair(1,1);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,1)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getPrevRt(i));

	//peak three
	i = make_pair(1,2);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,2)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getPrevRt(i));

	//peak four
	i = make_pair(1,3);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,2)
	TEST_EXCEPTION(FeatureFinderDefs::NoSuccessor, t.getPrevRt(i));

#ifdef OPENMS_DEBUG
	//test for corrupt index
	i = make_pair(5,0);
	TEST_PRECONDITION_VIOLATED(t.getPrevRt(i));
	i = make_pair(1,5);
	TEST_PRECONDITION_VIOLATED(t.getPrevRt(i));
#endif
END_SECTION

START_SECTION(void addConvexHull(const FeatureFinderDefs::IndexSet& set, Feature& feature) const)
	Peak2D p;
	std::vector<Peak2D> peak_array;
	p.getPosition()[0] = 1240.54;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1241.81;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 692.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 693.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;   p.getPosition()[1] = 695.2;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.07;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1247.09;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1248.41;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1249.76;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 689.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 692.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;   p.getPosition()[1] = 694.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1247.75;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1249.12;   p.getPosition()[1] = 688;       peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 689.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 691;       peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 692.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;   p.getPosition()[1] = 693.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;   p.getPosition()[1] = 695.2;     peak_array.push_back(p);
	p.getPosition()[0] = 1248.41;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1243.78;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1239.9;    p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1237.27;   p.getPosition()[1] = 692;       peak_array.push_back(p);
	p.getPosition()[0] = 1237.27;   p.getPosition()[1] = 691;       peak_array.push_back(p);
	p.getPosition()[0] = 1237.93;   p.getPosition()[1] = 688.4;     peak_array.push_back(p);

	std::sort(peak_array.begin(),peak_array.end(),Peak2D::PositionLess());
	PeakMap exp2;
	exp2.set2DData(peak_array);

	FeaFiModule<Peak1D> t(&exp2,nullptr,nullptr);

	FeatureFinderDefs::IndexSet set;
	for (Size i=0; i<exp2.size(); ++i)
	{
		for (Size j=0; j<exp2[i].size(); ++j)
		{
			set.insert(std::make_pair(i,j));
		}
	}

	Feature f;
	t.addConvexHull(set,f);
	ConvexHull2D& hull = f.getConvexHulls()[0];
	ConvexHull2D::PointArrayType hullpoints = hull.getHullPoints();
	TEST_EQUAL(hullpoints.size(), 30);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1237.27, 691)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1237.93, 688.4)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1240.54, 687.6)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1252.39, 687.6)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1252.39, 693.8)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1251.73, 695.2)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1251.07, 695.4)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1239.9, 695.4)) != hullpoints.end(), true);
	TEST_EQUAL(find(hullpoints.begin(), hullpoints.end(), ConvexHull2D::PointType(1237.27, 692)) != hullpoints.end(), true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
