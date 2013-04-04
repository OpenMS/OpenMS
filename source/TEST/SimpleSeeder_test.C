// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Clemens Groepl$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
///////////////////////////

START_TEST(SimpleSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
SimpleSeeder<Peak1D,Feature>* ptr = 0;
SimpleSeeder<Peak1D,Feature>* nullPointer = 0;
START_SECTION(SimpleSeeder(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff))
	MSExperiment<Peak1D> exp;
	ptr = new SimpleSeeder<Peak1D,Feature>(&exp,0,0);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~SimpleSeeder())
	delete ptr;
END_SECTION


START_SECTION(IndexPair nextSeed())

	//create map
	MSExperiment<Peak1D> exp;
	MSExperiment<Peak1D>::PeakType p;
	
	exp.resize(1);
	exp.back().setRT(1.0);
	p.setMZ(500.0);
	p.setIntensity(10.0);
	exp.back().push_back(p);
	p.setMZ(600.0);
	p.setIntensity(20.0);
	exp.back().push_back(p);
	p.setMZ(800.0);
	p.setIntensity(30.0);
	exp.back().push_back(p);
	p.setMZ(1000.0);
	p.setIntensity(40.0);
	exp.back().push_back(p);
	p.setMZ(1200.0);
	p.setIntensity(110.0);
	exp.back().push_back(p);
	
	exp.resize(2);
	exp.back().setRT(2.0);
	p.setMZ(500.0);
	p.setIntensity(101.0);
	exp.back().push_back(p);
	p.setMZ(600.0);
	p.setIntensity(81.0);
	exp.back().push_back(p);
	p.setMZ(800.0);
	p.setIntensity(31.0);
	exp.back().push_back(p);
	p.setMZ(1000.0);
	p.setIntensity(11.0);
	exp.back().push_back(p);
	p.setMZ(1200.0);
	p.setIntensity(111.0);
	exp.back().push_back(p);
	
	exp.updateRanges();
	
	//create dummy FF instance (for flags)
  FeatureFinder ff;
  FeatureMap<Feature> features;
  ff.run("none", exp,features, Param(), FeatureMap<>());
  
  //make two peaks used
	ff.getPeakFlag(make_pair(0,4)) = FeatureFinderDefs::USED;
	ff.getPeakFlag(make_pair(1,4)) = FeatureFinderDefs::USED;
	
	//First test (3 unused peaks with intensity greater than 35)
	{
		SimpleSeeder<Peak1D,Feature> seeder(&exp, &features, &ff);
		Param param;
		param.setValue("min_intensity",35.0);
		param.setValue("signal_to_noise",0.0);
		
		seeder.setParameters(param);
		
		FeatureFinderDefs::IndexPair peak;
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,1);
		TEST_EQUAL(peak.second,0);
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,1);
		TEST_EQUAL(peak.second,1);
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,0);
		TEST_EQUAL(peak.second,3);
	
		TEST_EXCEPTION( FeatureFinderDefs::NoSuccessor , seeder.nextSeed() )
	}
	
	
	// Second test (2 unused peaks with intensity greater than 27,75)
	{
		SimpleSeeder<Peak1D,Feature> seeder(&exp, &features, &ff);
		Param param;
		param.setValue("min_intensity",0.0);
		param.setValue("signal_to_noise",0.1);
    param.setValue("SignalToNoiseEstimationParameter:win_len",1.0);
    param.setValue("SignalToNoiseEstimationParameter:bin_count",4);
    param.setValue("SignalToNoiseEstimationParameter:min_required_elements",1);
		
		seeder.setParameters(param);
		
		FeatureFinderDefs::IndexPair peak;
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,1);
		TEST_EQUAL(peak.second,2);
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,0);
		TEST_EQUAL(peak.second,2);
			
	//	TEST_EXCEPTION( FeatureFinderDefs::NoSuccessor , seeder.nextSeed() )		
	}
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


