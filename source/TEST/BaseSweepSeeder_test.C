// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>


#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/MSExperimentExtern.h>

///////////////////////////

START_TEST(BaseSweepSeeder, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

typedef MSExperimentExtern< >::SpectrumType SpectrumType;

class TestSweepSeeder : public BaseSweepSeeder
{
  public:
	
	/// charge state estimate with associated score
	typedef BaseSweepSeeder::ScoredChargeType ScoredChargeType;
	/// m/z position in spectrum with charge estimate and score
	typedef BaseSweepSeeder::ScoredMZType ScoredMZType;
	/// container of scored m/z positions
	typedef BaseSweepSeeder::ScoredMZVector ScoredMZVector;
	
	TestSweepSeeder()
		: BaseSweepSeeder()
	{
		setName(getProductName());
		
		check_defaults_ = false;
		
		defaultsToParam_();
	}

	TestSweepSeeder(const TestSweepSeeder& source)
		: BaseSweepSeeder(source)
	{
		updateMembers_();
	}
	
	virtual TestSweepSeeder& operator = (const TestSweepSeeder& source)
	{
		if (&source == this) return *this;
		
		BaseSweepSeeder::operator = (source);
		updateMembers_();
		
		return *this;
	}
	
	void updateMembers_()
	{
		BaseSweepSeeder::updateMembers_();
	}

	// dummy implementation 
	ScoredMZVector detectIsotopicPattern_(SpectrumType& scan)
	{
		ScoredMZVector scmzvec;
		
		Int count = 0;
		
		for (SpectrumType::const_iterator citer = scan.begin();
					citer != scan.end();
					++citer)
		{
			// charge estimate 2 and score 3.0 for everyone
			ScoredChargeType scc1;
			scc1.first = 2; scc1.second = 3.0;
			
			ScoredMZType scmz1;
			scmz1.first = count++; scmz1.second = scc1;				
			scmzvec.push_back( scmz1 );
		}	
	
		return scmzvec;
	}

	static TestSweepSeeder* create()
	{
		return new TestSweepSeeder();
	}
	
	static const String getProductName()
	{
		return "TestSweepSeeder";
	}

};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
TestSweepSeeder* ptr = 0;
CHECK((BaseSweepSeeder()))
	ptr = new TestSweepSeeder();
  TEST_EQUAL(ptr->getName(), "TestSweepSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~BaseSweepSeeder()))
	delete ptr;
RESULT

// assignment operator
CHECK((virtual BaseSweepSeeder& operator=(const BaseSweepSeeder &source)))
	TestSweepSeeder tss1;
  TestSweepSeeder tss2;
	
	Param p;
	p.setValue("scans_to_sumup",10);
  tss1.setParameters(p);
	
  tss1 = tss2;
	TEST_EQUAL(tss1 == tss2,true)
RESULT

// copy constructor
CHECK((BaseSweepSeeder(const BaseSweepSeeder &source)))
	TestSweepSeeder tss1;	
  
	Param p;
	p.setValue("scans_to_sumup",10);
  tss1.setParameters(p);
	
  TestSweepSeeder tss2(tss1);
	TEST_EQUAL(tss1 == tss2,true)
RESULT

CHECK((virtual IndexSet nextSeed()))
	TestSweepSeeder tss1;
	Param p;
	p.setValue("scans_to_sumup",0);
	p.setValue("min_number_scans",0);
	p.setValue("min_number_peaks",0);
	tss1.setParameters(p);
	
	FeaFiTraits traits;
	
	MSExperimentExtern< > exp;
	SpectrumType spec;
	Peak1D peak;
	
	for (UInt j=1;j<=5;++j)
	{
		for (UInt i=1;i<=5;++i)
		{
			peak.setMZ(i);
			peak.setIntensity(i);
		
			spec.push_back(peak);		
		}
		spec.setRT(j);
		exp.push_back(spec);		
		spec.clear();	
	}
	
	traits.setData(exp.begin(),exp.end(),10);
	
	tss1.setTraits( &traits );

	FeaFiModule::IndexSet set;
	
	for (UInt mz = 0;mz<5;++mz)
	{
		set = tss1.nextSeed();
		UInt rt = 0;
		for (FeaFiModule::IndexSet::const_iterator citer = set.begin();
				citer != set.end();
				++citer)
		{
			TEST_EQUAL(citer->first,rt)
			TEST_EQUAL(citer->second,mz)
			++rt;
		}	
	}
	
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
