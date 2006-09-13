// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DBAdapter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <qapplication.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DBAdapter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// check for credentials
// if they are not present, abort the test (successfully)
bool do_tests=true;
TextFile credentials;
try
{
	credentials.load("DB_credentials.txt",true);
}
catch(...)
{
	do_tests=false;
}

String db,host,user,password,port;

//read out connection data
for (TextFile::iterator it = credentials.begin(); it!= credentials.end(); ++it)
{
	//comments and empty lines
	if (it->hasPrefix('#') || *it == "")
	{
		continue;
	}
	
	//extract connection info
	if (it->hasPrefix("Host:")) host = it->suffix(':').trim();
	if (it->hasPrefix("Port:")) port = it->suffix(':').trim();
	if (it->hasPrefix("User:")) user = it->suffix(':').trim();
	if (it->hasPrefix("Password:")) password = it->suffix(':').trim();
	if (it->hasPrefix("DB:")) db = it->suffix(':').trim();
}

if (do_tests)
{
	// start QApplication. Otherwise the DB classes are not available
	QApplication qapp(argc,argv,false);
	
	//DB connection for DBAdapter
	DBConnection con;
	con.connect(db, user, password, host, port.toInt());
	
	DBAdapter* ptr = 0;

	CHECK(DBAdapter(DBConnection& db_con))
		ptr = new DBAdapter(con);
		TEST_NOT_EQUAL(ptr, 0)
	RESULT
		
	CHECK(~DBAdapter())
		delete ptr;
	RESULT
	
	
	
	//create test data
	MSExperiment<> exp_original;
	exp_original.setComment("bla");
	//MS spectrum
	MSExperiment<>::SpectrumType spec;
	MSExperiment<>::SpectrumType::PeakType p;
	p.setIntensity(565);
	p.getPosition()[0] = 600.1;
	spec.getContainer().push_back(p);
	p.setIntensity(620);
	p.getPosition()[0] = 700.1;
	spec.getContainer().push_back(p);
	p.setIntensity(701);
	p.getPosition()[0] = 800.1;
	spec.getContainer().push_back(p);
	spec.setRetentionTime(1.98);
	spec.setMSLevel(1);
	exp_original.push_back(spec);
	
	//MSMS spectrum
	spec.getContainer().clear();
	p.setIntensity(210);
	p.getPosition()[0] = 100.155;
	spec.getContainer().push_back(p);
	p.setIntensity(101);
	p.getPosition()[0] = 150.25;
	spec.getContainer().push_back(p);
	p.setIntensity(90);
	p.getPosition()[0] = 300.5;
	spec.getContainer().push_back(p);
	spec.setRetentionTime(3.96);
	spec.setMSLevel(2);
	spec.getPrecursorPeak().getPosition()[0] = 600.1;
	spec.getPrecursorPeak().setIntensity(4711);
	spec.getPrecursorPeak().setCharge(2);
	spec.getPrecursor().setMetaValue("icon",String("Precursor"));
	spec.setComment("bla");
	exp_original.push_back(spec);
	
	
	//meta info
	exp_original.setMetaValue("label",5.55);
	exp_original.setMetaValue("icon",String("MSExperiment"));
	exp_original.setMetaValue("color",5);
	exp_original[0].setMetaValue("icon",String("Spectrum1"));
	exp_original[1].setMetaValue("icon",String("Spectrum2"));
	
	// to store the id of reading and writing
	UID tmp_id,spec_tmp_id;
	
	
	
	CHECK(template <class ExperimentType> void DBAdapter::storeExperiment(ExperimentType& exp))
	  DBAdapter a(con);
	  a.storeExperiment(exp_original);
		tmp_id = exp_original.getPersistenceId();
		spec_tmp_id = exp_original[0].getPersistenceId();
	RESULT		

	CHECK(template <class SpectrumType> void DBAdapter::loadSpectrum(UID id, SpectrumType& spec))
  	DBAdapter a(con);
	  
		MSSpectrum<> spec;
		a.loadSpectrum(spec_tmp_id, spec);
					
	  TEST_EQUAL( spec.getRetentionTime() , exp_original.begin()->getRetentionTime() )
		TEST_EQUAL( spec.getMSLevel() , exp_original.begin()->getMSLevel() )
		TEST_EQUAL( spec.size() , exp_original.begin()->size() )
		for (UnsignedInt i=0; i<3; ++i)
		{
			TEST_REAL_EQUAL( spec.getContainer()[i].getIntensity() , exp_original.begin()->getContainer()[i].getIntensity() )
			TEST_REAL_EQUAL( spec.getContainer()[i].getPosition()[0] , exp_original.begin()->getContainer()[i].getPosition()[0] )
		}
	RESULT
	
  CHECK(template <class ExperimentType> void DBAdapter::loadExperiment(UID id, ExperimentType& exp))
	  DBAdapter a(con);
	  MSExperiment<> exp_new;
		
		a.loadExperiment(tmp_id, exp_new);
		TEST_EQUAL(exp_new.getPersistenceId(), tmp_id)
		TEST_EQUAL(exp_new.getComment() , "bla" )
		
		//------ test if values are correct ------
		
		//SPECTRUM 1
		MSExperiment<>::const_iterator itn(exp_new.begin());
		MSExperiment<>::const_iterator ito(exp_original.begin());
			
	  TEST_EQUAL( itn->getRetentionTime() , ito->getRetentionTime() )
		TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
		TEST_EQUAL( itn->size() , ito->size() )
		for (UnsignedInt i=0; i<3; ++i)
		{
			TEST_REAL_EQUAL( itn->getContainer()[i].getIntensity() , ito->getContainer()[i].getIntensity() )
			TEST_REAL_EQUAL( itn->getContainer()[i].getPosition()[0] , ito->getContainer()[i].getPosition()[0] )
		}
	
		//SPECTRUM 2
		++itn;
		++ito;
			
	  TEST_EQUAL( itn->getRetentionTime() , ito->getRetentionTime() )
		TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
		TEST_EQUAL( itn->getPrecursorPeak().getPosition()[0] , ito->getPrecursorPeak().getPosition()[0] )
		TEST_EQUAL( itn->getPrecursorPeak().getIntensity() , ito->getPrecursorPeak().getIntensity() )
		TEST_EQUAL( itn->getPrecursorPeak().getCharge() , ito->getPrecursorPeak().getCharge() )
		TEST_EQUAL( itn->getPrecursor().getMetaValue("icon") , "Precursor" )
		TEST_EQUAL( itn->getComment() , "bla" )
		TEST_EQUAL( itn->size() , ito->size() )
		for (UnsignedInt i=0; i<3; ++i)
		{
			TEST_REAL_EQUAL( itn->getContainer()[i].getIntensity() , ito->getContainer()[i].getIntensity() )
			TEST_REAL_EQUAL( itn->getContainer()[i].getPosition()[0] , ito->getContainer()[i].getPosition()[0] )
		}
		
		//META INFO
		TEST_REAL_EQUAL((double)exp_new.getMetaValue("label"),5.55)
		TEST_EQUAL((string)exp_new.getMetaValue("icon"),"MSExperiment")
		TEST_EQUAL((int)exp_new.getMetaValue("color"),5)
		TEST_EQUAL((string)exp_new[0].getMetaValue("icon"),"Spectrum1")
		TEST_EQUAL((string)exp_new[1].getMetaValue("icon"),"Spectrum2")
	RESULT

	CHECK([EXTRA] load and store of empty map)
  	DBAdapter a(con);
	  MSExperiment<> in, out;
	  a.storeExperiment(in);
		a.loadExperiment(in.getPersistenceId(),out);
		TEST_EQUAL(in==out, true)
	RESULT

}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



