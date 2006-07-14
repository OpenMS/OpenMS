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
	p.setIntensity(90);
	p.getPosition()[0] = 300.5;
	spec.getContainer().push_back(p);
	p.setIntensity(101);
	p.getPosition()[0] = 150.25;
	spec.getContainer().push_back(p);
	p.setIntensity(210);
	p.getPosition()[0] = 100.155;
	spec.getContainer().push_back(p);
	spec.setRetentionTime(3.96);
	spec.setMSLevel(2);
	spec.getPrecursorPeak().getPosition()[0] = 600.1;
	spec.getPrecursorPeak().setIntensity(4711);
	spec.getPrecursorPeak().setCharge(2);
	exp_original.push_back(spec);

	// to store the id of reading and writing
	UID tmp_id,tmp_id2,spec_tmp_id;
	
CHECK([EXTRA] operator<<)
	  DBAdapter a(con);
	  a << exp_original;
		tmp_id = exp_original.getPersistenceId();
		//test a second time to check if initialization and cleanup are complete
		a << exp_original;
		tmp_id2 = exp_original.getPersistenceId();
		spec_tmp_id = exp_original.begin()->getPersistenceId();
	RESULT
	
CHECK([EXTRA] operator>>)
	  DBAdapter a(con);
	  PersistentObject* exp_new(0);
		
		a.readId(tmp_id);
		a >> exp_new;
		TEST_NOT_EQUAL(exp_new,0)
		TEST_NOT_EQUAL(dynamic_cast< MSExperiment<>* >(exp_new),0)
		TEST_EQUAL(exp_new->getPersistenceId(), tmp_id)
		
		
		//------ test if values are correct ------
		
		//SPECTRUM 1
		MSExperiment<>::const_iterator itn((dynamic_cast< MSExperiment<>* >(exp_new))->begin());
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
		TEST_EQUAL( itn->size() , ito->size() )
		for (UnsignedInt i=0; i<3; ++i)
		{
			TEST_REAL_EQUAL( itn->getContainer()[i].getIntensity() , ito->getContainer()[i].getIntensity() )
			TEST_REAL_EQUAL( itn->getContainer()[i].getPosition()[0] , ito->getContainer()[i].getPosition()[0] )
		}
		
		
		//  ---- test a second time to check if initialization and cleanup are complete ------
		
		exp_new = 0;
		a.readId(tmp_id2);
		a >> exp_new;
		TEST_NOT_EQUAL(exp_new,0)
		TEST_NOT_EQUAL(dynamic_cast< MSExperiment<>* >(exp_new),0)
		TEST_EQUAL(exp_new->getPersistenceId(), tmp_id2)
		//-- test if values are correct --
		//SPECTRUM 1
		itn = (dynamic_cast< MSExperiment<>* >(exp_new))->begin();
		ito = exp_original.begin();
			
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
		TEST_EQUAL( itn->size() , ito->size() )
		for (UnsignedInt i=0; i<3; ++i)
		{
			TEST_REAL_EQUAL( itn->getContainer()[i].getIntensity() , ito->getContainer()[i].getIntensity() )
			TEST_REAL_EQUAL( itn->getContainer()[i].getPosition()[0] , ito->getContainer()[i].getPosition()[0] )
		}
	RESULT
	
CHECK(MSExperiment<>* loadMSExperiment(UID id))
		DBAdapter a(con);
	  
		MSExperiment<>* ptr;
		ptr = a.loadMSExperiment(tmp_id);
		
		TEST_EQUAL(ptr->size(),2)
		TEST_EQUAL(ptr->begin()->getMSLevel(),1)
		TEST_EQUAL((++(ptr->begin()))->getMSLevel(),2)
	RESULT
	
CHECK(MSSpectrum<>* loadSpectrum(UID id))
		DBAdapter a(con);
	  
		MSSpectrum<>* ptr;
		ptr = a.loadSpectrum(spec_tmp_id);
					
	  TEST_EQUAL( ptr->getRetentionTime() , exp_original.begin()->getRetentionTime() )
		TEST_EQUAL( ptr->getMSLevel() , exp_original.begin()->getMSLevel() )
		TEST_EQUAL( ptr->size() , exp_original.begin()->size() )
		for (UnsignedInt i=0; i<3; ++i)
		{
			TEST_REAL_EQUAL( ptr->getContainer()[i].getIntensity() , exp_original.begin()->getContainer()[i].getIntensity() )
			TEST_REAL_EQUAL( ptr->getContainer()[i].getPosition()[0] , exp_original.begin()->getContainer()[i].getPosition()[0] )
		}
	RESULT

}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



