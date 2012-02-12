// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/DataProcessing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DataProcessing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DateTime time;
time.set("2000-10-09 08:07:40");

DataProcessing* ptr = 0;
DataProcessing* nullPointer = 0;
START_SECTION(DataProcessing())
	ptr = new DataProcessing();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~DataProcessing())
	delete ptr;
END_SECTION

START_SECTION(const DateTime& getCompletionTime() const)
  DataProcessing tmp;
  TEST_EQUAL(tmp.getCompletionTime().get(),"0000-00-00 00:00:00");
END_SECTION

START_SECTION(void setCompletionTime(const DateTime& completion_time))
  DataProcessing tmp;
  tmp.setCompletionTime(time);
  TEST_EQUAL(tmp.getCompletionTime()==time,true);
END_SECTION

START_SECTION(Software& getSoftware())
  DataProcessing tmp;
  TEST_EQUAL(tmp.getSoftware()==Software(),true)
END_SECTION

START_SECTION(const Software& getSoftware() const)
  DataProcessing tmp;
  tmp.getSoftware().setName("name");
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"name")
END_SECTION

START_SECTION(void setSoftware(const Software& software))
  DataProcessing tmp;
  Software tmp2;
  tmp2.setName("name");
  tmp.setSoftware(tmp2);
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"name")
END_SECTION

START_SECTION(const std::set<ProcessingAction>& getProcessingActions() const)
  DataProcessing tmp;
  TEST_EQUAL(tmp.getProcessingActions().size(),0)
END_SECTION

START_SECTION(std::set<ProcessingAction>& getProcessingActions())
  DataProcessing tmp;
  tmp.getProcessingActions().insert(DataProcessing::DEISOTOPING);
  TEST_EQUAL(tmp.getProcessingActions().size(),1)
END_SECTION

START_SECTION(void setProcessingActions(const std::set<ProcessingAction>& actions))
  DataProcessing tmp;
  std::set<DataProcessing::ProcessingAction> tmp2;
  tmp2.insert(DataProcessing::DEISOTOPING);
  tmp2.insert(DataProcessing::CHARGE_DECONVOLUTION);
  tmp.setProcessingActions(tmp2);
  TEST_EQUAL(tmp.getProcessingActions().size(),2)
END_SECTION


START_SECTION(DataProcessing& operator= (const DataProcessing& source))
  DataProcessing tmp;
  tmp.setCompletionTime(time);
  tmp.getProcessingActions().insert(DataProcessing::DEISOTOPING);
  tmp.getSoftware().setName("name");
  tmp.setMetaValue("label",String("label"));

  DataProcessing tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getCompletionTime()==time,true);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getProcessingActions().size(),1)
  TEST_STRING_EQUAL(tmp2.getSoftware().getName(),"name")
END_SECTION

START_SECTION(DataProcessing(const DataProcessing& source))
  DataProcessing tmp;
  tmp.setCompletionTime(time);
  tmp.getProcessingActions().insert(DataProcessing::DEISOTOPING);
  tmp.getSoftware().setName("name");
  tmp.setMetaValue("label",String("label"));
  
  DataProcessing tmp2(tmp);
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getCompletionTime()==time,true);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getProcessingActions().size(),1)
  TEST_STRING_EQUAL(tmp2.getSoftware().getName(),"name")
END_SECTION

START_SECTION(bool operator== (const DataProcessing& rhs) const)
  DataProcessing edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
  edit.setCompletionTime(time);
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.getProcessingActions().insert(DataProcessing::DEISOTOPING);
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.getSoftware().setName("name");
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty, false);
END_SECTION

START_SECTION(bool operator!= (const DataProcessing& rhs) const)
  DataProcessing edit, empty;
  
  TEST_EQUAL(edit!=empty, false);
  
  edit.setCompletionTime(time);
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.getProcessingActions().insert(DataProcessing::DEISOTOPING);
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.getSoftware().setName("name");
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty, true);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



